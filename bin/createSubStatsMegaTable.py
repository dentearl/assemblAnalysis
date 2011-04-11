#!/usr/bin/env python
"""
"""
import createIndividualSection as cis
import glob
from optparse import OptionParser
import os
import re
import sys

class Assembly:
   """ Assembly objects are generated from lines 
   in the two substitution summary files, lower and upper
   """
   def __init__( self ):
      self.ID    = ''
      self.subStatsLower = {}
      self.subStatsUpper = {}

def initOptions( parser ):
   parser.add_option( '--subStatsDir', dest='subStatsDir',
                      type='string',
                      help=('Directory with subStats. Names: A1.subStats.upper.txt .'))
   parser.add_option( '--order', dest='order',
                      type='string',
                      help=('Order (left-right, top-bottom) of plots, comma '
                            'separated. Names must match file prefixes in the --dir.' ))

def checkOptions( args, options, parser ):
   dirs = { 'subStatsDir' : options.subStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('Error, specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('Error, --%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('Error, --%s %s is not a directory!\n' % (d, dirs[ d ]) )
   if options.order != None:
      options.order = options.order.split(',')
   else:
      options.order = []

def readSubStatsDir( assembliesDict, options ):
   lowerStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.lower.txt') )
   upperStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.upper.txt') )
   
   namereg = '^([A-Z0-9]{2,3})\.subStats.*'
   namepat = re.compile( namereg  )
   for l in lowerStatsFiles:
      m = re.match( namepat, os.path.basename( l ))
      if not m:
         sys.stderr.write('Error, unable to match regex "%s" against filename "%s"' % ( namereg, l ))
         sys.exit( 1 )
      ID = m.group( 1 )
      assembliesDict[ ID ] = Assembly()
      assembliesDict[ ID ].ID = ID
      f = open( l, 'r' )
      for line in f:
         line = line.strip()
         d = line.split('\t')
         assembliesDict[ ID ].subStatsLower[ d[0] ] = d[ 1 ]
      f.close()
   for u in upperStatsFiles:
      m = re.match( namepat, os.path.basename( u ))
      if not m:
         sys.stderr.write('Error, unable to match regex "%s" against filename "%s"' % ( namepat, u ))
         sys.exit( 1 )
      ID = m.group( 1 )
      if ID not in assembliesDict:
         sys.stderr.write('Error, unable to locate key %s in assembliesDict.\n')
         sys.exit( 1 )
      f = open( u, 'r' )
      for line in f:
         line = line.strip()
         d = line.split('\t')
         assembliesDict[ ID ].subStatsUpper[ d[0] ] = d[ 1 ]
      f.close()
   return assembliesDict

def prettyNumber( n ):
   pat = re.compile( '^[0-9]+\.[0-9]+$' )
   if re.match( pat, n ):
      return cis.prettyFloat( float(n), 2 )
   else:
      return cis.prettyInt( int(n) )

def printLine( a, kind ):
   sys.stdout.write( '%s' % a.ID )
   if kind == 'hom':
      for k in [ 'Total-calls-in-homozygous', 
                 'Total-correct-in-homozygous', 'Total-errors-in-homozygous' ]:
         sys.stdout.write( ' & %s -- %s' % ( prettyNumber( a.subStatsLower[ k ]), prettyNumber( a.subStatsUpper[ k ])))
      sys.stdout.write( ' \\\\\n' )
   elif kind == 'het':
      for k in [ 'Total-calls-in-heterozygous', 
                 'Total-correct-in-heterozygous', 'Total-errors-in-heterozygous' ]:
         sys.stdout.write( ' & %s -- %s' % ( prettyNumber( a.subStatsLower[ k ]), prettyNumber( a.subStatsUpper[ k ])))
      sys.stdout.write( ' \\\\\n' )
   elif kind == 'indel':
      for k in [ 'Total-calls-in-one-haplotype-only', 
                 'Total-correct-in-one-haplotype-only', 'Total-errors-in-one-haplotype-only' ]:
         sys.stdout.write( ' & %s -- %s' % ( prettyNumber( a.subStatsLower[ k ]), prettyNumber( a.subStatsUpper[ k ])))
      sys.stdout.write( ' \\\\\n' )

def printTables( assembliesDict, options ):
   tables = { 'Homozygous':'hom',
              'Heterozygous':'het',
              'Indel':'indel'}
   genericCaption = 'All cells in the table have an upper and a lower value, see main text section \\ref{sect:SubErrors} for details. The Calls column represents the total number of valid columns. The Correct column is a bit-score to score correct but ambiguous matches (IUPAC ambiguity characters) within valid columns. The Errors column represents number of positions where the assembly sequence\'s position has an IUPAC character that is not one of the two haplotypes\'s position\'s base pair(s).'
   captions = { 'Homozygous':'Although we do not allow structural rearrangements within blocks, blocks are tolerant of substitutions. Homozygous columns are those containing a member of both haplotypes, both of which have the same basepair. %s' % genericCaption,
                'Heterozygous':'Although we do not allow structural rearrangements within blocks, blocks are tolerant of substitutions. Heterozygous columns are those containing a member of both haplotypes, but which have distinct basepairs. %s' % genericCaption,
                'Indel':'Although we do not allow structural rearrangements within blocks, blocks are tolerant of substitutions. Indel columns are those containing a member of either haplotype, but not both. %s' % genericCaption }
   for t in [ 'Homozygous', 'Heterozygous', 'Indel' ]:
      printTable( t, captions[ t ], tables[ t ], assembliesDict, options )

def printTable( name, caption, key, assembliesDict, options ):
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{FPtable}
\\caption[Subsitution statistics table, %s columns.]{Subsitution statistics table, %s columns. %s}
\\scriptsize
\\centering
\\begin{tabular}{| r | c | c | c |}
\\hline
Assembly & Calls & Correct (bits) & Errors \\\\
\\hline
\\hline''' % ( name, name, caption )
   if len( options.order ) == 0:
      options.order = sorted( assembliesDict.keys(), key=lambda x: x[0], reverse=False )
      options.order = sorted( options.order, key=lambda x: int( x[1:] ), reverse=False )
   i = 0
   for a in options.order:
      i += 1
      printLine( assembliesDict[ a ], key )
      if not i % 10:
         print '\\hline'
   print '''\\hline
\\end{tabular}
\label{table:Subs%s}
\\end{FPtable}\par
\\normalsize
\\vspace{0.3in}''' % name
   
def main():
   usage = ( 'usage: %prog --subStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of substitution stats files ( --subStatsdir ) with\n'
             'names as NAME.subStats.[upper|lower].txt and writes to STDOUT a latex formated table.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   
   assembliesDict = {}
   assembliesDict = readSubStatsDir( assembliesDict, options )
   
   printTables( assembliesDict, options )

if __name__ == '__main__':
   main()
