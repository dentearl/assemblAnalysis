#!/usr/bin/env python
"""
"""
from libGeneral import prettyNumber
import glob
from optparse import OptionParser
import os
import re
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

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
                      help=('Directory with subStats. Names: A1.subStats.upper.xml .'))
   parser.add_option( '--order', dest='order',
                      type='string',
                      help=('Order (left-right, top-bottom) of plots, comma '
                            'separated. Names must match file prefixes in the --dir.' ))

def checkOptions( args, options, parser ):
   dirs = { 'subStatsDir' : options.subStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   if options.order != None:
      options.order = options.order.split(',')
   else:
      options.order = []

def readSubStatsDir( assembliesDict, options ):
   lowerStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.lower.xml') )
   upperStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.upper.xml') )
   
   namereg = '^([A-Z0-9]{2,3})\.subStats.*'
   namepat = re.compile( namereg  )
   for l in lowerStatsFiles:
      m = re.match( namepat, os.path.basename( l ))
      if not m:
         sys.stderr.write('unable to match regex "%s" against filename "%s"' % ( namereg, l ))
         sys.exit( 1 )
      ID = m.group( 1 )
      assembliesDict[ ID ] = Assembly()
      assembliesDict[ ID ].ID = ID
      try:
         xmlTree = ET.parse( l )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( l )
      root=xmlTree.getroot()
      assembliesDict[ ID ] = Assembly()
      assembliesDict[ ID ].ID = ID
      for elm in root.attrib.keys():
         assembliesDict[ ID ].subStatsLower[ elm ] = int(float( root.attrib[ elm ]))
   for u in upperStatsFiles:
      m = re.match( namepat, os.path.basename( u ))
      if not m:
         sys.stderr.write('unable to match regex "%s" against filename "%s"' % ( namepat, u ))
         sys.exit( 1 )
      ID = m.group( 1 )
      if ID not in assembliesDict:
         sys.stderr.write('unable to locate key %s in assembliesDict.\n')
         sys.exit( 1 )
      try:
         xmlTree = ET.parse( u )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( u )
      root=xmlTree.getroot()
      for elm in root.attrib.keys():
         assembliesDict[ ID ].subStatsUpper[ elm ] = int(float( root.attrib[ elm ]))
   return assembliesDict



def printLine( a, kind ):
   sys.stdout.write( '%s' % a.ID )
   if kind == 'hom':
      for k in [ 'totalCallsInHomozygous', 
                 'totalCorrectInHomozygous', 'totalErrorsInHomozygous' ]:
         sys.stdout.write( ' & %s -- %s' % ( prettyNumber( a.subStatsLower[ k ]), prettyNumber( a.subStatsUpper[ k ])))
      sys.stdout.write( ' \\\\\n' )
   elif kind == 'het':
      for k in [ 'totalCallsInHeterozygous', 
                 'totalCorrectInHeterozygous', 'totalErrorsInHeterozygous' ]:
         sys.stdout.write( ' & %s -- %s' % ( prettyNumber( a.subStatsLower[ k ]), prettyNumber( a.subStatsUpper[ k ])))
      sys.stdout.write( ' \\\\\n' )
   elif kind == 'indel':
      for k in [ 'totalCallsInOneHaplotypeOnly', 
                 'totalCorrectInOneHaplotypeOnly', 'totalErrorsInOneHaplotypeOnly' ]:
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
             'names as NAME.subStats.[upper|lower].xml and writes to STDOUT a latex formated table.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )
   
   assembliesDict = {}
   assembliesDict = readSubStatsDir( assembliesDict, options )
   
   printTables( assembliesDict, options )

if __name__ == '__main__':
   main()
