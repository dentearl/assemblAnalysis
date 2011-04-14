#!/usr/bin/env python
"""
createContigPathStatsTable.py
11 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the contig path stats table.

output is latex

"""
import createIndividualSection as cis
import glob
import libAssemblySubset as las
from optparse import OptionParser
import os
import re
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

class Assembly:
   def __init__( self ):
      self.ID          = ''
      self.valuesDict  = {}
      self.totalErrors = 0

def initOptions( parser ):
   parser.add_option( '--contigPathStatsDir', dest='contigPathStatsDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.contigPathStats.xml .'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='totalErrors',
                      help=('Column to sort the table on. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks',
                      action='store_true', default=False,
                      help=('Prints rankings as tab delimited data to STDOUT.'))
   parser.add_option( '--printAllowedKeys', dest='printAllowedKeys',
                      action='store_true', default=False,
                      help=('Prints out the allowed keys for --sortOn and exits.'))

def checkOptions( args, options, parser ):
   allowedKeys = set([ 'errorsPerContig', 'errorsPerMappedBase',
                       'coverage', 'totalHomoToHeteroSwitches',
                       'totalScaffoldGaps+totalBleedingHeartScaffoldGaps',
                       'totalContigEnds+totalContigEndsWithNs',
                       'totalContigEndsWithInsert',
                       'totalErrorsHaplotypeToHaplotype',
                       'totalErrorsHaplotypeToInsert',
                       'totalErrorsHaplotypeToDeletion',
                       'totalErrorsNonSpecific',
                       'totalErrors',
                       'scaffolds' ])
   if options.printAllowedKeys:
      for k in allowedKeys:
         print k
      sys.exit( 0 )
   dirs = { 'contigPathStatsDir'   : options.contigPathStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('Error, specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('Error, --%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('Error, --%s %s is not a directory!\n' % ( d, dirs[ d ]) )
   if options.sortOn not in allowedKeys:
      parser.error( 'Error, --sortOn %s is not in dict of allowed keys:\n %s.' % ( options.sortOn, allowedKeys ))

def readDir( options ):
   aFiles = glob.glob( os.path.join( options.contigPathStatsDir, '*.contigPathStats.xml'))
   namepat = re.compile( '^(\S{2,3})\.contigPathStats.xml' )
   assembliesList = []
   for f in aFiles:
      name = re.match( namepat, os.path.basename( f )).group( 1 )
      if options.subsetFile:
         if name not in options.assemblySubset:
            continue
      try:
         xmlTree = ET.parse( f )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( f )
      root=xmlTree.getroot()
      a = Assembly()
      a.ID = name
      for elm in root.attrib.keys():
         if elm in ('errorsPerContig', 'errorsPerMappedBase', 'coverage'):
            a.valuesDict[ elm ] = float( root.attrib[ elm ] )            
         elif elm in ( 'insertionErrorSizeDistribution', 'deletionErrorSizeDistribution' ):
            a.valuesDict[ elm ] = root.attrib[ elm ].split()
         else:
            a.valuesDict[ elm ] = int( root.attrib[ elm ] )
      assembliesList.append( a )
   return assembliesList

def prettyNumber( n ):
   if isinstance( n, float ):
      return cis.prettyFloat( float(n), 2 )
   elif isinstance( n, int ):
      return cis.prettyInt( int(n) )

def printTable( assembliesList, caption, options ):
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{FPtable}
\caption[Scaffold path statistics.]{Scaffold path statistics. %s}
\\tiny
\\centering
\\begin{tabular}{| r | p{.75in} | p{.75in} | p{.4in} | p{.4in} | p{.5in} | p{.5in} || c |}
\\hline
ID & Haplotype to haplotype, intra chromosomal & Haplotype to haplotype inter chromosomal & Haplotype to insert & haplotype to deletion & haplotype to insert and deletion & contig ends in insert & \(\sum\) errors \\\\
\\hline
\\hline''' % ( caption ) # hom/het sw & scf & ctgE+ctgN 
   i = 0
   for row in assembliesList:
      i += 1
      sys.stdout.write( '%s' % ( row.ID ))
      # sys.stdout.write( ' & %s & %s' % ( prettyNumber(row.valuesDict['totalHomoToHeteroSwitches']),
      #                                    prettyNumber( row.valuesDict['totalScaffoldGaps'] + 
      #                                                  row.valuesDict['totalBleedingHeartScaffoldGaps'] ) ))
      # sys.stdout.write( ' & %s' % ( prettyNumber( row.valuesDict['totalContigEnds'] + 
      #                                             row.valuesDict['totalContigEndsWithNs'] )))
      #sys.stdout.write( ' & %.4f' % ( row.valuesDict['errorsPerContig']))
      #sys.stdout.write( ' & %.2e' % ( row.valuesDict['errorsPerMappedBase']))
      for v in [ 'totalErrorsHaplotypeToHaplotypeSameChromosome',
                 'totalErrorsHaplotypeToHaplotypeDifferentChromosome',
                 'totalErrorsHaplotypeToInsertion',
                 'totalErrorsHaplotypeToDeletion',
                 'totalErrorsHaplotypeToInsertionAndDeletion',
                 'totalErrorsContigEndsWithInsert' ]:
         sys.stdout.write( ' & %s' % ( prettyNumber( row.valuesDict[v] )))
      #sys.stdout.write( ' & %s' % ( prettyNumber( row.valuesDict['totalErrorsNonSpecific'])))
      sys.stdout.write( ' & %s \\\\\n' % prettyNumber( row.totalErrors ))
      if not i % 10 and i != len( assembliesList ):
         print '\\hline'
   print '''\\hline
\\end{tabular}
\\label{table:contigPathStats}
\\end{FPtable}\par
\\normalsize
\\vspace{0.3in}'''

def printRanks( assembliesList, options ):
   print '#Assembly\tHap-Hap-sameChr\tHap-Hap-diffChr\tHap-Ins\tHap-Del\tHap-Ins&Del\tCtgEndsIns\tSumErr'
   for row in assembliesList:
      sys.stdout.write( '%s' % ( row.ID ))
      for v in [ 'totalErrorsHaplotypeToHaplotypeSameChromosome',
                 'totalErrorsHaplotypeToHaplotypeDifferentChromosome',
                 'totalErrorsHaplotypeToInsertion',
                 'totalErrorsHaplotypeToDeletion',
                 'totalErrorsHaplotypeToInsertionAndDeletion',
                 'totalErrorsContigEndsWithInsert' ]:
         sys.stdout.write( '\t%s' % ( row.valuesDict[v] ))
      #sys.stdout.write( ' & %s' % ( prettyNumber( row.valuesDict['totalErrorsNonSpecific'])))
      sys.stdout.write( '\t%s\n' % row.totalErrors )

def calculateErrors( assembliesList, options ):
   for a in assembliesList:
      for elm in [ 'totalErrorsHaplotypeToHaplotypeSameChromosome',
                   'totalErrorsHaplotypeToHaplotypeDifferentChromosome',
                   'totalErrorsHaplotypeToContamination',
                   'totalErrorsHaplotypeToInsertionToContamination',
                   'totalErrorsHaplotypeToInsertion',
                   'totalErrorsHaplotypeToDeletion',
                   'totalErrorsHaplotypeToInsertionAndDeletion',
                   'totalErrorsContigEndsWithInsert' ]:
         a.totalErrors += a.valuesDict[ elm ]
      if a.totalErrors != a.valuesDict[ 'totalErrors' ]:
         sys.stderr.write('Error, calculated sum of total errors not equal to value from xml! %s\n' % a.ID )

def performSort( assembliesList, options ):
   if options.sortOn == 'totalErrors':
      assembliesList = sorted( assembliesList, key=lambda x: x.totalErrors, reverse=False )
   elif ( options.sortOn == 'totalScaffoldGaps+totalBleedingHeartScaffoldGaps' or 
          options.sortOn == 'scaffolds' ):
      assembliesList = sorted( assembliesList, key=lambda x: x.valuesDict[totalScaffoldGaps] 
                               + x.valuesDict[totalBleedingHeartScaffoldGaps], reverse=False )
   elif options.sortOn == 'totalContigEnds+totalContigEndsWithNs':
      assembliesList = sorted( assembliesList, key=lambda x: x.valuesDict[totalContigEnds]
                               + x.valuesDict[totalContigEndsWithNs], reverse=False )
   else:
      assembliesList = sorted( assembliesList, key=lambda x: 
                               x.valuesDict[ options.sortOn ], reverse=False )
   return assembliesList

def main():
   usage = ( 'usage: %prog --contigPathStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes in the contig path statistics directory\n'
             '( --contigPathStatsDir ) and prints to STDOUT a latex formatted table.' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   las.checkOptions( options, parser )
   
   assembliesList = readDir( options )
   calculateErrors( assembliesList, options )

   assembliesList = performSort( assembliesList, options )
   
   sortString = { 'totalErrors':'the sum of the errors (\\(\\sum\\) err)',
                  'errorsPerMappedBase':'errors per mapped base (ePmb)',
                  'totalScaffoldGaps+totalBleedingHeartScaffoldGaps':'the sum of the total number of scaffold gaps and liberal scaffold gaps (scf)',
                  'scaffolds':'the sum of the total number of scaffold gaps and liberal scaffold gaps (scf)',
                  'totalContigEnds+totalContigEndsWithNs':'the sum of total contig ends and the total contigs ending with N\'s (ctgE+ctgN)',
                  'totalHomoToHeteroSwitches':'the total homozygous to heterozygous switches (hom/het sw)',
                  'errorsPerContig':'errors per contig (ePc)',
                  'totalContigEndsWithInsert':'total contigs ending with insert (e-i)',
                  'totalErrorsHaplotypeToHaplotype':'total haplotype to haplotype errors (h-h)',
                  'totalErrorsHaplotypeToInsert':'total haplotype to insert errors (h-i)',
                  'totalErrorsHaplotypeToDeletion':'total haplotype to deletion errors (h-d)',
                  'totalErrorsNonSpecific':'total non--specific errors (nonSpec)'
                  }

   caption = 'The table is sorted on %s. Column headers are total homozygous to heterozygous switches (hom/het sw), sum of the total scaffold gaps and liberal scaffold gaps (scf), the sum of total contig ends and the total contigs ending with N\'s (ctgE+ctgN), errors per contig (ePc), errors per mapped base (ePmb), total contigs ending with insert (e-i), total haplotype to haplotype errors (h-h), total haplotype to insert errors (h-i), total haplotype to deletion errors (h-d), total non--specific errors (nonSpec) and sum of the errors (\\(\\sum\\) err).' % sortString[ options.sortOn ]
   
   if options.outputRanks:
      printRanks( assembliesList, options )
   else:
      printTable( assembliesList, caption, options )

if __name__ == '__main__':
   main()
