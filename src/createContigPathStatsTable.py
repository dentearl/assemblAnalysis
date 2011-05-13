#!/usr/bin/env python
"""
createContigPathStatsTable.py
11 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the contig path stats table.

output is latex

"""
##############################
# Copyright (C) 2009-2011 by 
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedict.paten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
# ... and other members of the Reconstruction Team of David Haussler's 
# lab (BME Dept. UCSC).
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
import glob
import libAssemblySubset as las
import libGeneral as lgn
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
   parser.add_option( '--statsScaffoldsContigPathDir', dest='statsScaffoldsContigPathDir',
                      type='string',
                      help=('Directory with contigPathStats from Scaffolds alignment. '
                            'Names: A1.contigPathStats.xml .'))
   parser.add_option( '--statsContigsContigPathDir', dest='statsContigsContigPathDir',
                      type='string',
                      help=('Directory with contigPathStats from Contigs alignment. '
                            'Names: A1.contigPathStats.xml .'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='totalErrors',
                      help=('Column to sort the table on. default=%default'))
   parser.add_option( '--showAssemblyNumbers', dest='showAssemblyNumbers', default=False,
                      action='store_true',
                      help=('Shows the intra-team assembly number next to the name. default=%default'))
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
   dirs = { 'statsScaffoldsContigPathDir' : options.statsScaffoldsContigPathDir,
            'statsContigsContigPathDir'   : options.statsContigsContigPathDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % ( d, dirs[ d ]) )
   if options.sortOn not in allowedKeys:
      parser.error( '--sortOn %s is not in dict of allowed keys:\n %s.' % ( options.sortOn, allowedKeys ))

def readDirs( options ):
   """ readDirs reads two directories, the scaffolds contigPath dir and the contigns
   contigPath dir and creates a single dict that is keyed on assembly IDs. it combines
   the data from the two directories with a preference for the scaffolds data.
   The only data from the contigs dir is the contigN50 and contigNG50. Those same keys
   from the scaffolds dir are mapped to scaffoldN50 and scaffoldNG50, respectively.
   """
   # scaffolds
   assembliesDict = readDir( options.statsScaffoldsContigPathDir, options )
   # contigs
   assembliesDict = readDir( options.statsContigsContigPathDir, options, 
                             isScaffolds=False, assembliesDict=assembliesDict )
   return assembliesDict.values()

def readDir( theDir, options, isScaffolds=True, assembliesDict=None ):
   """ This is not a symmetric function! If isScaffolds is True, then everything
   is read out of the xml files, however, if isScaffolds is False, then only
   'contigN50' and 'contigNG50' are read out of the xml.
   """
   sFiles = glob.glob( os.path.join( theDir, '*.pathStats.xml'))
   namepat = re.compile( '^(\S{2,3})\.pathStats.xml' )
   if assembliesDict is None:
      assembliesDict = {}
   for f in sFiles:
      name = re.match( namepat, os.path.basename( f )).group( 1 )
      if 'subsetFile' in vars( options ):
         if options.subsetFile:
            if name not in options.assemblySubset:
               continue
      try:
         xmlTree = ET.parse( f )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( f )
      root=xmlTree.getroot()
      if name not in assembliesDict:
         a = Assembly()
         a.ID = name
      else:
         a = assembliesDict[ name ]
      if isScaffolds:
         for elm in root.attrib.keys():
            if elm in ('errorsPerContig', 'errorsPerMappedBase', 'coverage'):
               a.valuesDict[ elm ] = float( root.attrib[ elm ] )            
            elif elm in ( 'insertionErrorSizeDistribution', 'deletionErrorSizeDistribution' ):
               a.valuesDict[ elm ] = root.attrib[ elm ].split()
               for i in xrange( 0, len(a.valuesDict[ elm ])):
                  a.valuesDict[ elm ][i] = int(a.valuesDict[ elm ][i])
            else:
               a.valuesDict[ elm ] = int( root.attrib[ elm ] )
         a.valuesDict[ 'scaffoldN50' ] = int( root.attrib[ 'contigN50' ])
         a.valuesDict[ 'scaffoldNG50' ] = int( root.attrib[ 'contigNG50' ])
      else:
         a.valuesDict[ 'contigN50' ] = int( root.attrib[ 'contigN50' ])
         a.valuesDict[ 'contigNG50' ] = int( root.attrib[ 'contigNG50' ])
      if name not in assembliesDict:
         assembliesDict[ name ] = a
   return assembliesDict

def printTable( assembliesList, caption, options ):
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{FPtable}
\caption[Contig path statistics.]{Contig path statistics. %s}
\\tiny
\\centering
\\begin{tabular}{| r | p{.75in} | p{.75in} | p{.4in} | p{.4in} | p{.5in} | p{.5in} || c |}
\\hline
ID & Intra chromosomal joins & Inter chromosomal joins & Insertions & Deletions & Insertion and deletion & Insertion at ends & \(\sum\) errors \\\\
\\hline
\\hline''' % ( caption ) # hom/het sw & scf & ctgE+ctgN 
   i = 0
   for row in assembliesList:
      i += 1
      if options.showAssemblyNumbers:
         s = '%s.%s' % (lgn.idMap[row.ID[0]], row.ID[1:])
      else:
         s = '%s' % (lgn.idMap[row.ID[0]])
      sys.stdout.write( '%s' % s )
      # sys.stdout.write( ' & %s & %s' % ( lgn.prettyNumber(row.valuesDict['totalHomoToHeteroSwitches']),
      #                                    lgn.prettyNumber( row.valuesDict['totalScaffoldGaps'] + 
      #                                                  row.valuesDict['totalBleedingHeartScaffoldGaps'] ) ))
      # sys.stdout.write( ' & %s' % ( lgn.prettyNumber( row.valuesDict['totalContigEnds'] + 
      #                                             row.valuesDict['totalContigEndsWithNs'] )))
      #sys.stdout.write( ' & %.4f' % ( row.valuesDict['errorsPerContig']))
      #sys.stdout.write( ' & %.2e' % ( row.valuesDict['errorsPerMappedBase']))
      for v in [ 'totalErrorsHaplotypeToHaplotypeSameChromosome',
                 'totalErrorsHaplotypeToHaplotypeDifferentChromosome',
                 'totalErrorsHaplotypeToInsertion',
                 'totalErrorsHaplotypeToDeletion',
                 'totalErrorsHaplotypeToInsertionAndDeletion',
                 'totalErrorsContigEndsWithInsert' ]:
         sys.stdout.write( ' & %s' % ( lgn.prettyNumber( row.valuesDict[v] )))
      #sys.stdout.write( ' & %s' % ( lgn.prettyNumber( row.valuesDict['totalErrorsNonSpecific'])))
      sys.stdout.write( ' & %s \\\\\n' % lgn.prettyNumber( row.totalErrors ))
      if not i % 10 and i != len( assembliesList ):
         print '\\hline'
   print '''\\hline
\\end{tabular}
\\label{table:contigPathStats}
\\end{FPtable}\par
\\normalsize
\\vspace{0.3in}'''

def printRanks( assembliesList, options ):
   print '#Assembly\tSumErr\tHap-Hap-sameChr\tHap-Hap-diffChr\tHap-Ins\tHap-Del\tHap-Ins&Del\tCtgEndsIns'
   for row in assembliesList:
      sys.stdout.write( '%s' % ( row.ID ))
      sys.stdout.write( '\t%s' % row.totalErrors )
      for v in [ 'totalErrorsHaplotypeToHaplotypeSameChromosome',
                 'totalErrorsHaplotypeToHaplotypeDifferentChromosome',
                 'totalErrorsHaplotypeToInsertion',
                 'totalErrorsHaplotypeToDeletion',
                 'totalErrorsHaplotypeToInsertionAndDeletion',
                 'totalErrorsContigEndsWithInsert' ]:
         sys.stdout.write( '\t%s' % ( row.valuesDict[v] ))
      #sys.stdout.write( ' & %s' % ( lgn.prettyNumber( row.valuesDict['totalErrorsNonSpecific'])))
      sys.stdout.write( '\n' )

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
         sys.stderr.write('calculated sum of total errors not equal to value from xml! %s\n' % a.ID )

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
   usage = ( 'usage: %prog --statsScaffoldsContigPathDir=path/to/dir/ --statsContigsContigPathDir=path/to/dir/ [options]\n\n'
             '%prog takes in the contig path statistics directories from both the Scaffolds and Contigs alignments\n'
             '( --statsScaffoldsContigPathDir --statsContigsContigPathDir ) and prints to STDOUT a latex formatted table.' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )
   las.checkOptions( options, parser )
   
   assembliesList = readDirs( options )
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
      printTable( assembliesList, '', options )

if __name__ == '__main__':
   main()
