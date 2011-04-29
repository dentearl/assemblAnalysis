#!/usr/bin/env python
"""
createN50StatsTable.py
13 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the N50 stats table.

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
from libGeneral import prettyNumber
import createContigPathStatsTable as ccpst
import createN50StatsPlot as cnfsp
import xml.etree.ElementTree as ET
import glob
from optparse import OptionParser
import os
import re
import sys

def printTable( assembliesList, caption, maxes, options ):
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{table}
\caption[N50 statistics.]{N50 statistics. %s}
\\tiny
\\centering
\\begin{tabular}{| r | c | c | c | c | c | c | c |}
\\hline
ID & \# Contigs & N50 & NA50 & SPA50 & HPA50 & BNA50 & \(\sum\) Errors\\\\
\\hline
\\hline''' % caption
   i = 0
   eMin = calculateMinError( assembliesList )
   for row in assembliesList:
      i += 1
      sys.stdout.write( '%s' % ( row.ID ))
      for n in [  'totalContigNumber','contigN50', 'contigNA50',
                  'scaffoldPathN50' , 'haplotypePathN50', 'blockN50' ]:
         sys.stdout.write( ' & %s' % ( isMaxFormat( row.valuesDict[ n ], maxes[ n ] ) ))
      sys.stdout.write( ' & %s' % isMaxFormat( row.totalErrors, eMin ))
      sys.stdout.write( ' \\\\\n' ) 
      if not i % 10:
         print '\\hline'

   print '''\\hline
\\end{tabular}
\\label{table:N50Stats}
\\end{table}\par
\\normalsize
\\vspace{0.3in}'''

def isMaxFormat( n, m ):
   if n == m:
      return '\\textbf{ %s }' % prettyNumber( n )
   else:
      return '%s' % prettyNumber( n )

def calculateMaxesDict( assembliesList ):
   maxesDict = { 'totalContigNumber':0,
                 'contigN50':0,
                 'contigNA50':0,
                 'scaffoldPathN50':0,
                 'haplotypePathN50':0,
                 'blockN50':0 }
   for a in assembliesList:
      for m in maxesDict:
         if maxesDict[ m ] < a.valuesDict[ m ]:
            maxesDict[ m ] = a.valuesDict[ m ]
   return maxesDict

def calculateMinError( assembliesList ):
   eMin = sys.maxint
   for a in assembliesList:
      if eMin > a.totalErrors:
         eMin = a.totalErrors
   return eMin

def main():
   usage = ('usage: %prog --contigPathStatsDir=path/to/dir/ [options]\n\n'
            '%prog takes the contig path stats directory\n'
            '( --contigPathStatsDir ) with names as NAME.contigPathStats.xml and then\n'
            'writes to STDOUT a latex formatted table.')
   parser = OptionParser( usage=usage )
   cnfsp.initOptions( parser )
   options, args = parser.parse_args()
   cnfsp.checkOptions( options, parser )
   
   assembliesList = ccpst.readDirs( options )
   ccpst.calculateErrors( assembliesList, options )
   assembliesList = sorted( assembliesList, key=lambda x: x.valuesDict[ options.sortOn ], reverse=True )
   maxesDict = calculateMaxesDict( assembliesList )
   
   caption = 'Columns are the total number of contigs in the assembly, N50, N50 relative to the number of columns in the alignment (NA50) as defined in the main text section \\ref{sect:NA50}, the scaffold path 50 (SPA50), contig path (HPA50), block path 50 (BA50), and the sum of the sum of the total number of errors present in the assembly (\\(\\sum\\) Errors).'
   printTable( assembliesList, caption, maxesDict, options )

if __name__ == '__main__':
   main()
