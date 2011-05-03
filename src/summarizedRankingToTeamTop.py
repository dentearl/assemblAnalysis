#!/usr/bin/env python
"""
summarizedRankingsToTeamTop.py
dent earl, dearl (a) soe ucsc edu
27 April 2011

Script to take the rankings of all assemblies each teams and pull out 
only the top assembly from each team

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
import libAssemblySubset as las
from optparse import OptionParser
import re
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Assembly:
   def __init__(self):
      self.name   = ''
      self.ranks  = {}
      self.values = {}

def initOptions( parser ):
   parser.add_option( '--reRank', dest='reRank', 
                      action='store_true', default=False,
                      help=('After pulling out the top assembly from each team '
                            'the ranks will be recalculated relative to just the '
                            'top assemblies.'))

def checkOptions( options, parser ):
   pass

def simpleProcessStream( options ):
   """ simpleProcessStream() just looks for line breaks since the 
   assemblies will be grouped by team and sorted, so the top assembly
   in a group is the top for that team.
   """
   prevLineEmpty = False
   for line in sys.stdin:
      line = line.strip()
      if line.startswith('#'):
         print line
      if line == '':
         prevLineEmpty = True
         continue
      if prevLineEmpty == True:
         print line
         prevLineEmpty = False

def fullProcessStream( options ):
   """ fullProcessStream needs to read in everything from stdin and then 
   parse the input into objects, rank those objects, then recomput the total
   (or 'overall') column, and output the results.
   """
   if options.assemblySubset == set():
      sys.stderr.write('Error, you must specify a valid assembly subset file (--subsetFile) to use --reRank.\n')
      sys.exit(1)
   headerList = []
   assembliesList = []
   pat = re.compile('\s+\((\S+)\)')
   for line in sys.stdin:
      line = line.strip()
      if line.startswith('#'):
         headerList = line.split('\t')
         continue
      if line == '':
         continue
      a = Assembly()
      d = line.split('\t')
      a.name = d[0]
      if a.name not in options.assemblySubset:
         continue
      # [0] is the name, [1] is the overall rank
      for i in xrange(2, len(headerList)):
         m = re.search(pat, d[i])
         if not m:
            sys.stderr.write('Error, unable to locate pattern in %s for assembly %s. '
                             'Be sure to use --retainValue in summarizeRankings.py\n' % 
                             (d[i], a.name))
            sys.exit(1)
         a.values[ headerList[i] ] = float( m.group(1) )
      assembliesList.append( a )
   rankAssemblies( assembliesList, headerList[2:] )
   ranked = sorted( assembliesList, key=lambda x: sum(x.ranks.values()))
   
   print '\t'.join( headerList )
   for r in ranked:
      sys.stdout.write('%s\t%d' % (r.name, sum(r.ranks.values())))
      for j in xrange(2, len(headerList)):
         sys.stdout.write('\t%d (%s)' % (r.ranks[headerList[j]], 
                                         formattedValue( headerList[j], r.values[headerList[j]] )))
      sys.stdout.write('\n')

def formattedValue( header, value ):
   if header == 'N50_CPNG50 (value)':
      return '%.2e' % value
   elif header =='N50_SPNG50 (value)':
      return '%.2e' % value
   elif header == 'contiguousRanks (value)':
      return '%.2e' % value
   elif header == 'copyNumberErrors (value)':
      return '%.2e' % value
   elif header == 'coverageCDS (value)':
      return '%.1f' % value
   elif header == 'coverageTotal (value)':
      return '%.1f' % value
   elif header == 'structuralContigPathErrors (value)':
      return '%d' % value
   elif header == 'substitutionErrors (value)':
      return '%.2e' % value

def rankAssemblies( assembliesList, columns ):
   # for sortDir, key is the header, value is the boolean "Larger is better"
   sortDir = { 'N50_CPNG50 (value)':True, 'N50_SPNG50 (value)':True,
               'contiguousRanks (value)':True, 'copyNumberErrors (value)':False,
               'coverageCDS (value)':True, 'coverageTotal (value)':True,
               'structuralContigPathErrors (value)':False, 'substitutionErrors (value)':False }
   for c in columns:
      rankedCol = sorted( assembliesList, key = lambda x: x.values[c], reverse=sortDir[c] )
      rank  = 0
      count = 0
      prevValue = - sys.maxint
      for r in rankedCol:
         if r.values[c] == prevValue:
            count += 1
         else:
            count += 1
            rank = count
         r.ranks[c] = rank
         prevValue = r.values[c]

def main():
   usage = ( 'usage: %prog [options] < summarizeRankings.out\n\n'
             '%prog takes in the output of summarizeRankings.py and reports only the.\n'
             'top assembly per team. Use --reRank to rerank the top assembly per team\n'
             'the tops from all the other teams.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   las.checkOptions( options, parser )
   
   if not options.reRank:
      simpleProcessStream( options )
   else:
      fullProcessStream( options )
      
if __name__ == '__main__':
   main()
