#!/usr/bin/env python
"""
summarizeRankings.py
dent earl, dearl (a) soe ucsc edu
27 April 2011

Python script to take a set of rankings files and aggregate their results.

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
from optparse import OptionParser
import os
import sys
import signal # deal with broken pipes
signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Assembly:
   def __init__(self):
      self.teamName = ''
      self.name     = ''
      self.ranks    = []
      self.values   = []

class Tab:
   def __init__(self):
      self.name  = ''
      self.rank  = 0
      self.value = -1.0

def initOptions( parser ):
   parser.add_option( '--overall', dest='overall', default=False, action='store_true',
                      help=('Instead of ranking based on within Team, '
                            'ranks instead are global. default=%default'))
   parser.add_option( '--retainValues', dest='retainValues', default=False,
                      action='store_true',
                      help=('Stores the second column of the tab and outputs the value in '
                            'parenthesis following the ranking.'))
   parser.add_option( '--mode', dest='mode', default='s',
                      type='string',
                      help=('Mode can be "s" for a sum or "h" for harmonic mean.'))

def checkOptions( args, options, parser ):
   if len( args ) < 1:
      parser.error('no files in positional arguments. Input a filename.')
   for f in args:
      if not os.path.exists( f ):
         parser.error('file "%s" does not exist.\n' % f )
      if not f.endswith('.tab'):
         parser.error('file "%s" does not end in ".tab".\n' % f )
   options.mode = options.mode.lower()
   if options.mode not in ('s', 'h'):
      parser.error( 'Unrecognized --mode %s. Choose either "s" for '
                    'sum or "h" for harmonic.' % options.mode )

def readFiles( args, options ):
   assemblies = {}
   options.fileNames = []
   for aFile in args:
      options.fileNames.append( os.path.basename( aFile ) )
      f = open( aFile, 'r' )
      rank  = 0
      count = 0
      prevValue = - sys.maxint
      for line in f:
         line = line.strip()
         if line.startswith('#'):
            continue
         a = Tab()
         d = line.split()
         a.name  = d[0]
         a.value = d[1]
         if a.value == prevValue:
            count += 1
         else:
            count += 1
            rank = count
         a.rank  = rank
         prevValue = a.value
         if a.name not in assemblies:
            assemblies[ a.name ] = []
         assemblies[ a.name ].append( a )
   return assemblies

def printHeader( options ):
   sys.stdout.write( '#Assembly\tOverall' )
   for f in options.fileNames:
      if os.path.basename( f ).endswith('.tab'):
         name = os.path.basename( f )[:-4]
      else:
         name = os.path.basename( f )
      if options.retainValues:
         name += ' (value)'
      sys.stdout.write('\t%s' % name )
   sys.stdout.write('\n')

def reportRank( assemblies, options ):
   printHeader( options )
   if options.overall:
      # ranking is between all assemblies
      overallDict = {}
      for t in assemblies:
         overallDict[t] = 0
         for a in t:
            overallDict[t] += a.rank
      ranked = sorted( overallDict, key=lambda x: overallDict[x] )
      for r in ranked:
         print '%s\t%d' % ( r, overallDict[r] )
         for z in assemblies[r]:
            if options.retainValue:
               sys.stdout.write('\t%d (%s)' % (z.rank, z.value))
            else:
               sys.stdout.write('\t%d' % z.rank )
         sys.stdout.write('\n')
   else:
      # ranking is within teams
      teams = {}
      for a in assemblies:
         if a[0] not in teams:
            teams[ a[0] ] = []
         newAssemb = Assembly()
         newAssemb.teamName = a[0]
         newAssemb.name = a
         for z in assemblies[ a ]:
            newAssemb.ranks.append( z.rank )
            newAssemb.values.append( z.value )
         teams[ a[0] ].append( newAssemb )
      # sort alphabetically
      aemst = teams.keys()
      aemst.sort()
      for t in aemst:
         print ''
         if options.mode == 's':
            ranked = sorted( teams[ t ], key=lambda x: sum( x.ranks ))
         else:
            ranked = sorted( teams[ t ], key=lambda x: harmonic( x.ranks ))
         for r in ranked:
            if options.mode == 's':
               sys.stdout.write( '%s\t%d' % ( r.name, sum( r.ranks ) ))
            else:
               sys.stdout.write( '%s\t%.2f' % ( r.name, harmonic( r.ranks ) ))
            for j in xrange( 0, len(r.ranks)):
               if options.retainValues:
                  sys.stdout.write('\t%d (%s)' % (r.ranks[j], r.values[j] ))
               else:
                  sys.stdout.write('\t%d' % r.ranks[j] )
            sys.stdout.write('\n')

def harmonic( ranks ):
   """ returns the harmonic mean of the assembly's ranks
   """ 
   h = 0
   for r in ranks:
      h += 1.0 / r
   h = 1.0 / h
   h *= len( ranks )
   return h
   

def main():
   usage = ( 'usage: %prog [options] rankingFile1.tab rankingFile2.tab rankingFile3.tab\n\n'
             '%prog takes in ranking files and reports the best overall or best within a Team.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )
   
   assemblies = readFiles( args, options )
   reportRank( assemblies, options )

if __name__ == '__main__':
   main()

