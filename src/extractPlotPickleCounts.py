#!/usr/bin/env python
"""
extractPlotPickleCounts.py
4 March 2011
dent earl, dearl (a) soe ucsc edu

simple script to open up plot pickles and pull out all of the counts
stored in the numpy array under a particular key value.
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
from libMafGffPlot import Data
from libMafGffPlot import unpackData
import numpy
from optparse import OptionParser
import os
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '--key', dest='key',
                      type='string',
                      help='Key to extract from supplied maf plot pickles.' )
   parser.add_option( '--chr', dest='chr',
                      type='string',
                      help='Used to limit the output to just one or more chromosomes. Comma separated list.' )
   parser.add_option( '--printChromosomes', dest='printChroms',
                      default=False, action='store_true',
                      help='Prints the chromosome keys in a given pickle. default=%default' )
   parser.add_option( '--printAllowedKeys', dest='printAllowedKeys',               
                      action='store_true', default=False,                          
                      help=('Prints out the allowed keys for --sortOn and exits. default=%default'))
   parser.add_option( '--verify', dest='verify',               
                      action='store_true', default=False,                          
                      help=('Enables extra checks to verify the data structure is accurate. '
                            'Not necessary unless the output plots look odd. default=%default' ))
   
def checkOptions( args, options, parser, data ):
   allowedKeys = set( [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
                        'maf1e5', 'maf1e6', 'maf1e7', 'blockEdgeCount',
                        'mafCpl1e2', 'mafCpl1e3', 'mafCpl1e4', 
                        'mafCpl1e5', 'mafCpl1e6', 'mafCpl1e7', 
                        'mafCtg1e2', 'mafCtg1e3', 'mafCtg1e4', 
                        'mafCtg1e5', 'mafCtg1e6', 'mafCtg1e7', 
                        'mafSpl1e2', 'mafSpl1e3', 'mafSpl1e4', 
                        'mafSpl1e5', 'mafSpl1e6', 'mafSpl1e7', 
                        'mafCpEdgeCount', 'mafCpErrorCount','mafCpScafGapCount',
                        'CDSCount', 'UTRCount',
                        'NXECount', 'NGECount', 'islandCount', 'tandemCount', 
                        'repeatCount', 'CDSMax', 'UTRMax', 'NXEMax', 'NGEMax',
                        'islandMax', 'tandemMax', 'repeatMax', 'xAxis', 
                        'columnsInBlocks'
                        ])
   if options.printAllowedKeys:
      for k in allowedKeys:
         print k
      sys.exit( 0 )
   if options.key == '':
      parser.error('please specify --key\n')
   if options.key not in allowedKeys:
      parser.error('please specify --key from one of the %s\n' % allowedKeys )
   if len( args ) < 1:
      parser.error('please specify at least one pickle to inspect as a positional argument.\n' )
   for f in args:
      if not os.path.exists( f ):
         parser.error('file "%s" does not exist.\n' % f )
      if not f.endswith('.pickle'):
         parser.error('file "%s" does not end in ".pickle".\n' % f )
   if options.printChroms:
      for f in args:
         d = unpackData( f, options, {} )
         sys.stdout.write('%s\t' % f )
         for c in d:
            sys.stdout.write('\t%s' % c)
         sys.stdout.write('\n')
      sys.exit(0)
   if options.chr is not None:
      options.chrSet = set( options.chr.split(',') )
      for f in args:
         d = unpackData( f, options, {} )
         for c in options.chrSet:
            if c not in d:
               sys.stderr.write('chromosome %s is not present in %s\n' % ( c, f ))
   else:
      options.chrSet = set()
      for f in args:
         d = unpackData( f, options, {} )
         for c in d:
            if c not in options.chrSet:
               options.chrSet.add( c )

def printData( valuesDict, options, data ):
   for c in options.chrSet:
      checkKey( options.key, c, valuesDict )
      if ( isinstance( valuesDict[c][ options.key ], float ) or
           isinstance( valuesDict[c][ options.key ], int )):
         print valuesDict[c][ options.key ]
      elif isinstance( valuesDict[c][ options.key ], numpy.ndarray ):
         for i in range( 0, len( valuesDict[c][ options.key ]) ):
            print valuesDict[c][ options.key ][i]
      else:
         sys.stderr.write( 'unexpected object type '
                           'for valuesDict[ %s ][ %s ]: %s\n' % ( c, options.key, 
                                                                  valuesDict[c][ options.key ].__class__ ))
         sys.exit( 1 )

def checkKey( key, chr, dictionary ):
   if chr not in dictionary:
      sys.stderr.write( 'chromosome %s not in this dictionary.\n' % chr )
      sys.exit( 1 )
   if key not in dictionary[ chr ]:
      sys.stderr.write( 'key %s not in this dictionary for chr %s.\n' % ( key, chr ))
      sys.exit( 1 )

def verify( valuesDict, options, data ):
   for c in valuesDict:
      checkKey( options.key, c, valuesDict )
      if ( isinstance( valuesDict[c][ options.key ], float ) or
           isinstance( valuesDict[c][ options.key ], int )):
         pass
      elif isinstance( valuesDict[c][ options.key ], numpy.ndarray ):
         if sum( valuesDict[c][ options.key ] > 1.0 ) > 0:
            sys.stderr.write('elements greater than 1.0 detected.\n')
      else:
         sys.stderr.write( 'unexpected object type '
                           'for valuesDict[ %s ][ %s ]: %s\n' % ( options.key, c,
                                                                  valuesDict[c][ options.key ].__class__ ))
      sys.exit( 1 )

def loadPickles( args, options, data ):
   for f in args:
      valuesDict = unpackData( f, options, data )
      if options.verify:
         verify( valuesDict, options, data )
      printData( valuesDict, options, data )

def main():
   usage = ( 'usage: %prog --key=KEY file1.maf.pickle file2.maf.pickle ...\n\n'
             '%prog takes a valid key ( --key ) and one or more maf wiggle pickle(s)\n'
             'and then pulls out all the counts stored in the numpy array for the key.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser, data )
   
   loadPickles( args, options, data )
   

if __name__ == '__main__':
   main()
