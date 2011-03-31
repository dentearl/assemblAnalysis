#!/usr/bin/env python
"""
extractPlotPickleCounts.py
4 March 2011
dent earl, dearl (a) soe ucsc edu

simple script to open up plot pickles and pull out all of the counts
stored in the numpy array under a particular key value.
"""
import cPickle
from libMafGffPlot import Data
import numpy
from optparse import OptionParser
import os
import sys

def initOptions( parser ):
   parser.add_option( '--key', dest='key',
                      type='string',
                      help='Key to extract from supplied maf plot pickles.' )
   parser.add_option( '--printAllowedKeys', dest='printAllowedKeys',               
                      action='store_true', default=False,                          
                      help=('Prints out the allowed keys for --sortOn and exits.'))
   
def checkOptions( args, options, parser, data ):
   allowedKeys = set( [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
                        'maf1e5', 'maf1e6', 'maf1e7', 'blockEdgeCount',
                        'mafHpl1e2', 'mafHpl1e3', 'mafHpl1e4', 
                        'mafHpl1e5', 'mafHpl1e6', 'mafHpl1e7', 
                        'mafCtg1e2', 'mafCtg1e3', 'mafCtg1e4', 
                        'mafCtg1e5', 'mafCtg1e6', 'mafCtg1e7', 
                        'mafSpl1e2', 'mafSpl1e3', 'mafSpl1e4', 
                        'mafSpl1e5', 'mafSpl1e6', 'mafSpl1e7', 
                        'mafHpEdgeCount', 'mafHpErrorCount','mafHpScafGapCount',
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
      parser.error('Error, please specify --key\n')
   if options.key not in allowedKeys:
      parser.error('Error, please specify --key from one of the %s\n' % allowedKeys )
   if len( args ) < 1:
      parser.error('Error, please specify a list of pickles to inspect.\n' )
   for f in args:
      if not os.path.exists( f ):
         parser.error('Error, file "%s" does not exist.\n' % f )
      if f[-7:] != '.pickle':
         parser.error('Error, file "%s" does not end in ".pickle", aborting.\n' % f )

def printData( valuesDict, options, data ):
   if options.key not in valuesDict:
      sys.stderr.write( 'Error, key %s not in this dictionary.\n' % options.key )
      sys.exit( 1 )
   if ( isinstance( valuesDict[ options.key ], float ) or
        isinstance( valuesDict[ options.key ], int )):
      print valuesDict[ options.key ]
   elif isinstance( valuesDict[ options.key ], numpy.ndarray ):
      for i in range( 0, len( valuesDict[ options.key ]) ):
         print valuesDict[ options.key ][i]
   else:
      sys.stderr.write( 'Error, unexpected object type '
                        'for valuesDict[ %s ]: %s\n' % ( options.key, 
                                                         valuesDict[ options.key ].__class__ ))
      sys.exit( 1 )

def unpackData( filename, options, data ):
   if not os.path.exists( filename ):
      sys.stderr.write( 'Error, %s does not exist.\n' % filename)
      sys.exit( 1 )
   f = open( filename, 'rb' )
   d = cPickle.load( f )
   f.close()
   return d

def loadPickles( args, options, data ):
   for f in args:
      printData( unpackData( f, options, data ), options, data )

def main():
   usage = ( 'usage: %prog --key=KEY file1.maf.pickle file2.maf.pickle ...\n\n'
             '%prog takes a valid key ( --key ) and one or more maf wiggle pickle(s)\n'
             'and then pulls out all the counts stored in the numpy array for the key.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser, data )
   loadPickles( args, options, data )

if __name__ == '__main__':
   main()
