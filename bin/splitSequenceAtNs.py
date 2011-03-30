#!/usr/bin/env python
"""
Changes all headers to contain '.split%03d', even
if a split did not occur.
"""
from optparse import OptionParser
import signal # deal with broken pipes
import standardizeNumNs as snn
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '-n', '--splitAt', dest='n',
                      type='int',
                      help='Stretches of n or more will be split into two sequences.' )
   parser.add_option( '-l', '--lineLength', dest='lineLength',
                      type='int', default=50,
                      help='Changes the length of the output lines.' )
   parser.add_option( '--label', dest='label',
                       type='string', default='contig',
                       help='Will result in headers like: >PREFIX.label001')

def checkOptions( options, parser ):
   if options.n == None:
      parser.error( 'Error, specify --splitAt.\n' )
   if options.n < 1:
      parser.error( 'Error, --splitAt must be greater than 0.\n' )

def processStream( options ):
   register = 0
   header=''
   for line in sys.stdin:
      line = line.strip()
      if line == '':
         continue
      if line[0] == '>':
         if header != '' and register != 0:
            sys.stdout.write('\n')
         header = line
         splitCount = 1
         sys.stdout.write('%s.%s%03d\n' % ( header, options.label, splitCount ))
         register = 0
         nCount = 0
         continue
      for c in line:
         if c == 'N' or c == 'n':
            nCount += 1
            # we buffer the Ns in case we need to split the sequence
         else:
            if nCount < options.n:
               for j in range( 0, nCount ):
                  register = snn.myPrint( 'N', register, options )
            else:
               splitCount += 1
               if register !=0:
                  sys.stdout.write('\n')
               print '%s.%s%03d' % ( header, options.label, splitCount )
               register = 0
            nCount = 0
            register = snn.myPrint( c, register, options )
   if 0 < nCount < options.n:
      for j in range( 0, nCount ):
         register = snn.myPrint( 'N', register, options )
   if register != 0:
      sys.stdout.write( '\n' )

def main():
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   
   processStream( options )
   
if __name__ == '__main__':
   main()
