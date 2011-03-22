#!/usr/bin/env python
"""
standardizeNumNs.py
dent earl, dearl (a) soe ucsc edu
20 March 2011

standardize the number of Ns based on
--expandAt. --expandAt is treated as a
greater than or equal to, up to 25.

For example if --expandAt is 
at 3, then 
ACGTACGTNNNNNNNNNNACGTACGT
becomes
ACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGT

ACGTACGTNNACGTACGT
becomes
ACGTACGTNNACGTACGT

ACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGT
becomes
ACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGT

"""
from optparse import OptionParser
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '-n', '--expandAt', dest='n',
                      type='int',
                      help='Stretches of n or more (>=) will be expanded to 25Ns.' )
   parser.add_option( '-l', '--lineLength', dest='lineLength',
                      type='int', default=50,
                      help='Changes the length of the output lines.' )

def checkOptions( options, parser ):
   if options.n == None:
      parser.error( 'Error, specify --expandAt.\n' )
   if options.n < 1:
      parser.error( 'Error, --expandAt must be greater than 0.\n' )
   if options.n > 25:
      parser.error( 'Error, --expandAt must be less than 26.\n' )

def myPrint( s, register, options ):
   if register == ( options.lineLength - 1 ):
      sys.stdout.write('%s\n' % s )
      register = 0
   else:
      sys.stdout.write('%s' % s)
      register += 1
   return register

def processStream( options ):
   count = 0
   i = 0
   lineNumber = -1
   expanded = False
   for line in sys.stdin:
      lineNumber += 1
      line = line.strip()
      if line == '':
         continue
      if line[0] == '>':
         if i == 0:
            print '%s' % line
         else:
            print '\n%s' % line
         i = 0
         continue
      for c in line:
         if c == 'n' or c == 'N':
            count += 1
            if count == options.n:
               expanded = True
               for j in range( count, 26 ):
                  # expand to 25 Ns.
                  i = myPrint( 'N', i, options )
            elif not expanded:
               i = myPrint( 'N', i, options )
         else:
            count = 0
            expanded = False
            i = myPrint( c, i, options )
   if i != 0:
      sys.stdout.write('\n')

def main():
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   processStream( options )

if __name__ == '__main__':
   main()
