#!/usr/bin/env python
"""
Supply a character and then check to see the distribution
of repeats of that character in the STDIN file (fasta).
"""

from optparse import OptionParser
import splitSequenceAtNs as ssan
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '-c', '--char', dest='char',
                      type='string',
                      help='Character to look out for (case insensitive).' )

def checkOptions( options, parser ):
   if options.char == None:
      parser.error( 'Error, specify --char.\n' )
   if len( options.char ) > 1:
      parser.error( 'Error, --char should only be one character.\n')
   options.uChar = options.char.upper()
   options.lChar = options.char.lower()

def processStream( options ):
   count = 0
   for line in sys.stdin:
      line = line.strip()
      if line == '':
         continue
      if line[0] == '>':
         continue
      for c in line:
         if c == options.uChar or c == options.lChar:
            count += 1
         else:
            if count > 0:
               print count
            count = 0
   if count > 0:
      print count

def main():
   usage = ( 'usage: %prog [options] \n'
             'Supply a character (ex: --char=N ) and then see the distribution\n'
             'of repeats of that character in the STDIN file (fasta format).' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   
   processStream( options )

if __name__ == '__main__':
   main()

