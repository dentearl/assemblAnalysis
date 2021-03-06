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
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '-n', '--expandAt', dest='n',
                      type='int',
                      help='Stretches of N or more (>=) will be expanded to 25Ns.' )
   parser.add_option( '-l', '--lineLength', dest='lineLength',
                      type='int', default=50,
                      help='Changes the length of the output lines. default=%default' )

def checkOptions( options, parser ):
   if options.n is None:
      parser.error( 'specify --expandAt.\n' )
   if options.n < 1:
      parser.error( '--expandAt must be greater than 0.\n' )
   if options.n > 25:
      parser.error( '--expandAt must be less than 26.\n' )

def myPrint( c, register, options ):
   """ my print takes a charater, c and a 
   register position, and based on options.lineLength
   decides whether or not to simply print the character
   or the character and a new line.
   """
   if register == ( options.lineLength - 1 ):
      sys.stdout.write( '%s\n' % c )
      register = 0
   else:
      sys.stdout.write( '%s' % c )
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
      if line.startswith('>'):
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
               for j in xrange( count, 26 ):
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
   usage = ( 'usage: %prog --expandAt=N [options] < fasta.fa\n\n'
             '%prog takes in a fasta formated file in STDIN and a number\n'
             '( --expandAt ) and scans through the fasta looking for runs equal to or greater than \n'
             '--expandAt N (or n) characters and alters the sequence to contain exactly 25 N\'s instead.\n'
             'Runs of 26 N\'s or more are truncated to be 25. Output written to STDOUT.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   processStream( options )

if __name__ == '__main__':
   main()
