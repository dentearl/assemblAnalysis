#!/usr/bin/env python
"""
Supply a character and then check to see the distribution
of repeats of that character in the STDIN file (fasta).
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
import splitSequenceAtNs as ssan
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '-c', '--char', dest='char',
                      type='string',
                      help='Character to look out for (case insensitive).' )

def checkOptions( options, parser ):
   if options.char is None:
      parser.error( 'specify --char.\n' )
   if len( options.char ) > 1:
      parser.error( '--char should only be one character.\n')
   options.uChar = options.char.upper()
   options.lChar = options.char.lower()

def processStream( options ):
   count = 0
   for line in sys.stdin:
      line = line.strip()
      if line == '':
         continue
      if line.startswith('>'):
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
   options, args = parser.parse_args()
   checkOptions( options, parser )
   
   processStream( options )

if __name__ == '__main__':
   main()

