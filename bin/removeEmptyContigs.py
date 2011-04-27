#!/usr/bin/env python
"""
removeEmptyContigs.py
March? 2011
dent earl, dearl (a) soe ucsc edu

simple script to walk a fasta file
and print out all the non-empty sequences.

So ill-formed contigs such as myBrokenContig
below:
...
ACGTACGTACGTACGTACGTACGT
>myBrokenContig

>myOKAYContig
ACGTACGTACGTACGTACGTACGT
...

are removed from the output:
...
ACGTACGTACGTACGTACGTACGT
>myOKAYContig
ACGTACGTACGTACGTACGTACGT
...
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
import sys

def main():
   usage = ( 'usage: %prog < fasta.fa\n\n'
             '%prog takes in via STDIN a fasta formated file writes to\n'
             'STDOUT all of the non-empty sequences.' )
   parser = OptionParser( usage=usage )
   options, args = parser.parse_args()
   
   header = ''
   for line in sys.stdin:
      line = line.strip()
      if line != '':
         if line[0] == '>':
            header = line
            continue
         if header != '':
            print header
            print line
            header = ''
         else:
            print line
      else:
         header = ''

if __name__ == '__main__':
    main()
