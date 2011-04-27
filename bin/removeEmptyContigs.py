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
from optparse import OptionParser
import sys

def initOptions( parser ):
   pass

def checkOptions( options, parser ):
   pass

def main():
   usage = ( 'usage: %prog < fasta.fa\n\n'
             '%prog takes in via STDIN a fasta formated file writes to\n'
             'STDOUT all of the non-empty sequences.' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   
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
