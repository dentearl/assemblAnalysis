#!/usr/bin/env python
"""
simple script to walk a fasta file
and print out all the non-empty contigs.
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

import sys

def main():
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
