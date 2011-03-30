#!/usr/bin/env python
"""pipe in input. prints out the lengths of all sequences
in the fasta. --names option prints names and lengths, one 
per line.
"""
import sys
import os
from optparse import OptionParser
import signal # deal with broken pipes
signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
    parser.add_option( '--names', dest='names',
                       action='store_true', default=False,
                       help='Prints out sequence names in addition to their lengths [default False].' )

def checkOptions( options ):
    pass

def reportSeq( s, l, options ):
    if l > 0:
        if options.names:
            print '%s, [ %d ]' % ( s, l )
        else:
            print '%d' % l

def main():
    parser = OptionParser()
    initOptions( parser )
    ( options, args ) = parser.parse_args()
    checkOptions( options )
    curSeq = ''
    curLen = 0
    if options.names:
        print 'Name, [ length ]'
    for line in sys.stdin:
        line = line.strip()
        if line == '':
            continue
        if line[0] == '>':
            reportSeq( curSeq, curLen, options )
            curLen = 0
            curSeq = line
        else:
            curLen += len( line )
    reportSeq( curSeq, curLen, options )


if __name__ == '__main__':
    main()
