#!/usr/bin/env python
"""pipe in input. prints out the lengths of all sequences
in the fasta. --names option prints names and lengths, one 
per line.
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
import sys
import os
from optparse import OptionParser
import signal # deal with broken pipes
signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
    parser.add_option( '--names', dest='names',
                       action='store_true', default=False,
                       help='Prints out sequence names in addition to their lengths [default %default].' )

def checkOptions( options ):
    pass

def reportSeq( s, l, options ):
    if l > 0:
        if options.names:
            print '%s, [ %d ]' % ( s, l )
        else:
            print '%d' % l

def main():
    usage = ( 'usage: %prog [--names]\n\n'
              '%prog takes in a fasta file via STDIN and \n'
              'then prints out the lengths of sequences contained in the \n'
              'fasta, one length per line of STDOUT. If --names is \n'
              'specified, sequence names are printed next to their length.')
    parser = OptionParser( usage=usage )
    initOptions( parser )
    options, args = parser.parse_args()
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
