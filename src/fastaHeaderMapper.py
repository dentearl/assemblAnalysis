#!/usr/bin/env python
"""
fastaHeaderMapper.py 
14 Feb 2011
dent earl, dearl (a) soe ucsc edu

Simple script to both create mappings between a fasta
file's unique header IDs that are verbose and unqiue
header IDs that are short. This is mainly a way to get
around RepeatMasker's header line length limit of 50 chars.

Input and output are via stdin/stdout.

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
import cPickle
import os
from optparse import OptionParser
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
    parser.add_option( '--createMap', dest='createMap',
                       type='string',
                       help='Create the map file that contains the mapping between IDs.')
    parser.add_option( '--prefix', dest='prefix',
                       type='string',
                       help='When using --createMap one can specify a header prefix that will '
                       'result in headers like: >PREFIX.label000001')
    parser.add_option( '--label', dest='label',
                       type='string', default='contig',
                       help='When using --createMap one can specify a header label that will '
                       'result in headers like: >prefix.LABEL000001 . default=%default')
    parser.add_option( '--map', dest='map',
                       type='string',
                       help='Specify the map file that contains the mapping between IDs.')
    parser.add_option( '--goForward', dest='goForward',
                       action='store_true', default=False,
                       help='Use the supplied map file to move from column 1 to column 2. default=%default')
    parser.add_option( '--goBackward', dest='goBackward',
                       action='store_true', default=False,
                       help='Use the supplied map file to move from column 2 to column 1. default=%default')

def checkOptions( parser, options ):
    if options.createMap is not None:
        return
    if options.prefix is not None:
        parser.error('--prefix may only be used in conjunction with --createMap.')
    if options.map is None:
        parser.error('You must specify the map you wish to use with --map.')
    if not os.path.exists( options.map ):
            parser.error('%s Does not exist.' % options.map )
    if  not options.goForward and not options.goBackward:
        parser.error('You must specify the direction to map with either --goForward or --goBackward.')

def createMap( options ):
    faMap = {}
    num = 1
    if options.prefix is not None:
        prefix = '%s.' % options.prefix
    else:
        prefix = ''
    for line in sys.stdin:
        line=line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if line in faMap:
                sys.stderr.write( 'duplicate header found: %s.\n' % line )
                sys.exit( 1 )
            faMap[ line ] = '>%s%s%06d' % ( prefix, options.label, num )
            num += 1
            
    FILE = open( options.createMap, 'w' )
    cPickle.dump( faMap, FILE, 0) # 0=ASCII, 1=old binary, 2=python 2.3 binary
    FILE.close()

def readMap( options ):
    FILE = open( options.map )
    data = cPickle.load( FILE )
    FILE.close()
    return data

def translate( faMap, options ):
    for line in sys.stdin:
        line=line.strip()
        if line == '':
            print line
            continue
        if line[0] == '>':
            if line in faMap:
                print faMap[ line ]
            else:
                sys.stderr.write( 'unable to find header "%s" in map file "%s"\n' % 
                                  ( line, options.map) )
                sys.exit( 1 )
        else:
            print line

def main():
    usage = ( 'usage: %prog [options]\n\n'
              '%prog creates, reads and uses fasta sequence header mappings to\n'
              'map lengthy sequence headers to shorter standardized headers.')
    parser = OptionParser( usage=usage )
    initOptions( parser )
    options, args = parser.parse_args()
    checkOptions( parser, options )
    if options.createMap:
        createMap( options )
        return
    faMap = readMap( options )
    if options.goBackward:
        faMap = dict( (v,k) for k, v in faMap.iteritems() )
    translate( faMap, options )

if __name__ == '__main__':
   main()
