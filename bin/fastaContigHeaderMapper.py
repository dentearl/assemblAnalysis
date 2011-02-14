#!/usr/bin/env python
"""
fastaContigHeaderMapper.py 
14 Feb 2011
dent earl, dearl (a) soe ucsc edu

Simple script to both create mappings between a fasta
file's unique contig IDs that are verbose and unqiue
contig IDs that are short. This is mainly a way to get
around RepeatMasker's header line length limit of 50 chars.

Input and output are via stdin/stdout.

"""

import cPickle
import os
from optparse import OptionParser
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def usage():
    print 'USAGE: '+sys.argv[0]+' --map myMap.map --goForward [--goBackward] < inputFile.fa > outputFile.fa'
    print __doc__
    sys.exit( 2 )

def initOptions( parser ):
    parser.add_option( '--createMap', dest='createMap',
                       type='string',
                       help='Create the map file that contains the mapping between IDs.')
    parser.add_option( '--map', dest='map',
                       type='string',
                       help='Specify the map file that contains the mapping between IDs.')
    parser.add_option( '--goForward', dest='goForward',
                       action='store_true', default=False,
                       help='Use the supplied map file to move from column 1 to column 2.')
    parser.add_option( '--goBackward', dest='goBackward',
                       action='store_true', default=False,
                       help='Use the supplied map file to move from column 2 to column 1.')

def checkOptions( parser, options ):
    if options.createMap != None:
        return
    if options.map == None:
        parser.error('You must specify the map you wish to use with --map')
    if not os.path.exists( options.map ):
            parser.error('%s Does not exist.' % options.map )
    if ( not options.goForward ) and ( not options.goBackward ):
        parser.error('You must specify the direction to map with either --goForward or --goBackward.')

def createMap( options ):
    faMap = {}
    num = 1
    for line in sys.stdin:
        line=line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if line in faMap:
                sys.stderr.write( 'Error, duplicate contig header found: %s.\n' % line )
                sys.exit( 1 )
            faMap[ line ] = '>contig%06d' % num
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
                sys.stderr.write( 'Error, unable to find contig header "%s" in map file "%s"\n' % 
                                  ( line, options.map) )
                sys.exit( 1 )
        else:
            print line

def main():
    parser = OptionParser()
    initOptions( parser )
    ( options, args ) = parser.parse_args()
    checkOptions( parser, options )
    if options.createMap:
        createMap( options )
        return
    faMap = readMap( options )
    if options.goBackward:
        faMap = dict((v,k) for k, v in faMap.iteritems())
    translate( faMap, options )

if __name__ == '__main__':
   main()
