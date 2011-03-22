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
                       'result in headers like: >prefix.LABEL000001')
    parser.add_option( '--label', dest='label',
                       type='string', default='contig',
                       help='When using --createMap one can specify a header label that will '
                       'result in headers like: >PREFIX.label000001')
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
    if options.prefix != None:
        parser.error('--prefix may only be used in conjunction with --createMap.')
    if options.map == None:
        parser.error('You must specify the map you wish to use with --map.')
    if not os.path.exists( options.map ):
            parser.error('%s Does not exist.' % options.map )
    if ( not options.goForward ) and ( not options.goBackward ):
        parser.error('You must specify the direction to map with either --goForward or --goBackward.')

def createMap( options ):
    faMap = {}
    num = 1
    if options.prefix != None:
        prefix = '%s.' % options.prefix
    else:
        prefix = ''
    for line in sys.stdin:
        line=line.strip()
        if line == '':
            continue
        if line[0] == '>':
            if line in faMap:
                sys.stderr.write( 'Error, duplicate header found: %s.\n' % line )
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
                sys.stderr.write( 'Error, unable to find header "%s" in map file "%s"\n' % 
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
        faMap = dict( (v,k) for k, v in faMap.iteritems() )
    translate( faMap, options )

if __name__ == '__main__':
   main()
