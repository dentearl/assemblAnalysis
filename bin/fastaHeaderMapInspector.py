#!/usr/bin/env python
"""
fastaHeaderMapInspector.py 
22 March 2011
dent earl, dearl (a) soe ucsc edu

Script that uses the .map file that is created
by fastaHeaderMapper.py and pulls out the headers
based on the standardized headers.

"""
import cPickle
import fastaHeaderMapper as fhm
import os
from optparse import OptionParser
import sys

def initOptions( parser ):
    parser.add_option( '--map', dest='map',
                       type='string',
                       help='Specify the map file that contains the mapping between IDs.')
    parser.add_option( '--id', '--key', dest='id',
                       type='string',
                       help='Specify the mapped ID you want to look up, i.e. >Z1.scaffold000001')

def checkOptions( parser, options ):
    if options.map == None:
        parser.error('You must specify the map you wish to inspect with --map.')
    if not os.path.exists( options.map ):
            parser.error('%s Does not exist.' % options.map )
    if options.id == None:
        parser.error('You must specify the id you wish to inspect with --id.')
    if options.id[0] != '>':
      options.id = '>%s' % options.id

def idLookup( faMap, options ):
   if options.id not in faMap:
      sys.stderr.write( 'Error, cannot find %s in map file.\n' % options.id )
      sys.exit(1)
   print faMap[ options.id ]

def main():
    parser = OptionParser()
    initOptions( parser )
    ( options, args ) = parser.parse_args()
    checkOptions( parser, options )
    faMap = fhm.readMap( options )
    faMap = dict( (v,k) for k, v in faMap.iteritems() )
    idLookup( faMap, options )

if __name__ == '__main__':
   main()