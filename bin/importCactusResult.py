#!/usr/bin/env python
""" 
import CactusResult.py
15 April 2011
dent earl, dearl (a) soe ucsc edu

Script to import a cactus result directory into an assemblathon
analysis project. See the usage string in the main() function
for details.

"""
from optparse import OptionParser
import os
import re
import string
import sys

def initOptions( parser ):
   parser.add_option( '-i', '--inDir', dest='inDir',
                      type='string',
                      help='Cactus result directory you would like to import.' )
   parser.add_option( '-o', '--outDir', dest='outDir',
                      type='string',
                      help=('Directory where you want to place the '
                            'data. If it does not exist, it will be created.' ))
   parser.add_option('--name', dest='name', default='R1',
                     help=('Sets the name of the input assembly. default=%default'))
   parser.add_option('--type', dest='type',
                     type='string',
                     help=('Type of the alignment, either scaffold or contig.'))

def checkOptions( options, parser ):
   opts = { 'inDir'  : options.inDir,
            'outDir' : options.outDir,
            'type'   : options.type }
   xhelp = { 'inDir'  : '',
            'outDir' : '',
            'type'   : ' [scaffold|contig]' }
   for o in opts:
      if not opts[ o ]:
         parser.error('specify --%s%s\n' % (o, xhelp[o] ))
   dirs = { 'inDir'  : options.inDir }
   for d in dirs:
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   
   if options.type not in ( 'scaffold', 'contig' ):
      parser.error('--type must either be scaffold or contig, %s is not recognized.' % options.type )
   
   regex = '([a-zA-Z]+)(\d+)'
   pat = re.compile(regex)
   m = re.match( pat, options.name )
   if not m:
      parser.error('unable to match --name=%s to regex "%s"' % (options.name, regex))
   if len(m.group(1)) > 1:
      sys.stderr.write('Warning, --name=%s has more than a 1 letter ID, legends, '
                       'tables, axis labels may all break in other scripts.\n' % m.group(1))
   if len( options.name ) > 3:
      sys.stderr.write('Warning, --name=%s is in total more than 3 characters long, '
                       'legends, tables, axis labels may all break in other scripts.\n' % options.name )
   if os.path.exists( options.outDir ) and not os.path.isdir( options.outDir ):
      parser.error('--outDir %s is not a directory!\n' % options.outDir )

def verifyNameIsUnique( options ):
   pass

def populateDirectoryStructure( options ):
   if not os.path.exists( options.outDir ):
      os.makedirs( options.outDir )
   for d in [ 'mafsContigs', 'mafsScaffolds', 'statsContigsAggregateColumns',
              'statsScaffoldsAggregateColumns', 'statsScaffoldsContigPath', 'statsScaffoldsCopyNumber',
              'statsScaffoldsLinkage', 'statsScaffoldsSubstitions' ]:
      d = os.path.join( options.outDir, d )
      if os.path.exists( d ) and not os.path.isdir( d ):
         sys.stderr.write('%s exists but is not a directory!\n' % d )
         sys.exit(1)
      if not os.path.exists( d ):
         os.makedirs( d )

def migrate( options ):
   pass

def main():
   usage = ( 'usage: %prog --inDir=path/to/input --outDir=path/to/out --type=[contig|scaffold] [options]\n\n'
             '%prog takes in a directory created by cactus ( --inDir )\n'
             'and a directory where you are staging the data for analysis ( --outDir )\n'
             'the type of alignment, either contig or scaffold, and then migrates the\n'
             'relevant files into the outDir. Use the --name flag to specify the name\n'
             '(e.g. --name R1) of the assembly in the analysis.' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )

   verifyNameIsUnique( options )
   populateDirectoryStructure( options )

   migrate( options )


if __name__ == '__main__':
   main()
