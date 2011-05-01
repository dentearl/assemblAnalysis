#!/usr/bin/env python
""" 
import CactusResult.py
15 April 2011
dent earl, dearl (a) soe ucsc edu

Script to import a cactus result directory into an assemblathon
analysis project. See the usage string in the main() function
for details.

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
import os
import re
import string
import shutil
import sys

def initOptions( parser ):
   parser.add_option( '-i', '--inDir', dest='inDir',
                      type='string',
                      help='Cactus result directory you would like to import.' )
   parser.add_option( '-o', '--outDir', dest='outDir',
                      type='string',
                      help=('Directory where you want to place the '
                            'data. If it does not exist, it will be created.' ))
   parser.add_option('--type', dest='type',
                     type='string',
                     help=('Type of the alignment, either scaffold or contig.'))
   parser.add_option('--name', dest='name', default='R1',
                     help=('Sets the name of the input assembly. default=%default'))


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
   typeMap = {'scaffold':'Scaffolds',
              'contig':'Contigs' }
   if os.path.exists( os.path.join( options.outDir, 'mafs%s' % typeMap[options.type], options.name )):
      sys.stderr.write('Error, %s already exists in the project dir.\n' % options.name)
      sys.exit(1)
   f = open( os.path.join( options.outDir, 'mafs%s' % typeMap[options.type], options.name+'.maf' ), 'w' )
   f.close()

def populateDirectoryStructure( options ):
   if not os.path.exists( options.outDir ):
      os.makedirs( options.outDir )
   for d in [ 'mafsContigs', 'mafsScaffolds', 'statsScaffoldsAggregateColumns', 
              'statsScaffoldsContigPath', 'statsScaffoldsCopyNumber',
              'statsScaffoldsContiguity', 'statsScaffoldsSubstitutions',
              'statsContigsContigPath']:
      d = os.path.join( options.outDir, d )
      if os.path.exists( d ) and not os.path.isdir( d ):
         sys.stderr.write('%s exists but is not a directory!\n' % d )
         sys.exit(1)
      if not os.path.exists( d ):
         os.makedirs( d )

def myCopy( src, dst ):
   if not os.path.exists( src ):
      sys.stderr.write('Error, %s does not exist!\n' % src)
      sys.exit(1)
   shutil.copy( src, dst )

def migrate( options ):
   if options.type == 'scaffold':
      template = { os.path.join( options.inDir, 'annotatedPaths.maf' ): os.path.join( options.outDir, 'mafsScaffolds', options.name+'.maf' ),
                   os.path.join( options.inDir, 'copyNumberStats_0.xml' ): os.path.join( options.outDir, 'statsScaffoldsCopyNumber', '%s.copyNumber_0.xml' % options.name ),
                   os.path.join( options.inDir, 'copyNumberStats_1000.xml' ): os.path.join( options.outDir, 'statsScaffoldsCopyNumber', '%s.copyNumber_1000.xml' % options.name ),
                   os.path.join( options.inDir, 'linkageStats.xml'): os.path.join( options.outDir, 'statsScaffoldsContiguity', '%s.contiguousStats.xml' % options.name),
                   os.path.join(options.inDir, 'substitutionStats_0_0_0.xml'): os.path.join( options.outDir, 'statsScaffoldsSubstitutions', '%s.subStats.upper.xml' % options.name),
                   os.path.join( options.inDir, 'substitutionStats_1000_98_5.xml'): os.path.join( options.outDir, 'statsScaffoldsSubstitutions', '%s.subStats.lower.xml' % options.name),
                   os.path.join( options.inDir, 'pathStats.xml'): os.path.join( options.outDir, 'statsScaffoldsContigPath', '%s.pathStats.xml' % options.name),
                   os.path.join( options.inDir, 'coveragePlots', 'blockLengthsVsCoverageOfAssemblyAndHaplotypes.txt'): os.path.join( options.outDir, 'statsScaffoldsAggregateColumns', '%s.blocks_haplotypes_agg.txt' % options.name ),
                   os.path.join( options.inDir, 'coveragePlots', 'contigPathLengthsVsCoverageOfAssemblyAndHaplotypes.txt'): os.path.join( options.outDir, 'statsScaffoldsAggregateColumns', '%s.contig_path_haplotypes_agg.txt' % options.name ),
                   os.path.join( options.inDir, 'coveragePlots', 'contigLengthsVsCoverageOfAssemblyAndHaplotypes.txt'): os.path.join( options.outDir, 'statsScaffoldsAggregateColumns', '%s.contigs_haplotypes_agg.txt' % options.name ),
                   os.path.join( options.inDir, 'coveragePlots', 'blockLengthsVsCoverageOfAssemblyAndHaplotypes.txt'): os.path.join( options.outDir, 'statsScaffoldsAggregateColumns', '%s.blocks_haplotypes_agg.txt' % options.name )
         }
   elif options.type == 'contig':
      template = { os.path.join( options.inDir, 'annotatedPaths.maf' ): os.path.join( options.outDir, 'mafsContigs', options.name+'.maf' ),
                   os.path.join( options.inDir, 'pathStats.xml'): os.path.join( options.outDir, 'statsContigsContigPath', '%s.pathStats.xml' % options.name)}

   for src in template:
      myCopy( src, template[src] )

def main():
   usage = ( 'usage: %prog --inDir=path/to/input --outDir=path/to/out --type=[contig|scaffold] [options]\n\n'
             '%prog takes in a directory created by cactus ( --inDir )\n'
             'and a directory where you are staging the data for analysis ( --outDir )\n'
             'the type of alignment ( --type ), either contig or scaffold, and then migrates the\n'
             'relevant files into the outDir. Use the --name flag to specify the name\n'
             '(e.g. --name R1) of the assembly in the analysis.' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )

   populateDirectoryStructure( options )
   verifyNameIsUnique( options )

   migrate( options )


if __name__ == '__main__':
   main()
