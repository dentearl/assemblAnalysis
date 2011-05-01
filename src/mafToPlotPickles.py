#!/usr/bin/env python
"""
mafToPlotPickles.py
17 February 2011
dent earl, dearl@soe.ucsc.edu

mafToPlotPickles.py is used to produce pickle
files from pairwise maf. The pickle files are used
to create a plot by plotPicklesToPlot.py.

mafToPlotPickles.py works by reading in a maf
and storing only the relevant information needed to
create a density plot of maf coverage relative to
the reference genome. Maf files whose '--reference'
genome contain multiple chromosomes are split into
a single pickle per chromosome.

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
from libMafGffPlot import Data
from libMafGffPlot import MafBlock
from libMafGffPlot import MafLine
from libMafGffPlot import newMafWigDict
from libMafGffPlot import objListToBinnedWiggle
from libMafGffPlot import objListUtility_xAxis
from libMafGffPlot import packData
import math
import numpy
from optparse import OptionParser
import os
import sys
import re

def initOptions( parser ):
   parser.add_option( '-a', '--referenceGenome', dest='ref',
                      type='string',
                      help='Establishes the genome in the maf that will be used as the reference.' )
   parser.add_option( '--maf', dest='maf',
                      type='string',
                      help='Establishes the maf that will be read.' )
   parser.add_option( '-b', '--comparisonGenome', dest='other',
                      type='string',
                      help='Establishes the genome in the maf that will be used as the comparison.' )
   parser.add_option( '--chrNames', dest='chrNames',
                      type='string',
                      help='comma separated list (no spaces) of chromosome names.' )
   parser.add_option( '--chrLengths', dest='chrLengths',
                      type='string',
                      help='comma separated list (no spaces) of chromosome lengths.' )
   parser.add_option( '--outDir', dest='outDir',
                      type='string',
                      help='Establishes where the pickles will be written.' )
   parser.add_option( '--name', dest='name',
                      help='changes prefix name from refGenome.chrN.pickle to NAME.chrN.pickle' )
   parser.add_option( '-n', '--numBins', dest='numBins', default=8*300,
                      type='int',
                      help='Number of bins to partion the the x axis into. default=%default' )
   parser.add_option( '--verify', dest='verify', default=False,
                      action='store_true',
                      help=('Enables extra checks to verify the data structure is accurate. '
                            'Not necessary unless the output plots look odd. default=%default' ))

def checkOptions( options, parser, data ):
   if options.maf is None:
      parser.error( 'specify --maf.\n' )
   if not os.path.exists( options.maf ):
      parser.error( '--maf %s does not exist.\n' % options.maf )
   if options.outDir is None:
      options.outDir = os.getcwd()
   if not os.path.exists( options.outDir ):
      parser.error( '--outDir %s does not exist.\n' % options.outDir )
   if not os.path.isdir( options.outDir ):
      parser.error( '--outDir %s is not a directory.\n' % options.outDir )
   if options.ref is None:
      parser.error( 'specify --referenceGenome.\n' )
   if options.other is None:
      parser.error( 'specify --comparisonGenome.\n' )
   data.genomesDict = { options.ref   : True,
                        options.other : True }
   if options.numBins < 1:
      parser.error('number of bins (%d) must be >= 1.' % options.numBins )
   opts = { 'chrLengths' : options.chrLengths,
            'chrNames'   : options.chrNames }
   for a in opts:
      if opts[ a ] is None:
         parser.error('specify --%s.\n' % a )
   data.chrNames   = options.chrNames.split(',')
   data.chrLengths = options.chrLengths.split(',')
   data.chrLengthsByChrom = {}
   if len( data.chrLengths ) != len( data.chrNames ):
      parser.error('number of elemnts in --chrLengths not equal to '
                   'number of elements in --chrNames.\n')
   for i in xrange( 0, len( data.chrLengths )):
      data.chrLengths[ i ] = int( data.chrLengths[ i ] )
      data.chrLengthsByChrom[ data.chrNames[ i ] ] = data.chrLengths[ i ]
   data.genomeLength = 0
   for c in data.chrLengths:
      data.genomeLength += c
   if options.numBins > data.genomeLength:
      parser.error('number of bins (%d) must be '
                   '< length of genome (%d).' % ( options.numBins, data.genomeLength ))
   # filename to store pickle
   if options.name:
      prefix = options.name
   else:
      prefix = options.ref + '.' + options.other
   options.filename = os.path.join( options.outDir, prefix + '.maf.pickle' )

def extractMafLine( line, order, pat, options, data ):
   """ parse a given line from a maf file into a 
   MafLine object.
   """
   m = re.match( pat, line )
   if m is None:
      return ( None, order )
   ml = MafLine()
   ml.genome = m.group( 1 )
   if ml.genome not in data.genomesDict:
      return ( 'notOurGenome', order )
   if ml.genome == options.other:
      order += 1
      ml.order = order
   ml.chr = m.group( 2 )
   if ml.genome == options.ref:
      if ml.chr not in data.chroms:
         data.chroms[ ml.chr ] = True
         data.mafBlocksByChrom[ ml.chr ] = []

   ml.start  = int( m.group( 3 ) )
   ml.length = int( m.group( 4 ) )
   ml.strand = int( m.group( 5 ) + '1')
   ml.totalLength = int( m.group( 6 ) )
   if ml.genome == options.ref:
      if ml.totalLength != data.chrLengthsByChrom[ ml.chr ]:
         sys.stderr.write( 'file %s: maf block on chromosome "%s" has sequence length (%d) '
                           'that does not equal the corresponding input from --chrLengths (%d). '
                           'Line below:\n%s\n' % ( options.maf, ml.chr, ml.totalLength,
                                                   data.chrLengthsByChrom[ ml.chr ], line ))
         sys.exit( 1 )
   ml.sequence = m.group( 7 )
   for b in ml.sequence:
      if b == '-':
         sys.stderr.write( 'file %s: maf line contains \'-\' character. Mafs are assumed '
                           'to be gapless. Bad line:\n%s\n' % ( options.maf, line ) )
         sys.exit( 1 )
   return ( ml, order )

def createMafBlockFromPair( iLine, jLine, hplList, options, data ):
   """for a pair of maf line objects, create a mafBlock and store all relevant
   information.
   We change the start position from being 0 based to being 1 based here.
   Recall that the .start is always the lower position on the positive strand,
   and the .end is always the higher position on the positive strand
   """
   mb = MafBlock()
   if iLine.genome != options.ref:
      iLine, jLine = jLine, iLine
   mb.refGenome  = iLine.genome
   mb.refChr     = iLine.chr
   mb.refStrand  = iLine.strand
   if mb.refStrand == 1:
      mb.refStart  = iLine.start + 1
      mb.refEnd    = iLine.start + iLine.length
   else:
      mb.refStart  = iLine.totalLength - iLine.start - iLine.length + 1
      mb.refEnd    = iLine.totalLength - iLine.start
      mb.refStrand = 1
   mb.refTotalLength = iLine.totalLength
   mb.pairGenome = jLine.genome
   mb.pairChr    = jLine.chr
   mb.pairStrand = jLine.strand
   if mb.pairStrand == 1:
      mb.pairStart = jLine.start + 1
      mb.pairEnd   = jLine.start + jLine.length
   else:
      mb.pairStart = jLine.totalLength - jLine.start - jLine.length + 1
      mb.pairEnd   = jLine.totalLength - jLine.start
      mb.pairStrand= 1
   mb.pairTotalLength = jLine.totalLength
   if len( hplList ) > 0:
      if jLine.order < len( hplList ):
         mb.hpl        = hplList[ jLine.order ][ 'hpl' ]
         mb.hplStart   = hplList[ jLine.order ][ 'hFive' ]
         mb.hplEnd     = hplList[ jLine.order ][ 'hThree' ]
         mb.spl        = hplList[ jLine.order ][ 'spl' ]
      else:
         sys.stderr.write( 'creating mafBlock but jLine.order (%d) is '
                           'greating than the length of the hpl list (%d))\n' % ( jLine.order, len( hplList ) ))
         sys.exit(1)
   if mb.refEnd > mb.refTotalLength:
      sys.stderr.write( 'creating mafBlock but reference end is greater than total length! %d > %d %s\n' 
                        % ( mb.refEnd, mb.refTotalLength, mb.refChr ))
      sys.exit(1)
   if mb.refStart > mb.refTotalLength:
      sys.stderr.write( 'creating mafBlock but reference start is greater than total length! %d > %d %s\n' 
                        % ( mb.refStart, mb.refTotalLength, mb.refChr ))
      sys.exit(1)
   if mb.pairEnd > mb.pairTotalLength:
      sys.stderr.write( 'creating mafBlock but pair end is greater than total length! %d > %d %s\n' 
                        % ( mb.pairEnd, mb.pairTotalLength, mb.pairChr ))
      sys.exit(1)
   if mb.pairStart > mb.pairTotalLength:
      sys.stderr.write( 'creating mafBlock but pair start is greater than total length! %d > %d %s\n' 
                        % ( mb.pairStart, mb.pairTotalLength, mb.pairChr ))
      sys.exit(1)
   data.mafBlocksByChrom[ mb.refChr ].append( mb )

def extractBlockPairs( mafLineList, hplList, options, data ):
   """ loop through all pairs in the block list, store
   the discovered blocks.
   Elements of blockList are MafLine() objects.
   """
   for i in xrange( 0, len( mafLineList )):
      if mafLineList[ i ].genome not in data.genomesDict:
         continue
      for j in xrange( i + 1, len( mafLineList )):
         if mafLineList[ j ].genome not in data.genomesDict:
            continue
         if mafLineList[ i ].genome == mafLineList[ j ].genome:
            continue
         createMafBlockFromPair( mafLineList[ i ], mafLineList[ j ], hplList, options, data )

def readMaf( options, data ):
   """ read the maf, populate the mafLineList, then
   eventually things get stuffed in data.mafBlocksByChrom
   by way of extractBlockPairs() -> createMafBlockFromPair()
   """
   regex = 's\s+([\w\d\-]+?)\.([\w\d\.\+\-]+?)\s+(\d+)\s+(\d+)\s+([-+])\s+(\d+)\s+([\-actgurykmswbdhvnACTGURYKMSWBDHVN]+)'
   pat = re.compile( regex )
   mf = open( options.maf )
   mafLineList = []
   order = -1
   hplList = []
   hpl   = ''
   five  = ''
   three = ''
   for line in mf:
      if line.startswith('#HPL'):
         d = line.split(' ')
         # example line: "#HPL=12049 5=1 3=1 SPL=123412 S5=0 S3=12"
         # there will be one hpl line per options.other line
         # in blocks that contain the options.ref
         hpl    = int( d[0][5:] ) # comment at start of this field
         hFive  = int( d[1][2] )
         hThree = int( d[2][2] )
         spl    = int( d[3][4:] ) # no comment at start of this field
         hplList.append( { 'hpl': hpl, 'hFive': hFive, 
                           'hThree': hThree, 'spl': spl } )
         continue
      if line.startswith('s'):
         line = line.strip()
         ml, order = extractMafLine( line, order, pat, options, data )
         if ml is None:
            sys.stderr.write( 'regexp fail on file %s line: \'%s\'\n'
                              'Regex: \'%s\'\n' % ( options.maf, line, regex ) )
            sys.exit( 1 )
         if ml == 'notOurGenome':
            continue
         if ml.length != len( ml.sequence ):
            sys.stderr.write( 'Error while working on file %s :\n   '
                              'printed sequence length (%d) not equal to actual sequence '
                              'length (%d) ref genome:%s other genome:%s line below:\n%s\n' % 
                              ( options.maf, ml.length, len( ml.sequence ), options.ref, options.other, line ) )
            sys.exit( 1 )
         mafLineList.append( ml )
      else:
         # end of the block
         if len( mafLineList ) > 0:
            extractBlockPairs( mafLineList, hplList, options, data )
            mafLineList = []
            order = -1
            hplList = []
            hpl   = ''
            five  = ''
            three = ''
   if len( mafLineList ) > 0:
      extractBlockPairs( mafLineList, hplList, options, data )

def mafDataOrNone( mafBlocksByChrom, c ):
   """ return d which is either the list of MafBlocks
   for a given chromosome, c, or return None.
   """
   if c not in mafBlocksByChrom:
      return None
   if mafBlocksByChrom[ c ] is None:
      return None
   if len( mafBlocksByChrom[ c ] ) < 1:
      return None
   return mafBlocksByChrom[ c ]

def convertDataToWiggle( options, data ):
   """ the mafWigDict is keyed on chromosome and then
   on either maf or xAxis. The dict will be pulled apart
   by chromosome in the packData() function.
   """
   mafWigDict = {}
   for c in data.chrNames:
      thisChrNumBins = int( math.floor( ( float( data.chrLengthsByChrom[ c ] ) / 
                                          data.genomeLength ) * 
                                        options.numBins ))
      mafWigDict[ c ] = {}
      d = mafDataOrNone( data.mafBlocksByChrom, c )
      if d is None:
         mafWigDict[ c ] = newMafWigDict( thisChrNumBins )
         mafWigDict[ c ] ['xAxis'] = objListUtility_xAxis( data.chrLengthsByChrom[c], thisChrNumBins )
      else:
         mafWigDict[ c ] = objListToBinnedWiggle( d, data.chrLengthsByChrom[ c ], 
                                                  thisChrNumBins, options.maf )
   data.mafWigDict = mafWigDict

def switchToPositiveStrandCoordinates( options, data ): 
   """ The refStart and refEnd in the maf are not necessarily
   in a set order. If the strand is negative then the refEnd will be
   a smaller number than the refStart. 
   This is only confusing for our purposes, so we explicity move to
   a positive strand only coordinate system.
   """
   for c in data.mafBlocksByChrom:
      for m in data.mafBlocksByChrom[ c ]:
         if m.refStart > m.refEnd:
            m.refStart, m.refEnd = m.refEnd, m.refStart
            m.refStrand *= -1
            m.hplStart, m.hplEnd = m.hplStart, m.hplEnd # this is now left-right draw order
         # sanity check
         if m.refStart > data.chrLengthsByChrom[ c ] or m.refEnd > data.chrLengthsByChrom[ c ]:
            sys.stderr.write( 'file %s has maf block on chr %s with '
                              'bounds [%d - %d] which are beyond featLen (%d)\n' %
                              ( options.maf, m.refChr, m.refStart, m.refEnd, data.chrLengthsByChrom[ c ] ))
            sys.exit( 1 )
      
def trimDups( options, data ):
   """ Walk the data.mafBlockByChrom structure and 
   look for any mafBlocks that overlap on the reference.
   if we find overlaps, either trim them back or remove 
   them.
   When looking for coverage, as we are here, having more
   than one base in the assembly map to a base in the 
   reference shouldn't get you a higher score. Otherwise
   you could just have the same maf block in your maf 10
   times and you'd get a higher score than just having it
   in once.
   """
   for c in data.chrNames:
      prevBlock = MafBlock()
      replacement = []
      if c not in data.mafBlocksByChrom:
         data.mafBlocksByChrom[ c ] = replacement
         continue
      for m in data.mafBlocksByChrom[ c ]:
         if m.refStart <= prevBlock.refEnd:
            if m.refEnd > prevBlock.refEnd:
               # only add in the new, distinct, bases
               m.refStart = prevBlock.refEnd + 1
            else:
               # this block is totally covered by the previous block
               continue
         replacement.append( m )
         prevBlock = m
      data.mafBlocksByChrom[ c ] = replacement

def recordCoverage( options, data ):
   """ walks the mafWigDict and figures out out of the 
   total length of the genome, how much of it is covered
   by alignments. Stores this under the key 'coverage' as
   a float.
   """
   for c in data.chrNames:
      data.mafWigDict[ c ]['columnsInBlocks'] = 0
      for m in data.mafBlocksByChrom[ c ]:
         if m.refEnd > m.refStart:
            data.mafWigDict[ c ]['columnsInBlocks'] += ( m.refEnd + 1 ) - m.refStart
         else:
            data.mafWigDict[ c ]['columnsInBlocks'] += ( m.refStart + 1 ) - m.refEnd

def verifyStacks( options, data ):
   """ For both blocks and paths, the stacked data structure must be a proper stack.
   That is, there must be a monotonically decreasing count in the categories.
   Logically this should be impossible to fail, but I want to be able to point to this test 
   whenever we see a weird spike on the plot and say 'this is due to a visual artifact
   of the plotting library we're using.'
   """
   tests = { 'blocks' : [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4',
                         'maf1e5', 'maf1e6', 'maf1e7' ],
             'paths'  : [ 'maf', 'mafCpl1e2', 'mafCpl1e3', 'mafCpl1e4',
                         'mafCpl1e5', 'mafCpl1e6', 'mafCpl1e7' ],
             'contigs': [ 'maf', 'mafCtg1e2', 'mafCtg1e3', 'mafCtg1e4',
                          'mafCtg1e5', 'mafCtg1e6', 'mafCtg1e7' ],
             'scaffolds':[ 'maf', 'mafSpl1e2', 'mafSpl1e3', 'mafSpl1e4',
                           'mafSpl1e5', 'mafSpl1e6', 'mafSpl1e7' ] }
   for c in data.chrNames:
      for t in tests:
         for i in xrange( 0, len( tests[t] ) - 1):
            if sum( data.mafWigDict[ c ][ tests[t][ i ]] < 
                    data.mafWigDict[ c ][ tests[t][ i + 1]] ):
               sys.stderr.write('stack validation failed on file:%s '
                                'chr:%s %s %s > %s !\n' 
                                % ( options.maf, c, tests[ t ][ i + 1 ], tests[ t ][ i ] ))
               sys.exit( 1 )
   sys.stderr.write('Verify monotonically decreasing property for all categories: OK.\n')

def verifyElements( options, data ):
   """ The largest value in any of the arrays should never be greater than 1.0.
   """
   for c in data.chrNames:
      for t in [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4',
                 'maf1e5', 'maf1e6', 'maf1e7', 'mafCpl1e2', 
                 'mafCpl1e3', 'mafCpl1e4', 'mafCpl1e5', 
                 'mafCpl1e6', 'mafCpl1e7', 'mafCtg1e2', 
                 'mafCtg1e3', 'mafCtg1e4', 'mafCtg1e5', 
                 'mafCtg1e6', 'mafCtg1e7', 'mafSpl1e2', 
                 'mafSpl1e3', 'mafSpl1e4', 'mafSpl1e5', 
                 'mafSpl1e6', 'mafSpl1e7' ]:
         if sum( data.mafWigDict[ c ][ t ] > 1.0 ) > 0:
            sys.stderr.write('element validation failed on file:%s '
                             'chr:%s type:%s contains values '
                             'greater than 1.0\n' % ( options.maf, c, t ))
            sys.exit( 1 )
   sys.stderr.write('Verify Elements <= 1.0: OK\n')

def verifyDistinct( options, data ):
   """ There should be no duplicate bases in any of the maf blocks.
   They should have all been trimed out by trimDups().
   """
   tot = 0
   for c in data.chrNames:
      s = set()
      d = mafDataOrNone( data.mafBlocksByChrom, c )
      if d is None:
         continue
      for mb in d:
         for i in xrange( mb.refStart, mb.refEnd + 1):
            if i in s:
               sys.stderr.write('duplicate base found! %s %d [%d-%d], %s [%d-%d]\n'
                                % (mb.refChr, i, mb.refStart, mb.refEnd, 
                                   mb.pairChr, mb.pairStart, mb.pairEnd ))
               sys.exit( 1 )
            else:
               s.add( i )
      tot += len( s )
   sys.stderr.write( 'Verify all bases sent to be binned are distinct: Found %s distinct bases in the alignment to the reference genome, no duplicates, OK.\n' % tot)

def verifyLengths( options, data ):
   """ The lengths of the arrays should all be the same.
   """
   types = [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4',
             'maf1e5', 'maf1e6', 'maf1e7', 'mafCpl1e2', 
             'mafCpl1e3', 'mafCpl1e4', 'mafCpl1e5', 
             'mafCpl1e6', 'mafCpl1e7', 'mafCtg1e2', 
             'mafCtg1e3', 'mafCtg1e4', 'mafCtg1e5', 
             'mafCtg1e6', 'mafCtg1e7', 'mafSpl1e2', 
             'mafSpl1e3', 'mafSpl1e4', 'mafSpl1e5', 
             'mafSpl1e6', 'mafSpl1e7', 'xAxis',
             'mafCpEdgeCount', 'mafCpErrorCount', 
             'mafCpScafGapCount', 'blockEdgeCount' ]
   if len( data.chrNames ) != len( data.mafWigDict ): 
      sys.stderr.write('the expected length of the data wig '
                       'dictionary is %d (i.e. number of chromosomes), but actual is %d\n' 
                       % ( len( data.chrNames ), len( data.mafWigDict )))
      sys.exit( 1 )
   for c in data.chrNames:
      if len( types ) + 5 != len( data.mafWigDict[c] ): # extra 5 are from the *Max records
         sys.stderr.write('the expected length of the data wig '
                          'dictionary for %s is %d, but actual is %d\n' 
                          % ( c, len( types ) + 5, len( data.mafWigDict[c] )))
         sys.stderr.write( '%s\n' % str( data.mafWigDict[ c ].keys() ))
         sys.exit( 1 )
   sys.stderr.write('Verify number of records in data structure = %d, OK.\n' % (len(types) + 4))
   for c in data.chrNames:
      for i in xrange(0, len( types ) - 1):
         if len( data.mafWigDict[c][ types[i] ] ) != len( data.mafWigDict[c][ types[i+1] ]):
            sys.stderr.write('the lengths of all vectors must the '
                             'same for a given chromosome. %s, %s (%d) != %s (%d)\n' 
                             % ( c, types[i], len(data.mafWigDict[c][types[i]]), 
                                 types[i+1], len(data.mafWigDict[c][types[i+1]]) ))
            sys.exit( 1 )
      sys.stderr.write('Verify length of records in data structure for chr %s are all %d, OK.\n' 
                       % ( c, len(data.mafWigDict[c][ types[0] ])))
   sys.stderr.write('Verify lengths of arrays inside data structure, OK.\n')

def main():
   usage = ( 'usage: %prog --maf=file.maf --referenceGenome=A --comparisonGenome=B --chrNames=c0,c1,... --chrLengths=N1,N2,... --outDir=path/to/dir/\n\n'
             '%prog takes in a maf filename ( --maf ), a reference genome name as\n'
             'it appears in the maf ( --referenceGenome ), a genome name to compare against\n'
             '( --comparisonGenome ), a paired list of chromosome names ( --chrNames comma\n'
             'separated ) and chromosome lengths ( --chrLengths comma separated ) and a path\n'
             'to a directory where the maf pickles will be written ( --outDir ), one pickle\n'
             'per chromosome.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser, data )

   readMaf( options, data )
   switchToPositiveStrandCoordinates( options, data )
   
   for c in data.chroms:
      data.mafBlocksByChrom[ c ].sort( key = lambda x: x.refStart, reverse=False )
   trimDups( options, data )
   if options.verify:
      verifyDistinct( options, data )
   
   convertDataToWiggle( options, data )
   recordCoverage( options, data )
   
   if options.verify:
      verifyStacks( options, data )
      verifyElements( options, data )
      verifyLengths( options, data )
   
   packData( data.mafWigDict, options.filename, options )

if __name__ == '__main__':
   main()
