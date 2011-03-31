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
import cPickle
from libMafGffPlot import Data
from libMafGffPlot import MafBlock
from libMafGffPlot import MafLine
from libMafGffPlot import newMafWigDict
from libMafGffPlot import objListToBinnedWiggle
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

def checkOptions( options, parser, data ):
   if options.maf == None:
      parser.error( 'Error, specify --maf.\n' )
   if not os.path.exists( options.maf ):
      parser.error( 'Error, --maf %s does not exist.\n' % options.maf )
   if options.outDir == None:
      options.outDir = os.getcwd()
   if not os.path.exists( options.outDir ):
      parser.error( 'Error, --outDir %s does not exist.\n' % options.outDir )
   if not os.path.isdir( options.outDir ):
      parser.error( 'Error, --outDir %s is not a directory.\n' % options.outDir )
   if options.ref == None:
      parser.error( 'Error, specify --referenceGenome.\n' )
   if options.other == None:
      parser.error( 'Error, specify --comparisonGenome.\n' )
   data.genomesDict = { options.ref   : True,
                        options.other : True }
   if options.numBins < 1:
      parser.error('Error, number of bins (%d) must be >= 1.' % options.numBins )
   opts = { 'chrLengths' : options.chrLengths,
            'chrNames'   : options.chrNames }
   for a in opts:
      if opts[ a ] == None:
         parser.error('Error, specify --%s.\n' % a )
   data.chrNames   = options.chrNames.split(',')
   data.chrLengths = options.chrLengths.split(',')
   data.chrLengthsByChrom = {}
   if len( data.chrLengths ) != len( data.chrNames ):
      parser.error('Error, number of elemnts in --chrLengths not equal to '
                   'number of elements in --chrNames.\n')
   for i in range( 0, len( data.chrLengths )):
      data.chrLengths[ i ] = int( data.chrLengths[ i ] )
      data.chrLengthsByChrom[ data.chrNames[ i ] ] = data.chrLengths[ i ]
   data.genomeLength = 0
   for c in data.chrLengths:
      data.genomeLength += c
   if options.numBins > data.genomeLength:
      parser.error('Error, number of bins (%d) must be '
                   '< length of genome (%d).' % ( options.numBins, data.genomeLength ))

def packData( options, data, prot='py23Bin' ):
   """ prot refers to the protocol to use.
   """
   protocols = { 'ASCII' : 0,
                 'pre23Bin' : 1,
                 'py23Bin'  : 2 }
   if options.name:
      prefix = options.name
   else:
      prefix = options.ref + '.' + options.other
   for c in data.chroms:
      f = open( os.path.join( options.outDir, prefix +
                              '.maf.' + c + '.pickle' ), 'wb' )
      cPickle.dump( data.mafWigDict[ c ], f, protocol=protocols[ prot ] )
      f.close()

def extractMafLine( line, order, pat, options, data ):
   """ parse a given line from a maf file into a 
   MafLine object.
   """
   m = re.match( pat, line )
   if m == None:
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
         sys.stderr.write( 'Error, file %s: maf block on chromosome "%s" has sequence length (%d) '
                           'that does not equal the corresponding input from --chrLengths (%d). '
                           'Line below:\n%s\n' % ( options.maf, ml.chr, ml.totalLength,
                                                   data.chrLengthsByChrom[ ml.chr ], line ))
         sys.exit( 1 )
   if ml.strand == -1:
      ml.start = ml.totalLength - ml.start + 1
   ml.sequence = m.group( 7 )
   for b in ml.sequence:
      if b == '-':
         sys.stderr.write( 'Error, file %s: maf line contains \'-\' character. Mafs are assumed '
                           'to be gapless. Bad line:\n%s\n' % ( options.maf, line ) )
         sys.exit( 1 )
   return ( ml, order )

def createMafBlockFromPair( iLine, jLine, hplList, options, data ):
   """for a pair of maf line objects, create a mafBlock and store all relevant
   information.
   """
   mb = MafBlock()
   if iLine.genome != options.ref:
      iLine, jLine = jLine, iLine
   mb.refGenome  = iLine.genome
   mb.refChr     = iLine.chr
   mb.refStart   = iLine.start + 1
   mb.refEnd     = iLine.start + iLine.strand * iLine.length
   mb.refStrand  = iLine.strand
   mb.refTotalLength = iLine.totalLength
   mb.pairGenome = jLine.genome
   mb.pairChr    = jLine.chr
   mb.pairStart  = jLine.start + 1
   mb.pairEnd    = jLine.start - 1 + jLine.strand * jLine.length
   mb.pairStrand = jLine.strand
   mb.pairTotalLength = jLine.totalLength
   if len( hplList ) > 0:
      if jLine.order < len( hplList ):
         mb.hpl        = hplList[ jLine.order ][ 'hpl' ]
         mb.hplStart   = hplList[ jLine.order ][ 'hFive' ]
         mb.hplEnd     = hplList[ jLine.order ][ 'hThree' ]
         mb.spl        = hplList[ jLine.order ][ 'spl' ]
      else:
         sys.stderr.write( 'Error, creating mafBlock but jLine.order (%d) is '
                           'greating than the length of the hpl list (%d))\n' % ( jLine.order, 
                                                                                  len( hplList ) ))
   data.mafBlocksByChrom[ mb.refChr ].append( mb )

def extractBlockPairs( mafLineList, hplList, options, data ):
   """ loop through all pairs in the block list, store
   the discovered blocks.
   Elements of blockList are MafLine() objects.
   """
   for i in range( 0, len( mafLineList )):
      if mafLineList[ i ].genome not in data.genomesDict:
         continue
      for j in range( i + 1, len( mafLineList )):
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
      if line[ 0 ] == '#':
         if line[ :4 ] == '#HPL':
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
      if line[ 0 ] == 's':
         line = line.strip()
         ( ml, order ) = extractMafLine( line, order, pat, options, data )
         if ml == None:
            sys.stderr.write( 'Error, regexp fail on file %s line: \'%s\'\n'
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
   if mafBlocksByChrom[ c ] == None:
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
      thisChrNumBins = int( ( float( data.chrLengthsByChrom[ c ] ) / data.genomeLength ) * options.numBins )
      mafWigDict[ c ] = {}
      d = mafDataOrNone( data.mafBlocksByChrom, c )
      if d == None:
         mafWigDict[ c ] = newMafWigDict( thisChrNumBins )
         for i in range( 0, thisChrNumBins ):
            mafWigDict[c]['xAxis'][ i ] = ((float( i ) / ( thisChrNumBins - 1.0 )) * 
                                           float( data.chrLengthsByChrom[ c ] ) )
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
            sys.stderr.write( 'Error, file %s has maf block on chr %s with '
                              'bounds [%d - %d] which are beyond featLen (%d)\n' %
                              ( options.maf, m.refChr, m.refStart, m.refEnd, data.chrLengthsByChrom[ c ] ))
            sys.exit( 1 )
      
def trimDups( options, data ):
   """ Walk the data.mafBlockByChrom structure and 
   look for any mafBlocks that overlap on the reference.
   if we find overlaps, either trim them back or remove 
   them.
   When looking for coverage as we are here, having more
   than one base in the assembly map to a base in the 
   reference shouldn't get you a higher score. Otherwise
   you could just have the same maf block in your .maf 10
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
         if m.refStart < prevBlock.refEnd:
            if m.refEnd > prevBlock.refEnd:
               m.refStart = prevBlock.refEnd + 1
            else:
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
   Logically this should be impossible, but I want to be able to point to this test 
   whenever we see a weird spike on the plot and say 'this is due to a visual artifact
   of the plotting library we're using.'
   """
   tests = { 'blocks': [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4',
                         'maf1e5', 'maf1e6', 'maf1e7' ],
             'paths' : [ 'maf', 'mafHpl1e2', 'mafHpl1e3', 'mafHpl1e4',
                         'mafHpl1e5', 'mafHpl1e6', 'mafHpl1e7' ] }
   for c in data.chrNames:
      for i in range( 0, len( data.mafWigDict[ c ][ 'xAxis' ] )):
         for t in tests:
            for j in range( 0, len( tests[ t ] ) - 1):
               if ( data.mafWigDict[ c ][ tests[ t ][ j ] ][ i ] < 
                    data.mafWigDict[ c ][ tests[ t ][ j + 1 ] ][ i ] ):
                  sys.stderr.write('Error, validation failed on file:%s '
                                   'chr:%s %s (%d) > %s (%d)!\n' % ( options.maf, c,
                                                                     tests[ t ][ j + 1],
                                                                     data.mafWigDict[ c ][ tests[ t ][ j + 1 ] ][i],
                                                                     tests[ t ][ j ],
                                                                     data.mafWigDict[ c ][ tests[ t ][ j ] ][i] ))
                  sys.exit( 1 )

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
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser, data )

   readMaf( options, data )
   switchToPositiveStrandCoordinates( options, data )
   for c in data.chroms:
      data.mafBlocksByChrom[ c ].sort( key = lambda x: x.refStart, reverse=False )
   trimDups( options, data )
   convertDataToWiggle( options, data )
   recordCoverage( options, data )
   
   verifyStacks( options, data )
   packData( options, data )

if __name__ == '__main__':
   main()
