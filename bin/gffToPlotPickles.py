#!/usr/bin/env python
"""
gffToPlotPickles.py
17 February 2011
dent earl, dearl@soe.ucsc.edu

gffToPlotPickles.py is used to produce pickle
files from a gff annotation file. The pickle files 
are used to create a plot by plotPicklesToPlot.py.

gffToPlotPickles.py works by reading in a gff
and storing only the relevant information needed to
create a density plot of annotations relative to
the reference genome. GFF files which multiple 
chromosomes are split into a single pickle per 
chromosome.

"""
import cPickle
from libMafGffPlot import Data
from libMafGffPlot import GffRecord
from libMafGffPlot import objListToBinnedWiggle
import numpy
from optparse import OptionParser
import os
import sys
import re

def initOptions( parser ):
   parser.add_option( '--gff', dest='gff',
                      type='string',
                      help='Establishes the gff that will be read.' )
   parser.add_option( '--outDir', dest='outDir',
                      type='string',
                      help='Establishes where the pickles will be written.' )
   parser.add_option( '--prefix', dest='prefix',
                      type='string',
                      help='Allows you to specify the prefix on the file names, '
                      'i.e., prefix.annots.chrX.pickle ')
   parser.add_option( '-n', '--numBins', dest='numBins', default=8*300,
                      type='int',
                      help='Number of bins to partion the the x axis into.' )
   parser.add_option( '--chrNames', dest='chrNames',
                      type='string',
                      help='comma separated list (no spaces) of chromosome names, as you want them '
                      'to appear in l-r order in the figure.' )
   parser.add_option( '--chrLengths', dest='chrLengths',
                      type='string',
                      help='comma separated list (no spaces) of chromosome lengths.' )
   parser.add_option( '--verbose', dest='isVerbose', default=False,
                      action='store_true',
                      help='Turns on verbose output.' )

def checkOptions( options, parser, data ):
   if options.gff == None:
      parser.error( 'Error, specify --gff.\n' )
   if not os.path.exists( options.gff ):
      parser.error( 'Error, --gff %s does not exist.\n' % options.gff )
   if options.outDir == None:
      options.outDir = os.getcwd()
   if not os.path.exists( options.outDir ):
      parser.error( 'Error, --outDir %s does not exist.\n' % options.outDir )
   if not os.path.isdir( options.outDir ):
      parser.error( 'Error, --outDir %s is not a directory.\n' % options.outDir )
   if options.prefix != None:
      options.prefix = options.prefix + '.'
   else:
      options.prefix = ''
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
      parser.error('Error, number of elemnts in --chrLengths not equal to number of elements in --chrNames.\n')
   for i in range( 0, len( data.chrLengths )):
      data.chrLengths[ i ] = int( data.chrLengths[ i ] )
      data.chrLengthsByChrom[ data.chrNames[ i ] ] = data.chrLengths[ i ]
   data.genomeLength = 0
   for c in data.chrLengths:
      data.genomeLength += c
   if options.numBins > data.genomeLength:
      parser.error('Error, number of bins (%d) must be < length of genome (%d).' % ( options.numBins, data.genomeLength ))

def packData( options, data, prot='py23Bin' ):
   """ prot refers to the protocol to use.
   """
   protocols = { 'ASCII' : 0,
                 'pre23Bin' : 1,
                 'py23Bin'  : 2 }
   for c in data.chrNames:
      f = open( os.path.join( options.outDir, options.prefix + 'annots.' + c + 
                              '.pickle' ), 'wb' )
      cPickle.dump( data.annotWigDict[ c ], f, protocol = protocols[ prot ] )
      f.close()

def readGff( options, data ):
   """ read the gff file, parse each line into
   a gff record, store the record first by chromosome,
   then by annotation type.
   """
   data.chromsInGff = {}
   for c in data.chrNames:
      data.gffRecordsByChrom[ c ] = {}
   gf = open( options.gff )
   for line in gf:
      line = line.strip()
      if line[0] == '#':
         # comments
         continue
      t = line.split('\t')
      if len( t ) != 9 and len( t ) != 8:
         continue
      r = GffRecord()
      r.chr = t[ 0 ]
      if r.chr not in data.chromsInGff:
         data.chromsInGff[ r.chr ] = True
         data.gffRecordsByChrom[ r.chr ] = {}
      r.source = t[ 1 ]
      r.type   = t[ 2 ]
      if r.type not in data.gffRecordsByChrom[ r.chr ]:
         data.gffRecordsByChrom[ r.chr ][ r.type ] = []
      r.start  = int( t[ 3 ] )
      r.end    = int( t[ 4 ] )
      if t[ 5 ] != '.':
         r.score  = int( t[ 5 ] )
      else:
         r.score = '.'
      r.strand = t[ 6 ]
      r.frame  = t[ 7 ]
      # we don't need this information, it just takes up space
      #if len( t ) == 9:
      #   r.group  = t[ 8 ]
      data.gffRecordsByChrom[ r.chr ][ r.type ].append( r )

def annotDataOrNone( gffRecordsByChrom, c, a ):
   """ return d which is either the list of GffRecords
   for a given chromosome, c, and a given annotation type,
   a, or return None.
   """
   if c not in gffRecordsByChrom:
      d = None
   else:
      if gffRecordsByChrom[ c ] == None:
         d = None
      else:
         if a not in gffRecordsByChrom[ c ]:
            d = None
         else:
            d = gffRecordsByChrom[ c ][ a ]
   return d

def convertDataToWiggle( options, data ):
   """ the annotWigDict is keyed first on chromosome, then on annotation or
   xAxis. This will be pulled apart by chr in the packData() function.
   """
   annotWigDict = {}
   annotOrder = [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'tandem', 'repeat']
   annotOrder.reverse() # reversing it will allow it to appear in the 'correct' top-bot order
   for c in data.chrNames:
      thisChrNumBins = int( ( float( data.chrLengthsByChrom[ c ] ) / data.genomeLength ) * options.numBins )
      annotWigDict[ c ] = {}
      annotWigDict[ c ]['xAxis'] = numpy.zeros( shape = thisChrNumBins )
      for a in annotOrder:
         d = annotDataOrNone( data.gffRecordsByChrom, c, a )
         annotWigDict[ c ][ a ] = objListToBinnedWiggle( d, data.chrLengthsByChrom[ c ], thisChrNumBins, options.gff )
      for j in range( 0, thisChrNumBins ):
         annotWigDict[ c ]['xAxis'][ j ] = ( float( j ) / ( thisChrNumBins - 1 )) * data.chrLengthsByChrom[ c ]
   data.annotWigDict = annotWigDict

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser, data )
   readGff( options, data )
   for c in data.chrNames:
      for a in data.gffRecordsByChrom[ c ]:
         data.gffRecordsByChrom[ c ][ a ].sort( key = lambda x: x.start, reverse=False )
   convertDataToWiggle( options, data )
   packData( options, data )

if __name__ == '__main__':
   main()
