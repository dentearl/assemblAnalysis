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
from libMafGffPlot import Data
from libMafGffPlot import GffRecord
from libMafGffPlot import objListToBinnedWiggle
from libMafGffPlot import packData
import numpy
from optparse import OptionParser
import os
import sys

def initOptions( parser ):
   parser.add_option( '--gff', dest='gff',
                      type='string',
                      help='Establishes the gff that will be read.' )
   parser.add_option( '--outDir', dest='outDir',
                      type='string',
                      help='Establishes where the pickles will be written.' )
   parser.add_option( '--chrNames', dest='chrNames',
                      type='string',
                      help='comma separated list (no spaces) of chromosome names. ' )
   parser.add_option( '--chrLengths', dest='chrLengths',
                      type='string',
                      help='comma separated list (no spaces) of chromosome lengths.' )
   parser.add_option( '--prefix', dest='prefix',
                      type='string',
                      help='Allows you to specify the prefix on the file names, '
                      'i.e., prefix.annots.chrX.pickle ')
   parser.add_option( '-n', '--numBins', dest='numBins', default=8*300,
                      type='int',
                      help='Number of bins to partion the the x axis into. default=%default' )
   parser.add_option( '--verbose', dest='isVerbose', default=False,
                      action='store_true',
                      help='Turns on verbose output. default=%default' )

def checkOptions( options, parser, data ):
   if options.gff is None:
      parser.error( 'specify --gff.\n' )
   if not os.path.exists( options.gff ):
      parser.error( '--gff %s does not exist.\n' % options.gff )
   if options.outDir is None:
      options.outDir = os.getcwd()
   if not os.path.exists( options.outDir ):
      parser.error( '--outDir %s does not exist.\n' % options.outDir )
   if not os.path.isdir( options.outDir ):
      parser.error( '--outDir %s is not a directory.\n' % options.outDir )
   if options.prefix != None:
      options.prefix = options.prefix + '.'
   else:
      options.prefix = ''
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
      parser.error('number of elemnts in --chrLengths not equal to number of elements in --chrNames.\n')
   for i in range( 0, len( data.chrLengths )):
      data.chrLengths[ i ] = int( data.chrLengths[ i ] )
      data.chrLengthsByChrom[ data.chrNames[ i ] ] = data.chrLengths[ i ]
   data.genomeLength = 0
   for c in data.chrLengths:
      data.genomeLength += c
   if options.numBins > data.genomeLength:
      parser.error('number of bins (%d) must be < length of genome (%d).' % ( options.numBins, data.genomeLength ))
   options.filename = os.path.join( options.outDir, options.prefix + 'annots.pickle' )

def readGff( options, data ):
   """ read the gff file, parse each line into
   a gff record, store the record first by chromosome,
   then by annotation type.
   """
   data.chromsInGff = {}
   for c in data.chrNames:
      data.gffRecordsByChrom[ c ] = []
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
         data.gffRecordsByChrom[ r.chr ] = []
      r.source = t[ 1 ]
      r.type   = t[ 2 ]
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
      data.gffRecordsByChrom[ r.chr ].append( r )

def annotDataOrNone( gffRecordsByChrom, c ):
   """ return d which is either the list of GffRecords
   for a given chromosome, c, or return None.
   """
   if c not in gffRecordsByChrom:
      return None
   if gffRecordsByChrom[ c ] is None:
      return None
   if len( gffRecordsByChrom[ c ] ) < 1:
      return None
   return gffRecordsByChrom[ c ]

def convertDataToWiggle( options, data ):
   """ the annotWigDict is keyed on chromosome, then on annotationCount,
   annotationMax and xAxis. This will be pulled apart by chr keys in the 
   packData() function.
   """
   annotWigDict = {}
   for c in data.chrNames:
      thisChrNumBins = int( ( float( data.chrLengthsByChrom[ c ] ) / data.genomeLength ) * options.numBins )
      annotWigDict[ c ] = {}
      d = annotDataOrNone( data.gffRecordsByChrom, c )
      if d is None:
         annotWigDict[ c ]['xAxis'] = numpy.zeros( shape = ( thisChrNumBins ))
         for t in [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'tandem', 'repeat' ]:
            annotWigDict[ c ][ t + 'Count' ] = numpy.zeros( shape = ( thisChrNumBins ))
            annotWigDict[ c ][ t + 'Max' ]   = 0

         for i in range( 0, thisChrNumBins ):
            annotWigDict[c]['xAxis'][ i ] = ((float( i ) / ( thisChrNumBins - 1.0 )) * 
                                             float( data.chrLengthsByChrom[ c ] ) )
      else:
         annotWigDict[ c ] = objListToBinnedWiggle( d, data.chrLengthsByChrom[ c ], thisChrNumBins, options.gff )
   data.annotWigDict = annotWigDict

def verbosePrint( s, options, data ):
   if options.isVerbose:
      print s

def main():
   usage = ( 'usage: %prog --gff=file.gff --outDir=path/to/dir/ --chrLengths=N1,N2,... --chrNames=A,B,...\n\n'
             '%prog takes in a gff file ( --gff ), an output directory ( --outDir ), and\n'
             'pairs of chromosome names ( --chrNames comma separated ) and chromosome \n'
             'lengths ( --chrLengths comma separated ) and returns one annotation wig\n'
             'pickle per chromosome in the output directory.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser, data )
   readGff( options, data )
   for c in data.chrNames:
      data.gffRecordsByChrom[ c ].sort( key = lambda x: x.start, reverse=False )
   convertDataToWiggle( options, data )
   
   packData( data.annotWigDict, options.filename, options )

if __name__ == '__main__':
   main()
