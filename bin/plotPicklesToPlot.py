#!/usr/bin/env python
"""
plotPicklesToPlot.py
18 February 2011
dent earl, dearl@soe.ucsc.edu


"""
import cPickle
import glob
from libMafGffPlot import Data
from libMafGffPlot import MafBlock
from libMafGffPlot import GffRecord

import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
import numpy
from optparse import OptionParser
import os
import sys
import re
import time

def initOptions( parser ):
   parser.add_option( '--annotDir', dest='annotDir',
                      type='string',
                      help='Directory where annotation pickles will be read from.' )
   parser.add_option( '--mafDir', dest='mafDir',
                      type='string',
                      help='Directory where maf pickles will be read from.' )
   parser.add_option( '--out', dest='out', default='myPlot',
                       type='string',
                       help='output pdf where figure will be created. No extension.' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|both]' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                       type='int',
                       help='Dots per inch of the output.')
   parser.add_option( '--chrLengths', dest='chrLengths',
                      type='string',
                      help='comma separated list (no spaces) of chromosome lengths.' )
   parser.add_option( '--chrNames', dest='chrNames',
                      type='string',
                      help='comma separated list (no spaces) of chromosome names, as you want them '
                      'to appear in l-r order in the figure.' )
   parser.add_option( '--chrLabels', dest='chrLabels', default='',
                       type='string',
                       help='comma separated list (no spaces) of chromosome labels, as the will appear '
                      'in the plot.')
   parser.add_option( '--gridLinesMajor', dest='gridLinesMajor', default=0,
                       type='int',
                       help='Place thick grid lines on the plot every X many bases.' )
   parser.add_option( '--verbose', dest='isVerbose', default=False,
                      action='store_true',
                      help='Turns on verbose output.' )

def checkOptions( options, parser, data ):
   dirs = { 'annotDir' : options.annotDir,
            'mafDir'   : options.mafDir }
   for d in dirs:
      if dirs[ d ] == None:
         parser.error( 'Error, specify --%s.' % d)
      if not os.path.exists( dirs[ d ] ):
         parser.error( 'Error, --%s %s does not exist.\n' % ( d, dirs[ d ]))
      if not os.path.isdir( dirs[ d ] ):
         parser.error( 'Error, --%s %s is not a directory.\n' % ( d, dirs[ d ] ))
   opts = { 'chrLengths' : options.chrLengths,
            'chrNames'   : options.chrNames }
   for a in opts:
      if opts[ a ] == None:
         parser.error('Error, specify --%s.\n' % a )

   data.chrLengths = options.chrLengths.split(',')
   data.chrLengthsByChrom = {}
   data.chrLabelsByChrom  = {}
   data.chrNames   = options.chrNames.split(',')
   if len( data.chrLengths ) != len( data.chrNames ):
      parser.error('Error, number of elemnts in --chrLengths not equal to number of elements in --chrNames.\n')
   if options.chrLabels == '':
      data.chrLabels = data.chrNames
   else:
      data.chrLabels = options.chrLabels.split(',')
   if len( data.chrNames ) != len( data.chrLabels ):
      parser.error('Error, number of elemnts in --chrLabels not equal to number of elements in --chrNames.\n')
   
   for i in range( 0, len( data.chrLengths )):
      data.chrLengths[ i ] = int( data.chrLengths[ i ] )
      data.chrLengthsByChrom[ data.chrNames[ i ] ] = data.chrLengths[ i ]
      data.chrLabelsByChrom[ data.chrNames[ i ] ] = data.chrLabels[ i ]

   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.' % options.dpi )
   data.genomeLength = 0
   for c in data.chrLengths:
      data.genomeLength += c
   options.annotColors = { 'CDS':'#1f77b4', 'UTR':'#aec7e8',
                           'NXE':'#ff600e', 'NGE':'#ffbb78',
                           'island':'#00662c', 'tandem':'#00e32c',
                           'repeat':'#662D91' }
   
def unpackData( filename, options, data ):
   t0 = time.time()
   if not os.path.exists( filename ):
      sys.stderr.write( 'Error, %s does not exist.\n' % filename)
      sys.exit( 1 )
   f = open( filename, 'rb' )
   d = cPickle.load( f )
   f.close()
   t1 = time.time()
   return d

def loadAnnots( options, data ):
   data.annotWigDict = {}
   for c in data.chrNames:
      f = glob.glob( os.path.join( options.annotDir, '*annots.%s.pickle' % c ))
      if len( f ) == 1:
         f = f[0]
      else:
         data.gffRecordsByChrom[ c ] = None
         continue
      if not os.path.exists( f ):
         sys.stderr.write('Error, unable to locate annot file for chr %s' % c )
         sys.exit( 1 )
      data.annotWigDict[ c ] = unpackData( f, options, data )
   
def loadMafs( options, data ):
   # sort of like loadAnnots, but needs an added loop that pulls 
   # from a glob of all the mafs in the maf directory.
   data.mafWigDict = {}
   data.mafNamesDict = {}
   patStr = '\S+\.(\S+)\.maf.*'
   pat = re.compile( patStr )
   for c in data.chrNames:
      mafFiles = glob.glob( os.path.join( options.mafDir, '*maf.%s.pickle' % c ))
      data.mafWigDict[ c ] = {}
      for f in mafFiles:
         m = re.search( pat, f )
         if m == None:
            sys.stderr.write('Error, unable to find genome name in %s using regex %s\n' % f, patStr )
            sys.exit( 1 )
         name = m.group(1)
         if name not in data.mafNamesDict:
            data.mafNamesDict[ name ] = 0
         data.mafWigDict[ c ][ name ] = unpackData( f, options, data )
   for c in data.chrNames:
      for n in data.mafNamesDict:
         data.mafNamesDict[ n ] += data.mafWigDict[ c ][ n ]['columnsInBlocks']
   data.orderedMafs = sorted( data.mafNamesDict, key=lambda key: data.mafNamesDict[ key ], reverse=True)
   for n in data.orderedMafs:
      print '%s %d / %d = %.3f' % ( n, data.mafNamesDict[n], data.genomeLength, 
                                    float(data.mafNamesDict[n])/data.genomeLength )

def initImage( options ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'both':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 10), dpi=options.dpi, facecolor='w' )
   return ( fig, pdf )

def establishAxes( fig, options, data ):
   """ create one axes per chromosome
   """
   axDict = {}
   options.axLeft = 0.1
   options.axWidth = 0.85
   options.axTop = 0.95
   options.axBottom = 0.05
   options.axHeight = options.axTop - options.axBottom
   options.chrMargin = 0.02
   curXPos = options.axLeft
   for c in data.chrNames:
      w = (( data.chrLengthsByChrom[ c ] / float( data.genomeLength ) ) * 
            ( options.axWidth - ( options.chrMargin * float( len( data.chrNames ) - 1) )))
      axDict[ c ] = fig.add_axes( [ curXPos, options.axBottom,
                                    w , options.axHeight] )
      curXPos += w + options.chrMargin
      plt.box( on=False )
   setAxisLimits( axDict, options, data )
   return ( axDict )

def setAxisLimits( axDict, options, data ):
   for c in axDict:
      axDict[ c ].set_ylim( 0, 1.01 )
      axDict[ c ].set_xlim( 0, data.chrLengthsByChrom[ c ] )
      axDict[ c ].xaxis.set_major_locator( pylab.NullLocator() )
      axDict[ c ].yaxis.set_major_locator( pylab.NullLocator() )

def writeImage( fig, pdf, options, data ):
   if options.outFormat == 'pdf':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
   elif options.outFormat == 'png':
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'both':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )

def drawChrLines( ax, options, data):
   for c in data.chrNames:
      if data.chrOffsets[ c ] != 0 and data.chrOffsets[ c ] != data.genomeLength:
         ax.add_line( lines.Line2D( xdata=[ data.chrOffsets[ c ], data.chrOffsets[ c ] ], ydata=[ 0,1 ], c='#87B6F9'))

def labelAxes( fig, axDict, options, data ):
   for c in data.chrNames:
      # chromosome names
      fs = scaleFont(c, data.chrLengthsByChrom[ c ], data.genomeLength, options.axWidth, options, data )
      xPos = ((( axDict[ c ].get_position().get_points()[1][0] - 
                 axDict[ c ].get_position().get_points()[0][0] ) / 2.0 ) + 
              axDict[ c ].get_position().get_points()[0][0])
      fig.text( x=xPos, y=0.96, s= data.chrLabelsByChrom[ c ], horizontalalignment='center',
                  verticalalignment='bottom', fontsize=fs )
   print options.axHeight
   increment = options.axHeight / 40.0
   j = 0.0
   for a in [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'tandem', 'repeat']:
      yPos = options.axTop - 0.018 - j
      fig.text( x= options.axLeft - 0.02, y= yPos, s = a, 
                horizontalalignment='right', verticalalignment='center', fontsize=8 )
      j +=  increment
   j += increment
   for n in data.orderedMafs:
      yPos = options.axTop - 0.018 - j
      fig.text( x= options.axLeft - 0.02, y= yPos, s = n, 
                horizontalalignment='right', verticalalignment='center', fontsize=8 )
      j += increment
         

def scaleFont( c, chrLen, genLen, axLen, options, data):
   """ find the approximate font size that will allow a string, c, to 
   fit in a given space.
   """
   fs = 9.0
   if len(c) * float(fs)/900.0 <= ( float(chrLen)/ genLen) * axLen:
      return fs
   while fs > 1:
      if len(c) * float(fs)/900.0 <= ( float(chrLen)/ genLen) * axLen:
         return fs
      fs -= .1
   return fs

def drawAnnotations( axDict, options, data ):
   annotOrder = [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'tandem', 'repeat']
   annotOrder.reverse()
   for c in data.chrNames:
      for i in range( 0, len( data.annotWigDict[ c ]['xAxis']) ):
         j = 33.0 / 40.0
         for a in annotOrder:
            data.annotWigDict[ c ][ a ][ i ] = j + float( data.annotWigDict[ c ][ a ][ i ] ) / 40.0
            j +=  1.0 / 40.0
      for a in annotOrder:
         axDict[ c ].add_line( lines.Line2D( xdata = data.annotWigDict[ c ]['xAxis'], 
                                             ydata = data.annotWigDict[ c ][ a ], 
                                             c = options.annotColors[ a ], linewidth = 0.5 ))

def drawMafs( axDict, options, data ):
   colors = { True: ( 0.1, 0.1, 0.1 ),
              False: ( 0.4, 0.4, 0.4 ) }
   for c in data.chrNames:
      col = True
      j = 31.0 / 40.0
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis']) ):
            data.mafWigDict[ c ][ n ][ 'maf' ][ i ] = j + float( data.mafWigDict[ c ][ n ][ 'maf' ][ i ] ) / 40.0
         j -=  1.0 / 40.0
         print 'draw %s %s' % ( n, c )
         axDict[ c ].add_line( lines.Line2D( xdata = data.mafWigDict[ c ][ n ]['xAxis'], 
                                             ydata = data.mafWigDict[ c ][ n ]['maf'], 
                                             c = colors[ col ], linewidth = 0.5 ))
         col = not col

def drawGridLines( mafAx, annotAx, options, data ):
   if options.gridLinesMajor < 1:
      return
   for ax in [ mafAx, annotAx ]:
      j = -1
      for c in data.chrNames:
         j += 1
         for i in range( data.chrOffsets[ c ], data.chrOffsets[ c ] + data.chrLengths[ j ], options.gridLinesMajor ):
            ax.add_line( lines.Line2D( xdata=[i, i], ydata=[0, 1], c='pink', linewidth=0.15 ))

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser, data )
   loadAnnots( options, data )
   loadMafs( options, data )
   ( fig, pdf ) = initImage( options )
   axDict = establishAxes( fig, options, data )
   labelAxes( fig, axDict, options, data )
   #drawGridLines( mafAx, annotAx, options, data )
   drawAnnotations( axDict, options, data )
   drawMafs( axDict, options, data )

   setAxisLimits( axDict, options, data )
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
