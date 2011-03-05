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
import math
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
   parser.add_option( '-a', '--referenceGenome', dest='ref',
                      type='string',
                      help='Establishes the genome in the maf that will be used as the reference.' )
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
   parser.add_option( '--forceOrder', dest='forceOrder', 
                      help='Specify either the complete ordering of the assemblies or '
                      'a partial ordering. In the case of a partial ordering, the remaining '
                      'assemblies will be listed in name sorted order.' )
   parser.add_option( '--fill', dest='fill', default=False,
                      action='store_true',
                      help='Turns on the fill color for the coverage wiggles. Useful for viewing '
                      'highly variable coverage alignments.')
   parser.add_option( '--stackFillBlocks', dest='stackFillBlocks', default=False,
                      action='store_true',
                      help='Turns on the fill color for the block wiggles. Shows different coverage '
                      'thresholds in different colors. Thresholds: 0, 1e2, 1e3,...,1e7.')
   parser.add_option( '--stackFillHapPaths', dest='stackFillHapPaths', default=False,
                      action='store_true',
                      help='Turns on the fill color for the hap path wiggles. Shows different coverage '
                      'thresholds in different colors. Thresholds: 0, 1e2, 1e3,...,1e7.')
   parser.add_option( '--blockEdgeDensity', dest='blockEdgeDensity', default=False,
                      action='store_true',
                      help='Turns on the wiggle track that shows relative density of block edges. ' )
   parser.add_option( '--hapPathEdgeDensity', dest='hapPathEdgeDensity', default=False,
                      action='store_true',
                      help='Turns on the wiggle track that shows relative density of haplotype path edges. ' )
   parser.add_option( '--hapPathErrorDensity', dest='hapPathErrorDensity', default=False,
                      action='store_true',
                      help=( 'Turns on the wiggle track that shows relative density of haplotype '
                             'path segment adjacency errors. ' ))
   parser.add_option( '--relative', dest='relative', default=False,
                      action='store_true',
                      help='Plots errors as relative to the genome max. Otherwise is to global genomes max.' )
   parser.add_option( '--transform', dest='transform', default=False,
                      action='store_true',
                      help=('Transform the block and haplotype errors by y^(1/4) to '
                            'reduce spikes from extreme values.' ))
   parser.add_option( '--zerosToNan', dest='zerosToNan', default=False,
                      action='store_true',
                      help=( 'Transform the errors that are 0 to NaNs, so they are not plotted.'
                             'Creates discontinuous plots.' ))
   parser.add_option( '--edgeErrorCeiling', dest='edgeErrorCeiling',
                      type='int',
                      help=( 'Changes the edge and error scales to go from 0 to --edgeErrorCeiling '
                             'and *clip* values greater.' ))
   parser.add_option( '--frames', dest='frames', default=False,
                      action='store_true',
                      help='Debug option, turns on the plotting of all axes frame boxes.' )

def checkOptions( options, parser, data ):
   if options.ref == None:
      parser.error( 'Error, specify --referenceGenome.\n' )
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
   if ( options.stackFillBlocks or options.stackFillHapPaths ) and options.fill:
      parser.error('Error, specify either --stackFillBlocks or --stackFillHapPaths or --fill, not more than one.\n')
   if options.stackFillBlocks and options.stackFillHapPaths:
      parser.error('Error, specify either --stackFillBlocks or --stackFillHapPaths not more than one.\n')
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
                           'island':'#00662c', 'repeat':'#00e32c',
                           'tandem':'#662D91' }
   if options.out[-4:] == '.png' or options.out[-4:] == '.pdf':
      options.out = options.out[:-4]
   data.stackFillColors = [ ( '#17becf' ), # dark blue
                            ( '#9edae5' ), # light blue
                            ( '#9467bd' ), # dark purple
                            ( '#c5b0d5' ), # light purple
                            ( '#7f7f7f' ), # dark gray
                            ( '#c7c7c7' ), # light gray
                            ( '#ff7f0e' ), # bright orange
                            ( '#ffbb78' )  # light orange
                            ]
   data.stackFillColors = [ ( '#4B4C5E' ), # dark slate gray
                            ( '#9edae5' ), # light blue 
                            ( '#7F80AB' ), # purple-ish slate blue
                            ( '#c7c7c7' ), # light gray
                            ( '#ff7f0e' ), # bright orange
                            ( '#ffbb78' ), # light orange
                            ( '#9467bd' ), # dark purple
                            ( '#c5b0d5' )  # light purple
                            ]
   
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
      f = glob.glob( os.path.join( options.annotDir, '%s.annots.%s.pickle' % ( options.ref, c )))
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
      mafFiles = glob.glob( os.path.join( options.mafDir, '%s*maf.%s.pickle' % ( options.ref, c )))
      data.mafWigDict[ c ] = {}
      for f in mafFiles:
         m = re.search( pat, f )
         if m == None:
            sys.stderr.write('Error, unable to find genome name in filename %s using regex %s\n' % f, patStr )
            sys.exit( 1 )
         name = m.group(1)
         if name not in data.mafNamesDict:
            data.mafNamesDict[ name ] = 0 # this serves the duel purpose of storing 
                                          # all seen names and the count of bases aligned
         data.mafWigDict[ c ][ name ] = unpackData( f, options, data )
   for c in data.chrNames:
      for n in data.mafNamesDict:
         data.mafNamesDict[ n ] += data.mafWigDict[ c ][ n ]['columnsInBlocks']
   if not options.forceOrder:
      data.orderedMafs = sorted( data.mafNamesDict, key=lambda key: data.mafNamesDict[ key ], reverse=True )
   else:
      spokenFor = {}
      data.orderedMafs = options.forceOrder.split(',')
      for n in data.orderedMafs:
         spokenFor[ n ] = True
      sortNames = sorted( data.mafNamesDict, key=lambda key: key )
      for n in sortNames:
         if n not in spokenFor:
            data.orderedMafs.append( n )
   data.numberOfMafs = len( data.mafNamesDict )
   data.numRows = data.numberOfMafs + 7.0 # 12.0 + 10.0 # 55.0 # data.numberOfMafs + 10 # number of total rows in the figure
   
   # discover which size categories are absent from all datasets... used in legend plotting
   labs = [ '1e2', '1e3', '1e4',
         '1e5', '1e6', '1e7' ]
   data.lengthThresholdPresent = {}
   for c in data.chrNames:
      for n in data.mafWigDict[ c ]:
         for l in labs:
            if l in data.lengthThresholdPresent:
               continue
            key = 'mafHpl' + l
            if key not in data.mafWigDict[ c ][ n ]:
               print 'thats weird, this key: %s is not in file: %s chr: %s' % ( key, n, c )
               continue
            if sum( data.mafWigDict[ c ][ n ][ key ] ) > 0:
               data.lengthThresholdPresent[ l ] = True

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'both':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   figHeight = ( data.numberOfMafs + 6.5 ) / 4.0 
   fig = plt.figure( figsize=( 8, figHeight ), dpi=options.dpi, facecolor='w' )
   return ( fig, pdf )

def establishAxes( fig, options, data ):
   """ create one axes per chromosome
   """
   axDict = {}
   options.axLeft = 0.1
   options.axWidth = 0.88
   options.axTop = 0.98
   options.axBottom = 0.08
   options.axHeight = options.axTop - options.axBottom
   options.chrMargin = 0.02
   data.footerAx = fig.add_axes( [ 0.02, 0.01, 0.96, options.axBottom - 0.02] )
   if not options.frames:
      plt.box( on=False )
   curXPos = options.axLeft
   data.labelAx = fig.add_axes( [ 0.02, options.axBottom, 0.08, options.axHeight] )
   if not options.frames:
      plt.box( on=False )
   for c in data.chrNames:
      w = (( data.chrLengthsByChrom[ c ] / float( data.genomeLength ) ) * 
            ( options.axWidth - ( options.chrMargin * float( len( data.chrNames ) - 1) )))
      axDict[ c ] = fig.add_axes( [ curXPos, options.axBottom,
                                    w , options.axHeight] )
      curXPos += w + options.chrMargin
      if not options.frames:
         plt.box( on=False )
   setAxisLimits( axDict, options, data )
   return ( axDict )

def setAxisLimits( axDict, options, data ):
   for ax in [ data.labelAx, data.footerAx ]:
      ax.set_ylim( 0.0, 1.01 )
      ax.set_xlim( 0.0, 1.0 )
      ax.xaxis.set_major_locator( pylab.NullLocator() )
      ax.yaxis.set_major_locator( pylab.NullLocator() )
   for c in axDict:
      axDict[ c ].set_ylim( 0.0, 1.01 )
      axDict[ c ].set_xlim( 0.0, data.chrLengthsByChrom[ c ] )
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
   data.annotYPos = []
   data.mafYPos = []
   for c in data.chrNames:
      # chromosome names
      fs = scaleFont(c, data.chrLengthsByChrom[ c ], data.genomeLength, options.axWidth, options, data )
      # get_position() returns a BBox object, and we can get out the bounding points with get_points()
      xPos = ((( axDict[ c ].get_position().get_points()[1][0] - 
                 axDict[ c ].get_position().get_points()[0][0] ) / 2.0 ) + 
              axDict[ c ].get_position().get_points()[0][0])
      fig.text( x=xPos, y=options.axTop + 0.005, s = data.chrLabelsByChrom[ c ], horizontalalignment='center',
                  verticalalignment='bottom', fontsize=fs )
      fig.text( x=xPos + 0.025, y=options.axTop + 0.005, 
                s = '(%s)' % prettyPrintLength( data.chrLengthsByChrom[ c ] ), 
                horizontalalignment='left',
                verticalalignment='bottom', 
                color= (0.5, 0.5, 0.5,), fontsize=6 )
   data.increment = 1.0 / float( data.numRows )
   # value of 1.0 will plot track tops and bottoms on top of each other
   # value of 0.9 will have a small amount of margin between tracks
   # value of 1.1 will have overlap of tracks onto one another
   j = 0.02
   for a in [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'repeat']: # 'tandem'
      yPos =  1.0 - j
      data.labelAx.text( x= 0.8, y= yPos, s = a, 
                         horizontalalignment='right', verticalalignment='bottom', fontsize=8 )
      data.annotYPos.append( yPos )
      j += data.increment
   j += data.increment / 2.0
   for n in data.orderedMafs:
      yPos =  1.0 - j
      data.labelAx.text( x= 0.8, y= yPos + data.increment/3.0, s = '%s' % n, 
                         horizontalalignment='right', verticalalignment='bottom', fontsize=7 )
      data.labelAx.text( x= 0.45, y= yPos + data.increment/3.0, s = '%.4f' % ( float( data.mafNamesDict[ n ]) / data.genomeLength ), 
                         horizontalalignment='right', verticalalignment='bottom', fontsize=7, 
                         color=(0.5, 0.5, 0.5) )
      data.mafYPos.append( yPos )
      j += data.increment

def drawLegend( options, data ):
   # LEGEND
   if options.stackFillBlocks or options.stackFillHapPaths:
      data.footerAx.text( x=0.22, y = 0.75, horizontalalignment='right',
                          verticalalignment = 'center',
                          s = 'Fill Color Key', fontsize = 8)
      data.footerAx.text( x=0.22, y = 0.56, horizontalalignment='right',
                          verticalalignment = 'top',
                          s = 'Item >=', fontsize = 7)
      xPos = 0.2
      xunit = 0.05
      labs = [ '0', '1e2', '1e3', '1e4',
               '1e5', '1e6', '1e7', '1e8' ]
      i = -1
      for col in data.stackFillColors:
         xPos += xunit
         i += 1
         if labs[ i ] not in data.lengthThresholdPresent and labs[ i ] != '0':
            # only print colored boxes and names for data present in the figure.
            continue
         data.footerAx.add_patch( patches.Rectangle( xy = ( xPos, 0.6 ), width = 0.05, height = 0.3,
                                                     color= col, edgecolor=None, linewidth=0.0 ))
         data.footerAx.text( x=xPos + xunit/2.0, y=0.56, s= labs[i], horizontalalignment='center',
                             verticalalignment='top', fontsize = 7 )
   if not options.stackFillBlocks and not options.fill and not options.stackFillHapPaths:
      data.footerAx.text( x=0.9, y = 0.5, horizontalalignment='right',
                          verticalalignment = 'center',
                          s = 'Coverage', fontsize = 8 )
      # baseline
      data.footerAx.add_line( lines.Line2D( xdata=[0.91, 0.97],
                                            ydata=[ 0.45, 0.45],
                                            color= ( 0.8, 0.8, 0.8),
                                            linewidth=0.3) )
      # maf demo line
      data.footerAx.add_line( lines.Line2D( xdata = [ 0.91, 0.92, 0.93, 0.94, 0.949, 0.95, 0.951,  0.959,.96, 0.97 ],
                                            ydata = [ 0.55, 0.55, 0.55, 0.55, 0.52,  0.505, 0.45, 0.45, 0.54, 0.55 ],
                                            c = (0.2, 0.2, 0.2), linewidth = 0.3 ))
      
   if options.blockEdgeDensity:
      data.footerAx.text( x=0.9, y = 0.25, horizontalalignment='right',
                          verticalalignment = 'center',
                          s = 'Gapless block edge density', fontsize = 8 )
      # block edge denisty demo line
      data.footerAx.add_line( lines.Line2D( xdata = [ 0.91, 0.92, 0.93, 0.94, 0.949, 0.95, 0.951,  0.959, 0.96, 0.961, 0.97 ],
                                            ydata = [ 0.25, 0.24, 0.25, 0.24, 0.24,  0.35, 0.23,  0.23, 0.36, 0.25, 0.245 ],
                                            c = '#FA698D', linewidth = 0.3 ))
         

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
   annotOrder = [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'repeat'] # 'tandem'
   for c in data.chrNames:
      if c not in data.annotWigDict:
         sys.stderr.write('Error, unable to locate chromosome %s in annotWigDict!\n' % c )
         sys.exit( 1 )
      if 'xAxis' not in data.annotWigDict[ c ]:
         sys.stderr.write('Error, unable to locate xAxis in annotWigDict[ %s ]!\n' % c )
         sys.exit( 1 )
      for i in range( 0, len( data.annotWigDict[ c ]['xAxis']) ):
         j = 0
         for a in annotOrder:
            
            data.annotWigDict[ c ][ a ][ i ] =  ( data.annotYPos[ j ] + 
                                                  float( data.annotWigDict[ c ][ a ][ i ] ) / 
                                                  data.numRows )
            j += 1
      j = 0
      for a in annotOrder:
         axDict[ c ].add_line( lines.Line2D( xdata=[0, data.chrLengthsByChrom[ c ]],
                                                ydata=[data.annotYPos[ j ], data.annotYPos[ j ]],
                                                color= options.annotColors[ a ],
                                                linewidth=0.3))
         axDict[ c ].fill_between( x=data.annotWigDict[ c ]['xAxis'], 
                                   y1=data.annotWigDict[ c ][ a ], 
                                   y2=data.annotYPos[ j ],
                                   facecolor = options.annotColors[ a ], 
                                   linewidth = 0.0 )
         j += 1

def drawMafs( axDict, options, data ):
   alternatingColors = { True: ( 0.2, 0.2, 0.2 ),
                         False: ( 0.2, 0.2, 0.2 ) }
   myGray = ( 0.8, 0.8, 0.8 )

   for c in data.chrNames:
      col = True
      j = 0
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis']) ):
            for r in [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
                       'maf1e5', 'maf1e6', 'maf1e7', 'blockEdgeCount',
                       'mafHpl1e2', 'mafHpl1e3', 'mafHpl1e4', 
                       'mafHpl1e5', 'mafHpl1e6', 'mafHpl1e7', 'mafHpEdgeCount',
                       'mafHpErrorCount' ]:
               # adjust the height and the position of the track to fit in the plot
               data.mafWigDict[ c ][ n ][ r ][ i ] = ( data.mafYPos[j] + 
                                                       float( data.mafWigDict[ c ][ n ][ r ][ i ] ) *
                                                       ( data.increment * 0.92 ) )
         # draw the baseline
         axDict[ c ].add_line( lines.Line2D( xdata=[0, data.chrLengthsByChrom[ c ]],
                                             ydata=[data.mafYPos[ j ], data.mafYPos[ j ]],
                                             color= myGray,
                                             linewidth=0.3))
         # Basic fills
         if options.fill:
            axDict[ c ].add_line( lines.Line2D( xdata=[0, data.chrLengthsByChrom[ c ]],
                                                ydata=[data.mafYPos[ j ], data.mafYPos[ j ]],
                                                color= myGray,
                                                linewidth=0.3))
            
            axDict[ c ].fill_between( x=data.mafWigDict[ c ][ n ]['xAxis'], 
                                      y1=data.mafWigDict[ c ][ n ]['maf'],
                                      y2=data.mafYPos[ j ],
                                      facecolor = myGray,
                                      linewidth = 0.0 )
            # Stack Fills
         elif options.stackFillBlocks:
            k = -1
            for r in [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
                       'maf1e5', 'maf1e6', 'maf1e7' ]:
               k += 1
               axDict[ c ].fill_between( x=data.mafWigDict[ c ][ n ]['xAxis'], 
                                         y1=data.mafWigDict[ c ][ n ][ r ],
                                         y2=data.mafYPos[ j ],
                                         facecolor = data.stackFillColors[ k ],
                                         linewidth = 0.0 )
         elif options.stackFillHapPaths:
            k = -1
            for r in [ 'maf', 'mafHpl1e2', 'mafHpl1e3', 'mafHpl1e4', 
                       'mafHpl1e5', 'mafHpl1e6', 'mafHpl1e7' ]:
               k += 1
               axDict[ c ].fill_between( x=data.mafWigDict[ c ][ n ]['xAxis'], 
                                         y1=data.mafWigDict[ c ][ n ][ r ],
                                         y2=data.mafYPos[ j ],
                                         facecolor = data.stackFillColors[ k ],
                                         linewidth = 0.0 )
            
         # old red = "#FA698D"
         # new red = "#FA9AAB"
         myRed  = '#FA9AAB'
         myBlue = '#C5C3E2'
         # --blockEdgeDensity track
         if options.blockEdgeDensity:
            axDict[ c ].add_line( lines.Line2D( xdata = data.mafWigDict[ c ][ n ]['xAxis'], 
                                                ydata = data.mafWigDict[ c ][ n ]['blockEdgeCount'], 
                                                c = myRed, linewidth=0.3)) # , linestyle='None',
                                                #marker='o', markerfacecolor=myRed, mec='None',
                                                #markersize=0.3) )
         # --hapPathEdgeDensity track
         if options.hapPathEdgeDensity:
            axDict[ c ].add_line( lines.Line2D( xdata = data.mafWigDict[ c ][ n ]['xAxis'], 
                                                ydata = data.mafWigDict[ c ][ n ]['mafHpEdgeCount'], 
                                                c = myBlue, linewidth=0.3)) #linestyle='None', linewidth=0.0, 
                                                #marker='.', markerfacecolor='b', markersize=1.0) )
         # --hapPathErrorDensity track
         if options.hapPathErrorDensity:
            axDict[ c ].add_line( lines.Line2D( xdata = data.mafWigDict[ c ][ n ]['xAxis'], 
                                                ydata = data.mafWigDict[ c ][ n ]['mafHpErrorCount'], 
                                                c = myRed, linewidth=0.3))# ,linestyle='None',
                                                #marker='o', markerfacecolor=myRed, mec='None', 
                                                #markersize=0.3) )

         if (not options.stackFillBlocks) and (not options.stackFillHapPaths) and (not options.fill):
            # No Fills, basic wiggle
            axDict[ c ].add_line( lines.Line2D( xdata = data.mafWigDict[ c ][ n ]['xAxis'], 
                                                ydata = data.mafWigDict[ c ][ n ]['maf'], 
                                                c = alternatingColors[ col ], linewidth = 0.3 ))

         j +=1
         col = not col

def prettyPrintLength( n ):
    """ takes an integer with a number of bases,
    returns a string with the number displayed nicely
    with units.
    """
    if float( n ) / 1000000000.0 >= 1:
        v = '%.2f' % ( float(n) / 1000000000.0 )
        v = v.strip('0').strip('.')
        units = 'Gb'
    elif float( n ) / 1000000.0 >= 1:
        v = '%.2f' % ( float(n) / 1000000.0 )
        v = v.strip('0').strip('.')
        units = 'Mb'
    elif float( n ) / 1000.0 >= 1:
        v = '%.2f' % ( float(n) / 1000.0 )
        v = v.strip('0').strip('.')
        units = 'Kb'
    elif n == 1:
        v = '%d' % n
        units = 'base'
    else:
        v = '%d' % n
        units = 'bases'
    return '%s %s' % ( v, units )

def normalizeErrorDensities( options, data ):
   if options.relative:
      localErrorMaxes = {}
      localEdgeMaxes  = {}
      localMaxes      = {}
      for n in data.orderedMafs:
         localErrorMaxes[ n ] = 0
         localEdgeMaxes[ n ]  = 0
         localMaxes[ n ]      = 0
         for c in data.chrNames:
            if localErrorMaxes[ n ] < data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ]:
               localErrorMaxes[ n ] = data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ]
            if localEdgeMaxes[ n ] < data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ]:
               localEdgeMaxes[ n ] = data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ]
            if localMaxes[ n ] < data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ]:
               localMaxes[ n ] = data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ]
            if localMaxes[ n ] < data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ]:
               localMaxes[ n ] = data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ]
   else:
      globalErrorMax = 0
      globalEdgeMax = 0
      for c in data.chrNames:
         for n in data.orderedMafs:
            if data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ] > globalErrorMax:
               globalErrorMax = data.mafWigDict[ c ][ n ][ 'mafHpErrorMax' ]
            if data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ] > globalEdgeMax:
               globalEdgeMax = data.mafWigDict[ c ][ n ][ 'mafHpEdgeMax' ]
            
   for c in data.chrNames:
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis']) ):
            if options.zerosToNan and data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] == 0:
               data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] = float( 'nan' )
            if options.edgeErrorCeiling:
               if data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] >= float( options.edgeErrorCeiling ):
                  data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] = 1.0
               else:
                  data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] /= float( options.edgeErrorCeiling )
                  
               if data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] > float( options.edgeErrorCeiling ):
                  data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] = 1.0
               else:
                  data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] /= float( options.edgeErrorCeiling )
            else:
               if options.relative:
                  data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] /= float( localMaxes[ n ] )
                  data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ]  /= float( localMaxes[ n ] )
               else:
                  data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] /= float( globalErrorMax )
                  data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ]  /= float( globalEdgeMax  )
                  

def normalizeBlockEdgeDensities( options, data ):
   if not options.relative:
      ultimateMax = 0
      for c in data.chrNames:
         for n in data.orderedMafs:
            if data.mafWigDict[ c ][ n ][ 'blockEdgeMax' ]> ultimateMax:
               ultimateMax = data.mafWigDict[ c ][ n ][ 'blockEdgeMax' ]
   for c in data.chrNames:
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis'])):
            if not options.relative:
               data.mafWigDict[ c ][ n ]['blockEdgeCount'][ i ] = data.mafWigDict[ c ][ n ]['blockEdgeCount'][ i ] / float( ultimateMax )
            else:
               data.mafWigDict[ c ][ n ]['blockEdgeCount'][ i ] = ( data.mafWigDict[ c ][ n ]['blockEdgeCount'][ i ] / float( data.mafWigDict[ c ][ n ][ 'blockEdgeMax' ] ))

def normalizeAnnotatinos( options, data ):
   pass

def transformErrorDensities( options, data ):
   for c in data.chrNames:
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis'])):
            data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] = data.mafWigDict[ c ][ n ]['mafHpErrorCount'][ i ] ** 0.25
            data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] = data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] ** 0.25

def transformBlockEdgeDensities( options, data ):
   for c in data.chrNames:
      for n in data.orderedMafs:
         for i in range( 0, len( data.mafWigDict[ c ][ n ]['xAxis'])):
            data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] = data.mafWigDict[ c ][ n ]['mafHpEdgeCount'][ i ] ** 0.25

def normalizeData( options, data ):
   normalizeErrorDensities( options, data )
   normalizeBlockEdgeDensities( options, data )
   normalizeAnnotations( options, data )

def transformData( options, data ):
   if options.transform:
      transformErrorDensities( options, data )
      transformBlockEdgeDensities( options, data )

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser, data )
   loadAnnots( options, data )
   loadMafs( options, data )

   normalizeData( options, data )
   transformData( options, data )

   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxes( fig, options, data )
   labelAxes( fig, axDict, options, data )
   drawAnnotations( axDict, options, data )
   drawMafs( axDict, options, data )
   drawLegend( options, data )

   setAxisLimits( axDict, options, data )
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
