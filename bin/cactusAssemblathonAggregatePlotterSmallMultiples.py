#!/usr/bin/env python
"""
cactusAssemblathonAggregatePlotterSmallMultiples.py
9 March 2011
dent earl, dearl@soe.ucsc.edu

This script takes a set of aggregate text files and produces
a pretty picture. Files look like:

[dearl@hgwdev demo]$ cat agg.A1.txt 
block/contig_lengths/haplotype_path_lengths	0_0_0	100_0_0	1000_0_0	10000_0_0	100000_0_0	1000000_0_0	10000000_0_0	100000000_0_0
hapA1/hapA2/assembly	109803876	109803876	109688632	108193470	74325404	0	0	0
hapA1/hapA2/!assembly	813068	813068	928312	2423474	36291540	110616944	110616944	110616944
hapA1/!hapA2/assembly	814571	814571	812441	784434	490059	0	0	0
hapA1/!hapA2/!assembly	761105	761105	763235	791242	1085617	1575676	1575676	1575676
!hapA1/hapA2/assembly	795714	795714	791657	765705	487629	0	0	0
!hapA1/hapA2/!assembly	834204	834204	838261	864213	1142289	1629918	1629918	1629918
!hapA1/!hapA2/assembly	1001654	1001338	923700	823749	508292	437	0	0
"""
from libMafGffPlot import Data
import cactusAssemblathonAggregatePlotter as cacPlot

import glob
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import sys
import re

def initOptions( parser ):
   parser.add_option( '--dir', dest='dir',
                      type='string',
                      help='Directory of Aggregate files to be read.' )
   parser.add_option( '--crazyMax', dest='crazyMax',
                      type='int',
                      help='Sets the height of the crazy bar plot axis.' )
   parser.add_option( '--out', dest='out', default='myAggPlot',
                      type='string',
                      help='output pdf where figure will be created. No extension.' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]' )
   parser.add_option( '--order', dest='order',
                      type='string',
                      help=('Order (left-right, top-bottom) of plots, comma '
                            'separated. Names must match file prefixes in the --dir.' ))
   parser.add_option( '--mode', dest='mode',
                      type='string', default='',
                      help='Plotting mode [contigs|hapPaths|blocks|contamination].' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output.')
   parser.add_option( '--frames', dest='frames', default=False,
                      action='store_true',
                      help='Debug option, turns on the printing of frames around axes.' )

def checkOptions( options, parser ):
   if options.dir == None:
      parser.error( 'Error, specify --dir.\n' )
   if not os.path.exists( options.dir ):
      parser.error( 'Error, --dir %s does not exist.\n' % options.dir )
   if not os.path.isdir( options.dir ):
      parser.error( 'Error, --dir %s is not a directory.\n' % options.dir )
   options.dir = os.path.abspath( options.dir )
   if ( options.out[-4:] == '.png' or options.out[-4:] == '.pdf' or 
        options.out[-4:] == '.eps' ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )
   if ( options.mode != 'contigs' and options.mode != 'contamination' and
        options.mode != 'blocks' and options.mode != 'hapPaths' ):
      parser.error('Error, you must specify one of the modes listed under --mode in --help.\n')
   if options.mode == 'blocks' or options.mode == 'hapPaths' or options.mode == 'contigs':
      options.topBotOrder = [ 'hapA1/hapA2/!assembly', 'hapA1ORhapA2/!assembly',
                              'hapA1ORhapA2/assembly','hapA1/hapA2/assembly' ]
   elif options.mode == 'contamination':
      options.topBotOrder = [ 'ecoli/!assembly', 'ecoli/assembly' ]
   if options.order != None:
      options.order = options.order.split(',')
   else:
      options.order = []
   options.SMM = True

def readFiles( options, data ):
   matches = glob.glob( os.path.join( options.dir, '*.' + options.mode+ '.*'))
   if len( matches ) < 1:
      sys.stderr.write( 'Error, unable to locate any %s files in %s' % ( options.mode, options.dir ))
      sys.exit( 1 )
   data.recordsDict = {} # keyed by assemblyID
   regEx = '([A-Z]\d{1,2})\.' + options.mode + '.*'
   pat = re.compile( regEx )
   for m in matches:
      aName = re.match( pat, os.path.basename( m ) ).group( 1 )
      data.recordsDict[ aName ] = cacPlot.readFile( m, options )

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(7, 8), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   data.pdf = pdf

def establishAxis( options, data ):
   options.axLeft  = 0.05
   options.axWidth = 0.9
   options.axTop = 0.95
   options.axHeight = 0.9
   options.axRight = options.axLeft + options.axWidth
   options.axBottom = options.axTop - options.axHeight
   options.margins = 0.015
   data.ax = data.fig.add_axes( [ options.axLeft, options.axBottom, 
                                  options.axWidth, options.axHeight ] )
   data.ax.yaxis.set_major_locator( pylab.NullLocator() )
   data.ax.xaxis.set_major_locator( pylab.NullLocator() )
   if not options.frames:
      plt.box( on=False )

def drawPlaceHolder( i, left, top, width, height, options, data ):
   data.ax.add_patch( patches.Rectangle( xy=(left, top - height ), width = width,
                                         height = height, color = (0.2, 0.2, 0.2),
                                         alpha=0.5, edgecolor='None'))
   data.ax.text( x=left + width/2.0, y=top - height / 2.0, s = str(i),
                 fontsize = 14, horizontalalignment='center',
                 verticalalignment = 'center', family='Helvetica',
                 color='w' )

def createAxes( i, left, top, width, height, options, data ):
   # transform coordinates
   figLeft   = options.axLeft + left * options.axWidth
   figTop    = options.axBottom + options.axHeight * top
   figWidth  = width * options.axWidth
   figHeight = height * options.axHeight
   figBottom = figTop - figHeight
   axMain  = data.fig.add_axes( [ figLeft, figBottom,
                                  figWidth, figHeight * 0.65 ] )
   axMain.yaxis.set_major_locator( pylab.NullLocator() )
   axMain.xaxis.set_major_locator( pylab.NullLocator() )
   axMain.text( x=0.5, y=0.5, s = str(i),
                fontsize = 14, horizontalalignment='center',
                verticalalignment = 'center', family='Helvetica',
                color=(0.7, 0.7, 0.7) )
   
   if not options.frames:
      plt.box( on=False )
   axCrazy = data.fig.add_axes( [ figLeft, figBottom + figHeight * 0.68,
                                  figWidth, figHeight * 0.04 ] )
   axCrazy.yaxis.set_major_locator( pylab.NullLocator() )
   axCrazy.xaxis.set_major_locator( pylab.NullLocator() )

   if not options.frames:
      plt.box( on=False )
   axBlowUp = data.fig.add_axes( [ figLeft, figBottom + figHeight * 0.75,
                                   figWidth, figHeight * 0.25 ] )
   axBlowUp.yaxis.set_major_locator( pylab.NullLocator() )
   axBlowUp.xaxis.set_major_locator( pylab.NullLocator() )

   if not options.frames:
      plt.box( on=False )
   return ( axMain, axCrazy, axBlowUp )

def drawPlots( options, data ):
   row = -1
   numRows = 7.0
   numCols = 9.0
   plotHeight = ( 1.0 - ( numRows - 1.0 ) * options.margins ) / numRows
   plotWidth  = ( 1.0 - ( numCols - 1.0 ) * options.margins ) / numCols
   for i in range( 0, 63 ):
      if not i % numCols:
         row += 1
      top  = 1.0 - row * float( plotHeight + options.margins )
      left = ( i % numCols ) * float( plotWidth + options.margins )
      ( axMain, axCrazy, axBlowUp ) = createAxes( i, left, top, plotWidth, 
                                                  plotHeight, options, data )
      #drawPlaceHolder( i, left, top, plotWidth, plotHeight, options, data )
   
def writeImage( options, data ):
   if options.outFormat == 'pdf':
      data.fig.savefig( data.pdf, format='pdf' )
      data.pdf.close()
   elif options.outFormat == 'png':
      data.fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'all':
      data.fig.savefig( data.pdf, format='pdf' )
      data.pdf.close()
      data.fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
      data.fig.savefig( options.out + '.eps', format='eps' )
   elif options.outFormat == 'eps':
      data.fig.savefig( options.out + '.eps', format='eps' )   

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   
   readFiles( options, data )
   initImage( options, data )
   establishAxis( options, data )
   
   drawPlots( options, data )
   
   writeImage( options, data )

if __name__ == '__main__':
   main()