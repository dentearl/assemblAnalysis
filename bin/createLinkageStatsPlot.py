#!/usr/bin/env python
"""
createLinkageStatsPlot.py
31 March 2011
dent earl dearl (a) soe ucsc edu

used in the assemblathon report project to
create a plot from a single linkage
stats xml file.

"""
from libMafGffPlot import Data
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import sys
import xml.etree.ElementTree as ET

class Bucket:
   def __init__( self ):
      self.start   = -1
      self.end     = -1
      self.mid     = -1.0
      self.correct = -1
      self.samples = -1
      self.cumCorrect = -1
      self.cumSamples = -1

def initOptions( parser ):
   parser.add_option( '--title', dest='title',
                      type='string', default='Cumulative Linkage Statistics',
                      help='Title placed at the top of the plot. default=%default' )
   parser.add_option( '--legendElements', dest='legendElements',
                      type='string', help=('Specify the legend text. Comma separated list.'))
   parser.add_option( '--out', dest='out', default='myLinkageStatsPlot',
                      type='string',
                      help='filename where figure will be created. No extension. default=%default' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps] default=%default' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output. default=%default' )

def checkOptions( args, options, parser ):
   if len( args ) < 1:
      parser.error('Error, please specify at least one linkage file to inspect as a positional argument.\n' )
   options.files = []
   for f in args:
      if not os.path.exists( f ):
         parser.error('Error, %s does not exist!\n' % ( f ))
      if not f.endswith('.xml'):
         parser.error('Error, file "%s" does not end in ".xml".\n' % f )
      options.files.append( os.path.abspath( f ) )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )
   if options.legendElements != None:
      options.legendElements = options.legendElements.split(',')

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 10), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   return ( fig, pdf )

def writeImage( fig, pdf, options, data ):
   if options.outFormat == 'pdf':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
   elif options.outFormat == 'png':
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'all':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
      fig.savefig( options.out + '.eps', format='eps' )
   elif options.outFormat == 'eps':
      fig.savefig( options.out + '.eps', format='eps' )

def establishAxes( fig, options, data ):
   axDict = {}
   options.axLeft = 0.12
   options.axWidth = 0.83
   axDict[ 'main' ] = fig.add_axes( [ options.axLeft, 0.1,
                                      options.axWidth , 0.8 ] )
   #plt.box( on=False )
   data.axDict = axDict
   return ( axDict )

def readFiles( options ):
   # these lists have one item per input file
   statsList = []
   xData = []
   for f in options.files:
      xmlTree = ET.parse( f )
      root=xmlTree.getroot()
      fileBucketList = []
      fileXData = []
      for elm in root.findall( 'bucket' ):
         b = Bucket()
         b.start = int( elm.attrib[ 'from' ] )
         b.end   = int( elm.attrib[ 'to' ] )
         b.mid   = (b.end - b.start) / 2.0 + b.start
         b.correct = int( elm.attrib[ 'correct' ] )
         b.samples = int( elm.attrib[ 'samples' ] )
         b.cumCorrect = int( elm.attrib[ 'cumulative_correct' ] )
         b.cumSamples = int( elm.attrib[ 'cumulative_samples' ] )
         fileBucketList.append( b )
         fileXData.append( b.mid )
      statsList.append( fileBucketList )
      xData.append( fileXData )
   return ( statsList, xData )

def setAxisLimits( ax, xData, options, data ):
   ax.set_xscale('log')
   ax.set_ylim( 0.2, 1.001 )
   #ax.set_xlim( 1, xData[ -1 ] )

def establishTicks( ax, xData, options, data ):
   # turn off ticks
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   minorLocator = LogLocator( base=10, subs = range(1,10) )
   ax.xaxis.set_minor_locator( minorLocator )
   
def drawLegend( options, data ):
   if len( options.files ) < 2:
      return
   if options.legendElements == None:
      pltListLabels = []
      for f in options.files:
         pltListLabels.append( os.path.basename( f ) )
   elif len( options.legendElements ) == len( options.files ):
      pltListLabels = options.legendElements
   else:
      sys.stderr.write('Error, length of items in --legendElements not equal to number of linkage xml files.\n')
      sys.exit( 1 )
   leg = plt.legend( data.pltList, pltListLabels, 'lower left' )
   leg._drawFrame=False

def drawAxisLabels( fig, options, data ):
   data.axDict['main'].set_title( options.title )
   plt.xlabel('Distance between points' )
   plt.ylabel('Proportion')

def drawData( ax, xData, sList, options, data ):
   colors = [ "#1f77b4", "#aec7e8", # blues 
              "#ff7f0e", "#ffbb78", # oranges
              "#2ca02c", "#98df8a", # greens
              "#d62728", "#ff9896", # reds
              "#9467bd", "#c5b0d5" ] # lavenders
   styles = { 0:'-', 1:'--' }
   data.pltList = [] # used for legends
   # grey 0.50 horizontal line:
   ax.add_line( lines.Line2D( xdata=[1, xData[0][-1]],
                              ydata=[0.5, 0.5],
                              linewidth=0.25,
                              color=(0.8, 0.8, 0.8)
                              ))
   for i in range( 0, len( sList )):
      yData = []
      for b in sList[ i ]:
         yData.append( float( b.correct ) / b.samples )
      p = ax.plot( xData[ i ], yData, 
                   color=colors[ i % len( colors ) ], 
                   linestyle=styles[ i >= len( colors )] )
      data.pltList.append( p )
   for loc, spine in ax.spines.iteritems():
      if loc in ['left','bottom']:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in ['right','top']:
         spine.set_color('none') # don't draw spine               
      else:
         raise ValueError('unknown spine location: %s' % loc )

def main():
   usage = ( 'usage: %prog [options] file1.xml file2.xml\n\n'
             '%prog takes in linkage path statistics file(s)\n'
             'and creates an image file.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxes( fig, options, data )
   
   ( data.statsList, data.xData ) = readFiles( options )
   for i in range(0, len( data.statsList )):
      data.statsList[i] = sorted( data.statsList[i], key=lambda x: x.mid, reverse=False )
   
   drawData( axDict['main'], data.xData, data.statsList, options, data )
   drawLegend( options, data )
   drawAxisLabels( fig, options, data )
   setAxisLimits( axDict['main'], data.xData, options, data )
   establishTicks( axDict['main'], data.xData, options, data )
   
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
