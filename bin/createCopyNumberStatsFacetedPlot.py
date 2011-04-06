#!/usr/bin/env python
"""
createCopyNumberStatsFacetedPlot.py
5 April 2011
dent earl dearl (a) soe ucsc edu

used in the assemblathon report project to
create a plot of excess, deficient and total copy number bases
from a single copy number stats xml file.

"""
import createCopyNumberStatsPlot as ccnsp
import glob
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

class CopyNumberStat:
   def __init__( self ):
      self.name     = ''
      self.defUpper = -1.0
      self.defLower = -1.0
      self.excUpper = -1.0
      self.excLower = -1.0
      self.sumUpper = -1.0
      self.sumLower = -1.0

def initOptions( parser ):
   parser.add_option( '--dir', dest='dir',
                      type='string',
                      help=('Location of all upper (_0.xml) and lower (_1000.xml) files.'))
   parser.add_option( '--title', dest='title',
                      type='string', default='Copy Number Statistics',
                      help='Title placed at the top of the plot. default=%default' )
   parser.add_option( '--out', dest='out', default='myCopyNumberStatsPlot',
                      type='string',
                      help='filename where figure will be created. No extension. default=%default' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps] default=%default' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output. default=%default' )

def checkOptions( args, options, parser ):
   if len( args ) > 0:
      parser.error('Error, unanticipated arguments: %s.\n' % args )
   if options.dir == None:
      parser.error( 'Error, specify --dir\n' )
   if not os.path.exists( options.dir ):
      parser.error('Error, --dir %s does not exist!\n' % ( options.dir ))
   if not os.path.isdir( options.dir ):
      parser.error('Error, --dir %s is not a directory!\n' % ( options.dir ))
   options.dir = os.path.abspath( options.dir )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=( 8, 10 ), dpi=options.dpi, facecolor='w' )
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
   options.axLeft   = 0.09
   options.axRight  = 0.98
   options.axWidth    = options.axRight - options.axLeft 
   options.axBottom = 0.05
   options.axTop    = 0.98
   options.axHeight = options.axTop - options.axBottom
   margin = 0.05
   facetHeight = ( options.axHeight - 2.0 * margin) / 3.0
   yPos = 0.0
   for ax in [ 'def', 'exc', 'sum' ]:
      axDict[ ax ] = fig.add_axes( [ options.axLeft, options.axBottom + yPos, 
                                     options.axWidth, facetHeight ] )
      axDict[ ax ].yaxis.set_major_locator( pylab.NullLocator() )
      axDict[ ax ].xaxis.set_major_locator( pylab.NullLocator() )
      yPos += facetHeight + margin
      #plt.box( on=False )
   for ax in axDict:
      for loc, spine in axDict[ ax ].spines.iteritems():
         if loc in ['left','bottom']:
            spine.set_position(('outward',10)) # outward by 10 points
         elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine               
         else:
            raise ValueError('unknown spine location: %s' % loc )
      # turn off ticks where there is no spine
      axDict[ ax ].xaxis.set_ticks_position('bottom')
      axDict[ ax ].yaxis.set_ticks_position('left')
   data.axDict = axDict
   return ( axDict )

def drawLegend( options, data ):
   pass

def drawAxisLabels( axDict, cDict, options, data ):
   pass

def setAxisLimits( axDict, options, data ):
   pass

def drawData( axDict, sList, options, data ):
   xNames = []
   yMax = 0
   yMin = sys.maxint
   lGray = ( 0.9, 0.9, 0.9 )
   
   axDict[ 'sum' ].set_yscale('log')
   for s in sList:
      if yMax < float( s.sumUpper ):
         yMax = float( s.sumUpper )
      if yMin > float( s.sumLower ): 
         yMin = float( s.sumLower )
   for i in xrange( 1, len( sList ) ):
      if not i % 10:
         axDict[ 'sum' ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ yMin, yMax ],
                                                 color=lGray))
   i=0
   for s in sList:
      i += 1
      axDict[ 'sum' ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ s.sumLower, s.sumUpper ],
                                              color='#1f77b4'))
      xNames.append( s.name )
   axDict[ 'sum' ].set_xticks( range( 1, len(xNames) + 1 ))
   axDict[ 'sum' ].set_xticklabels( xNames )
   for tick in axDict[ 'sum' ].xaxis.get_major_ticks():
      tick.label1.set_fontsize( 6 )
   for label in axDict[ 'sum' ].xaxis.get_ticklabels():
      label.set_rotation( 90 )
   axDict[ 'sum' ].set_ylim( [ yMin * 0.9, yMax * 1.1] )
   #plt.ylabel( 'log proportion ' )
   
   axDict[ 'exc' ].set_yscale('log')
   for s in sList:
      if yMax < float( s.excUpper ):
         yMax = float( s.excUpper )
      if yMin > float( s.excLower ): 
         yMin = float( s.excLower )
   for i in xrange( 1, len( sList ) ):
      if not i % 10:
         axDict[ 'exc' ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ yMin, yMax ],
                                                 color=lGray))
   i=0
   for s in sList:
      i += 1
      axDict[ 'exc' ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ s.excLower, s.excUpper ],
                                              color='#1f77b4'))
   axDict[ 'exc' ].set_xticks( range( 1, len(xNames) + 1 ))
   axDict[ 'exc' ].set_xticklabels( xNames )
   for tick in axDict[ 'exc' ].xaxis.get_major_ticks():
      tick.label1.set_fontsize( 6 )
   for label in axDict[ 'exc' ].xaxis.get_ticklabels():
      label.set_rotation( 90 )
   axDict[ 'exc' ].set_ylim( [ yMin * 0.9, yMax * 1.1] )
   #plt.ylabel( 'log proportion ' )
   
   axDict[ 'def' ].set_yscale('log')
   for s in sList:
      if yMax < float( s.defUpper ):
         yMax = float( s.defUpper )
      if yMin > float( s.defLower ): 
         yMin = float( s.defLower )
   for i in xrange( 1, len( sList ) ):
      if not i % 10:
         axDict[ 'def' ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ yMin, yMax ],
                                                 color=lGray))
   i=0
   for s in sList:
      i += 1
      axDict[ 'def' ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ s.defLower, s.defUpper ],
                                              color='#1f77b4'))
   axDict[ 'def' ].set_xticks( range( 1, len(xNames) + 1 ))
   axDict[ 'def' ].set_xticklabels( xNames )
   for tick in axDict[ 'def' ].xaxis.get_major_ticks():
      tick.label1.set_fontsize( 6 )
   for label in axDict[ 'def' ].xaxis.get_ticklabels():
      label.set_rotation( 90 )
   axDict[ 'def' ].set_ylim( [ yMin * 0.9, yMax * 1.1] )
   #plt.ylabel( 'log proportion ' )
   
   axDict[ 'sum' ].text( x=0.01, y=0.98, s = 'Sum of Proportional Copy Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['sum'].transAxes )
   axDict[ 'exc' ].text( x=0.01, y=0.98, s = 'Proportional Excess Copy Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['exc'].transAxes )
   axDict[ 'def' ].text( x=0.01, y=0.98, s = 'Proportional Deficient Copy Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['def'].transAxes )

def readFiles( options ):
   ups = glob.glob( os.path.join( options.dir, '*_0.xml'))
   los = glob.glob( os.path.join( options.dir, '*_1000.xml'))
   stats = {}
   for u in ups:
      c = CopyNumberStat()
      c.name = os.path.basename( u ).split('.')[0]
      xmlTree = ET.parse( u )
      root=xmlTree.getroot()
      elm = root.find( 'excessCopyNumberCounts' )
      c.excUpper  = float( elm.attrib['totalProportionOfColumns'] )
      elm = root.find( 'deficientCopyNumberCounts' )
      c.defUpper = float( elm.attrib['totalProportionOfColumns'] )
      c.sumUpper = c.excUpper + c.defUpper
      stats[ c.name ] = c
   for l in los:
      name = os.path.basename( l ).split('.')[0]
      if name not in stats:
         continue
      xmlTree = ET.parse( l )
      root=xmlTree.getroot()
      elm = root.find( 'excessCopyNumberCounts' )
      stats[name].excLower  = float( elm.attrib['totalProportionOfColumns'] )
      elm = root.find( 'deficientCopyNumberCounts' )
      stats[name].defLower = float( elm.attrib['totalProportionOfColumns'] )
      stats[name].sumLower = stats[name].excLower + stats[name].defLower
   validStats = {}
   for s in stats:
      if stats[s].excLower == -1.0 or stats[s].defLower == -1 or stats[s].sumLower == -1.0:
         continue
      validStats[s] = stats[s]
      
   return validStats

def main():
   usage = ( 'usage: %prog [options] --dir=path/to/dir/\n\n'
             '%prog takes in a copy number statistics file\n'
             'and creates an image file.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   ( fig, pdf ) = initImage( options, data )

   stats = readFiles( options )
   sortedOrder = sorted( stats.values(), key=lambda x: x.sumLower, reverse=False )

   #establishGlobalMinMax( storedCategories, options, data )

   axDict = establishAxes( fig, options, data )
   
   drawData( axDict, sortedOrder, options, data )
   drawLegend( options, data )
   drawAxisLabels( axDict, stats, options, data )
   setAxisLimits( axDict, options, data )
   
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
