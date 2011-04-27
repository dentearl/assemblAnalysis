#!/usr/bin/env python
"""
createCopyNumberStatsFacetedPlot.py
5 April 2011
dent earl dearl (a) soe ucsc edu

used in the assemblathon report project to
create a plot of excess, deficient and total copy bases
from a single copy stats xml file.

"""
import glob
import libAssemblySubset as las
from libMafGffPlot import Data
import libPlotting as lpt
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

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
                      type='string', default='Copy Statistics',
                      help='Title placed at the top of the plot. default=%default' )
   parser.add_option( '--out', dest='out', default='myCopyNumberStatsPlot',
                      type='string',
                      help='filename where figure will be created. No extension. default=%default' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps] default=%default' )
   parser.add_option( '--log', dest='log', default=False, action='store_true',
                      help='Turns on log scale y axes. default=%default')
   parser.add_option( '--outputRanks', dest='outputRanks', default=False, action='store_true',
                      help='Prints out rankings in tab delimited format. default=%default')
   parser.add_option( '--markers', dest='markers', default=False, action='store_true',
                      help='Turns on filled markers for lower values, open markers for uppers. default=%default')
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output. default=%default' )

def checkOptions( args, options, parser ):
   if len( args ) > 0:
      parser.error('unanticipated arguments: %s.\n' % args )
   if options.dir == None:
      parser.error( 'specify --dir\n' )
   if not os.path.exists( options.dir ):
      parser.error('--dir %s does not exist!\n' % ( options.dir ))
   if not os.path.isdir( options.dir ):
      parser.error('--dir %s is not a directory!\n' % ( options.dir ))
   options.dir = os.path.abspath( options.dir )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )

def establishAxes( fig, options, data ):
   axDict = {}
   options.axLeft   = 0.09
   options.axRight  = 0.97
   options.axWidth    = options.axRight - options.axLeft 
   options.axBottom = 0.06
   options.axTop    = 0.95
   options.axHeight = options.axTop - options.axBottom
   margin = 0.07
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
         if loc in ['left', 'bottom']:
            spine.set_position(('outward',10)) # outward by 10 points
         elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine               
         else:
            raise ValueError('unknown spine location: %s' % loc )
      # turn off ticks where there is no spine
      axDict[ ax ].xaxis.set_ticks_position('bottom')
      axDict[ ax ].yaxis.set_ticks_position('both')
   data.axDict = axDict
   return ( axDict )

def drawLegend( options, data ):
   pass

def drawAxisLabels( axDict, cDict, options, data ):
   pass

def setAxisLimits( axDict, options, data ):
   pass

def value( s, v ):
   """ returns the correct attribute of s by v=string
   """
   if v == 'sumUpper':
      return s.sumUpper
   elif v == 'sumLower':
      return s.sumLower
   elif v == 'excUpper':
      return s.excUpper
   elif v == 'excLower':
      return s.excLower
   elif v == 'defUpper':
      return s.defUpper
   elif v == 'defLower':
      return s.defLower
   else:
      sys.stderr.write('Unrecognized value: %s\n' % v)
      sys.exit( 1 )

def drawData( axDict, sList, options, data ):
   lGray = ( 0.8, 0.8, 0.8 )
   
   for ax in [ 'sum', 'exc', 'def' ]:
      xNames = []
      yMax = 0
      yMin = sys.maxint
      if options.log:
         axDict[ ax ].set_yscale('log')
      for s in sList:
         if yMax < float( value(s, '%sUpper'%ax) ):
            yMax = float( value(s, '%sUpper'%ax) )
         if yMin > float( value(s, '%sLower'%ax) ): 
            yMin = float( value(s, '%sLower'%ax) )
      for i in xrange( 1, len( sList ) + 1 ):
         if not i % 5:
            axDict[ ax ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ yMin, yMax * 1.1 ],
                                                 color=lGray,
                                                 linestyle='dotted'))
      i=0
      for s in sList:
         i += 1
         if options.markers:
            axDict[ ax ].add_line( lines.Line2D( xdata=[ i ],
                                                 ydata=[ value(s, '%sLower' % ax) ],
                                                 marker='o',
                                                 markerfacecolor='#1f77b4',
                                                 markeredgecolor='#1f77b4',
                                                 markersize=4.0))
            axDict[ ax ].add_line( lines.Line2D( xdata=[ i ],
                                                 ydata=[ value(s, '%sUpper' % ax) ],
                                                 markeredgecolor='#1f77b4',
                                                 marker='o',
                                                 markerfacecolor='none',
                                                 markersize=4.0))
         axDict[ ax ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ value(s, '%sLower' % ax), value(s, '%sUpper' % ax) ],
                                              color='#1f77b4', linewidth=4.0, solid_capstyle='round'))
         xNames.append( s.name )
      axDict[ ax ].set_xlim( 0, len( xNames ) + 1 )
      if ax != 'exc':
         axDict[ ax ].set_xticks( range( 1, len(xNames) + 1 ))
         axDict[ ax ].set_xticklabels( xNames )
         for tick in axDict[ ax ].xaxis.get_major_ticks():
            if options.subsetFile:
               tick.label1.set_fontsize( 12 )
            else:
               tick.label1.set_fontsize( 6 )
         #for label in axDict[ ax ].xaxis.get_ticklabels():
         #   label.set_rotation( 90 )
      else:
         axDict[ ax ].set_xticks( range( 1, len(xNames) + 1 ))
         axDict[ ax ].set_xticklabels( [] )
      axDict[ ax ].set_ylim( [ yMin * 0.9, yMax * 1.1] )
      # grid
      mts = axDict[ax].yaxis.get_majorticklocs()
      for m in mts:
         axDict[ax].add_line( lines.Line2D( xdata=[ 1, len( sList ) ],
                                            ydata=[ m, m ],
                                            linewidth=1,
                                            color=lGray,
                                            linestyle='dotted'))
   
   axDict['sum'].set_title('Sum of Proportional Copy Errors')
   axDict['exc'].set_title('Proportional Excess Copy Errors')
   axDict['def'].set_title('Proportional Deficient Copy Errors')

def readFiles( options ):
   ups = glob.glob( os.path.join( options.dir, '*_0.xml'))
   los = glob.glob( os.path.join( options.dir, '*_1000.xml'))
   stats = {}
   for u in ups:
      c = CopyNumberStat()
      c.name = os.path.basename( u ).split('.')[0]
      if options.subsetFile:
         if c.name not in options.assemblySubset:
            continue
      try:
         xmlTree = ET.parse( u )
      except expat.ExpatError: # empty xml file
         continue
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
      try:
         xmlTree = ET.parse( l )
      except expat.ExpatError: # empty xml file
         continue
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

def rankings( sortedOrder, options, data ):
   print ('#Assembly\tSum Errors Lower\tSum Errors Upper\t'
          'Excess Errors Lower\tExcess Errors Upper\t'
          'Deficient Errors Lower\tDeficient Errors Upper')
   for s in sortedOrder:
      sys.stdout.write('%s' % s.name )
      for v in [ s.sumLower, s.sumUpper, s.excLower, s.excUpper,
                 s.defLower, s.defUpper, ]:
         sys.stdout.write('\t%s' % v )
      sys.stdout.write('\n')

def main():
   usage = ( 'usage: %prog [options] --dir=path/to/dir/\n\n'
             '%prog takes in a copy statistics file\n'
             'and creates an image file.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )
   las.checkOptions( options, parser )
   if not options.outputRanks:
      fig, pdf = lpt.initImage( 8.0, 10.0, options, data )

   stats = readFiles( options )
   sortedOrder = sorted( stats.values(), key=lambda x: x.sumLower, reverse=False )
   
   if options.outputRanks:
      rankings( sortedOrder, options, data )
      return

   axDict = establishAxes( fig, options, data )
   drawData( axDict, sortedOrder, options, data )
   drawLegend( options, data )
   drawAxisLabels( axDict, stats, options, data )
   setAxisLimits( axDict, options, data )
   
   lpt.writeImage( fig, pdf, options )

if __name__ == '__main__':
   main()
