#!/usr/bin/env python
"""
cacutsAssemblerSnpErrorPlotter.py
10 March 2011
dent earl, dearl(a) soe ucsc edu



"""
from libMafGffPlot import Data

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

class Assembly:
   """ Assembly objects are generated from lines 
   in the two snp summary files, lower and upper
   """
   def __init__( self ):
      self.ID    = ''
      self.snpStatsLower = {}
      self.snpStatsUpper = {}
      self.allUp = -1
      self.allLo = -1

def initOptions( parser ):
   parser.add_option( '--snpStatsDir', dest='snpStatsDir',
                      type='string',
                      help=('Directory with snpStats. Names: A1.snpStats.upper.txt .'))
   parser.add_option( '--out', dest='out', default='mySnpStatsPlot',
                      type='string',
                      help='output pdf where figure will be created. No extension.' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output.')

def checkOptions( options, parser ):
   dirs = { 'snpStatsDir' : options.snpStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('Error, specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('Error, --%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('Error, --%s %s is not a directory!\n' % (d, dirs[ d ]) )
   if ( options.out[-4:] == '.png' or options.out[-4:] == '.pdf' or 
        options.out[-4:] == '.eps' ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )

def readSnpStatsDir( assembliesDict, options ):
   lowerStatsFiles = glob.glob( os.path.join( options.snpStatsDir, '*.snpStats.lower.txt') )
   upperStatsFiles = glob.glob( os.path.join( options.snpStatsDir, '*.snpStats.upper.txt') )
   
   namereg = '^([A-Z0-9]{2,3})\.snpStats.*'
   namepat = re.compile( namereg  )
   for l in lowerStatsFiles:
      m = re.match( namepat, os.path.basename( l ))
      if not m:
         sys.stderr.write('Error, unable to match regex "%s" against filename "%s"' % ( namereg, l ))
         sys.exit( 1 )
      ID = m.group( 1 )
      assembliesDict[ ID ] = Assembly()
      assembliesDict[ ID ].ID = ID
      f = open( l, 'r' )
      for line in f:
         line = line.strip()
         d = line.split('\t')
         assembliesDict[ ID ].snpStatsLower[ d[0] ] = d[ 1 ]
      f.close()
   for u in upperStatsFiles:
      m = re.match( namepat, os.path.basename( u ))
      if not m:
         sys.stderr.write('Error, unable to match regex "%s" against filename "%s"' % ( namepat, u ))
         sys.exit( 1 )
      ID = m.group( 1 )
      if ID not in assembliesDict:
         sys.stderr.write('Error, unable to locate key %s in assembliesDict.\n')
         sys.exit( 1 )
      f = open( u, 'r' )
      for line in f:
         line = line.strip()
         d = line.split('\t')
         assembliesDict[ ID ].snpStatsUpper[ d[0] ] = d[ 1 ]
      f.close()
   return assembliesDict

def sumErrors( assembliesDict ):
   for a in assembliesDict:
      assembliesDict[ a ].allUp = 0
      assembliesDict[ a ].allLo = 0
      for e in [ 'Total-errors-in-homozygous', 'Total-errors-in-heterozygous', 
                 'Total-errors-in-one-haplotype-only' ]:
         assembliesDict[ a ].allLo += int( assembliesDict[ a ].snpStatsLower[ e ] )
         assembliesDict[ a ].allUp += int( assembliesDict[ a ].snpStatsUpper[ e ] )

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(9, 11), dpi=options.dpi, facecolor='w' )
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
   """ create one axes per chromosome
   """
   axDict = {}
   options.axLeft   = 0.11
   options.axRight  = 0.95
   options.axWidth  = options.axRight - options.axLeft
   options.margin   = 0.05
   options.axTop    = 0.95
   options.axBot    = 0.05
   options.axHeight = options.axTop - options.axBot
   options.indHeight = float( options.axHeight - 3.0 * options.margin ) / 4.0
   axDict[ 'all' ] = fig.add_axes( [ options.axLeft, options.axTop - options.indHeight,
                                     options.axWidth , options.indHeight ] )
   axDict[ 'hom' ] = fig.add_axes( [ options.axLeft, options.axTop - options.indHeight * 2.0 - options.margin,
                                     options.axWidth , options.indHeight ] )
   axDict[ 'het' ] = fig.add_axes( [ options.axLeft, options.axTop - options.indHeight * 
                                     3.0 - options.margin * 2.0,
                                     options.axWidth , options.indHeight ] )
   axDict[ 'indel' ] = fig.add_axes( [ options.axLeft, options.axTop - options.indHeight * 
                                       4.0 - options.margin * 3.0,
                                       options.axWidth , options.indHeight ] )
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
      # if ax != 'all':
         # axDict[ ax ].set_xticks( [] )
   return axDict

def drawData( assembliesDict, sortOrder, axDict, options, data ):
   lGray = ( 0.9, 0.9, 0.9 )
   for ax in axDict:
      axDict[ ax ].set_xlim( 0, len(assembliesDict) + 3 )

   # all plot
   yMax = 0
   yMin = sys.maxint
   xNames = []
   i = 0
   for aName in sortOrder:
      i += 1
      a = assembliesDict[ aName ]
      if yMax < int( a.allUp ):
         yMax = int( a.allUp )
      if yMin > int( a.allLo ): 
         yMin = int( a.allLo )
      
   yMin = logLower( yMin )
   for i in range( 1, len( assembliesDict ) ):
      if not i % 10:
         axDict[ 'all' ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ 1, yMax ],
                                                 color=lGray))
   i = 0
   for aName in sortOrder:
      a = assembliesDict[ aName ]
      i += 1
      axDict[ 'all' ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ a.allLo, a.allUp ],
                                              color='#1f77b4'))
      xNames.append( aName )
   axDict[ 'all' ].set_yscale('log')
   axDict[ 'all' ].set_ylim( [ yMin, yMax] )
   axDict[ 'all' ].set_xticks( range( 1, len(xNames) + 1 ))
   axDict[ 'all' ].set_xticklabels( xNames )
   for tick in axDict[ 'all' ].xaxis.get_major_ticks():
      tick.label1.set_fontsize( 6 )
   for label in axDict[ 'all' ].xaxis.get_ticklabels():
      label.set_rotation( 90 )

   

   # other plots
   plotAxes = { 'Total-errors-in-homozygous' : 'hom',
                'Total-errors-in-heterozygous': 'het',
                'Total-errors-in-one-haplotype-only':'indel' }
   for key in plotAxes:
      yMax = 0
      yMin = sys.maxint
      xNames = []
      i = 0
      for aName in sortOrder:
         i += 1
         a = assembliesDict[ aName ]
         if yMax < int( a.snpStatsUpper[ key ] ):
            yMax = int( a.snpStatsUpper[ key ] )
         if yMin > int( a.snpStatsLower[ key ] ):
            yMin = int( a.snpStatsLower[ key ] )
      yMin = logLower( yMin )
      for i in range( 1, len( assembliesDict ) ):
         if not i % 10:
            axDict[ plotAxes[ key ] ].add_line( lines.Line2D( xdata=[ i, i ],
                                                    ydata=[ 1, yMax ],
                                                    color=lGray))
      i = 0
      for aName in sortOrder:
         i += 1
         a = assembliesDict[ aName ]
         axDict[ plotAxes[ key ] ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ a.snpStatsLower[ key ],
                                                         a.snpStatsUpper[ key ]],
                                                 color='#1f77b4'))
         xNames.append( aName )
      axDict[ plotAxes[ key ] ].set_yscale('log')
      axDict[ plotAxes[ key ] ].set_ylim( [yMin, yMax] )
      axDict[ plotAxes[ key ] ].set_xticks( range( 1, len(xNames) + 1 ))
      axDict[ plotAxes[ key ] ].set_xticklabels( xNames )
      for tick in axDict[ plotAxes[ key ] ].xaxis.get_major_ticks():
         tick.label1.set_fontsize( 6 )
      for label in axDict[ plotAxes[ key ] ].xaxis.get_ticklabels():
         label.set_rotation( 90 )
   
   axDict[ 'all' ].text( x=0.01, y=0.98, s = 'All Snp Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['all'].transAxes )
   axDict[ 'hom' ].text( x=0.01, y=0.98, s = 'Homozygous Snp Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['hom'].transAxes )
   axDict[ 'het' ].text( x=0.01, y=0.98, s = 'Heterozygous Snp Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict['het'].transAxes )
   axDict[ 'indel' ].text( x=0.01, y=0.98, s = 'Indel Snp Errors',
                  fontsize = 12, horizontalalignment='left',
                  verticalalignment = 'top', family='Helvetica',
                  color=( 0.3, 0.3, 0.3 ),
                  transform=axDict[ 'indel' ].transAxes )
   
   
def logLower( y ):
   """ find the approprate lower bound for y
   if y is going to be displayed on a log plot
   """
   for i in range( 1, 8 ):
      if y == ( y % float( 10.0 ** i)):
         return ( 10.0 ** ( i - 1 ) )

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxes( fig, options, data )
   
   assembliesDict = {}
   assembliesDict = readSnpStatsDir( assembliesDict, options )
   sumErrors( assembliesDict )
   sortOrder = sorted( assembliesDict, key=lambda key: assembliesDict[ key ].allLo, reverse=False )

   drawData( assembliesDict, sortOrder, axDict, options, data )
   
   
   writeImage( fig, pdf, options, data )
   

if __name__ == '__main__':
   main()
