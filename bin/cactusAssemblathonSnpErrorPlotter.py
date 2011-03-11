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

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(6, 8), dpi=options.dpi, facecolor='w' )
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
      axDict[ ax ].set_xticks( [] )

def drawData( assembliesDict, axDict, options, data ):
   pass

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
   
   drawData( assembliesDict, axDict, options, data )
   
   writeImage( fig, pdf, options, data )
   

if __name__ == '__main__':
   main()
