#!/usr/bin/env python
"""
createN50StatsPlot.py
14 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the N50 stats table.

output is latex

"""
from createHapPathStatsTable import readDir
import createIndividualSection as cis
import createN50StatsTable as cnfst
import createSortedCoveragesPlot as cscp
import glob
import libAssemblySubset as las
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import re
import sys
import xml.etree.ElementTree as ET

class Data:
   """Dummy class to hold data to 
   pass between functions
   """
   pass

def initOptions( parser ):
   parser.add_option( '--hapPathStatsDir', dest='hapPathStatsDir',
                      type='string',
                      help=('Directory with hapPathStats. Names: A1.hapPathStats.xml .'))
   parser.add_option('--labels', dest='labels',
                     type='string', default=( 'Block N50,Contig Path N50,Scaffold Path N50,Contig NA50'),
                     help=('Lables of the facets, comma separated, from top '
                           'to bottom: default=%default'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='contigNA50',
                      help=('Allows a different sort order. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks', action='store_true',
                      default=False, help=('Outputs rankings as tab delimited '
                                           'stream to STDOUT. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'hapPathStatsDir'   : options.hapPathStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('Error, specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('Error, --%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('Error, --%s %s is not a directory!\n' % (d, dirs[ d ]) )
   if options.outputRanks:
      return
   allowedKeys = set(['contigNA50', 'scaffoldPathN50',
                      'haplotypePathN50', 'blockN50'])
   if options.sortOn not in allowedKeys:
      parser.error('Error, --sortOn %s is not in the dict of allowed keys: %s' % 
                   ( options.sortOn, allowedKeys ))
   labelList = options.labels.split(',')
   if len( labelList ) != 4:
      parser.error('Error, --label needs to contain 4 items, comma separated.')
   k = [ 'blockN50', 'haplotypePathN50', 'scaffoldPathN50',  'contigNA50']
   options.labelDict = dict( zip( k, labelList ))

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 8), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   return ( fig, pdf )

def establishAxis( fig, options, data ):
   """ 
   """
   axDict = {}
   columns = [ 'blockN50', 'haplotypePathN50', 'scaffoldPathN50','contigNA50' ] 
   # 'totalContigNumber', 'contigN50'
   options.axLeft    = 0.12
   options.axWidth   = 0.85
   options.axBottom  = 0.05
   options.axHeight  = 0.9
   options.axTop     = options.axBottom + options.axHeight
   options.margin    = 0.08
   indvHeight = ( options.axHeight - (len( columns) - 1.0) * options.margin ) / float( len( columns ))
   prevY = options.axTop
   i = -1
   for n in columns:
      i += 1
      axDict[ n ] = fig.add_axes( [options.axLeft, prevY - indvHeight,
                                   options.axWidth, indvHeight ] )
      prevY = ( prevY - indvHeight ) - options.margin
   return axDict

def getVals( assembliesList, key ):
   v = []
   for a in assembliesList:
      v.append( a.valuesDict[ key ] )
   return v

def getIDs( assembliesList ):
   v = []
   for a in assembliesList:
      v.append( a.ID )
   return v

def drawData( assembliesList, maxesDict, axDict, options ):
   
   for ax in axDict:
      for i in range(1, len( assembliesList )):
         if not i % 10:
            axDict[ax].add_line( lines.Line2D( xdata=[ i-1, i-1 ],
                                               ydata=[ 1, maxesDict[ ax ] ],
                                               linewidth=0.5,
                                               color=(0.8, 0.8, 0.8) ))
      p1 = axDict[ax].plot( range( 0,len( assembliesList ) ), getVals( assembliesList, ax ), 
                            '.', color='#1f77b4', markersize=10.0)
      for loc, spine in axDict[ax].spines.iteritems():
         if loc in [ 'left'  ]:
            spine.set_position(('outward',10)) # outward by 10 points
         elif loc in [ 'top',  'right' ]:
            spine.set_color('none') # don't draw spine               
         elif loc in [ 'bottom' ]:
            pass
         else:
            raise ValueError('unknown spine location: %s' % loc )
      if ax != 'blockN50':
         if ax == 'contigNA50':
            axDict[ax].set_xticks( range( 0, len( assembliesList ) ))
            axDict[ax].set_xticklabels( getIDs( assembliesList )  )
            for tick in axDict[ax].xaxis.get_major_ticks():
               if options.subsetFile:
                  tick.label1.set_fontsize( 12 )
               else:
                  tick.label1.set_fontsize( 6 )
            #for label in axDict[ax].xaxis.get_ticklabels():
            #   label.set_rotation( 90 )
            axDict[ax].xaxis.set_ticks_position('bottom')
         else:
            axDict[ ax ].set_xticks( [] )
         axDict[ ax ].set_yscale( 'log' )
         axDict[ax].set_ylim( [ min( getVals( assembliesList, ax ))*.6, 
                                maxesDict[ ax ]* 1.6] )
         
      else:
         axDict[ax].set_xticks( range( 0, len( assembliesList ) ))
         axDict[ax].set_xticklabels( getIDs( assembliesList )  )
         for tick in axDict[ax].xaxis.get_major_ticks():
            if options.subsetFile:
               tick.label1.set_fontsize( 12 )
            else:
               tick.label1.set_fontsize( 6 )
         #for label in axDict[ax].xaxis.get_ticklabels():
         #      label.set_rotation( 90 )
         axDict[ax].xaxis.set_ticks_position('bottom')
         axDict[ax].set_ylim( [ min( getVals( assembliesList, ax ))*.9, 
                                maxesDict[ ax ]* 1.05] )
      axDict[ax].set_xlim( [ -0.5, len( assembliesList )] )
      if options.sortOn == ax:
         axDict[ax].set_title( options.labelDict[ ax ] + ' (sorted)' )
      else:
         axDict[ax].set_title( options.labelDict[ ax ] )

def rankings( assembliesList, options ):
   print '#Assembly\tNA50\tContig Path N50\tScaffold Path N50'
   for a in assembliesList:
      sys.stdout.write('%s' % a.ID )
      for e in [ 'contigNA50', 'haplotypePathN50', 'scaffoldPathN50' ]:
         sys.stdout.write('\t%s' % a.valuesDict[ e ])
      sys.stdout.write('\n')

def main():
   usage = ( 'usage: %prog --hapPathStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of haplotype path stats xml files\n'
             '( --hapPathStatsDir ) named as NAME.hapPathStats.xml and creates a plot.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   cscp.initOptions( parser )
   las.initOptions( parser )
   ( options, args ) = parser.parse_args()
   cscp.checkOptions( options, parser )
   las.checkOptions( options, parser )
   checkOptions( options, parser )
   
   assembliesList = readDir( options )
   assembliesList = sorted( assembliesList, key=lambda x: x.valuesDict[ options.sortOn ], 
                            reverse=True )
   maxesDict = cnfst.calculateMaxesDict( assembliesList )

   if options.outputRanks:
      rankings( assembliesList, options )
      return

   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxis( fig, options, data )

   drawData( assembliesList, maxesDict, axDict, options )
   
   cscp.writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
