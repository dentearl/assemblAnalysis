#!/usr/bin/env python
"""
createN50StatsPlot.py
14 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the N50 stats table.

output is latex

"""
from createContigPathStatsTable import readDir
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
   parser.add_option( '--contigPathStatsDir', dest='contigPathStatsDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.contigPathStats.xml .'))
   parser.add_option('--labels', dest='labels',
                     type='string',
                     default=( 'Block N50,Contig Path N50,Scaffold Path N50,Contig NA50, N50 Statistics'),
                     help=('Lables of the facets, comma separated, from top '
                           'to bottom: default=%default'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='contigNA50',
                      help=('Allows a different sort order. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks', action='store_true',
                      default=False, help=('Outputs rankings as tab delimited '
                                           'stream to STDOUT. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'contigPathStatsDir'   : options.contigPathStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   if options.outputRanks:
      return
   allowedKeys = set(['contigNA50', 'scaffoldPathN50', 'haplotypePathN50', 'blockN50'])
   if options.sortOn not in allowedKeys:
      parser.error('--sortOn %s is not in the dict of allowed keys: %s' % 
                   ( options.sortOn, allowedKeys ))
   labelList = options.labels.split(',')
   if len( labelList ) != 5:
      parser.error('--label needs to contain 5 items, comma separated.')
   k = [ 'blockN50', 'haplotypePathN50', 'scaffoldPathN50',  'contigNA50', 'main']
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
   options.axLeft    = 0.09
   options.axRight   = 0.98
   options.axWidth   = options.axRight - options.axLeft
   options.axBottom  = 0.05
   options.axHeight  = 0.9
   options.axTop     = options.axBottom + options.axHeight
   options.margin    = 0.08
   indvHeight = ( options.axHeight - (len( columns) - 1.0) * options.margin ) / float( len( columns ))
   prevY = options.axTop
   i = -1
   axDict['main'] = fig.add_axes( [options.axLeft, options.axBottom,
                                   options.axWidth, options.axHeight ] )
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

def drawData( assembliesList, maxesMax, minsMin, axDict, options ):
   # partition
   for i in xrange(1, len( assembliesList )+1):
      if not i % 5:
         axDict['main'].add_line( lines.Line2D( xdata=[ i-1, i-1 ],
                                                ydata=[ 1, maxesMax * 1.6 ],
                                                linewidth=1.0,
                                                linestyle='dotted',
                                                color=(0.8, 0.8, 0.8) ))
   columns = ['contigNA50', 'contigN50', 'scaffoldPathN50', 'haplotypePathN50','blockN50' ]
   columnLabels = { 'blockN50':'Block N50', 
                    'haplotypePathN50':'Contig Path N50',
                    'scaffoldPathN50':'Scaffold Path N50',
                    'contigN50':'Contig N50',
                    'contigNA50':'Contig NG50' }
   colors = { 'blockN50':'#d62728', # darker red
              'haplotypePathN50':'#ffbb78', # light orange
              'scaffoldPathN50':'#ff7f0e', # darker orange
              'contigN50':'#aec7e8', # lighter blue
              'contigNA50':'#1f77b4' # darker blue }
              }
   shapes = { 'blockN50':'.', 
              'haplotypePathN50':'.',
              'scaffoldPathN50':'.',
              'contigN50':'.',
              'contigNA50':'.' }
   plots = []
   # jitter the odd numbered points 0.1 to the 
   # left and the even numbered points 0.1 to the right
   jitter = 0.1
   side = 1
   for c in columns:
      side = -side
      plots.append( axDict['main'].plot( numpy.arange(0, len( assembliesList )) + jitter*side,
                                         getVals( assembliesList, c ), 
                                         marker=shapes[c], color=colors[c], markersize=18.0,
                                         linestyle='none', markeredgecolor='w'))
   for loc, spine in axDict['main'].spines.iteritems():
      if loc in [ 'left'  ]:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in [ 'top',  'right' ]:
         spine.set_color('none') # don't draw spine               
      elif loc in [ 'bottom' ]:
         pass
      else:
         raise ValueError('unknown spine location: %s' % loc )
   axDict['main'].set_xticks( range( 0, len( assembliesList ) ))
   axDict['main'].set_xticklabels( getIDs( assembliesList )  )
   for tick in axDict['main'].xaxis.get_major_ticks():
      if options.subsetFile:
         tick.label1.set_fontsize( 12 )
      else:
         tick.label1.set_fontsize( 6 )
      if len(assembliesList) > 25:
         for label in axDict['main'].xaxis.get_ticklabels():
            label.set_rotation( 90 )
      axDict['main'].xaxis.set_ticks_position('bottom')
   axDict['main'].set_yscale( 'log' )
   axDict['main'].set_ylim( [ minsMin *.6, 
                              maxesMax * 1.6] )
  # grid
   mts = axDict['main'].yaxis.get_majorticklocs()
   for m in mts:
      axDict['main'].add_line( lines.Line2D( xdata=[ 0, len(assembliesList) - 1 ],
                                             ydata=[ m, m ],
                                             linewidth=1,
                                             color=(0.8, 0.8, 0.8),
                                             linestyle='dotted'))
   axDict['main'].set_xlim( [ -0.5, len( assembliesList )] )
   axDict['main'].set_title( options.labelDict[ 'main' ] )
   legendLabels = []
   for c in columns:
      legendLabels.append( columnLabels[c] )
   leg = plt.legend( plots, legendLabels, 'upper right', numpoints=1 )
   leg._drawFrame=False

def rankings( assembliesList, options ):
   print '#Assembly\tNA50\tContig Path N50\tScaffold Path N50'
   for a in assembliesList:
      sys.stdout.write('%s' % a.ID )
      for e in [ 'contigNA50', 'haplotypePathN50', 'scaffoldPathN50' ]:
         sys.stdout.write('\t%s' % a.valuesDict[ e ])
      sys.stdout.write('\n')

def findMaxMin( assembliesList, options ):
   theMax = 0
   theMin = sys.maxint
   for a in assembliesList:
      for c in ['contigN50', 'contigNA50', 'scaffoldPathN50', 'haplotypePathN50','blockN50' ]:
         if theMax < a.valuesDict[ c ]:
            theMax = a.valuesDict[ c ]
         if theMin > a.valuesDict[ c ]:
            theMin = a.valuesDict[ c ]
   return ( theMax, theMin )

def main():
   usage = ( 'usage: %prog --contigPathStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of contig path stats xml files\n'
             '( --contigPathStatsDir ) named as NAME.contigPathStats.xml and creates a plot.')
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
   ( maxesMax, minsMin ) = findMaxMin( assembliesList, options )

   if options.outputRanks:
      rankings( assembliesList, options )
      return

   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxis( fig, options, data )

   drawData( assembliesList, maxesMax, minsMin, axDict, options )
   
   cscp.writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
