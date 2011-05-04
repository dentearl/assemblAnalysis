#!/usr/bin/env python
"""
createN50StatsPlot.py
14 March 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the N50 stats plots.

"""
##############################
# Copyright (C) 2009-2011 by 
# Dent Earl (dearl@soe.ucsc.edu, dent.earl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedict.paten@gmail.com)
# Mark Diekhans (markd@soe.ucsc.edu)
# ... and other members of the Reconstruction Team of David Haussler's 
# lab (BME Dept. UCSC).
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
##############################
from createContigPathStatsTable import readDirs
import createSortedCoveragesPlot as cscp
import libAssemblySubset as las
import libGeneral as lgn
import libPlotting as lpt
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

class Data:
   """Dummy class to hold data to 
   pass between functions
   """
   pass

def initOptions( parser ):
   parser.add_option( '--statsScaffoldsContigPathDir', dest='statsScaffoldsContigPathDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.contigPathStats.xml .'))
   parser.add_option( '--statsContigsContigPathDir', dest='statsContigsContigPathDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.contigPathStats.xml .'))
   parser.add_option('--title', dest='title',
                     type='string',
                     default=( 'N50 Statistics'),
                     help=('Lables of the facets, comma separated, from top '
                           'to bottom: default=%default'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='contigNG50',
                      help=('Allows a different sort order. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks', action='store_true',
                      default=False, help=('Outputs rankings as tab delimited '
                                           'stream to STDOUT. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'statsScaffoldsContigPathDir' : options.statsScaffoldsContigPathDir,
            'statsContigsContigPathDir'   : options.statsContigsContigPathDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   options.columns = [ 'blockNG50', 'contigPathNG50', 'scaffoldPathNG50', 
                       'contigN50', 'contigNG50', 'scaffoldN50', 'scaffoldNG50' ]
   if options.outputRanks:
      return
   allowedKeys = set([ 'blockNG50', 'contigPathNG50', 'scaffoldPathNG50', 
                       'scaffoldN50', 'scaffoldNG50', 'contigNG50', 'contigN50' ])
   if options.sortOn not in allowedKeys:
      parser.error('--sortOn %s is not in the dict of allowed keys: %s' % 
                   ( options.sortOn, allowedKeys ))
   options.columnLabels = { 'blockNG50':'Block NG50', 
                            'contigPathNG50':'Contig Path NG50',
                            'scaffoldPathNG50':'Scaffold Path NG50',
                            'contigN50':'Contig N50',
                            'contigNG50':'Contig NG50',
                            'scaffoldN50':'Scaffold N50',
                            'scaffoldNG50':'Scaffold NG50'
                            }
   options.colors = { 'blockNG50':'#CCFFCC',        # light greenish
                      'contigPathNG50':'#ffbb78',   # light orange
                      'scaffoldPathNG50':'#ff7f0e', # darker orange
                      'contigN50':'#aec7e8',        # lighter blue
                      'contigNG50':'#1f77b4',       # darker blue
                      'scaffoldN50':(0.8, 0.8, 0.8), # light gray
                      'scaffoldNG50':(0.3, 0.3, 0.3) # dark gray
                      }
   options.shapes = { 'blockNG50':'.', 
                      'contigPathNG50':'.',
                      'scaffoldPathNG50':'.',
                      'contigN50':'.',
                      'contigNG50':'.',
                      'scaffoldN50':'.',
                      'scaffoldNG50':'.'
                      }

def establishAxis( fig, options, data ):
   """ 
   """
   axDict = {}
   # 'totalContigNumber', 'contigN50'
   options.axLeft    = 0.075
   options.axRight   = 0.98
   options.axWidth   = options.axRight - options.axLeft
   options.axBottom  = 0.08
   options.axTop     = 0.95
   options.axHeight  = options.axTop - options.axBottom
   options.margin    = 0.08
   indvHeight = ( options.axHeight - (len( options.columns) - 1.0) * options.margin ) / float( len( options.columns ))
   prevY = options.axTop
   i = -1
   axDict['main'] = fig.add_axes( [options.axLeft, options.axBottom,
                                   options.axWidth, options.axHeight ] )
   return axDict

def getVals( assembliesList, key ):
   """ returns a list of values given a key and the list of assemblies.
   """
   v = []
   for a in assembliesList:
      v.append( a.valuesDict[ key ] )
   return v

def getIDs( assembliesList ):
   """ returns a list of the names of the assemblies
   """ 
   v = []
   for a in assembliesList:
      v.append( lgn.idMap[ a.ID[0] ] )
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
   plots = []
   # jitter the odd numbered points 0.1 to the 
   # left and the even numbered points 0.1 to the right
   jitter = 0.1
   side = 1
   for c in options.columns:
      side = -side
      plots.append( axDict['main'].plot( numpy.arange(0, len( assembliesList )) + jitter*side,
                                         getVals( assembliesList, c ), 
                                         marker=options.shapes[c], color=options.colors[c], markersize=18.0,
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
         for label in axDict['main'].xaxis.get_ticklabels():
            label.set_rotation( 45 )
      else:
         tick.label1.set_fontsize( 6 )
      # if len(assembliesList) > 25:
      #    for label in axDict['main'].xaxis.get_ticklabels():
      #       label.set_rotation( 45 )
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
   axDict['main'].set_title( options.title )
   legendLabels = []
   for c in options.columns:
      legendLabels.append( options.columnLabels[c] )
   # I want the legend to print in the same order as the data appears
   legendLabels.reverse()
   plots.reverse()
   leg = plt.legend( plots, legendLabels, 'upper right', numpoints=1 )
   leg._drawFrame=False

def rankings( assembliesList, options ):
   print '#Assembly\tNG50\tContig Path NG50\tScaffold Path NG50\tContig N50\tScaffold N50'
   for a in assembliesList:
      sys.stdout.write('%s' % a.ID )
      for e in [ 'scaffoldNG50', 'contigPathNG50', 'scaffoldPathNG50', 'contigN50', 'scaffoldN50' ]:
         sys.stdout.write('\t%s' % a.valuesDict[ e ])
      sys.stdout.write('\n')

def findMaxMin( assembliesList, options ):
   theMax = 0
   theMin = sys.maxint
   for a in assembliesList:
      for c in options.columns:
         if theMax < a.valuesDict[ c ]:
            theMax = a.valuesDict[ c ]
         if theMin > a.valuesDict[ c ]:
            theMin = a.valuesDict[ c ]
   return ( theMax, theMin )

def main():
   usage = ( 'usage: %prog --statsScaffoldsContigPathDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of contig path stats xml files\n'
             '( --statsScaffoldsContigPathDir ) named as NAME.contigPathStats.xml and creates a plot.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   cscp.initOptions( parser )
   las.initOptions( parser )
   lpt.initOptions( parser )
   options, args = parser.parse_args()
   cscp.checkOptions( options, parser )
   las.checkOptions( options, parser )
   lpt.checkOptions( options, parser )
   checkOptions( options, parser )
   
   assembliesList = readDirs( options )
   assembliesList = sorted( assembliesList, key=lambda x: x.valuesDict[ options.sortOn ], 
                            reverse=True )

   maxesMax, minsMin = findMaxMin( assembliesList, options )
   if options.outputRanks:
      rankings( assembliesList, options )
      return

   fig, pdf = lpt.initImage( 8.0, 8.0, options, data )
   axDict = establishAxis( fig, options, data )

   drawData( assembliesList, maxesMax, minsMin, axDict, options )
   
   lpt.writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
