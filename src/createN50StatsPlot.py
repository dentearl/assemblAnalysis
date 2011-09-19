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
                      help=('Directory with contigPathStats. Names: A1.pathStats.xml .'))
   parser.add_option( '--statsContigsContigPathDir', dest='statsContigsContigPathDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.pathStats.xml .'))
   parser.add_option('--title', dest='title',
                     type='string',
                     default=( 'N50 Statistics'),
                     help=('Title of the plot. default=%default'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='contigNG50',
                      help=('Allows a different sort order. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks', action='store_true',
                      default=False, help=('Outputs rankings as tab delimited '
                                           'stream to STDOUT. default=%default'))
   parser.add_option('--cheapskates', dest = 'cheapskateMode', default = False, 
                     action = 'store_true',
                     help = 'Turns on garbage mode.')

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
   if options.cheapskateMode:
      options.colors = { 'blockNG50':(0.3, 0.3, 0.3), # dark gray
                         'contigPathNG50':(0.8, 0.8, 0.8), # light gray
                         'scaffoldPathNG50':(0.3, 0.3, 0.3), # dark gray
                         'contigN50':(0.8, 0.8, 0.8), # light gray
                         'contigNG50':(0.3, 0.3, 0.3), # dark gray
                         'scaffoldN50':(0.8, 0.8, 0.8), # light gray
                         'scaffoldNG50':(0.3, 0.3, 0.3) # dark gray
                         }
      options.shapes = { 'blockNG50':'v', 
                         'contigPathNG50':'^',
                         'scaffoldPathNG50':'^',
                         'contigN50':'s',
                         'contigNG50':'s',
                         'scaffoldN50':'.',
                         'scaffoldNG50':'.'
                         }
      options.sizes = { 'blockNG50':10., 
                        'contigPathNG50':10.,
                        'scaffoldPathNG50':10.,
                        'contigN50':8.5,
                        'contigNG50':8.5,
                        'scaffoldN50':18.,
                        'scaffoldNG50':18.
                        }
   else:
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
      options.sizes = { 'blockNG50':18., 
                        'contigPathNG50':18.,
                        'scaffoldPathNG50':18.,
                        'contigN50':18.,
                        'contigNG50':18.,
                        'scaffoldN50':18.,
                        'scaffoldNG50':18.
                        }

def establishAxis( fig, options, data ):
   """ 
   """
   axDict = {}
   options.axLeft    = 0.09
   options.axRight   = 0.98
   options.axWidth   = options.axRight - options.axLeft
   options.axBottom  = 0.08
   options.axTop     = 0.95
   options.axHeight  = options.axTop - options.axBottom
   options.margin    = 0.08
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

def getIDs( assembliesList, options ):
   """ returns a list of the names of the assemblies
   """ 
   v = []
   for a in assembliesList:
      if options.subsetFile:
         v.append( lgn.idMap[ a.ID[0] ] )
      else:
         v.append( lgn.idMap[ a.ID[0] ]+'.'+a.ID[1:] )
   return v

def drawData( assembliesList, maxesMax, minsMin, axDict, options ):
   ax = axDict['main']
   # partition
   for i in xrange(1, len( assembliesList )+1):
      if not i % 5:
         ax.add_line( lines.Line2D( xdata=[ i-1, i-1 ],
                                                ydata=[ 1, maxesMax * 1.6 ],
                                                linewidth=1.0,
                                                linestyle='dotted',
                                                color=(0.8, 0.8, 0.8) ))
   plots = []
   # nudge the odd numbered points 0.1 to the 
   # left and the even numbered points 0.1 to the right
   nudge = 0.1
   side = 1
   for c in options.columns:
      side = -side
      plots.append( ax.plot( numpy.arange(0, len( assembliesList )) + nudge*side,
                             getVals( assembliesList, c ), 
                             marker=options.shapes[c], color=options.colors[c], 
                             markersize=options.sizes[c],
                             linestyle='none', markeredgecolor='w'))
   for loc, spine in ax.spines.iteritems():
      if loc in [ 'left'  ]:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in [ 'top',  'right' ]:
         spine.set_color('none') # don't draw spine               
      elif loc in [ 'bottom' ]:
         pass
      else:
         raise ValueError('unknown spine location: %s' % loc )
   ax.set_xticks( range( 0, len( assembliesList ) ))
   ax.set_xticklabels( getIDs( assembliesList, options )  )
   for tick in ax.xaxis.get_major_ticks():
      if options.subsetFile:
         tick.label1.set_fontsize( 12 )
      else:
         tick.label1.set_fontsize( 6 )
      for label in ax.xaxis.get_ticklabels():
            label.set_rotation( 45 )
      ax.xaxis.set_ticks_position('bottom')
   ax.set_yscale( 'log' )
   ax.set_ylim( [ minsMin *.6, 
                              maxesMax * 1.6] )
  # grid
   mts = ax.yaxis.get_majorticklocs()
   for m in mts:
      ax.add_line( lines.Line2D( xdata=[ 0, len(assembliesList) - 1 ],
                                             ydata=[ m, m ],
                                             linewidth=1,
                                             color=(0.8, 0.8, 0.8),
                                             linestyle='dotted'))
   ax.set_xlim( [ -0.5, len( assembliesList )] )
   ax.set_title( options.title )
   legendLabels = []
   for c in options.columns:
      legendLabels.append( options.columnLabels[c] )
   # I want the legend to print in the same order as the data appears
   legendLabels.reverse()
   plots.reverse()
   leg = plt.legend( plots, legendLabels, 'upper right', numpoints=1 )
   leg._drawFrame=False
   plt.ylabel('Bases')

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
   usage = ( 'usage: %prog --statsScaffoldsContigPathDir=path/to/dir/ '
             '--statsContigssContigPathDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of scaffold-alignment contig path stats xml files\n'
             '( --statsScaffoldsContigPathDir ) named as NAME.pathStats.xml, contig-alignment '
             'contig path stats xml files ( --statsContigsContigPathDir ) named as NAME.pathStats.xml,'
             ' and creates a plot.\n')
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

   fig, pdf = lpt.initImage( 10.0, 8.0, options, data )
   axDict = establishAxis( fig, options, data )

   drawData( assembliesList, maxesMax, minsMin, axDict, options )
   
   lpt.writeImage( fig, pdf, options )

if __name__ == '__main__':
   main()
