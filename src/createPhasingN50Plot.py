#!/usr/bin/env python
"""
createPhasingN50Plot.py
18 May 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the phasing N50 stats plots.

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
from createContigPathStatsTable import readDirs, Assembly
from createN50StatsPlot import getVals, getIDs
import createSortedCoveragesPlot as cscp
import glob
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
import re
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
   parser.add_option( '--statsScaffoldsContigPathPhasingDir', dest='statsScaffoldsContigPathPhasingDir',
                      type='string',
                      help=('Directory with contigPathStats. Names: A1.contigPathStats.xml .'))
   parser.add_option('--title', dest='title',
                     type='string',
                     default=( 'Phasing N50 Statistics'),
                     help=('Title of the plot. default=%default'))
   parser.add_option( '--sortOn', dest='sortOn',
                      type='string', default='scaffoldPathNG50',
                      help=('Allows a different sort order. default=%default'))
   parser.add_option( '--outputRanks', dest='outputRanks', action='store_true',
                      default=False, help=('Outputs rankings as tab delimited '
                                           'stream to STDOUT. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'statsScaffoldsContigPathDir'        : options.statsScaffoldsContigPathDir,
            'statsContigsContigPathDir'          : options.statsContigsContigPathDir,
            'statsScaffoldsContigPathPhasingDir' : options.statsScaffoldsContigPathPhasingDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   options.columns = [ 'scaffoldPathNG50', 'contigPathNG50', 
                       'hap1ScaffoldPathN50', 'hap2ScaffoldPathN50',
                       'hap1ContigPathN50', 'hap2ContigPathN50'  ]
   if options.outputRanks:
      return
   allowedKeys = set([ 'contigPathNG50', 'scaffoldPathNG50', 
                       'hap1ContigPathN50', 'hap1ScaffoldPathN50',
                       'hap2ContigPathN50', 'hap2ScaffoldPathN50' ])
   if options.sortOn not in allowedKeys:
      parser.error('--sortOn %s is not in the dict of allowed keys: %s' % 
                   ( options.sortOn, allowedKeys ))
   options.columnLabels = { 
                            'scaffoldPathNG50':'Scaffold Path NG50',
                            'contigPathNG50':'Contig Path NG50',
                            'hap1ScaffoldPathN50':'Hap1 Scaffold Path N50',
                            'hap2ScaffoldPathN50':'Hap2 Scaffold Path N50',
                            'hap2ContigPathN50':'Hap2 Contig Path N50', 
                            'hap1ContigPathN50':'Hap1 Contig Path N50'
                            }
   options.colors = { 'hap1ContigPathN50':'#1f77b4',  # darker blue
                      'hap1ScaffoldPathN50':'#ff7f0e',  # darker orange
                      'hap2ContigPathN50':'#aec7e8',    # lighter blue
                      'hap2ScaffoldPathN50':'#ffbb78',    # light orange
                      'contigPathNG50':(0.8, 0.8, 0.8),  # light gray
                      'scaffoldPathNG50':(0.3, 0.3, 0.3) # dark gray
                      }
   options.shapes = { 'hap1ContigPathN50':'.',
                      'hap1ScaffoldPathN50':'.',
                      'hap2ContigPathN50':'.',
                      'hap2ScaffoldPathN50':'.',
                      'contigPathNG50':'.',
                      'scaffoldPathNG50':'.'
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

def readPhasingDir( hapNum, options ):
   files = glob.glob( os.path.join( options.statsScaffoldsContigPathPhasingDir, 
                                     '*.hap%d.pathStats.xml' % hapNum))
   namepat = re.compile( r'^(\S{2,3})\.hap\d\.pathStats\.xml' )
   assembliesDict = {}
   for f in files:
      name = re.match( namepat, os.path.basename( f )).group( 1 )
      if 'subsetFile' in vars( options ):
         if options.subsetFile:
            if name not in options.assemblySubset:
               continue
      try:
         xmlTree = ET.parse( f )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( f )
      root=xmlTree.getroot()
      if name not in assembliesDict:
         a = Assembly()
         a.ID = name
      else:
         a = assembliesDict[ name ]
      
      a.valuesDict[ 'contigPathNG50' ] = int( root.attrib[ 'contigPathNG50' ] )
      a.valuesDict[ 'scaffoldPathNG50' ] = int( root.attrib[ 'scaffoldPathNG50' ] )
      if name not in assembliesDict:
         assembliesDict[ name ] = a
   if len(assembliesDict) == 0:
      sys.stderr.write('Error, no phasing information.\n')
      sys.exit(1)
   return assembliesDict.values()

def mergeLists( aList, h1List, h2List ):
   aDict = {}
   for a in aList:
      aDict[a.ID] = a
   i = 0
   for hapList in [ h1List, h2List ]:
      i += 1
      for h in hapList:
         if h.ID not in aDict:
            sys.stderr.write('Error, ID %s found in haplotype %d path '
                             'file but not in standard path file.\n' % ( h.ID, i ) )
            sys.exit(1)
         aDict[h.ID].valuesDict[ 'hap%dContigPathN50'%i ] = h.valuesDict['contigPathNG50']
         aDict[h.ID].valuesDict[ 'hap%dScaffoldPathN50'%i ] = h.valuesDict['scaffoldPathNG50']
   return aDict.values()

def readData( options ):
   aList = readDirs( options )
   hap1List = readPhasingDir( 1, options )
   hap2List = readPhasingDir( 2, options )
   assembliesList = mergeLists( aList, hap1List, hap2List )
   return assembliesList

def findMaxMin( assembliesList, options ):
   theMax = 0
   theMin = sys.maxint
   for a in assembliesList:
      for c in options.columns:
         if c not in a.valuesDict:
            sys.stderr.write('Error, can not find %s in %s.\n' % (c, a.ID))
            sys.exit(1)
         if theMax < a.valuesDict[ c ]:
            theMax = a.valuesDict[ c ]
         if theMin > a.valuesDict[ c ]:
            theMin = a.valuesDict[ c ]
   return ( theMax, theMin )

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
   nudge = -0.25
   span = abs(2.0 * nudge)
   for c in options.columns:
      plots.append( ax.plot( numpy.arange(0, len( assembliesList )) + nudge,
                                         getVals( assembliesList, c ), 
                                         marker=options.shapes[c], color=options.colors[c], markersize=18.0,
                                         linestyle='none', markeredgecolor='w'))
      if c == 'hap1ScaffoldPathN50':
         vals = getVals( assembliesList, c )
         nextVals = getVals( assembliesList, 'hap2ScaffoldPathN50' )
         i = -1
         for v in vals:
            i += 1
            ax.add_line( lines.Line2D( xdata=[i + nudge, i + nudge + (span / float(len(options.columns)))],
                                        ydata=[v, nextVals[i]],
                                        linewidth=.75,
                                        color='red'))
      elif c == 'hap1ContigPathN50':
         vals = getVals( assembliesList, c )
         nextVals = getVals( assembliesList, 'hap2ContigPathN50' )
         i = -1
         for v in vals:
            i += 1
            ax.add_line( lines.Line2D( xdata=[i + nudge, i + nudge + (span / float(len(options.columns)))],
                                        ydata=[v, nextVals[i]],
                                        linewidth=.75,
                                        color='red'))
      nudge +=  span / float(len(options.columns))
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
   leg = plt.legend( plots, legendLabels, 'upper right', numpoints=1 )
   leg._drawFrame=False
   plt.ylabel('Bases')

def main():
   usage = ( 'usage: %prog --statsScaffoldsContigPathDir=path/to/dir/ '
             '--statsContigssContigPathDir=path/to/dir/ '
             '--statsScaffoldsContigPathPhasingDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of scaffold-alignment contig path stats xml files\n'
             '( --statsScaffoldsContigPathDir ) named as NAME.pathStats.xml, contig-alignment '
             'contig path stats xml files ( --statsContigsContigPathDir ) named as NAME.pathStats.xml,'
             'scaffold-alignment contig path phasing stats xml files ( --statsScaffoldsContigPathPhasingDir )'
             ' named as NAME.hap%d.pathStats.xml and creates a  plot.\n' )
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
   
   assembliesList = readData( options )
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
