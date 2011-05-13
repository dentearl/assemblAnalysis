#!/usr/bin/env python
""" createIndelistributionPlot.py
dent earl, dearl (a) soe ucsc edu
13 May 2011

Script to look at the distribution of insertions and deletions
detected in each assembly.
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
from createContigPathStatsTable import readDir
import libAssemblySubset as las
import libGeneral as lgn
import libPlotting as lpt
import math
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
   parser.add_option('--title', dest='title',
                     type='string',
                     default=( 'Indel Distribution'),
                     help=('Title of the plot. default=%default'))
   parser.add_option( '--normalize', dest='normalize', default='self', 
                      help=('Normalization method. May either be self, global, max, log2.'
                            'Global sums all assemblies, max takes the max. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'statsScaffoldsContigPathDir' : options.statsScaffoldsContigPathDir}
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   options.normalize = options.normalize.lower()
   if options.normalize not in ['self', 'global', 'max', 'log2']:
      parser.error( '--normalize %s is not recognized. Must be '
                    'either self, global, or max' % options.normalize )

def establishAxis( numAssemblies, fig, options, data ):
   """ 
   """
   axDict = {}
   if numAssemblies <= 20:
      options.numFacets = 2
      options.axLeft    = 0.07
      options.axRight   = 0.98
      options.axWidth   = options.axRight - options.axLeft
      options.axBottom  = 0.08
      options.axTop     = 0.95
      options.axHeight  = options.axTop - options.axBottom
      options.margin    = 0.08
      indvHeight = ( options.axHeight - (options.numFacets - 1.0) * options.margin ) / float( options.numFacets )
      prevY = options.axTop
      i = -1
      for n in ['insertions', 'deletions']:
         i += 1
         axDict[n] = fig.add_axes( [options.axLeft, prevY - ( indvHeight + (i * options.margin)),
                                    options.axWidth, indvHeight ] )
         prevY = prevY - ( indvHeight + (i * options.margin))
   else:
      options.numFacets = int(math.ceil(2.0 * numAssemblies / 20.0))
      options.axLeft    = 0.05
      options.axRight   = 0.98
      options.axWidth   = options.axRight - options.axLeft
      options.axBottom  = 0.03
      options.axTop     = 0.98
      options.axHeight  = options.axTop - options.axBottom
      options.margin    = 0.03
      indvHeight = ( options.axHeight - (options.numFacets - 1.0) * options.margin ) / float( options.numFacets )
      prevY = options.axTop
      for n in ['insertions', 'deletions']:
         for j in xrange(0, int(options.numFacets/2.0)):
            axDict['%s%d' % (n, j)] = fig.add_axes( [options.axLeft, prevY - indvHeight,
                                                     options.axWidth, indvHeight ] )
            prevY = prevY - ( indvHeight + options.margin)
      
   return axDict

def createXYData( aList, options ):
   """ Take the data from aList, which is in a raw state,
   and collapse it into counts which will later be normalized
   """
   for t in ['insertionErrorSizeDistribution', 'deletionErrorSizeDistribution']:
      for a in aList:
         counts = {}
         for i in a.valuesDict[t]:
            if i not in counts:
               counts[i] = 0
            counts[i] += 1
         xd = counts.keys()
         xd.sort()
         yd = []
         for x in xd:
            yd.append( counts[x] )
         if 'xData' not in vars(a):
            a.xData = { t: numpy.array(xd) }
            a.yData = { t: yd }
         else:
            a.xData[t] = numpy.array(xd)
            a.yData[t] = yd
   return aList

def normalizeDist( aList, options ):
   """ takes the sum of the indel dist list, divides each
   member by that sum to normalize the dist.
   """
   for t in ['insertionErrorSizeDistribution', 'deletionErrorSizeDistribution']:
      normBy = 0.0
      if options.normalize == 'global':
         for a in aList:
            a.yData[t] = numpy.array( a.yData[t], dtype='float' )
            normBy += sum( a.yData[t] )
         for a in aList:
            a.yData[t] /= normBy
      elif options.normalize == 'max':
         for a in aList:
            a.yData[t] = numpy.array( a.yData[t], dtype='float' )
            if normBy < sum( a.yData[t] ):
               normBy = sum( a.yData[t] )
         for a in aList:
            a.yData[t] /= normBy
      elif options.normalize == 'self':
         for a in aList:
            a.yData[t] = numpy.array( a.yData[t], dtype='float' )
            a.yData[t] /= sum( a.yData[t] )
      elif options.normalize == 'log2':
         for a in aList:
            a.yData[t] = numpy.array( a.yData[t], dtype='float' )
            a.yData[t] = numpy.log2( a.yData[t] )
            a.xData[t] = numpy.log2( a.xData[t] )
      else:
         raise RuntimeError('Error, unexpected value for options.normalize: %s\n' % options.normalize)
   return aList

def drawData( assembliesList, axDict, options, data ):
   colors = [ "#1f77b4", "#aec7e8", # blues 
              "#ff7f0e", "#ffbb78", # oranges
              "#2ca02c", "#98df8a", # greens
              "#d62728", "#ff9896", # reds
              "#9467bd", "#c5b0d5", # lavenders
              "#8c564b", "#c49c94", # browns
              "#e377c2", "#f7b6d2", # strawberry pinks
              "#7f7f7f", "#c7c7c7"  # greys
              ]
   styles = { 0:'solid', 1:'dashed', 2:'dashdot', 3:'dotted' }
   data.pltList = [] # used for legends
   if len(assembliesList) <= 20:
      axToData = { 'insertions': 'insertionErrorSizeDistribution',
                   'deletions': 'deletionErrorSizeDistribution' }
      for n in ['insertions', 'deletions']:
         styleIndex = -1
         colorIndex = -1
         for a in assembliesList:
            styleIndex = ( styleIndex + 1 ) % len( styles )
            if not styleIndex: colorIndex += 1
            p = axDict[n].plot( a.xData[axToData[n]], a.yData[axToData[n]], 
                                color=colors[ colorIndex % len( colors ) ], 
                                linestyle=styles[styleIndex],
                                linewidth=2.0)
            if n == 'insertions':
               # it's symmetrical, we only need to record one or the other
               data.pltList.append( p )
         for loc, spine in axDict[n].spines.iteritems():
            if loc in ['left','bottom']:
               spine.set_position(('outward',10)) # outward by 10 points
            elif loc in ['right','top']:
               spine.set_color('none') # don't draw spine               
            else:
               raise ValueError('unknown spine location: %s' % loc )
         if options.normalize != 'log2':
            axDict[n].set_xscale( 'log' )
         axDict[n].set_title( n )
      maxX = max( axDict['insertions'].axis()[1], axDict['deletions'].axis()[1])
      for n in ['insertions', 'deletions']:
         axDict[n].set_xlim( 1, maxX )
      axDict['deletions'].set_xlabel('Length')
      if options.normalize == 'global':
         axDict['deletions'].set_ylabel('Global proportion')
      elif options.normalize == 'self':
         axDict['deletions'].set_ylabel('Per assembly proportion')
      elif options.normalize == 'max':
         axDict['deletions'].set_ylabel('Proportion relative to max assembly')
      elif options.normalize == 'log2':
         axDict['deletions'].set_ylabel(r'$\log_2$ Count')
         axDict['deletions'].set_xlabel(r'$\log_2$ Length')
   else:
      axToData = { 'insertions': 'insertionErrorSizeDistribution',
                   'deletions': 'deletionErrorSizeDistribution' }
      for n in ['insertions', 'deletions']:
         styleIndex = -1
         colorIndex = -1
         plotIndex  = -1
         for a in assembliesList:
            styleIndex = ( styleIndex + 1 ) % len( styles )
            plotIndex += 1
            if not styleIndex: colorIndex += 1
            facetIndex = int(math.floor( plotIndex / (len(assembliesList)/float(options.numFacets / 2.0))))
            key = '%s%d' % (n, facetIndex)
            p = axDict[key].plot( a.xData[axToData[n]], a.yData[axToData[n]], color=colors[ colorIndex % len( colors ) ], linestyle=styles[styleIndex], linewidth=2.0)
            if n == 'insertions':
               # it's symmetrical, we only need to record one or the other
               data.pltList.append( p )
         for j in xrange(0, int(options.numFacets/2.0)):
            key = '%s%d' % (n, j)
            for loc, spine in axDict[key].spines.iteritems():
               if loc in ['left','bottom']:
                  spine.set_position(('outward',10)) # outward by 10 points
               elif loc in ['right','top']:
                  spine.set_color('none') # don't draw spine               
               else:
                  raise ValueError('unknown spine location: %s' % loc )
            if options.normalize != 'log2':
               axDict[key].set_xscale( 'log' )
            axDict[key].set_title( key )
      maxX = 0.0
      maxY = 0.0
      for j in xrange(0, int(options.numFacets/2.0)):
         maxX = max( maxX, axDict['insertions%d' % j].axis()[1])
         maxX = max( maxX, axDict['deletions%d' % j].axis()[1])
         maxY = max( maxY, axDict['insertions%d' % j].axis()[3])
         maxY = max( maxY, axDict['deletions%d' % j].axis()[3])
      for n in ['insertions', 'deletions']:
         for j in xrange(0, int(options.numFacets/2.0)):
            axDict['%s%d' % (n, j)].set_xlim( 1, maxX )
            axDict['%s%d' % (n, j)].set_ylim( 0, maxY )
      key = int((options.numFacets/2.0) - 1)
      axDict['deletions%d' % key ].set_xlabel('Length')
      if options.normalize == 'global':
         axDict['deletions%d' % key].set_ylabel('Global proportion')
      elif options.normalize == 'self':
         axDict['deletions%d' % key ].set_ylabel('Per assembly proportion')
      elif options.normalize == 'max':
         axDict['deletions%d' % key ].set_ylabel('Proportion relative to max assembly')
      elif options.normalize == 'log2':
         axDict['deletions%d' % key ].set_ylabel(r'$\log_2$ Count')
         axDict['deletions%d' % key ].set_xlabel(r'$\log_2$ Length')

def drawLegend( aList, axDict, options, data ):
   pltListLabels = []
   for a in aList:
      if options.subsetFile:
         pltListLabels.append( lgn.idMap[a.ID[0]] )
      else:
         pltListLabels.append( lgn.idMap[a.ID[0]]+ '.'+ a.ID[1:] )
   if len(aList) > 20:
      prev=-1
      for j in xrange(0, int(options.numFacets/2.0)):
         start  = prev+1
         finish = int( (j+1) * math.floor(len(aList) / float(options.numFacets/2.0)) )
         if j == int(options.numFacets/2.0):
            finish = int(len(aList))
         prev = finish
         leg = axDict['insertions%d' % j ].legend( data.pltList[start:finish],
                                                   pltListLabels[start:finish], 'upper right', ncol=3 )
         for t in leg.get_texts():
            t.set_fontsize('x-small')    # the legend text fontsize      
         leg._drawFrame=False
      prev = -1
      for j in xrange(0, int(options.numFacets/2.0)):
         start  = prev+1
         finish = int( (j+1) * math.floor(len(aList) / float(options.numFacets/2.0)) )
         if j == int(options.numFacets/2.0):
            finish = int(len(aList))
         prev = finish
         leg = axDict['deletions%d' % j].legend( data.pltList[start:finish], 
                                                 pltListLabels[start:finish], 'upper right', ncol=3 )
         for t in leg.get_texts():
            t.set_fontsize('x-small')
         leg._drawFrame=False
   else:
      leg = axDict['insertions'].legend( data.pltList, pltListLabels, 'upper right', ncol=2 )
      for t in leg.get_texts():
         t.set_fontsize('small')    # the legend text fontsize
      leg._drawFrame=False

def main():
   usage = ( 'usage: %prog --statsScaffoldsContigPathDir=path/to/dir/ [options]\n\n'
             '%prog takes a directory of contig path stats xml files\n'
             '( --statsScaffoldsContigPathDir ) named as NAME.contigPathStats.xml and creates a plot.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   lpt.initOptions( parser )
   options, args = parser.parse_args()
   las.checkOptions( options, parser )
   lpt.checkOptions( options, parser )
   checkOptions( options, parser )
   
   assembliesDict = readDir( options.statsScaffoldsContigPathDir, options )
   assembliesList = assembliesDict.values()
   
   if len(assembliesList) <= 20:
      fig, pdf = lpt.initImage( 14.0, 8.0, options, data )
   else:
      fig, pdf = lpt.initImage( 14.0, 24.0, options, data )
   axDict = establishAxis( len(assembliesList), fig, options, data )
   
   assembliesList = createXYData( assembliesList, options )

   assembliesList = normalizeDist( assembliesList, options )
   assembliesList = sorted( assembliesList, 
                            key=lambda x: max(x.yData['insertionErrorSizeDistribution']), 
                            reverse=True )

   drawData( assembliesList, axDict, options, data )
   drawLegend( assembliesList, axDict, options, data )
   
   lpt.writeImage( fig, pdf, options )
   
if __name__ == '__main__':
   main()
