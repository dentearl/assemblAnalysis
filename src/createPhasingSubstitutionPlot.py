#!/usr/bin/env python
"""
createPhasingSubstitutionPlot.py
18 May 2011
dent earl, dearl(a) soe ucsc edu

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
import glob
import createSubStatsPlot as cssp
import libAssemblySubset as las
import libGeneral as lgn
import libPlotting as lpt
import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter # minor tick marks
import numpy
from optparse import OptionParser
import os
import signal # deal with broken pipes
import sys
import re
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Data:
   """ Dummy class
   """

def initOptions( parser ):
   parser.add_option( '--subStatsDir', dest='subStatsDir',
                      type='string',
                      help=('Directory with subStats. Names: A1.subStats.upper.xml .'))

def checkOptions( options, parser ):
   dirs = { 'subStatsDir' : options.subStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )
   options.columns = [ 'totalCallsInHeterozygous'
                       'totalCorrectHap1InHeterozygous',
                       'totalCorrectHap2InHeterozygous' ]
   options.columnLabels = { 'totalCallsInHeterozygous':'Total Correct Hets / 2.0',
                            'totalCorrectHap1InHeterozygous':'Total Correct Hets Hap 1',
                            'totalCorrectHap2InHeterozygous':'Total Correct Hets Hap 2'
                            }
   options.colors = { 'totalCorrectHap1InHeterozygous':'#1f77b4',       # darker blue
                      'totalCorrectHap2InHeterozygous':'#aec7e8',       # lighter blue
                      'totalCallsInHeterozygous':(0.8, 0.8, 0.8) # light gray
                      }
   # 'hap2ScaffoldPathN50':'#ffbb78',    # light orange
   # '':'#ff7f0e',  # darker orange
   # 'scaffoldPathNG50':(0.3, 0.3, 0.3) # dark gray
   options.shapes = { 'totalCallsInHeterozygous':'.',
                      'totalCorrectHap1InHeterozygous':'.',
                      'totalCorrectHap2InHeterozygous':'.',
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

def createValsList( aDict, options ):
   vals = []
   for a in aDict:
      vals.append((a.ID, 
                    a.subStatsLower['totalCorrectHap1InHeterozygous'],
                    a.subStatsLower['totalCorrectHap2InHeterozygous'],
                    a.subStatsUpper['totalCorrectHap1InHeterozygous'],
                    a.subStatsUpper['totalCorrectHap2InHeterozygous'],
                    )
                   )
   return vals

def drawData( aList, axDict, options, data ):
   ax = axDict['main']
   
   ax.add_line( lines.Line2D( xdata=[ -0.5, len(aList) ],
                              ydata=[ 1, 1],
                              linewidth=1.0,
                              color=(0.8, 0.8, 0.8) ))
   
   labelList = []
   i = -1
   plots = ['', '']
   for a in aList:
      i += 1
      if options.subsetFile:
         labelList.append( lgn.idMap[ a[0][0] ] )
      else:
         labelList.append( lgn.idMap[ a[0][0] ]+'.'+a[0][1:] )
      plots[0] = ax.plot( [ i-.1 ], [ float(a[1]) / a[2] ],
                          color='#1f77b4', marker='.',
                          markersize=18.0, linestyle='none')
      plots[1] = ax.plot( [ i+.1 ], [ float(a[3]) / a[4] ],
                          color='#ff7f0e', marker='.',
                          markersize=18.0, linestyle='none')

   for loc, spine in ax.spines.iteritems():
      if loc in [ 'left'  ]:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in [ 'top',  'right' ]:
         spine.set_color('none') # don't draw spine               
      elif loc in [ 'bottom' ]:
         pass
      else:
         raise ValueError('unknown spine location: %s' % loc )
   ax.set_xticks( range( 0, len( aList ) ))
   ax.set_xticklabels( labelList )
   for label in ax.xaxis.get_ticklabels():
         label.set_rotation( 45 )
   if not options.subsetFile:
      for tick in ax.xaxis.get_major_ticks():
         tick.label1.set_fontsize( 6 )
   ax.set_yscale('log')
   ax.set_xlim( [ -0.5, len( aList )] )
   # partition
   yMin = ax.axis()[2]
   yMax = ax.axis()[3]
   for i in xrange(1, len( aList )+1):
      if not i % 5:
         ax.add_line( lines.Line2D( xdata=[ i-1, i-1 ],
                                    ydata=[ yMin, yMax ],
                                    linewidth=1.0,
                                    linestyle='dotted',
                                    color=(0.8, 0.8, 0.8) ))
   plt.ylabel('log Ratio Haplotype 1 : Haplotype 2')
   plt.title('Ratios of substitutions in haplotypes')
   leg = plt.legend( plots, ['Lower', 'Upper'], 'upper left', numpoints=1 )
   leg._drawFrame=False
 
def main():
   usage = ( 'usage: %prog --subStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes in a directory of substitution stats files ( --subStatsDir )\n'
             'with filenames as NAME.subStats.[upper|lower].xml and produces a plot showing\n'
             'the difference in subs in hap1 versus hap2.\n')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   lpt.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   las.checkOptions( options, parser )
   lpt.checkOptions( options, parser )
   
   fig, pdf = lpt.initImage( 11., 8., options, data )
   axDict = establishAxis( fig, options, data )
   
   assembliesDict = {}
   assembliesDict = cssp.readSubStatsDir( assembliesDict, options )
   valuesList = createValsList( assembliesDict.values(), options )
   valuesList = sorted( valuesList, key = lambda x: float(x[1])/x[2], reverse=False)
   drawData( valuesList, axDict, options, data )
   
   lpt.writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
