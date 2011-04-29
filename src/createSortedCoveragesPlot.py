#!/usr/bin/env python
"""
createSortedCoveragesPlot.py
13 March 2011
dent earl, dearl(a)soe ucsc edu

This script takes a list of numbers on STDIN and produces
a pretty picture. 

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
import libPlotting as lpt
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, LogFormatter # minor tick marks
from optparse import OptionParser
import os
import sys

class Data:
   """Dummy class to hold data to 
   pass between functions
   """
   pass

class Assembly:
   def __init__( self ):
      self.name = ''
      self.tot  = -1
      self.hap1 = -1
      self.hap2 = -1
      self.bac  = -1

def initOptions( parser ):
   pass

def checkOptions( options, parser ):
   pass

def readStream( options ):
   values = []
   for line in sys.stdin:
      line = line.strip()
      d = line.split()
      v = Assembly()
      v.name = d[0]
      v.tot  = float( d[ 1 ] )
      v.hap1 = float( d[ 2 ] )
      v.hap2 = float( d[ 3 ] )
      v.bac  = float( d[ 5 ] )
      values.append( v )
   return values

def establishAxis( fig, options, data ):
   options.axLeft  = 0.12
   options.axWidth = 0.85
   options.axBottom  = 0.1
   options.axHeight  = 0.8
   ax = fig.add_axes( [options.axLeft, options.axBottom,
                       options.axWidth, options.axHeight ] )
   return ax

def extractTots( values ):
   a = []
   for v in values:
      a.append( v.tot )
   return a

def extractBacs( values ):
   a = []
   for v in values:
      a.append( v.bac )
   return a

def drawData( values, ax, options ):
   ax.set_title( 'Total coverage for all assemblies' )
   for i in xrange(1, len(values)):
      if not i % 10:
         ax.add_line( lines.Line2D( xdata=[ i-1, i-1 ],
                                    ydata=[ 0, 1 ],
                                    linewidth=0.5,
                                    color=(0.8, 0.8, 0.8) ))
   ax.add_line( lines.Line2D( xdata=[ -2, len(values) + 1 ],
                              ydata=[ 0.95, 0.95 ],
                              linewidth=0.5,
                              color=(0.8, 0.8, 0.8) ))
   bData = extractBacs( values )
   p1 = ax.plot( range(0,len(values)), bData, '.', color=(0.7, 0.7, 0.7))
   yData = extractTots( values )
   p1 = ax.plot( range(0,len(values)), yData, '.', color='#1f77b4')
   for loc, spine in ax.spines.iteritems():
      if loc in [ 'left', 'right' ]:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in [ 'top', 'bottom' ]:
         spine.set_color('none') # don't draw spine               
      else:
         raise ValueError('unknown spine location: %s' % loc )
   ax.set_xticks( [ ] )
   allValues = extractAllValues( values )
   ax.set_xlim( [ -0.5, len( yData )] )
   rng = 1.0 - min( allValues )
   lo = min( allValues ) - rng * 0.1
   ax.set_ylim( [ -0.02 , 1.01 ] )
   # turn off ticks where there is no spine
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('both') # left
   plt.ylabel('Coverage')
   ax.text(x=0.5, y=-0.06, s='Assemblies sorted by coverage',
           horizontalalignment='center', verticalalignment='top',
           transform= ax.transAxes )
   ax.set_yticks( [ 0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0 ] )
   ax.set_xticks( range( 0, len( yData ) ))
   ax.set_xticklabels( extractNames( values )  )
   for tick in ax.xaxis.get_major_ticks():    
      tick.label1.set_fontsize( 6 )                        
   for label in ax.xaxis.get_ticklabels():    
      label.set_rotation( 90 )

def extractNames( values ):
   a = []
   for v in values:
      a.append( v.name )
   return a

def extractAllValues( values ):
   a = []
   for v in values:
      a.append( v.tot )
      a.append( v.hap1 )
      a.append( v.hap2 )
      a.append( v.bac )
   return a

def main():
   usage = ( 'usage: %prog [options] < rankedAssemblies.txt\n\n'
             '%prog takes via STDIN a list of coverage values, each line formatted as:\n'
             '#assembly ave.CovBothGenomes hap1 hap2 delta bac\n'
             'P1 .9885178 .9888077 .9882265 5.782e-04 0\n'
             'B1 .9869388 .9871948 .9866798 5.257e-04 .9978954\n'
             'F5 .9869094 .9872685 .9865338 7.295e-04 .9993371\n'
             '...\n'
             'And produces a plot showing the Total and Bacterial coverages for\n'
             'all assemblies in the input.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   lpt.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   lpt.checkOptions( options, parser )
   
   valuesList = readStream( options )
   valuesList = sorted( valuesList, key=lambda key: key.tot, reverse=True )
   
   fig, pdf = lpt.initImage( 8.0, 6.0, options, data )
   ax = establishAxis( fig, options, data )
   
   drawData( valuesList, ax, options )

   lpt.writeImage( fig, pdf, options )

if __name__ == '__main__':
   main()
