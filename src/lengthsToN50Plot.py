#!/usr/bin/env python
"""
lengthsToN50Plot.py
9 March 2011
dent earl, dearl(a)soe ucsc edu

This script takes two files which consit of lengths, one per line,
and a genome length, it produces a figure showing the cumulative
plot of the N statistic for both files.

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
import numpy
from optparse import OptionParser
import os
import sys
import re

class Data:
   """Dummy class to hold data to 
   pass between functions
   """
   pass

def initOptions(parser):
   parser.add_option('--scaffoldsFile', dest='scaffoldsFile',
                      type='string',
                      help='First size file.')
   parser.add_option('--contigsFile', dest='contigsFile',
                      type='string',
                      help='Second size file.')
   parser.add_option('--size', dest='size',
                      type='float',
                      help='Total size of the genome.')
   parser.add_option('--title', dest='title',
                      type='string',
                      help='Title of the plot.')
   parser.add_option('--log', dest='log', default=False,
                      action='store_true',
                      help='Puts y axis into log scale. default=%default')
   parser.add_option('--n50Line', dest='n50Line', default=False,
                      action='store_true',
                      help=('Adds straight lines from the y-axis to the curves. default=%default'))
   parser.add_option('--xlabel', dest='xlabel',
                      type='string', default='Cumulative length proportional to genome size',
                      help='Label on the x-axis. default=%default')
   
   
def checkOptions(options, parser):
   if options.scaffoldsFile is None:
      parser.error('specify --scaffoldsFile.\n')
   if not os.path.exists(options.scaffoldsFile):
      parser.error('--scaffoldsFile %s does not exist.\n' % options.scaffoldsFile)
   if options.contigsFile is None:
      parser.error('specify --contigsFile.\n')
   if not os.path.exists(options.contigsFile):
      parser.error('--contigsfile %s does not exist.\n' % options.contigsFile)
   # if options.size is None:
   #    parser.error('specify --size\n')
   if options.title is None:
      parser.error('specify --title.\n')

def readFile(filename):
   d = []
   f = open(filename, 'r')
   for line in f:
      line = line.strip()
      d.append(int(line))
   f.close()
   return d

def establishAxis(fig, options):
   """ create one axes per chromosome
   """
   options.axLeft  = 0.1
   options.axWidth = 0.85
   options.axBottom  = 0.15
   options.axHeight  = 0.75
   ax = fig.add_axes([options.axLeft, options.axBottom,
                       options.axWidth, options.axHeight])
   return ax

def drawData(scaffolds, contigs, ax, options):
   ax.set_title(options.title + ' N Stats')
   # create the N50 line
   globalMin = min(min(scaffolds['values']), min(contigs['values']))
   if options.n50Line:
      color50 = (0.4, 0.4, 0.4)
      # vertical line
      # ax.add_line(lines.Line2D(xdata=[0.5, 0.5],
      #                            ydata=[globalMin,
      #                                    scaffolds['values'][- sum(numpy.array(scaffolds['xData']) > 0.5)]],
      #                            color=color50,
      #                            linewidth= 0.75,
      #                            linestyle= ':'))
      for d in [scaffolds, contigs]:
         # horizontal lines
         ax.add_line(lines.Line2D(xdata=[0.0, 0.5],
                                    ydata=[d['values'][- sum(numpy.array(d['xData']) > 0.5)],
                                            d['values'][- sum(numpy.array(d['xData']) > 0.5)]],
                                    color=color50,
                                    linewidth= 0.75,
                                    linestyle= ':'))
   
   p1 = ax.plot(scaffolds['xData'], scaffolds['values'], color='#1f77b4')
   p2 = ax.plot(contigs['xData'], contigs['values'], color='#aec7e8')
   for loc, spine in ax.spines.iteritems():
      if loc in ['left','bottom']:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in ['right','top']:
         spine.set_color('none') # don't draw spine               
      else:
         raise ValueError('unknown spine location: %s' % loc)

   if options.log:
      ax.set_yscale('log')
      ax.yaxis.set_minor_locator(LogLocator(base=10, subs = range(1,10)))
   plt.ylabel('Length')

   ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
   ax.xaxis.set_ticklabels([0, '', '', '', '', 0.5, '', '', '', '', 1.0])
   # turn off ticks where there is no spine
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   plt.xlabel(options.xlabel)
   
   leg = plt.legend([p1, p2], ['Scaffolds', 'Contigs'])
   leg._drawFrame=False

def processData(scaffs, contigs, options):
   scaffs.sort(reverse = True)
   contigs.sort(reverse = True)
   pScaffs = { 'values': [],
               'xData' : [] }
   pContigs = { 'values': [],
                'xData' : [] }

   cum = 0
   for i in xrange(0, len(scaffs)):
      cum += scaffs[i]
      pScaffs['values'].append(scaffs[i])
      if options.size is not None:
         pScaffs['xData'].append(float(cum) / float(options.size))
      else:
         pScaffs['xData'].append(float(cum))
   cum = 0
   for i in xrange(0, len(contigs)):
      cum += contigs[i]
      pContigs['values'].append(contigs[i])
      if options.size is not None:
         pContigs['xData'].append(float(cum) / float(options.size))
      else:
         pContigs['xData'].append(float(cum))
   if options.size is None:
      options.size = max(pContigs['xData'][-1], pScaffs['xData'][-1])
      for d in (pContigs['xData'], pScaffs['xData']):
         for i, e in enumerate(d):
            d[i] = e / float(options.size)
   return (pScaffs, pContigs)

def main():
   usage = ('usage: %prog --scaffoldsFile=sFile.txt --contigsFile=cFile.txt --size=N --title=TITLE\n\n'
             '%prog takes in a scaffolds file (--scaffoldsFile), a contigs\n'
             'file (--contigs), the size of the genome (--size) and a title (--title)\n'
             'and then produces an N50 style figure.')
   data = Data()
   parser = OptionParser(usage=usage)
   initOptions(parser)
   lpt.initOptions(parser)
   options, args = parser.parse_args()
   checkOptions(options, parser)
   lpt.checkOptions(options, parser)
   
   scaffolds = readFile(options.scaffoldsFile)
   contigs   = readFile(options.contigsFile)
   
   pScaffs, pContigs = processData(scaffolds, contigs, options)
   fig, pdf = lpt.initImage(8.0, 5.0, options, data)
   ax = establishAxis(fig, options)
   
   drawData(pScaffs, pContigs, ax, options)

   lpt.writeImage(fig, pdf, options)
   

if __name__ == '__main__':
   main()
