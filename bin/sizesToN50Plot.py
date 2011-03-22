#!/usr/bin/env python
"""
sizesToN50Plot.py
9 March 2011
dent earl, dearl(a)soe ucsc edu

This script takes two files which consit of sizes, one per line,
and a genome length, it produces a figure showing the cumulative
plot of the N statistic for both files.

"""

import matplotlib.backends.backend_pdf as pltBack
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

def initOptions( parser ):
   parser.add_option( '--scaffoldsFile', dest='scaffoldsFile',
                      type='string',
                      help='First size file.' )
   parser.add_option( '--contigsFile', dest='contigsFile',
                      type='string',
                      help='Second size file.' )
   parser.add_option( '--size', dest='size',
                      type='int',
                      help='Total size of the genome.' )
   parser.add_option( '--out', dest='out', default='myAggPlot',
                      type='string',
                      help='output pdf where figure will be created. No extension.' )
   parser.add_option( '--title', dest='title',
                      type='string',
                      help='Title of the plot.' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps]' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output.')
   parser.add_option( '--log', dest='log', default=False,
                      action='store_true',
                      help='Puts y axis into log scale.')
   parser.add_option( '--n50Line', dest='n50Line', default=False,
                      action='store_true',
                      help=('Adds straight lines from y axis and x axis to the curves.'))
   
   
def checkOptions( options, parser ):
   if options.scaffoldsFile == None:
      parser.error( 'Error, specify --scaffoldsFile.\n' )
   if not os.path.exists( options.scaffoldsFile ):
      parser.error( 'Error, --scaffoldsFile %s does not exist.\n' % options.scaffoldsFile )
   if options.contigsFile == None:
      parser.error( 'Error, specify --contigsFile.\n' )
   if not os.path.exists( options.contigsFile ):
      parser.error( 'Error, --contigsfile %s does not exist.\n' % options.contigsFile )
   if options.size == None:
      parser.error('Error, specify --size\n')
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )
   if options.title == None:
      parser.error('Error, specify --title.\n')
   if ( options.out[-4:] == '.png' or options.out[-4:] == '.pdf' or 
        options.out[-4:] == '.eps' ):
      options.out = options.out[:-4]
   

def readFile( filename ):
   d = []
   f = open( filename, 'r' )
   for line in f:
      line = line.strip()
      d.append( int( line ))
   f.close()
   return d

def initImage( options ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 5), dpi=options.dpi, facecolor='w' )
   return ( fig, pdf )

def establishAxis( fig, options ):
   """ create one axes per chromosome
   """
   options.axLeft  = 0.1
   options.axWidth = 0.85
   options.axBottom  = 0.15
   options.axHeight  = 0.75
   ax = fig.add_axes( [options.axLeft, options.axBottom,
                       options.axWidth, options.axHeight ] )
   #ax = fig.add_subplot(111)
   #plt.box( on= False )
   return ax

def drawData( scaffolds, contigs, ax, options ):
   ax.set_title( options.title + ' N Stats' )
   # create the N50 line
   globalMin = min( min( scaffolds['values']), min( contigs['values'] ))
   if options.n50Line:
      color50 = ( 0.4, 0.4, 0.4 )
      # vertical line
      # ax.add_line( lines.Line2D( xdata=[ 0.5, 0.5],
      #                            ydata=[ globalMin,
      #                                    scaffolds['values'][ - sum( numpy.array( scaffolds[ 'xData' ] ) > 0.5 ) ]],
      #                            color=color50,
      #                            linewidth= 0.75,
      #                            linestyle= ':'))
      for d in [ scaffolds, contigs ]:
         # horizontal lines
         ax.add_line( lines.Line2D( xdata=[ 0.0, 0.5],
                                    ydata=[ d['values'][ - sum( numpy.array( d[ 'xData' ] ) > 0.5 ) ],
                                            d['values'][ - sum( numpy.array( d[ 'xData' ] ) > 0.5 ) ]],
                                    color=color50,
                                    linewidth= 0.75,
                                    linestyle= ':'))
   
   p1 = ax.plot( scaffolds['xData'], scaffolds['values'], color='#1f77b4' )
   p2 = ax.plot( contigs['xData'], contigs['values'], color='#aec7e8' )
   for loc, spine in ax.spines.iteritems():
      if loc in ['left','bottom']:
         spine.set_position(('outward',10)) # outward by 10 points
      elif loc in ['right','top']:
         spine.set_color('none') # don't draw spine               
      else:
         raise ValueError('unknown spine location: %s' % loc )

   if options.log:
      ax.set_yscale('log')
      plt.ylabel('log Size')
      ax.yaxis.set_minor_locator( LogLocator( base=10, subs = range(1,10) ) )
   else:
      plt.ylabel('Size')

   ax.set_xticks( [ 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 ] )
   ax.xaxis.set_ticklabels( [ 0, '', '', '', '', 0.5, '', '', '', '', 1.0 ] )
   # turn off ticks where there is no spine
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   plt.xlabel('Cumulative length proportional to Haplotype 2')
   
   leg = plt.legend([p1, p2], ['Scaffolds', 'Contigs'])
   leg._drawFrame=False

def writeImage( fig, pdf, options ):
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

def processData( scaffs, contigs, options ):
   scaffs.sort( reverse = True )
   contigs.sort( reverse = True )
   pScaffs = { 'values': [],
               'xData' : [] }
   pContigs = { 'values': [],
                'xData' : [] }

   cum = 0
   for i in range(0, len(scaffs)):
      cum += scaffs[i]
      pScaffs[ 'values' ].append( scaffs[i] )
      pScaffs[ 'xData' ].append( float(cum) / float( options.size ))
   cum = 0
   for i in range(0, len(contigs)):
      cum += contigs[i]
      pContigs[ 'values' ].append( contigs[i] )
      pContigs[ 'xData' ].append( float(cum) / float( options.size ))
   return ( pScaffs, pContigs )

def main():
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   
   scaffolds = readFile( options.scaffoldsFile )
   contigs   = readFile( options.contigsFile )
   
   ( pScaffs, pContigs ) = processData( scaffolds, contigs, options )
   
   ( fig, pdf ) = initImage( options  )
   ax = establishAxis( fig, options )
   
   drawData( pScaffs, pContigs, ax, options )

   writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
