#!/usr/bin/env python
"""
cactusAssemblathonAggregatePlotter.py
22 February 2011
dent earl, dearl@soe.ucsc.edu


"""
from libMafGffPlot import Data

import matplotlib.backends.backend_pdf as pltBack
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pylab  as pylab
import matplotlib.pyplot as plt
import numpy
from optparse import OptionParser
import os
import sys
import re

def initOptions( parser ):
   parser.add_option( '--file', dest='file',
                      type='string',
                      help='Aggregate file to read.' )
   parser.add_option( '--crazyMax', dest='crazyMax',
                      type='int',
                      help='Sets the height of the crazy bar plot axis.' )
   parser.add_option( '--title', dest='title',
                      type='string',
                      help='Title placed at the top of the plot.' )
   parser.add_option( '--out', dest='out', default='myAggPlot',
                      type='string',
                      help='output pdf where figure will be created. No extension.' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|both]' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output.')

def checkOptions( options, parser ):
   if options.file == None:
      parser.error( 'Error, specify --file.\n' )
   if options.out[-4:] == '.png' or options.out[-4:] == '.pdf':
      options.out = options.out[:-4]
   if not os.path.exists( options.file ):
      parser.error( 'Error, --file %s does not exist.\n' % options.file )
   options.file = os.path.abspath( options.file )
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.' % options.dpi )
   options.topBotOrder = [ 'hapA1/hapA2/!assembly', 'hapA1ORhapA2/!assembly',
                           'hapA1ORhapA2/assembly','hapA1/hapA2/assembly' ]
   

def readFile( filename ):
   f = open( filename, 'r' )
   data = {}
   for line in f:
      line = line.strip()
      d = line.split('\t')
      data[ d[0] ] = d[ 1: ]
   for d in data:
      for i in range(0, len( data[d] )):
         data[d][i] = int( data[d][i] )
   f.close()
   return data

def initImage( options ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'both':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 10), dpi=options.dpi, facecolor='w' )
   return ( fig, pdf )

def setAxisLimits( axDict, options, data ):
   axDict[ 'main' ].set_ylim( 0.0, 1.01 )
   axDict[ 'main' ].set_xlim( 0.9, 8.1 )
   #axDict[ 'main' ].xaxis.set_major_locator( pylab.NullLocator() )
   axDict[ 'crazy' ].set_ylim( 0.0, 1.02 )
   axDict[ 'crazy' ].set_xlim( 0.9, 8.6 )
   axDict[ 'crazy' ].yaxis.set_major_locator( pylab.NullLocator() )
   axDict[ 'crazy' ].xaxis.set_major_locator( pylab.NullLocator() )
   axDict[ 'blowUp' ].set_ylim( 0.9, 1.01 )
   axDict[ 'blowUp' ].set_xlim( 0.9, 8.1 )
   axDict[ 'blowUp' ].xaxis.set_major_locator( pylab.NullLocator() )
   #axDict[ 'main' ].yaxis.set_major_locator( pylab.NullLocator() )

def establishAxes( fig, options, data ):
   """ create one axes per chromosome
   """
   axDict = {}
   options.axLeft = 0.1
   options.axWidth = 0.85
   axDict[ 'main' ] = fig.add_axes( [ options.axLeft, 0.07,
                                      options.axWidth , 0.58 ] )
   plt.box( on=False )
   axDict[ 'crazy' ] = fig.add_axes( [ options.axLeft, 0.655,
                                       options.axWidth , 0.085 ] )
   plt.box( on=False )
   axDict[ 'blowUp' ] = fig.add_axes( [ options.axLeft, 0.75,
                                       options.axWidth , 0.22 ] )
   plt.box( on=False )
   setAxisLimits( axDict, options, data )
   data.axDict = axDict
   return ( axDict )

def writeImage( fig, pdf, options, data ):
   if options.outFormat == 'pdf':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
   elif options.outFormat == 'png':
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'both':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )

def prettyList( uglyList ):
   """ takes a list of numbers in str format,
   shortens their names and returns a nicer str format
   """
   pl = []
   for l in uglyList:
      if l == '0':
         pl.append('0')
      else:
         pl.append( '%.0e' % int( l ))
   return pl

def vectorAddition( v1, v2 ):
   if len( v1 ) != len( v2 ):
      sys.stderr.write( 'Error, lists are not the same length, cannot add %s to %s\n' % (v1, v2))
      sys.exit(1)
   r = []
   for i in range( 0, len(v1)):
      r.append( v1[i] + v2[i] )
   return r

def normalizeData( options, data ):
   if options.crazyMax == None:
      data.crazyMax = 0
      for v in data.valuesDict[ '!hapA1/!hapA2/assembly' ]:
         if v > data.crazyMax:
            data.crazyMax = v
   else:
      data.crazyMax = options.crazyMax
   # normalize crazy data against itself
   for i in range( 0, len( data.valuesDict[ '!hapA1/!hapA2/assembly' ])):
      data.valuesDict[ '!hapA1/!hapA2/assembly' ][i] /= float( data.crazyMax )
      
   # verify the remaining columns all have the same sum
   colSum = [ 0 ] * 8
   for i in range( 0, 8 ):
      for j in ['hapA1/hapA2/assembly','hapA1/hapA2/!assembly','hapA1/!hapA2/assembly',
                'hapA1/!hapA2/!assembly','!hapA1/hapA2/assembly','!hapA1/hapA2/!assembly']:
         colSum[ i ] += data.valuesDict[ j ][ i ]
   # check column sums
   for i in range(1, len( colSum )):
      if colSum[ 0 ] != colSum[ i ]:
         sys.stderr.write('Error, column sums do not equal one another, col 0 != col %d\n' % i)
         sys.exit(1)
   # create the OR categories where we collapse hap1/!hap2 and !hap1/hap2 into hap1ORhap2
   data.valuesDict[ 'hapA1ORhapA2/assembly' ] = vectorAddition( data.valuesDict['!hapA1/hapA2/assembly'], 
                                                                data.valuesDict['hapA1/!hapA2/assembly'] )
   data.valuesDict[ 'hapA1ORhapA2/!assembly' ] = vectorAddition( data.valuesDict['!hapA1/hapA2/!assembly'], 
                                                                 data.valuesDict['hapA1/!hapA2/!assembly'] )
   # normalize the data
   for i in range( 0, 8 ):
      for j in ['hapA1/hapA2/assembly','hapA1/hapA2/!assembly','hapA1ORhapA2/assembly',
                'hapA1ORhapA2/!assembly' ]:
         data.valuesDict[ j ][ i ] /= float( colSum[ 0 ] )
   # stack the data
   options.topBotOrder.reverse()
   for i in range( 0, 8 ):
      cumSum = 0.0
      for j in options.topBotOrder:
         data.valuesDict[ j ][ i ] += cumSum
         cumSum = data.valuesDict[ j ][ i ]
   options.topBotOrder.reverse()

def drawData( fig, options, data ):
   data.colors = [ "#9467bd", "#c5b0d5", "#17becf", 
                   "#9edae5", "#ff7f0e", "#ffbb78" ]
   i = -1
   for n in options.topBotOrder:
      i += 1
      data.axDict['main'].fill_between( x=data.xData,
                                        y1=data.valuesDict[ n ],
                                        y2=[0]*8, 
                                        facecolor = data.colors[ i ],
                                        linewidth = 0.0)
      data.axDict['blowUp'].fill_between( x=data.xData,
                                          y1=data.valuesDict[ n ],
                                          y2=[0]*8, 
                                          facecolor = data.colors[ i ], 
                                          linewidth = 0.0)
   for i in range(1, 9):
      data.axDict['crazy'].add_patch( patches.Rectangle( xy=(i, 0), 
                                                         height = data.valuesDict['!hapA1/!hapA2/assembly'][i-1],
                                                         width=0.5,
                                                         color='r',
                                                         edgecolor=None) )
   #data.axDict['crazy'].bar( range(1, 9), data.valuesDict['!hapA1/!hapA2/assembly'],
   #                          width=0.5, color='r', linewidth = 0.0 )
   #fig.text( x = options.axLeft, y = 0.8, s=data.crazyMax )

def drawLegend( options, data ):
   width = 3.0
   left = 1.2
   right = 1.2 + width
   bottom = 0.05
   height = 0.3
   top = bottom + height
   data.axDict['main'].add_patch( patches.Rectangle( xy=(left, bottom), width = width,
                                                     height = height, color = (0.95, 0.95, 0.95)))
   data.axDict['main'].text( x = (left + right)/2.0, y = top - 0.04, s = 'Legend', 
                             color = (0.1, 0.1, 0.1), horizontalalignment='center' )
   data.axDict['main'].add_line( lines.Line2D( xdata=[left+0.1, right-0.1], ydata=[top - 0.05, top - 0.05],
                                               color=(0.8, 0.8, 0.8)))
   yPos = top - 0.085
   
   legendText = ['hap1 and hap2, no Assembly',
                 'hap1 xor hap2, no Assembly',
                 'hap1 xor hap2, Assembly',
                 'hap1 and hap2, Assembly']
   options.topBotOrder.reverse()
   i = -1
   for n in options.topBotOrder:
      i += 1
      data.axDict['main'].add_patch( patches.Rectangle( xy=(left + 0.1, yPos - 0.01), width = 0.2,
                                                     height = 0.025, color = data.colors[i]) )
      data.axDict['main'].text( x = left + .35, y = yPos, s = legendText[i], 
                                color = (0.1, 0.1, 0.1), horizontalalignment='left',
                                verticalalignment='center')
      yPos -= 0.045
   data.axDict['main'].add_patch( patches.Rectangle( xy=(left + 0.1, yPos - 0.01), width = 0.2,
                                                     height = 0.025, color = 'r') )
   data.axDict['main'].text( x = left + .35, y = yPos, s = 'no hap1, no hap2, Assembly', 
                                color = (0.1, 0.1, 0.1), horizontalalignment='left',
                                verticalalignment='center')
   options.topBotOrder.reverse()

def drawAxisLabels( fig, options, data ):
   if options.title != None:
      fig.text(x = 0.5, y = 0.96, s = options.title,
               fontsize = 18, horizontalalignment='center',
               verticalalignment='bottom')
   fig.text(x = 0.5, y = 0.02, s = 'Contig length',
            fontsize = 14, horizontalalignment='center',
            verticalalignment='bottom')
   fig.text(x = options.axLeft - 0.06, y = 0.28, s = 'Stacked Proportion',
            fontsize = 14, horizontalalignment='center',
            verticalalignment='bottom',
            rotation=90)

def main():
   data = Data()
   parser = OptionParser()
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   ( fig, pdf ) = initImage( options )
   axDict = establishAxes( fig, options, data )
   
   data.valuesDict = readFile( options.file )
   data.xData = range(1, 9)
   axDict['main'].set_xticks( data.xData )
   axDict['main'].set_xticklabels( prettyList( data.valuesDict['columnLength'] ))

   normalizeData( options, data )
   drawData( fig, options, data )
   drawLegend( options, data )
   drawAxisLabels( fig, options, data )
   
   setAxisLimits( axDict, options, data )
   writeImage( fig, pdf, options, data )
   

if __name__ == '__main__':
   main()
