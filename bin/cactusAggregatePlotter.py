#!/usr/bin/env python
"""
cactusAggregatePlotter.py
22 February 2011
dent earl, dearl(a)soe ucsc edu

This script takes an aggregate text file and produces
a pretty picture. Files look like:

[dearl@hgwdev demo]$ cat agg.A1.txt 
columnLength	0	100	1000	10000	100000	1000000	10000000	100000000
hapA1/hapA2/assembly	101277168	101277168	101180952	91579210	1828234	0	0	0
hapA1/hapA2/!assembly	9440056	9440056	9536272	19138014	108888990	110717224	110717224	110717224
hapA1/!hapA2/assembly	670370	670370	669794	594560	15208	0	0	0
hapA1/!hapA2/!assembly	805664	805664	806240	881474	1460826	1476034	1476034	1476034
!hapA1/hapA2/assembly	733168	733168	732480	664854	9869	0	0	0
!hapA1/hapA2/!assembly	799267	799267	799955	867581	1522566	1532435	1532435	1532435
!hapA1/!hapA2/assembly	195281	195281	192837	151029	694	0	0	0

"""
from libMafGffPlot import Data
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
   parser.add_option( '--file', dest='file',
                      type='string',
                      help='Aggregate file to read.' )
   parser.add_option( '--mode', dest='mode',
                      type='string', default='',
                      help='Plotting mode [scaffPaths|contigs|scaffolds|hapPaths|blocks|contamination].' )
   parser.add_option( '--crazyMax', dest='crazyMax',
                      type='int',
                      help='Sets the height of the crazy bar plot axis.' )
   parser.add_option( '--title', dest='title',
                      type='string',
                      help='Title placed at the top of the plot.' )
   parser.add_option( '--xLabel', dest='xLabel',
                      type='string', default='',
                      help='The x axis label placed at the bottom of the plot.' )
   parser.add_option( '--out', dest='out', default='myAggPlot',
                      type='string',
                      help='filename where figure will be created. No extension. default=%default' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps] default=%default' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output. default=%default' )
   parser.add_option( '--smallMultipleMode', dest='SMM',
                      action='store_true', default=False,
                      help=('Turns off the printing of the legend and other '
                            'details, turns on the printing of the --title as the ID. default=%default' ))

def checkOptions( options, parser ):
   if options.file == None:
      parser.error( 'Error, specify --file.\n' )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]
   if not os.path.exists( options.file ):
      parser.error( 'Error, --file %s does not exist.\n' % options.file )
   options.file = os.path.abspath( options.file )
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )
   modes = set(['contigs', 'contamination', 'blocks', 'hapPaths', 'scaffPaths', 'scaffolds'])
   if options.mode not in modes:
      parser.error('Error, you must specify one of the modes listed under --mode in --help.\n')
   if options.mode in set([ 'blocks', 'hapPaths', 'contigs', 'scaffPaths', 'scaffolds']):
      options.topBotOrder = [ 'hapA1/hapA2/!assembly', 'hapA1ORhapA2/!assembly',
                              'hapA1ORhapA2/assembly','hapA1/hapA2/assembly' ]
   elif options.mode == 'contamination':
      options.topBotOrder = [ 'ecoli/!assembly', 'ecoli/assembly' ]

def readFile( filename, options ):
   f = open( filename, 'r' )
   data = {}
   for line in f:
      line = line.strip()
      d = line.split('\t')
      data[ d[0] ] = d[ 1: ]
   
   #data[ 'columnLength' ] = data[ 'block/contig_lengths/haplotype_path_lengths' ]
   #del data[ 'block/contig_lengths/haplotype_path_lengths' ]
   data[ 'columnLength' ] = data[ 'category' ]
   del data[ 'category' ]
   
   redundantColumns = {}
   prev = -1
   i = -1
   for d in data[ 'columnLength' ]:
      i += 1
      if d == prev:
         redundantColumns[ i ] = True
      prev = d
   trimmedData = {}
   for d in data:
      trimmedData[ d ] = []
      for i in range(0, len( data[d])):
         if i not in redundantColumns:
            trimmedData[ d ].append( int( data[d][i] ) )
   f.close()
   return trimmedData

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=(8, 10), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   return ( fig, pdf )

def setAxisLimits( axMain, axCrazy, axBlowUp, xData, options, data ):
   axMain.set_xscale('log')
   axMain.set_xlim( 1, xData[ -1 ] )
   axMain.set_ylim( 0.0, 1.0 )
   #if options.SMM:
   #   axDict[ 'main' ].yaxis.set_major_locator( pylab.NullLocator() )
   if options.mode in set( [ 'blocks', 'contigs', 'hapPaths', 'scaffPaths','scaffolds' ]):
      if options.mode != 'hapPaths' and options.mode != 'scaffPaths':
         axCrazy.set_ylim( 0.0, 1.02 )
         axCrazy.set_xscale('log')
         axCrazy.set_xlim( 1, xData[ -1 ] )
         axCrazy.xaxis.set_ticklabels( [] )   
      
   axBlowUp.set_xscale('log')
   axBlowUp.set_xlim( 1, xData[ -1 ] )
   axBlowUp.set_ylim( 0.9, 1.0 )
   axBlowUp.xaxis.set_ticklabels( [] )   

   # turn off ticks
   for ax in [ axMain, axCrazy, axBlowUp ]:
      if not ax == None:
         ax.xaxis.set_ticks_position('bottom')
         ax.yaxis.set_ticks_position('left')

   if options.SMM:
      if options.mode != 'hapPaths':
         axCrazy.yaxis.set_major_locator( pylab.NullLocator() )
         axCrazy.xaxis.set_major_locator( pylab.NullLocator() )
         axCrazy.xaxis.set_minor_locator( pylab.NullLocator() )
      axMain.xaxis.set_major_locator( pylab.NullLocator() )
      axMain.xaxis.set_minor_locator( pylab.NullLocator() )
      axBlowUp.xaxis.set_minor_locator( pylab.NullLocator() )

   #if options.SMM:
   #   axDict[ 'blowUp' ].yaxis.set_major_locator( pylab.NullLocator() )

def establishAxes( fig, options, data ):
   axDict = {}
   options.axLeft = 0.11
   options.axWidth = 0.85
   if options.mode in ( [ 'blocks', 'contigs', 'hapPaths', 'scaffPaths', 'scaffolds' ]):
      if options.mode == 'hapPaths' or options.mode == 'scaffPaths':
         axDict[ 'crazy' ] = None
         axDict[ 'main' ] = fig.add_axes( [ options.axLeft, 0.07,
                                            options.axWidth , 0.66 ] )
         plt.box( on=False )
         axDict[ 'blowUp' ] = fig.add_axes( [ options.axLeft, 0.75,
                                              options.axWidth , 0.20 ] )
         plt.box( on=False )
      else:
         axDict[ 'main' ] = fig.add_axes( [ options.axLeft, 0.07,
                                            options.axWidth , 0.58 ] )
         plt.box( on=False )
         axDict[ 'crazy' ] = fig.add_axes( [ options.axLeft, 0.68, # 0.655
                                             options.axWidth , 0.06 ] ) # 0.085
         plt.box( on=False )
         axDict[ 'blowUp' ] = fig.add_axes( [ options.axLeft, 0.75,
                                              options.axWidth , 0.20 ] )
         plt.box( on=False )
   else:
      axDict[ 'main' ] = fig.add_axes( [ options.axLeft, 0.07,
                                         options.axWidth , 0.60 ] )
      plt.box( on=False )
      axDict[ 'blowUp' ] = fig.add_axes( [ options.axLeft, 0.71,
                                           options.axWidth , 0.27 ] )
      plt.box( on=False )
   data.axDict = axDict
   return ( axDict )

def establishTicks( axMain, axCrazy, axBlowUp, options, data ):
   #data.axDict['main'].set_xticks( data.xData )
   #data.axDict['main'].set_xticklabels( prettyList( data.valuesDict['columnLength'] ))
   if options.mode in set( [ 'blocks', 'contigs', 'hapPaths', 'scaffPaths', 'scaffolds' ]):
      if not options.SMM:
         if options.mode != 'hapPaths' and options.mode != 'scaffPaths':
            axCrazy.set_yticks( [0, 1 ] )
            axCrazy.set_yticklabels( [ 0, '%d' % data.crazyMax ] )
   minorLocator = MultipleLocator( 5 )
   minorLocator = LogLocator( base=10, subs = range(1,10) )
   
   if options.SMM:
      axBlowUp.set_yticklabels( [] )
      axBlowUp.set_xticklabels( [] )
      axMain.set_yticklabels( [] )
      axMain.set_xticklabels( [] )
      for i in range( 1, 8 ):
         axMain.add_line( lines.Line2D( xdata=[ 10**i, 10**i ],
                                        ydata=[ 0, 0.1 ],
                                        color=( 0.8, 0.8, 0.8 ),
                                        linewidth = 0.5 )
                          )
   else:
      axBlowUp.set_yticks( [ 0.9, 0.92, 0.94, 0.96, 0.98, 1.0 ], minor=False )
      axBlowUp.set_yticks( [ 0.91, 0.92, 0.93, 0.94, 0.95,
                             0.96, 0.97, 0.98, 0.99, 1.0 ], minor=True )
      axMain.xaxis.set_minor_locator( minorLocator )
      axBlowUp.xaxis.set_minor_locator( minorLocator )
      if axCrazy:
         axCrazy.xaxis.set_minor_locator( minorLocator )

def writeImage( fig, pdf, options, data ):
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
      sys.stderr.write( 'Error, lists are not the same length, cannot add %s to %s\n' % ( str(v1), str(v2) ))
      sys.exit(1)
   r = []
   for i in range( 0, len(v1)):
      r.append( v1[i] + v2[i] )
   return r

def normalizeDataNormalMode( valuesDict, options, data ):
   if options.crazyMax == None:
      data.crazyMax = 0
      for v in valuesDict[ '!hapA1/!hapA2/assembly' ]:
         if v > data.crazyMax:
            data.crazyMax = v
   else:
      data.crazyMax = options.crazyMax
   # normalize crazy data against itself
   for i in range( 0, len( valuesDict[ '!hapA1/!hapA2/assembly' ])):
      valuesDict[ '!hapA1/!hapA2/assembly' ][i] /= float( data.crazyMax )
      
   # collect column sums, they should all be the same
   upperlimit = len( valuesDict[ '!hapA1/!hapA2/assembly' ])
   colSum = [ 0 ] * upperlimit
   for i in range( 0, upperlimit ):
      for j in ['hapA1/hapA2/assembly','hapA1/hapA2/!assembly','hapA1/!hapA2/assembly',
                'hapA1/!hapA2/!assembly','!hapA1/hapA2/assembly','!hapA1/hapA2/!assembly' ]:
         colSum[ i ] += valuesDict[ j ][ i ]
   # verify the columns all have the same sum
   for i in range(1, len( colSum )):
      if colSum[ 0 ] != colSum[ i ]:
         sys.stderr.write('Error, column sums do not equal one another, col 0 != col %d\n' % i)
         sys.exit( 1 )
   # create the OR categories where we collapse hap1/!hap2 and !hap1/hap2 into hap1ORhap2
   valuesDict[ 'hapA1ORhapA2/assembly' ] = vectorAddition( valuesDict['!hapA1/hapA2/assembly'], 
                                                           valuesDict['hapA1/!hapA2/assembly'] )
   valuesDict[ 'hapA1ORhapA2/!assembly' ] = vectorAddition( valuesDict['!hapA1/hapA2/!assembly'], 
                                                            valuesDict['hapA1/!hapA2/!assembly'] )
   # normalize the data
   for i in range( 0, upperlimit ):
      for j in ['hapA1/hapA2/assembly','hapA1/hapA2/!assembly','hapA1ORhapA2/assembly',
                'hapA1ORhapA2/!assembly' ]:
         valuesDict[ j ][ i ] /= float( colSum[ 0 ] )
   # stack the data
   options.topBotOrder.reverse()
   for i in range( 0, upperlimit ):
      cumSum = 0.0
      for j in options.topBotOrder:
         valuesDict[ j ][ i ] += cumSum
         cumSum = valuesDict[ j ][ i ]
   options.topBotOrder.reverse()
   return valuesDict

def normalizeDataContaminationMode( options, data ):
   # collect column sums, they should all be the same
   colSum = [ 0 ] * 7
   for i in range( 0, 7 ):
      for j in [ 'ecoli/!assembly', 'ecoli/assembly' ]:
         colSum[ i ] += data.valuesDict[ j ][ i ]
   # verify the columns all have the same sum
   for i in range(1, len( colSum )):
      if colSum[ 0 ] != colSum[ i ]:
         sys.stderr.write( 'Error, column sums do not equal one another, '
                           'col 0 (%d) != col %d (%d)\n' % ( colSum[0], i, colSum[i] ))
         sys.exit( 1 )
   # normalize the data
   for i in range( 0, 7 ):
      for j in ['ecoli/!assembly', 'ecoli/assembly' ]:
         data.valuesDict[ j ][ i ] /= float( colSum[ 0 ] )
   # stack the data
   options.topBotOrder.reverse()
   for i in range( 0, 7 ):
      cumSum = 0.0
      for j in options.topBotOrder:
         data.valuesDict[ j ][ i ] += cumSum
         cumSum = data.valuesDict[ j ][ i ]
   options.topBotOrder.reverse()

def drawData( axMain, axCrazy, axBlowUp, xData, yData, options, data ):
   i = -1
   # colors order is top to bottom
   if options.mode == 'contigs':
      data.colors = [ "#9467bd", "#c5b0d5", "#17becf", 
                      "#9edae5", "#ff7f0e", "#ffbb78" ]
   elif options.mode == 'hapPaths':
      data.colors = [ '#dd791f', '#fd8913',
                      '#cbdb2a', '#fff200' ]
   elif options.mode == 'blocks':
      data.colors = [ '#a89e89', '#6e5d3a',
                      '#f2aad2', '#ba759e' ]
   elif options.mode == 'scaffPaths':
      data.colors = [ '#6FB586', '#F2DC9D',
                      '#1C4169', '#72929D' ]
   elif options.mode == 'scaffolds':
      data.colors = [ '#542D54', '#EDCB23',
                      '#636991', '#ADB8FF' ]
   if options.mode in set( [ 'blocks', 'contigs', 'hapPaths', 'scaffPaths', 'scaffolds' ]):
      for n in options.topBotOrder:
         i += 1
         axMain.fill_between( x=xData,
                              y1=yData[ n ],
                              y2= [0] * len( xData ), 
                              facecolor = data.colors[ i ],
                              linewidth = 0.0 )
         axBlowUp.fill_between( x=xData,
                                y1=yData[ n ],
                                y2= [0] * len( xData ), 
                                facecolor = data.colors[ i ], 
                                linewidth = 0.0 )
      if options.mode != 'hapPaths' and options.mode != 'scaffPaths':
         # add baseline for homespun bar plot:
         axCrazy.add_line( lines.Line2D( xdata=[1, xData[-1]],
                                         ydata=[0,0],
                                         color=( 0.6, 0.6, 0.6 ),
                                         linewidth=0.5 ))
         # Error fills
         axCrazy.fill_between( x=xData, 
                               y1=yData[ '!hapA1/!hapA2/assembly' ],
                               y2= [0] * len( xData ), 
                               facecolor = 'r', 
                               linewidth = 0.0 )
      
      # create the 50 line
      if not options.SMM:
         if options.mode =='hapPaths':
            color50 = ( 0.3, 0.3, 0.3 )
         else:
            color50 = 'w'
         axMain.add_line( lines.Line2D( xdata=[ xData[ sum( numpy.array( yData[ 'hapA1/hapA2/assembly' ] ) > 0.5 ) ],
                                                xData[ sum( numpy.array( yData[ 'hapA1/hapA2/assembly' ] ) > 0.5 ) ]],
                                        ydata=[ 0, 0.5],
                                        color=color50,
                                        linewidth= 0.5))
         axMain.add_line( lines.Line2D( xdata=[ 1,
                                                xData[ sum( numpy.array( yData[ 'hapA1/hapA2/assembly' ] ) > 0.5 ) ]],
                                        ydata=[ 0.5, 0.5],
                                        color=color50,
                                        linewidth= 0.5))
                                     
   else:
      data.colors = [ "#8ca252", "#b5cf6b" ]
      for n in options.topBotOrder:
         i += 1
         axMain.fill_between( x=xData,
                              y1=yData[ n ],
                              y2=[0]* len( xData ), 
                              facecolor = data.colors[ i ],
                              linewidth = 0.0)
         axBlowUp.fill_between( x=xData,
                                y1=yData[ n ],
                                y2=[0] * len( xData ), 
                                facecolor = data.colors[ i ], 
                                linewidth = 0.0)

def drawLegend( options, data ):
   if options.SMM:
      return
   if options.mode != 'contamination':
      
      left = 0.03
      if options.mode == 'scaffPaths' or options.mode == 'hapPaths':
         height = 0.25
      else:
         height = 0.3
      bottom = 0.05
      width = 0.4
      right = left + width
      top = bottom + height
      data.axDict['main'].add_patch( patches.Rectangle( xy=(left, bottom), width = width,
                                                        height = height, color = (0.95, 0.95, 0.95),
                                                        edgecolor='None' ,
                                                        transform=data.axDict['main'].transAxes ))
      data.axDict['main'].text( x = (left + right)/2.0, y = top - 0.04, s = 'Legend', 
                                color = (0.1, 0.1, 0.1), horizontalalignment='center',
                                transform = data.axDict['main'].transAxes )
      data.axDict['main'].add_line( lines.Line2D( xdata=[left+0.07, right-0.07], ydata=[top - 0.05, top - 0.05],
                                                  color=(0.8, 0.8, 0.8), 
                                                  transform = data.axDict['main'].transAxes))
      yPos = top - 0.085
      legendText = ['hap1 and hap2, no Assembly',
                    'hap1 xor hap2, no Assembly',
                    'hap1 xor hap2, Assembly',
                    'hap1 and hap2, Assembly']
      options.topBotOrder.reverse()
      i = -1
      for n in options.topBotOrder:
         i += 1
         data.axDict['main'].add_patch( patches.Rectangle( xy=(left + 0.03, yPos - 0.01), width = 0.02 ,
                                                           height = 0.025, color = data.colors[i], 
                                                           edgecolor = '#ffffff',
                                                           transform=data.axDict['main'].transAxes ))
         data.axDict['main'].text( x = left + .08, y = yPos, s = legendText[i], 
                                   color = (0.1, 0.1, 0.1), horizontalalignment='left',
                                   verticalalignment='center', transform=data.axDict['main'].transAxes,
                                   fontsize=9.5)
         yPos -= 0.045
      if (not options.mode == 'scaffPaths') and ( not options.mode == 'hapPaths' ):
         data.axDict['main'].add_patch( patches.Rectangle( xy=(left + 0.03, yPos - 0.01), width = 0.02,
                                                           height = 0.025, color = 'r',
                                                           transform=data.axDict['main'].transAxes ))
         data.axDict['main'].text( x = left + .08, y = yPos, s = 'no hap1, no hap2, Assembly', 
                                   color = (0.1, 0.1, 0.1), horizontalalignment='left',
                                   verticalalignment='center', transform=data.axDict['main'].transAxes,
                                   fontsize = 9.5)
      options.topBotOrder.reverse()
   elif options.mode == 'contamination':
      width = 3.0
      left = 1.2
      right = 1.2 + width
      bottom = 0.05
      height = 0.17
      top = bottom + height
      data.axDict['main'].add_patch( patches.Rectangle( xy=(left, bottom), width = width,
                                                        height = height, color = (0.95, 0.95, 0.95)))
      data.axDict['main'].text( x = (left + right)/2.0, y = top - 0.04, s = 'Legend', 
                                color = (0.1, 0.1, 0.1), horizontalalignment='center' )
      data.axDict['main'].add_line( lines.Line2D( xdata=[left+0.1, right-0.1], ydata=[top - 0.05, top - 0.05],
                                                  color=(0.8, 0.8, 0.8)))
      yPos = top - 0.085
      legendText = ['Contamination, no Assembly',
                    'Contamination, Assembly' ]
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
      options.topBotOrder.reverse()

def drawAxisLabels( fig, options, data ):
   if not options.SMM:
      if options.title != None:
         fig.text(x = 0.5, y = 0.96, s = options.title,
                  fontsize = 18, horizontalalignment='center',
                  verticalalignment='bottom')
      fig.text(x = 0.5, y = 0.02, s = options.xLabel,
               fontsize = 14, horizontalalignment='center',
               verticalalignment='bottom')
      fig.text(x = options.axLeft - 0.06, y = 0.28, s = 'Stacked Proportion',
               fontsize = 14, horizontalalignment='center',
               verticalalignment='bottom',
               rotation=90 )
   else:
      data.axDict['blowUp'].text( x=0.9, y=0.8, s = options.title,
                                  fontsize = 80, horizontalalignment='right',
                                  verticalalignment = 'top', family='Helvetica',
                                  color='w',
                                  transform=data.axDict['blowUp'].transAxes )

def main():
   usage = ( '%prog --file=file.txt --mode=[scaffPaths|contigs|hapPaths|blocks|contamination] [options]\n\n'
             '%prog takes an aggregate text file ( --file ) and a mode \n'
             '( --mode ) and then produces a pretty picture.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( options, parser )
   ( fig, pdf ) = initImage( options, data )
   axDict = establishAxes( fig, options, data )
   
   data.valuesDict = readFile( options.file, options )
   data.xData = data.valuesDict['columnLength']
   
   if options.mode != 'contamination':
      data.valuesDict = normalizeDataNormalMode( data.valuesDict, options, data )
   else:
      normalizeDataContaminationMode( options, data )

   setAxisLimits( axDict['main'], axDict['crazy'], 
                  axDict['blowUp'], data.xData, 
                  options, data )
   drawData( axDict['main'], axDict['crazy'], 
             axDict['blowUp'], data.xData, data.valuesDict, options, data )
   drawLegend( options, data )
   drawAxisLabels( fig, options, data )
   
   setAxisLimits( axDict['main'], axDict['crazy'], 
                  axDict['blowUp'], data.xData,
                  options, data )

   establishTicks( axDict['main'], axDict['crazy'], 
                   axDict['blowUp'], options, data )
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
