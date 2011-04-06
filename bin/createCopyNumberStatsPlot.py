#!/usr/bin/env python
"""
createStatsPlot.py
1 April 2011
dent earl dearl (a) soe ucsc edu

used in the assemblathon report project to
create a plot from a single copy number
stats xml file.

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
import xml.etree.ElementTree as ET

class copyNumberCategory:
   def __init__( self ):
      self.key = ( -1, -1 )
      self.assemblyColumnCounts = {}

def initOptions( parser ):
   parser.add_option( '--title', dest='title',
                      type='string', default='Copy Number Statistics',
                      help='Title placed at the top of the plot. default=%default' )
   parser.add_option( '--out', dest='out', default='myCopyNumberStatsPlot',
                      type='string',
                      help='filename where figure will be created. No extension. default=%default' )
   parser.add_option( '--outFormat', dest='outFormat', default='pdf',
                      type='string',
                      help='output format [pdf|png|all|eps] default=%default' )
   parser.add_option( '--dpi', dest='dpi', default=300,
                      type='int',
                      help='Dots per inch of the output. default=%default' )

def checkOptions( args, options, parser ):
   if len( args ) < 1:
      parser.error('Error, please specify at least one linkage file to inspect as a positional argument.\n' )
   options.files = []
   for f in args:
      if not os.path.exists( f ):
         parser.error('Error, %s does not exist!\n' % ( f ))
      if not f.endswith('.xml'):
         parser.error('Error, file "%s" does not end in ".xml".\n' % f )
      options.files.append( os.path.abspath( f ) )
   if ( options.out.endswith('.png') or options.out.endswith('.pdf') or 
        options.out.endswith('.eps') ):
      options.out = options.out[:-4]
   if options.dpi < 72:
      parser.error('Error, I refuse to have a dpi less than screen res, 72. (%d) must be >= 72.\n' % options.dpi )

def initImage( options, data ):
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=( 11, 3.25 ), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   return ( fig, pdf )

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

def establishAxes( fig, categories, options, data ):
   axDict = {}
   data.backgroundAx = fig.add_axes( [ 0.0, 0.0, 1.0, 1.0 ] )
   data.backgroundAx.yaxis.set_major_locator( pylab.NullLocator() )
   data.backgroundAx.xaxis.set_major_locator( pylab.NullLocator() )
   plt.box( on=False )
   options.axLeft   = 0.01
   options.axRight  = 0.99
   options.width    = options.axRight - options.axLeft 
   options.axBottom = 0.1
   options.axTop    = 0.85
   options.axHeight = options.axTop - options.axBottom
   margin = 0.015
   width  = ( options.width - ( len(categories) - 1 ) * margin ) / len(categories)
   xpos = options.axLeft
   sortedOrder = categories.keys()
   sortedOrder.sort()
   options.axDictSortedOrder = sortedOrder
   for c in sortedOrder:
      axDict[ c ] = fig.add_axes( [ xpos, options.axBottom, 
                                    width, options.axHeight ] )
      axDict[ c ].yaxis.set_major_locator( pylab.NullLocator() )
      axDict[ c ].xaxis.set_major_locator( pylab.NullLocator() )
      xpos += width + margin
      plt.box( on=False )
   data.axDict = axDict
   return ( axDict )

def readFiles( options ):
   # these lists have one item per input file
   storedCategories = {}
   for f in options.files:
      xmlTree = ET.parse( f )
      root=xmlTree.getroot()
      for elm in root.findall( 'copy_number_category' ):
         t = ( int( elm.attrib[ 'minimumHaplotypeCopyNumber' ] ),
               int( elm.attrib[ 'maximumHaplotypeCopyNumber' ] ))
         if t not in storedCategories:
            storedCategories[ t ] = copyNumberCategory()
            storedCategories[ t ].key = t
         c = storedCategories[ t ]
         acn = int( elm.attrib[ 'assemblyCopyNumber' ] )
         if acn in c.assemblyColumnCounts:
            sys.stderr.write( 'Error, dupilicate assemblyCopyNumber= %d found in %s' % ( acn, f ))
            sys.exit( 1 )
         c.assemblyColumnCounts[ acn ] = int( elm.attrib[ 'columnCount' ] )
   return storedCategories

def setAxisLimits( axDict, options, data ):
   for c in axDict:
      axDict[c].set_ylim( options.globalCopyMin, options.globalCopyMax + 1.0 )
      axDict[c].set_xlim( 0.0, 1.0 )

def drawLegend( options, data ):
   pass

def prettyInt( i ):
   s = ''
   r = '%d' % i
   for j in range(0, len( r )):
      if j > 0 and not j % 3:
         s = '%s,%s' % ( r[ (len(r) - 1) - j], s )
      else:
         s = '%s%s' % ( r[ (len(r) - 1) - j], s )
   return s

def drawAxisLabels( axDict, cDict, options, data ):
   i = 0
   for c in options.axDictSortedOrder:
      i = not i
      axDict[c].set_title( c, fontsize=9 )
      axDict[c].text( x=options.facetText, 
                      y = - options.globalCopyMax / ( 12.0 + 18.0 * i ),
                      s= prettyInt( max( cDict[ c ].assemblyColumnCounts.values() )),
                      horizontalalignment='left',
                      verticalalignment='top',
                      fontsize=8 )
      
   # grey line across top of plot
   data.backgroundAx.add_line( lines.Line2D( xdata=[ options.axLeft, options.axRight ],
                                             ydata=[ options.axTop, options.axTop ],
                                             color=( 0.8, 0.8, 0.8),
                                             linewidth=0.1 ))
   # plot title
   data.backgroundAx.text( x=0.01, y = .98,
                           s=options.title,
                           horizontalalignment='left',
                           verticalalignment='top',
                           fontsize=10 )
   
   

def drawOneDataAxis( ax, key, cDict, options, data ):
   height = 0.2
   left   = 0.2 # to print numbers
   options.facetText = left
   margin = 0.05
   
   # print text:
   cMin = min( [ key[0] ] + cDict[ key ].assemblyColumnCounts.keys() )
   cMax = max( [ key[1] ] + cDict[ key ].assemblyColumnCounts.keys() )
   vMax = float( max( cDict[ key ].assemblyColumnCounts.values() ) )
   
   for i in range( cMin, cMax + 1 ):
      if key[0] <= i <= key[1]:
         color  = 'k'
         weight = 'medium'
      else:
         color = (0.5, 0.5, 0.5)
         weight = 'normal'
      ax.text( x=left, y=i+0.5, s = str(i),
               horizontalalignment='right',
               verticalalignment='center',
               fontsize=9, color=color, weight=weight )
      if i in cDict[ key ].assemblyColumnCounts:
         if key[0] <= i <= key[1]:
            color = ( 0.3, 0.3, 0.3 )
         else:
            color = 'r'
         ax.add_patch( patches.Rectangle( xy=( left + margin, i + 0.5 - height/2.0 ),
                                          width = cDict[ key ].assemblyColumnCounts[i] / vMax,
                                          height = height,
                                          color = color,
                                          edgecolor=None,
                                          linewidth=0.0))
   

def drawData( axDict, cDict, options, data ):
   for c in cDict:
      drawOneDataAxis( axDict[ c ], c, cDict , options, data )

def establishGlobalMinMax( cDict, options, data ):
   options.globalCopyMin =  sys.maxint
   options.globalCopyMax = -sys.maxint
   for c in cDict:
      if options.globalCopyMin > min( [ c[0] ] + cDict[ c ].assemblyColumnCounts.keys() ):
         options.globalCopyMin = min( [ c[0] ] + cDict[ c ].assemblyColumnCounts.keys() )
      if options.globalCopyMax < max( [ c[1] ] + cDict[ c ].assemblyColumnCounts.keys() ):
         options.globalCopyMax = max( [ c[1] ] + cDict[ c ].assemblyColumnCounts.keys() )

def main():
   usage = ( 'usage: %prog [options] file1.xml\n\n'
             '%prog takes in a copy number statistics file\n'
             'and creates an image file.' )
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   ( fig, pdf ) = initImage( options, data )

   storedCategories = readFiles( options )
   
   establishGlobalMinMax( storedCategories, options, data )
   axDict = establishAxes( fig, storedCategories, options, data )
   
   drawData( axDict, storedCategories, options, data )
   drawLegend( options, data )
   drawAxisLabels( axDict, storedCategories, options, data )
   setAxisLimits( axDict, options, data )
   
   writeImage( fig, pdf, options, data )

if __name__ == '__main__':
   main()
