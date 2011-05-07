#!/usr/bin/env python
"""
createSubStatsMegaPlot.py
migrated from cactusAssemblathonSubErrorPlotter.py
16 March 2011
( 10 March 2011 )
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

class Assembly:
   """ Assembly objects are generated from lines 
   in the two sub summary files, lower and upper
   """
   def __init__( self ):
      self.ID    = ''
      self.subStatsLower = {}
      self.subStatsUpper = {}
      self.allUp = -1
      self.allLo = -1

def initOptions( parser ):
   parser.add_option( '--subStatsDir', dest='subStatsDir',
                      type='string',
                      help=('Directory with subStats. Names: A1.subStats.upper.xml .'))
   parser.add_option( '--outputRanks', dest='outputRanks', default=False,
                      action='store_true', 
                      help='Outputs tab delimited rankings. default=%default')
   parser.add_option( '--raw', dest='raw', default=False,
                      action='store_true', 
                      help=('Doesn\'t normalize errors by the "Correct (bits)" '
                            'field, print raw values. default=%default'))

def checkOptions( options, parser ):
   dirs = { 'subStatsDir' : options.subStatsDir }
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % (d, dirs[ d ]) )

def readSubStatsDir( assembliesDict, options ):
   lowerStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.lower.xml') )
   upperStatsFiles = glob.glob( os.path.join( options.subStatsDir, '*.subStats.upper.xml') )
   namereg = '^([A-Z0-9]{2,3})\.subStats.*'
   namepat = re.compile( namereg )
   for l in lowerStatsFiles:
      m = re.match( namepat, os.path.basename( l ))
      if not m:
         sys.stderr.write('unable to match regex "%s" against filename "%s"' % ( namereg, l ))
         sys.exit( 1 )
      ID = m.group( 1 )
      if options.subsetFile:
         if ID not in options.assemblySubset:
            continue
      try:
         xmlTree = ET.parse( l )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( l )
      root=xmlTree.getroot()
      assembliesDict[ ID ] = Assembly()
      assembliesDict[ ID ].ID = ID
      for elm in root.attrib.keys():
         assembliesDict[ ID ].subStatsLower[ elm ] = int(float( root.attrib[ elm ]))
   for u in upperStatsFiles:
      m = re.match( namepat, os.path.basename( u ))
      if not m:
         sys.stderr.write('unable to match regex "%s" against filename "%s"' % ( namepat, u ))
         sys.exit( 1 )
      ID = m.group( 1 )
      if options.subsetFile:
         if ID not in options.assemblySubset:
            continue
      if ID not in assembliesDict:
         sys.stderr.write('unable to locate key %s in assembliesDict.\n')
         sys.exit( 1 )
      try:
         xmlTree = ET.parse( u )
      except expat.ExpatError: # broken xml file
         continue
      xmlTree = ET.parse( u )
      root=xmlTree.getroot()
      for elm in root.attrib.keys():
         assembliesDict[ ID ].subStatsUpper[ elm ] = int(float( root.attrib[ elm ]))
   return assembliesDict

def establishAxes( fig, options, data ):
   """ create one axes per chromosome
   """
   axDict = {}
   options.axLeft   = 0.1
   options.axRight  = 0.95
   options.axWidth  = options.axRight - options.axLeft
   options.margin   = 0.10
   options.axTop    = 0.96
   options.axBot    = 0.08
   options.axHeight = options.axTop - options.axBot
   axesNames = [ 'all', 'hom', 'het' ]
   numberOfAxes = len( axesNames )
   options.indHeight = float( options.axHeight - ( numberOfAxes - 1.0) * options.margin ) / numberOfAxes
   yPos = options.indHeight
   for ax in axesNames:
      axDict[ ax ] = fig.add_axes( [ options.axLeft, options.axTop - yPos ,
                                        options.axWidth , options.indHeight ] )
      yPos += options.indHeight + options.margin
   for ax in axDict:
      for loc, spine in axDict[ ax ].spines.iteritems():
         if loc in ['left', 'bottom']:
            spine.set_position(('outward',10)) # outward by 10 points
         elif loc in ['right','top']:
            spine.set_color('none') # don't draw spine               
         else:
            raise ValueError('unknown spine location: %s' % loc )
      # turn off ticks where there is no spine
      axDict[ ax ].xaxis.set_ticks_position('bottom')
      axDict[ ax ].yaxis.set_ticks_position('both')
      # if ax != 'all':
         # axDict[ ax ].set_xticks( [] )
   return axDict

def drawData( assembliesDict, sortOrder, axDict, options, data ):
   lGray = ( 0.8, 0.8, 0.8 )
   # all plot
   yMax = 0
   yMin = sys.maxint
   xNames = []
   i = 0
   for aName in sortOrder:
      i += 1
      a = assembliesDict[ aName ]
      if yMax < float( a.allUp ):
         yMax = float( a.allUp )
      if yMin > float( a.allLo ): 
         yMin = float( a.allLo )
   if options.raw:
      yMin = logLower( yMin )
   # partitions
   for i in xrange( 1, len( assembliesDict ) + 1):
      if not i % 5:
         axDict[ 'all' ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ yMin, yMax * 1.1],
                                                 linestyle='dotted',
                                                 color=lGray))
   i = 0
   for aName in sortOrder:
      a = assembliesDict[ aName ]
      i += 1
      axDict[ 'all' ].add_line( lines.Line2D( xdata=[ i, i ],
                                              ydata=[ a.allLo, a.allUp ],
                                              color='#1f77b4',
                                              linewidth= 4.0,
                                              solid_capstyle='round'))
      if options.subsetFile:
         xNames.append( lgn.idMap[ aName[0] ] )
      else:
         xNames.append( lgn.idMap[ aName[0] ]+'.'+aName[1:] )
   #if not options.normalize:
   axDict[ 'all' ].set_yscale('log')
   if yMin > yMax:
      sys.stderr.write( 'Error, yMin > yMax: %f > %f\n' % ( yMin, yMax ))
      sys.exit( 1 )
   axDict[ 'all' ].set_ylim( [ yMin * 0.9, yMax * 1.1] )
   axDict[ 'all' ].set_xlim( 0, len(xNames) + 1 )
   axDict[ 'all' ].set_xticks( range( 1, len(xNames) + 1 ))
   axDict[ 'all' ].set_xticklabels( xNames )
   for label in axDict[ 'all' ].xaxis.get_ticklabels():
      label.set_rotation( 45 )
   if not options.subsetFile:
      for tick in axDict[ 'all' ].xaxis.get_major_ticks():
         tick.label1.set_fontsize( 6 )

   # all the other plots
   axNames = { 'hom':'totalErrorsInHomozygous', 
               'het':'totalErrorsInHeterozygous'}
   #            'indel':'Total-errors-in-one-haplotype-only' }
   for key in axNames:
      yMax = 0
      yMin = sys.maxint
      i = 0
      for aName in sortOrder:
         i += 1
         a = assembliesDict[ aName ]
         if yMax < float( a.subStatsUpper[ axNames[ key ] ]):
            yMax = float( a.subStatsUpper[ axNames[ key ] ])
         if a.subStatsLower[ axNames[ key ]] > 0:
            if yMin > float( a.subStatsLower[ axNames[ key ] ]): 
               yMin = float( a.subStatsLower[ axNames[ key ] ])
         if options.raw:
            yMin = logLower( yMin )
      # partitions
      for i in xrange( 1, len( assembliesDict ) + 1):
         if not i % 5:
            axDict[ key ].add_line( lines.Line2D( xdata=[ i, i ],
                                                  ydata=[ yMin, yMax * 1.1],
                                                  linestyle='dotted',
                                                  color=lGray))
      i = 0
      for aName in sortOrder:
         a = assembliesDict[ aName ]
         i += 1
         axDict[ key ].add_line( lines.Line2D( xdata=[ i, i ],
                                                 ydata=[ a.subStatsLower[ axNames[ key ]], 
                                                         a.subStatsUpper[ axNames[ key ]]],
                                                 color='#1f77b4',
                                               linewidth=4.0,
                                               solid_capstyle='round'))
      #if not options.normalize:
      axDict[ key ].set_yscale('log')
      axDict[ key ].set_ylim( [ yMin, yMax] )
      axDict[ key ].set_xlim( 0, len(xNames) + 1 )
      axDict[ key ].set_xticks( range( 1, len(xNames) + 1 ))

      # grid
      for ax in axDict:
         mts = axDict[ax].yaxis.get_majorticklocs()
         for m in mts:
            axDict[ax].add_line( lines.Line2D( xdata=[ 1, len(xNames) ],
                                               ydata=[ m, m ],
                                               linewidth=1,
                                               color=lGray,
                                               linestyle='dotted'))

      if key == 'het':
         axDict[ key ].set_xticklabels( xNames )
      else:
         axDict[ key ].set_xticklabels( [] )
      for label in axDict[ key ].xaxis.get_ticklabels():
         label.set_rotation( 45 )
      if not options.subsetFile:
         for tick in axDict[ key ].xaxis.get_major_ticks():
            tick.label1.set_fontsize( 6 )
         

   if not options.raw:
      suffix = ' / Correct (bits)'
   else:
      suffix = ''
   axDict['all'].set_title( 'Sum of Substitution Errors%s' % suffix )
   axDict['hom'].set_title( 'Homozygous Substitution Errors%s' % suffix )
   axDict['het'].set_title( 'Heterozygous Substitution Errors%s' % suffix )
   
def logLower( y ):
   """ find the approprate lower bound for y
   if y is going to be displayed on a log plot
   """
   for i in xrange( 1, 8 ):
      if y == ( y % (10.0 ** i)):
         return ( 10.0 ** ( i - 1 ) )

def sumErrors( assembliesDict, options ):
   for a in assembliesDict:
      assembliesDict[ a ].allUp = 0
      assembliesDict[ a ].allLo = 0
      for e in [ 'totalErrorsInHomozygous', 'totalErrorsInHeterozygous']: 
         # 'Total-errors-in-one-haplotype-only' ]:
         assembliesDict[ a ].allLo += float( assembliesDict[ a ].subStatsLower[ e ] )
         assembliesDict[ a ].allUp += float( assembliesDict[ a ].subStatsUpper[ e ] )
      assembliesDict[a].allLo /= ( assembliesDict[a].subStatsLower['totalCallsInHomozygous'] + 
                                   assembliesDict[a].subStatsLower['totalCallsInHeterozygous'])
      assembliesDict[a].allUp /= ( assembliesDict[a].subStatsUpper['totalCallsInHomozygous'] + 
                                   assembliesDict[a].subStatsUpper['totalCallsInHeterozygous'])

def normalizeData( assembliesDict, options ):
   # names is a dict keyed on numerators and valued by denominators
   names = { 'totalErrorsInHomozygous':'totalCorrectInHomozygous',
             'totalErrorsInHeterozygous':'totalCorrectInHeterozygous',
             'totalErrorsInOneHaplotypeOnly':'totalCorrectInOneHaplotypeOnly' }
   for a in assembliesDict:
      for key in names:
         if key not in assembliesDict[ a ].subStatsLower:
            sys.stderr.write('Error, assembly %s lower, is missing key %s!\n' % (a, key))
            sys.exit(1)
         if key not in assembliesDict[ a ].subStatsUpper:
            sys.stderr.write('Error, assembly %s upper, is missing key %s!\n' % (a, key))
            sys.exit(1)
         assembliesDict[ a ].subStatsLower[ key ] = ( float(assembliesDict[ a ].subStatsLower[ key ]) / 
                                                      float(assembliesDict[ a ].subStatsLower[ names[ key ] ]) )
         assembliesDict[ a ].subStatsUpper[ key ] = ( float(assembliesDict[ a ].subStatsUpper[ key ]) / 
                                                      float(assembliesDict[ a ].subStatsUpper[ names[ key ] ]) )
      
def rankings( assembliesDict, sortOrder, options, data ):
   print ('#Assembly\tAll Lower\tAll Upper\tHomozygous Lower\tHomozygous Upper\t'
          'Heterozygous Lower\tHeterozygous Upper') # \tIndel Lower\tIndel Upper')
   for aName in sortOrder:
      sys.stdout.write('%s' % aName )
      a = assembliesDict[ aName ]
      for v in [ a.allLo, a.allUp ]:
         sys.stdout.write('\t%s' % v )
      for key in [ 'totalErrorsInHomozygous', 'totalErrorsInHeterozygous',
                   'totalErrorsInOneHaplotypeOnly' ]:
         sys.stdout.write('\t%s\t%s' % ( a.subStatsLower[ key ], a.subStatsUpper[ key ] ))
      sys.stdout.write('\n')
      
def main():
   usage = ( 'usage: %prog --subStatsDir=path/to/dir/ [options]\n\n'
             '%prog takes in a directory of substitution stats files ( --subStatsDir )\n'
             'with filenames as NAME.subStats.[upper|lower].xml and produces a plot.')
   data = Data()
   parser = OptionParser( usage=usage )
   initOptions( parser )
   las.initOptions( parser )
   lpt.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   las.checkOptions( options, parser )
   lpt.checkOptions( options, parser )
   
   if not options.outputRanks:
      fig, pdf = lpt.initImage( 9., 11., options, data )
      axDict = establishAxes( fig, options, data )
   
   assembliesDict = {}
   assembliesDict = readSubStatsDir( assembliesDict, options )
   
   sumErrors( assembliesDict, options )
   normalizeData( assembliesDict, options )

   sortOrder = sorted( assembliesDict, key=lambda key: assembliesDict[ key ].allLo, reverse=False )

   if options.outputRanks:
      rankings( assembliesDict, sortOrder, options, data )
      return

   drawData( assembliesDict, sortOrder, axDict, options, data )
   
   lpt.writeImage( fig, pdf, options )
   

if __name__ == '__main__':
   main()
