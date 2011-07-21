#!/usr/bin/env python
"""
createFeatureOverlapTable.py
14 July 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the Feature (+Gene) Overlap table

output is tab delimited

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
from optparse import OptionParser
import os
import re
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

class Assembly:
   """ used to store information for a particular assembly.
   """
   def __init__( self ):
      self.id          = ''
      # these refer to the four xml files,
      # will be keyed based on the featureFile attribute of the <intervals> tag,
      self.valuesDict = {}

class IntervalsTag:
   """ used to store the results of an intervals tag from an xml file
   """
   def __init__( self ):
      self.file = ''
      self.annot = ''
      self.haplotype = 0
      self.complete = 0 # features all correct
      self.samples = 0
      self.baseLength = 0
      self.totalComplete = 0 # base length of complete

def initOptions( parser ):
   parser.add_option( '--statsScaffoldsFeatureOverlapDir', dest = 'statsScaffoldsFeatureOverlapDir',
                      type='string',
                      help=('Directory with feature overlap xmls from Scaffolds alignment. '
                            'Names: A1.scaffoldsOverlap.xml, etc .'))
   parser.add_option('--noContigs', dest = 'noContigs',
                     default=False, action='store_true',
                     help=('turns off the contig values'))
   parser.add_option('--scaffolds', dest = 'scaffolds',
                     default=False, action='store_true',
                     help=('turns on the scaffold values'))
   parser.add_option('--aveHaps', dest = 'aveHaps',
                     default=False, action='store_true',
                     help=('average the haplotype values'))
   parser.add_option('--bases', dest = 'bases',
                     default=False, action='store_true',
                     help=('show base numbers instead of feature counts'))
   parser.add_option( '--hideAssemblyNumbers', dest='hideAssemblyNumbers', default=False,
                      action='store_true',
                      help=('Hides the intra-team assembly number next to the name. default=%default'))

def checkOptions( args, options, parser ):
   if len( args ):
      parser.error('unexpected arguments: %s' % ' '.join(args))
   dirs = { 'statsScaffoldsFeatureOverlapDir' : options.statsScaffoldsFeatureOverlapDir}
   for d in dirs:
      if not dirs[ d ]:
         parser.error('specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('--%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('--%s %s is not a directory!\n' % ( d, dirs[ d ]) )
   if options.subsetFile:
      options.hideAssemblyNumbers = True

def processDirectory( options ):
   filenames = glob.glob(os.path.join(options.statsScaffoldsFeatureOverlapDir, '*.xml'))
   nameRegex = r'^(\w\d+)\.(.+)\.xml$'
   assemblyPat = re.compile( nameRegex )
   assemblyDict = {}
   for f in filenames:
      m = re.match(assemblyPat, os.path.basename( f ))
      assemblyName = m.group(1)
      if 'subsetFile' in vars( options ):
         if options.subsetFile:
            if assemblyName not in options.assemblySubset:
               continue
      filetype = m.group(2)
      # if m.group(2) != 'contigsOverlapGene':
      #    continue
      if m is None:
         raise RuntimeError( 'bad regex: %s for filename %s' % ( nameRegex, os.path.basename( f ) ))
      if assemblyName not in assemblyDict:
         a = Assembly()
         a.id = assemblyName
         assemblyDict[ a.id ] = a
      addData( assemblyDict[ assemblyName ], filetype, f )
   return assemblyDict

def addData( assembly, filetype, filename ):
   if filetype not in assembly.valuesDict:
      assembly.valuesDict[ filetype ] = {}
   try:
      xmlTree = ET.parse( filename )
   except expat.ExpatError: # broken xml file
      return
   xmlTree = ET.parse( filename )
   root = xmlTree.getroot()
   for elm in root.findall( 'intervals' ):
      annot = os.path.basename( elm.attrib['featureFile'] ).split('.')[1].lower()
      haplotype = int(os.path.split( os.path.dirname (elm.attrib['featureFile']))[1][-1])
      if annot == 'nxe' or annot == 'nge':
         annot = 'nxe+nge'
         if annot in assembly.valuesDict[ filetype ]:
            if haplotype in assembly.valuesDict[ filetype ][ annot ]:
               it = assembly.valuesDict[ filetype ][ annot ][ haplotype ]
            else:
               it = IntervalsTag()
         else:
            it = IntervalsTag()
      else:
         it = IntervalsTag()
      it.file = elm.attrib['featureFile']
      it.annot = annot
      it.haplotype = haplotype
      it.complete += int( elm.attrib['complete'] )
      it.samples += int( elm.attrib['samples'] )
      if 'baseLength' not in elm.attrib:
         raise RuntimeError('file %s %s lacks attrib baseLength' % (filename, it.file))
      it.baseLength += int( elm.attrib['baseLength'])
      it.totalComplete += int( elm.attrib['totalComplete'])
      
      if it.annot not in assembly.valuesDict[ filetype ]:
         assembly.valuesDict[ filetype ][ it.annot ] = { it.haplotype : it }
         continue
      if it.haplotype in assembly.valuesDict[ filetype ][ it.annot ] and annot != 'nxe+nge':
         raise RuntimeError('there should be only one entry of '
                            'haplotype %d for annot %s there are at '
                            'least two in file %s' % ( it.haplotype, it.annot, filename))
      assembly.valuesDict[ filetype ][ it.annot ][ it.haplotype ] = it

def outputDict( assembliesDict, options ):
   sortOrder = sorted( assembliesDict, key = lambda x: int(x[1:]) ) # number
   sortOrder = sorted( sortOrder, key = lambda x: assembliesDict[x].id[0] ) # letter
   fileOrder = []
   if not options.noContigs:
      fileOrder += ['contigsOverlapGene', 'contigsOverlap' ]
   # if options.scaffolds:
   #    fileOrder += ['scaffoldsOverlapGene', 'scaffoldsOverlap']

   annotOrderGenes = ['transcripts'] # 'cds'
   annotOrder = ['cds', 'utr', 'nxe+nge', 'repeat'] # 'island'
   fileNameMap = {'contigsOverlapGene':'COG', 'contigsOverlap':'CO'}
   annotNameMap = {'transcripts':'xcript', 'cds':'cds', 'utr':'utr',
                   'nxe+nge':'nxe+nge', 'repeat':'repeat'}
   
   # print the header
   sys.stdout.write('#ID')
   i = -1
   for f in fileOrder:
      i += 1
      if f.endswith('OverlapGene'):
         order = annotOrderGenes
      else:
         order = annotOrder
      for a in order:
         if options.bases:
            if options.aveHaps:
               totalStr = lgn.prettyNumber( assembliesDict[sortOrder[0]].valuesDict[f][a][1].baseLength +
                                            assembliesDict[sortOrder[0]].valuesDict[f][a][2].baseLength )
            else:
               totalStr = '%s, %s' % ( lgn.prettyNumber(assembliesDict[sortOrder[0]].valuesDict[f][a][1].baseLength),
                                       lgn.prettyNumber(assembliesDict[sortOrder[0]].valuesDict[f][a][2].baseLength ))
         else:
            if options.aveHaps:
               totalStr = lgn.prettyNumber( assembliesDict[sortOrder[0]].valuesDict[f][a][1].samples +
                                            assembliesDict[sortOrder[0]].valuesDict[f][a][2].samples )
            else:
               totalStr = '%s, %s' % ( lgn.prettyNumber(assembliesDict[sortOrder[0]].valuesDict[f][a][1].samples),
                                       lgn.prettyNumber(assembliesDict[sortOrder[0]].valuesDict[f][a][2].samples) )
         sys.stdout.write('\t%s-%s (%s)' % (fileNameMap[f], annotNameMap[a], totalStr))
   sys.stdout.write('\n')
   
   for assembly in sortOrder:
      if options.hideAssemblyNumbers:
         nameStr = '%s' % (lgn.idMap[assembliesDict[assembly].id[0]])
      else:
         nameStr = '%s.%s' % (lgn.idMap[assembliesDict[assembly].id[0]], assembliesDict[assembly].id[1:])
      sys.stdout.write('%s' % nameStr)
      for filetype in fileOrder:
         if filetype.endswith('OverlapGene'):
            order = annotOrderGenes
         else:
            order = annotOrder
         for annot in order:
            aveHap = IntervalsTag()
            for hap in assembliesDict[assembly].valuesDict[filetype][annot]:
               it = assembliesDict[assembly].valuesDict[filetype][annot][hap]
               if options.bases:
                  fracStr = '%.2f' % ( float( it.totalComplete ) / it.baseLength )
                  numer = it.complete
                  denom = it.samples
               else:
                  fracStr = '%.2f' % ( float( it.complete ) / it.samples )
                  numer = it.complete
                  denom = it.samples
               if fracStr == '1.00' and ( numer != denom ):
                  fracStr = '0.99'
               if not options.aveHaps:
                  sys.stdout.write('\thapA%d %s' 
                                   % (hap, fracStr))
               else:
                  aveHap.complete += it.complete
                  aveHap.samples += it.samples
                  aveHap.totalComplete += it.totalComplete
                  aveHap.baseLength += it.baseLength
            if options.aveHaps:
               if options.bases:
                  fracStr = '%.2f' % ( float( aveHap.totalComplete ) / aveHap.baseLength )
                  numer = aveHap.totalComplete
                  denom = aveHap.baseLength
               else:
                  fracStr = '%.2f' % ( float( aveHap.complete ) / aveHap.samples )
                  numer = aveHap.complete
                  denom = aveHap.samples
               if fracStr == '1.00' and ( numer != denom ):
                  fracStr = '0.99'
               sys.stdout.write('\t%s' % fracStr)
      sys.stdout.write('\n')

def main():
   usage = ('usage: %prog')
   parser = OptionParser( usage = usage )
   initOptions( parser )
   las.initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )
   las.checkOptions( options, parser )
   
   assembliesDict = processDirectory(options)
   
   outputDict(assembliesDict, options)

if __name__ == '__main__':
   main()
