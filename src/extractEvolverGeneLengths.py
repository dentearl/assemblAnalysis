#!/usr/bin/env python
"""
.py
18 July 2011
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
determine information about the gene structure of
an 'annots.gff' file produced for a simulated 
genome by Evolver

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
from optparse import OptionParser
import os
import re
import sys

class Gene:
   def __init__(self):
      # left and not start here to be pos/neg strand agnostic
      self.left = sys.maxint
      self.right = -sys.maxint
      self.length = 0
      self.numCds = 0
      self.numUtr = 0
   def updateLength(self):
      self.length = self.right - self.left

def initOptions(parser):
   pass

def checkOptions(options, args, parser):
   if len(args) != 1:
      parser.error('expected gff filename as only argument')
   if not os.path.exists(args[0]):
      parser.error('filename %s does not exist' % args[0])
   if not args[0].endswith('.gff'):
      parser.error('filename %s does not end in .gff' % args[0])

def processFile(filename):
   """ reads a gff file, stores only the UTR or CDS annotations.
   evolver records the gene_index value and we use that to determine 
   the region a particular gene covers.
   """
   geneDict = {}
   pat = re.compile(r'gene_index (\d+)')
   gf = open(filename)
   for line in gf:
      line = line.strip()
      if line.startswith('#'):
         continue
      t = line.split('\t')
      if len(t) != 9 and len(t) !=8:
         continue
      if t[2] not in ['CDS', 'UTR']:
         continue
      s = t[8].split(';')
      m = re.search(pat, s[1])
      if m is None:
         raise RuntimeError('bad regex, unable to pull out gene index for line %s' % line)
      geneIndex = m.group(1)
      if geneIndex not in geneDict:
         g = Gene()
         geneDict[geneIndex] = g
      else:
         g = geneDict[geneIndex]
      if t[2] == 'CDS':
         g.numCds += 1
      elif t[2] == 'UTR':
         g.numUtr += 1
      if g.left > int(t[3]):
         g.left = int(t[3])
         g.updateLength()
      if g.right < int(t[4]):
         g.right = int(t[4])
         g.updateLength()
   return geneDict

def summarizeResults(geneDict):
   sortOrder = sorted( geneDict, key = lambda x: geneDict[x].length, reverse = False)
   print '#gene_index\tLength\tNum CDS\tNum UTR'
   for g in sortOrder:
      print g, geneDict[g].length, geneDict[g].numCds, geneDict[g].numUtr

def main():
   usage = ('usage: %prog annots.gff [options]')
   parser = OptionParser(usage = usage)
   initOptions(parser)
   options, args = parser.parse_args()
   checkOptions(options, args, parser)
   
   geneDict = processFile(args[0])
   summarizeResults(geneDict)

if __name__ == '__main__':
   main()
