#!/usr/bin/env python
"""
removeTextFromTags.py
14 July 2011
dent earl dearl(a) soe ucsc edu

given an xml file, this walks through an xml and sets all of the .text
values to ''.

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
import sys
import xml.etree.ElementTree as ET
import xml.parsers.expat as expat # exception handling for empty xml

def initOptions( parser ):
   pass

def checkOptions( options, args, parser ):
   if len(args) > 1:
      parser.error('expected just one argument, the filename, saw: %s' % ' '.join(args))
   if len(args) != 1:
      parser.error('expected the filename as an argument')
   if not os.path.exists(args[0]):
      parser.error('filename %s does not exist' % args[0])
   if not args[0].endswith('.xml'):
      parser.error('filename %s does not end in ".xml"' % args[0])

def processFile( filename, options ):
   try:
      xmlTree = ET.parse( filename )
   except expat.ExpatError: # broken xml file
      raise RuntimeError('bad xml: %s' % filename)
   xmlTree = ET.parse( filename )
   root = xmlTree.getroot()
   
   killText( root )

   xmlTree.write( filename )

def killText( node ):
   node.text = None
   for n in node:
      killText(n)

def main():
   parser = OptionParser()
   initOptions(parser)
   options, args = parser.parse_args()
   checkOptions(options, args, parser)
   
   processFile(args[0], options)

if __name__ == '__main__':
   main()
