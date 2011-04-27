# libAssemblySubset.py
# dent earl, dearl (a) soe ucsc edu
# 8 April 2011
#
# Module to handle the creation of a set of assemblies
# that are to be used by the calling script.
# It is up to each script to make use of the set 
# options.assemblySubset in it's own way.
#
# One suggestion is to place code like the following in a readData()
# function:
#
# for input in inputAssemblies:
#    if options.subsetFile:
#       if input.name not in options.assemblySubset:
#          continue
#   ...
#
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
def initOptions( parser ):
   parser.add_option( '--subsetFile', dest='subsetFile',
                      type='string', help=('Subset file contains one assembly ID per line.'))

def checkOptions( options, parser ):
   import os
   options.assemblySubset = set()
   if options.subsetFile == None:
      return
   if not os.path.exists( options.subsetFile ):
      parser.error('%s does not exist!\n' % ( options.subsetFile ))
   f = open( options.subsetFile, 'r' )
   for line in f:
      line = line.strip()
      if line.startswith('#'):
         continue
      d = line.split()
      assembly = d[0]
      options.assemblySubset.add( assembly )
