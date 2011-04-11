#!/usr/bin/env python
"""
"""
import libAssemblySubset as las
from optparse import OptionParser
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def main():
   usage = ('usage: %prog --subsetFile=file [options] < table.tab\n\n')
   parser = OptionParser( usage=usage )
   las.initOptions( parser )
   ( options, args ) = parser.parse_args()
   las.checkOptions( options, parser )

   for line in sys.stdin:
      line = line.strip()
      if line == '':
         continue
      if line.startswith('#'):
         print line
         continue
      if line.split('\t')[0] not in options.assemblySubset:
         continue
      print line

if __name__ == '__main__':
   main()
