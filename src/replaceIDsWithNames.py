#!/usr/bin/env python
"""
replaceIDsWithNames.py
1 May 2011
dent earl, dearl (a) soe ucsc edu

Script takes text on STDIN, converts the assembly ID code into
the assembly name, as defined in libGeneral.py

"""
import libGeneral as lgn
from optparse import OptionParser
import re
import sys

def initOptions( parser ):
   parser.add_option( '--retainNumber', dest='retainNumber',
                      action='store_true', default=False,
                      help='Will print the assembly number after the name. default=%default')

def checkOptions( options, parser ):
   pass

def processStream( options ):
   for line in sys.stdin:
      line = line.strip()
      if line.startswith('#'):
         # respect comments
         print line
         continue
      for i in lgn.idMap:
         if options.retainNumber:
            line = re.sub( r'%s(\d+)' % i , r'%s \1' % lgn.idMap[i], line )
         else:
            line = re.sub( r'%s(\d+)' % i , r'%s' % lgn.idMap[i], line )
      print line

def main():
   usage = ('usage: %prog < coverageTotal.tab.tmp > coverageTotal.tab.tmptmp [options]\n\n'
            '%prog is a simple script to find and replace instances of assembly IDs with\n'
            'assembly names, as defined in the dict idMap contained in libGeneral.py.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   processStream( options )

if __name__ == '__main__':
   main()
