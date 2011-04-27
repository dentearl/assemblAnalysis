#!/usr/bin/env python
"""
tab2tex.py
dent earl, dearl (a) soe ucsc edu
27 april 2011

python script to take a tab delimited file and format it for latex tables.

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
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '--columnTitleCharLimit', dest='n',
                      type='int', default=9,
                      help='Limits each columns title to "n" characters. default=%default' )
   parser.add_option( '--tableStyle', dest='tableStyle',
                      default='FPtable', type='string',
                      help=('This string is inserted in \\begin{ STYLE }. default=%default'))   
   parser.add_option( '--partitionEvery', dest='mod',
                      type='int', default=5,
                      help=('A horizontal dividing line will be drawn every MOD rows. default=%default'))

def checkOptions( options, parser ):
   if options.n < 0:
      parser.error( '--columnTitleCharLimit must be >= 0, %d not allowed.\n' % options.n )
   if options.tableStyle not in ('FPtable', 'table'):
      parser.error('--tableStyle must be either FPtable or table, not %s' % options.tableStyle )
   if options.mod < 0:
      parser.error( '--partitionEvery must be >= 0, %d not allowed.\n' % options.mod )

def main():
   usage = ( 'usage: %prog [options] < fasta.fa\n\n'
             '%prog takes in a tab delimited table and produces a tex formated table.\n'
             '( --columnTitleCharLimit ) will limit the length of each column\'s title.\n' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( options, parser )
   linenumber = 0
   header = []
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{%s}
\caption[A table.]{A table.}
\\tiny
\\centering''' % options.tableStyle
   hline = False
   for line in sys.stdin:
      line = line.strip()
      if line == '':
         continue
      linenumber += 1
      if line.startswith('#'):
         if linenumber == 1:
            # header
            line = line[1:]
            header = line.split('\t')
            data = line.split('\t')
            sys.stdout.write('\\begin{tabular}{ | r |')
            for i in xrange(1, len(data) ):
               sys.stdout.write( ' c |' )
            sys.stdout.write( '}\n\\hline\n' )
            sys.stdout.write( '%s' % header[0][:options.n] )
            for i in xrange(1, len(header)):
               t = header[i][:options.n]
               t = t.replace('%', '\\%')
               sys.stdout.write(' & %s' % t )
            sys.stdout.write( ' \\\\\n\\hline\n\\hline\n' )
            continue
      data = line.split('\t')
      sys.stdout.write( '%s' % data[0] )
      for i in xrange( 1, len( data )):
         d = data[ i ]
         sys.stdout.write(' & %s' % d )
         
      sys.stdout.write(' \\\\\n')
      if not (linenumber - 1) % options.mod:
         print '\\hline'
         hline = True
      else:
         hline = False
   if not hline:
      print '\\hline'
   print '''\\end{tabular}
\\label{table:aTable}
\\end{%s}\par
\\normalsize
\\vspace{0.3in}''' % options.tableStyle

if __name__ == '__main__':
   main()
