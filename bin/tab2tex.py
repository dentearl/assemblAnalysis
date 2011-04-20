#!/usr/bin/env python
"""
"""
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

def checkOptions( options, parser ):
   if options.n < 0:
      parser.error( '--columnTitleCharLimit must be >= 0, %d not allowed.\n' % options.n )
   if options.tableStyle not in ('FPtable', 'table'):
      parser.error('--tableStyle must be either FPtable or table, not %s' % options.tableStyle )

def main():
   usage = ( 'usage: %prog [options] < fasta.fa\n\n'
             '%prog takes in a tab delimited table and produces a tex formated table.\n'
             '( --columnTitleCharLimit ) will limit the length of each column\'s title.\n' )
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
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
               t = header[i][:options.n].replace('_', ' ')
               t = t.replace('%', '\\%')
               sys.stdout.write(' & %s' % t )
            sys.stdout.write( ' \\\\\n\\hline\n\\hline\n' )
            continue
      data = line.split('\t')
      sys.stdout.write( '%s' % data[0].replace('_', ' ') )
      for i in xrange( 1, len( data ) ):
         d = data[ i ].replace('_', ' ')
         sys.stdout.write(' & %s' % d )
         
      sys.stdout.write(' \\\\\n')
      if not (linenumber - 1) % 10:
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
