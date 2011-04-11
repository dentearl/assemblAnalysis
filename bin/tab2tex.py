#!/usr/bin/env python
"""
"""
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def main():
   linenumber = 0
   header = []
   print '''
\\rowcolors{1}{tableShade}{white}
\\begin{FPtable}
\caption[A table.]{A table.}
\\tiny
\\centering'''
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
            sys.stdout.write( '%s' % header[0][:9] )
            for i in xrange(1, len(header)):
               sys.stdout.write(' & %s' % header[ i ][:9].replace('_', ' '))
            sys.stdout.write( ' \\\\\n\\hline\n\\hline\n' )
            continue
      data = line.split('\t')
      sys.stdout.write( '%s' % data[0] )
      for i in xrange( 1, len( data ) ):
         sys.stdout.write(' & %s' % data[ i ])
         
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
\\end{FPtable}\par
\\normalsize
\\vspace{0.3in}'''   

if __name__ == '__main__':
   main()
