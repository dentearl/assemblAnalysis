#!/usr/bin/env python
"""
"""
from optparse import OptionParser
import os
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Team:
   def __init__(self):
      self.name  = ''
      self.ranks = []

def initOptions( parser ):
   parser.add_option( '--overall', dest='overall', default=False, action='store_true',
                      help=('Instead of ranking based on within Team, '
                            'ranks instead are global. default=%default'))

def checkOptions( args, options, parser ):
   if len( args ) < 1:
      parser.error('Error, no files in positional arguments. Input a filename.')
   for f in args:
      if not os.path.exists( f ):
         parser.error('Error, file "%s" does not exist.\n' % f )
      if not f.endswith('.tab'):
         parser.error('Error, file "%s" does not end in ".tab".\n' % f )

def readFiles( args, options ):
   assemblies = {}
   options.fileNames = []
   for aFile in args:
      options.fileNames.append( os.path.basename( aFile ) )
      f = open( aFile, 'r' )
      rank = 0
      for line in f:
         line = line.strip()
         if line.startswith('#'):
            continue
         rank += 1
         d = line.split()
         if d[0] not in assemblies:
            assemblies[ d[0] ] = []
         assemblies[ d[0] ].append( rank )
   return assemblies

def reportRank( assemblies, options ):
   if options.overall:
      ranked = sorted( assemblies, key=lambda x: sum( assemblies[x] ) )
      print options.fileNames
      for r in ranked:
         print r, sum(assemblies[r]), assemblies[r]
   else:
      teams = {}
      for a in assemblies:
         if a[0] not in teams:
            teams[ a[0] ] = []
         t = Team()
         t.name = a
         t.ranks = assemblies[ a ]
         teams[ a[0] ].append( t )
      # sort alphabetically
      aemst = teams.keys()
      aemst.sort()
      for t in aemst:
         print ''
         ranked = sorted( teams[ t ], key=lambda x: sum( x.ranks ))
         for r in ranked:
            print r.name, sum( r.ranks ), r.ranks
      

def main():
   usage = ( 'usage: %prog [options] rankingFile1.tab rankingFile2.tab rankingFile3.tab\n\n'
             '%prog takes in ranking files and reports the best overall or best within a Team.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   
   assemblies = readFiles( args, options )
   reportRank( assemblies, options )

if __name__ == '__main__':
   main()

