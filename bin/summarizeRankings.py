#!/usr/bin/env python
"""
"""
from optparse import OptionParser
import os
import sys
import signal # deal with broken pipes

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

class Assembly:
   def __init__(self):
      self.teamName = ''
      self.name     = ''
      self.ranks    = []
      self.values   = []

class Tab:
   def __init__(self):
      self.name  = ''
      self.rank  = 0
      self.value = -1.0

def initOptions( parser ):
   parser.add_option( '--overall', dest='overall', default=False, action='store_true',
                      help=('Instead of ranking based on within Team, '
                            'ranks instead are global. default=%default'))
   parser.add_option( '--retainValues', dest='retainValues', default=False,
                      action='store_true',
                      help=('Stores the second column of the tab and outputs the value in '
                            'parenthesis following the ranking.'))

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
         a = Tab()
         d = line.split()
         a.name  = d[0]
         a.rank  = rank
         a.value = d[1]
         if a.name not in assemblies:
            assemblies[ a.name ] = []
         assemblies[ a.name ].append( a )
   return assemblies

def printHeader( options ):
   sys.stdout.write( '#Assembly\tOverall' )
   for f in options.fileNames:
      if os.path.basename( f ).endswith('.tab'):
         name = os.path.basename( f )[:-4]
      else:
         name = os.path.basename( f )
      if options.retainValues:
         name += ' (value)'
      sys.stdout.write('\t%s' % name )
   sys.stdout.write('\n')

def reportRank( assemblies, options ):
   printHeader( options )
   if options.overall:
      overallDict = {}
      for t in assemblies:
         overallDict[t] = 0
         for a in t:
            overallDict[t] += a.rank
      ranked = sorted( overallDict, key=lambda x: overallDict[x] )
      for r in ranked:
         print '%s\t%d' % ( r, overallDict[r] )
         for z in assemblies[r]:
            if options.retainValue:
               sys.stdout.write('\t%d (%s)' % (z.rank, z.value))
            else:
               sys.stdout.write('\t%d' % z.rank )
         sys.stdout.write('\n')
   else:
      teams = {}
      for a in assemblies:
         if a[0] not in teams:
            teams[ a[0] ] = []
         newAssemb = Assembly()
         newAssemb.teamName = a[0]
         newAssemb.name = a
         for z in assemblies[ a ]:
            newAssemb.ranks.append( z.rank )
            newAssemb.values.append( z.value )
         teams[ a[0] ].append( newAssemb )
      # sort alphabetically
      aemst = teams.keys()
      aemst.sort()
      for t in aemst:
         print ''
         ranked = sorted( teams[ t ], key=lambda x: sum( x.ranks ))
         for r in ranked:
            sys.stdout.write( '%s\t%d' % ( r.name, sum( r.ranks ) ))
            for j in xrange( 0, len(r.ranks)):
               if options.retainValues:
                  sys.stdout.write('\t%d (%s)' % (r.ranks[j], r.values[j] ))
               else:
                  sys.stdout.write('\t%d' % r.ranks[j] )
            sys.stdout.write('\n')

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

