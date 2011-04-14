#!/usr/bin/env python
"""
createIndividualSection.py
dent earl dearl(a) soe ucsc edu

used in the assemblathon report project to 
create the Individual Section part of the report.

output is latex

Creates a subsection for each assembly, shows info
about the assembly, creates links to images for the
assembly.

generally, this script is pretty cool.
"""
from optparse import OptionParser
import os
import sys

class Team:
   """ Team objects are generated from lines in the infoFile
   """
   def __init__( self ):
      self.ID           = 'empty'
      self.name         = 'empty'
      self.affiliations = 'empty'
      self.contact      = 'empty'
      self.numEntries   = -1
      self.software     = 'empty'
      self.entries      = []

class Assembly:
   """ Assembly objects are generated from lines in the rankFile
   """
   def __init__( self ):
      self.team  = None
      self.ID    = ''
      self.rank  = -1
      self.total = -1
      self.hap1  = -1
      self.hap2  = -1
      self.bac   = -1
      self.subStatsLower = {}
      self.subStatsUpper = {}
      self.sizeStatsScaffold = []
      self.sizeStatsContigs  = []

def initOptions( parser ):
   parser.add_option( '--rankFile', dest='rankFile',
                      type='string',
                      help=('File with columns (space delim) for assemblyID, total coverage, '
                            'hap1 coverage, hap2 coverage, delta which is abs(hap1 - hap2) and '
                            'bacterial coverage.' ))
   parser.add_option( '--infoFile', dest='infoFile',
                      type='string',
                      help=('File with columns (tab delim) for teamID, team name, '
                            'affiliations, contact, number of entries, software used.' ))
   parser.add_option( '--subStatsDir', dest='subStatsDir',
                      type='string',
                      help=('Directory with subStats. Names: A1.subStats.upper.txt .'))
   parser.add_option( '--submissionStatsDir', dest='submissionStatsDir',
                      type='string',
                      help=('Directory with submission stats. Names: A1.summary.txt .'))
   parser.add_option( '--submissionLengthsDir', dest='submissionLengthsDir',
                      type='string',
                      help=('Directory with submission lengths. Names: A1.contigs.txt or A1.scaffolds.txt'))
   parser.add_option( '--imagesDir', dest='imagesDir',
                      type='string',
                      help=('Directory where .eps images are already stored.'))
   parser.add_option( '--placeHolders', dest='placeHolders',
                      action='store_true', default=False,
                      help=('Creates frame boxes in place of missing images. default=%default'))

def checkOptions( args, options, parser ):
   if not options.rankFile:
      parser.error('Error, specify --rankFile\n')
   if not os.path.exists( options.rankFile ):
      parser.error('Error, --rankFile %s does not exist!\n' % options.rankFile )
   if not options.infoFile:
      parser.error('Error, specify --infoFile\n')
   if not os.path.exists( options.infoFile ):
      parser.error('Error, --infoFile %s does not exist!\n' % options.infoFile )
   dirs = { 'imagesDir'   : options.imagesDir,
            'subStatsDir' : options.subStatsDir,
            'submissionStatsDir' : options.submissionStatsDir,
            'submissionLengthsDir' : options.submissionLengthsDir}
   for d in dirs:
      if not dirs[ d ]:
         parser.error('Error, specify --%s\n' % d )
      if not os.path.exists( dirs[ d ] ):
         parser.error('Error, --%s %s does not exist!\n' % ( d, dirs[ d ] ))
      if not os.path.isdir( dirs[ d ] ):
         parser.error('Error, --%s %s is not a directory!\n' % (d, dirs[ d ]) )

def readInfoFile( options ):
   teamsDict = {}
   file = open( options.infoFile, 'r' )
   for line in file:
      line = line.strip()
      if line[0] == '#':
         continue
      data = line.split('\t')
      t = Team()
      t.ID = data[ 0 ]
      t.name = data[ 1 ]
      t.affiliations = data[ 2 ]
      t.contact = data[ 3 ]
      t.numEntries = int( data[ 4 ] )
      t.software = data[ 5 ]
      if t.name == '':
         t.name = t.contact
      teamsDict[ t.ID ] = t
   return teamsDict

def readRankFile( teamsDict, options ):
   assembliesList = []
   file = open( options.rankFile, 'r' )
   r = 0
   for line in file:
      line = line.strip()
      if line[0] == '#':
         continue
      r += 1
      data = line.split(' ')
      a = Assembly()
      a.ID = data[ 0 ]
      a.rank = r
      a.total = float( data[ 1 ] )
      a.hap1  = float( data[ 2 ] )
      a.hap2  = float( data[ 3 ] )
      a.bac   = float( data[ 5 ] )
      if a.ID[ 0 ] in teamsDict:
         a.team = teamsDict[ a.ID[ 0 ] ]
         teamsDict[ a.ID[ 0 ] ].entries.append( a )
      assembliesList.append( a )
   return assembliesList

def printRankRow( a, gray=False):
   if gray:
      print ('  \\textcolor{myGray40}{%s} '
             '& \\textcolor{myGray40}{%.5f} & \\textcolor{myGray40}{%.5f} '
             '& \\textcolor{myGray40}{%.5f} & \\textcolor{myGray40}{%.5f} \\\\' % ( a.ID, 
                                                                                    a.total, a.hap1, 
                                                                                    a.hap2, a.bac ))
   else:
      print '  %s & %.5f & %.5f & %.5f & %.5f \\\\' % ( a.ID, a.total, a.hap1, a.hap2, a.bac )

def prettyInt( i ):
   s = ''
   r = '%d' % i
   for j in range(0, len( r )):
      if j > 0 and not j % 3:
         s = '%s,%s' % ( r[ (len(r) - 1) - j], s )
      else:
         s = '%s%s' % ( r[ (len(r) - 1) - j], s )
   return s

def prettyFloat( f, n ):
   s = ''
   r = '%.*f' % ( n, f )
   for j in range(0, len( r ) ):
      if not (j - n - 1) % 3 and j > (n + 1):
         s = '%s,%s' % ( r[ (len(r) - 1) - j], s )
      else:
         s = '%s%s' % ( r[ (len(r) - 1) - j], s )
   return s

def printAssembly( a, assembliesList, captionsDict, options ):
   print '\\subsubsection{%s}' % a.ID
   #print '\\noindent Coverage: %d\par' % a.rank

   # print neighbor table
   print '\\noindent Coverage neighbor table:\par'
   print '\\vspace{0.25in}'
   print '\\rowcolors{1}{tableShade}{white}'
   print '\\begin{tabular}{ | r | c | c | c | c | }'
   print '\\hline'
   print 'ID & Total & Hap 1 & Hap 2 & Bac \\\\'
   print '\\hline \\hline'
   if  1 < a.rank < len( assembliesList ):
      printRankRow( assembliesList[ a.rank - 2 ], gray=True)
      printRankRow( assembliesList[ a.rank - 1 ])
      printRankRow( assembliesList[ a.rank ], gray=True)
   elif a.rank == 1:
      printRankRow( assembliesList[ a.rank - 1 ])
      printRankRow( assembliesList[ a.rank  ], gray=True)
      printRankRow( assembliesList[ a.rank + 1], gray=True)
   elif a.rank == len( assembliesList ):
      printRankRow( assembliesList[ a.rank - 3 ], gray=True)
      printRankRow( assembliesList[ a.rank - 2 ], gray=True)
      printRankRow( assembliesList[ a.rank - 1 ])
   print '\\hline'
   print '\\end{tabular}\par'
   print '\\vspace{0.3in}'

   # print n50 thing
   print '\\noindent Submitted assembly N stats plot\par'
   print '\\vspace{0.25in}'
   print '\\begin{center}'
   print '\\epsfig{file=images/n50.%s.eps, width=3.5in}' % a.ID
   print '\\end{center}'
   print '\\vspace{0.3in}'

   # print submission size stats table
   print '\\noindent Submitted assembly size stats table\par'
   print '\\vspace{0.25in}'
   print '\\rowcolors{1}{tableShade}{white}'
   print '\\tiny'
   print '\\begin{tabular}{ | r | c | c | c | c | c | c | c | c | c | }'
   print '\\hline'
   print 'Category & n & min & 1st Qu. & Median & Mean & 3rd Qu. & Max. & Stdev & Sum\\\\'
   print '\\hline \\hline'
   print ( 'Scaffolds & %s & %s & %s '
           '& %s & %s & %s & %s & %s '
           '& % s \\\\' % ( prettyInt( int( float(a.sizeStatsScaffold[ 0 ]))),
                            prettyInt(int( float(a.sizeStatsScaffold[ 1 ]))),
                            prettyFloat(float( a.sizeStatsScaffold[ 2 ]),2),
                            prettyInt(int( float(a.sizeStatsScaffold[ 3 ]))),
                            prettyFloat(float( a.sizeStatsScaffold[ 4 ]),2),
                            prettyFloat(float( a.sizeStatsScaffold[ 5 ]),2),
                            prettyInt(int( float(a.sizeStatsScaffold[ 6 ]))),
                            prettyFloat(float( a.sizeStatsScaffold[ 7 ]),2),
                            prettyInt(int( float(a.sizeStatsScaffold[ 8 ])))
                            ))
   
   print ( 'Contigs & %s & %s & %s '
           '& %s & %s & %s & %s & %s '
           '& %s \\\\' % ( prettyInt(int( float(a.sizeStatsContigs[ 0 ]))),
                           prettyInt(int( float(a.sizeStatsContigs[ 1 ]))),
                           prettyFloat(float( a.sizeStatsContigs[ 2 ]),2),
                           prettyInt(int( float(a.sizeStatsContigs[ 3 ]))),
                           prettyFloat(float( a.sizeStatsContigs[ 4 ]),2),
                           prettyFloat(float( a.sizeStatsContigs[ 5 ]),2),
                           prettyInt(int( float(a.sizeStatsContigs[ 6 ]))),
                           prettyFloat(float( a.sizeStatsContigs[ 7 ]),2),
                           prettyInt(int( float(a.sizeStatsContigs[ 8 ]))),
                           ))
   print '\\hline'
   print '\\end{tabular}\par'
   print '\\normalsize'
   print '\\vspace{0.3in}'
   
   # print subStats table
   print '\\noindent SNP stats table\par'
   print '\\vspace{0.25in}'
   print '\\rowcolors{1}{tableShade}{white}'
   print '\\tiny'
   print '\\begin{tabular}{ | r | c | c | c | c | }'
   print '\\hline'
   print 'Category & Total & Calls & Correct (bits) & Errors \\\\'
   print '\\hline \\hline'
   if len( a.subStatsLower ) > 0 and len( a.subStatsUpper ) > 0:
      lower = a.subStatsLower
      upper = a.subStatsUpper
      print 'Homozygous & %s -- %s & %s -- %s & %s -- %s & %s -- %s \\\\' % ( prettyInt(int( lower[ 'Total-homozygous' ])), prettyInt(int( upper[ 'Total-homozygous' ])),
                                                                                  prettyInt(int( lower[ 'Total-calls-in-homozygous' ])), prettyInt(int( upper[ 'Total-calls-in-homozygous' ])),
                                                                                  prettyFloat(float( lower[ 'Total-correct-in-homozygous' ]),1), prettyFloat(float( upper[ 'Total-correct-in-homozygous' ]),1),
                                                                                  prettyInt(int( lower[ 'Total-errors-in-homozygous' ])), prettyInt(int( upper[ 'Total-errors-in-homozygous' ])))
      print 'Heterozygous & %s -- %s & %s -- %s & %s -- %s & %s -- %s \\\\' % ( prettyInt(int(lower[ 'Total-heterozygous' ])), prettyInt(int(upper[ 'Total-heterozygous' ])),
                                                                                    prettyInt(int(lower[ 'Total-calls-in-heterozygous' ])), prettyInt(int(upper[ 'Total-calls-in-heterozygous' ])),
                                                                                    prettyFloat(float(lower[ 'Total-correct-in-heterozygous' ]),1), prettyFloat(float(upper[ 'Total-correct-in-heterozygous' ]),1),
                                                                                    prettyInt(int(lower[ 'Total-errors-in-heterozygous' ])), prettyInt(int(upper[ 'Total-errors-in-heterozygous' ])))
      print 'Indel & %s -- %s & %s -- %s & %s -- %s & %s -- %s \\\\' % ( prettyInt(int(lower[ 'Total-in-one-haplotype-only' ])), prettyInt(int(upper[ 'Total-in-one-haplotype-only' ])),
                                                                           prettyInt(int(lower[ 'Total-calls-in-one-haplotype-only' ])), prettyInt(int(upper[ 'Total-calls-in-one-haplotype-only' ])),
                                                                           prettyFloat(float(lower[ 'Total-correct-in-one-haplotype-only' ]),1), prettyFloat(float(upper[ 'Total-correct-in-one-haplotype-only' ]),1),
                                                                           prettyInt(int(lower[ 'Total-errors-in-one-haplotype-only' ])), prettyInt(int(upper[ 'Total-errors-in-one-haplotype-only' ])))
      
   print '\\hline'
   print '\\end{tabular}\par'
   print '\\normalsize'
   print '\\vspace{0.3in}'
   
   # print aggregate plots   
   if os.path.exists( os.path.join( options.imagesDir , a.ID+ '.contigs.eps') ):
      if options.placeHolders:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\fbox{\\begin{minipage}{6in} '
         print '\\hspace{1in}'
         print '\\vspace{6in}'
         print '\\end{minipage}}'
         print '\\caption[%s contig length cumulative plot.]{%s contig length cumulative plot. %s}' % ( a.ID, a.ID, captionsDict['contigs'] )
         print '\\label{fig:%sContigs}' % a.ID
         print '\\end{figure}'
      else:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\epsfig{file=images/%s.contigs.eps, height=6in}' % a.ID
         print '\\caption[%s contig length cumulative plot.]{%s contig length cumulative plot. %s}' % ( a.ID, a.ID, captionsDict['contigs'] )
         print '\\label{fig:%sContigs}' % a.ID
         print '\\end{figure}'
         print '\\clearpage'

   if os.path.exists( os.path.join( options.imagesDir , a.ID+ '.scaffPaths.eps') ):
      if options.placeHolders:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\fbox{\\begin{minipage}{6in} '
         print '\\hspace{1in}'
         print '\\vspace{6in}'
         print '\\end{minipage}}'
         print '\\caption[%s scaffold path length cumulative plot.]{%s scaffold path length cumulative plot. %s}' % ( a.ID, a.ID, captionsDict['scaffolds'] )
         print '\\label{fig:%sScaffolds}' % a.ID
         print '\\end{figure}'
      else:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\epsfig{file=images/%s.scaffPaths.eps, height=6in}' % a.ID
         print '\\caption[%s scaffold path length cumulative plot.]{%s scaffold path length cumulative plot. %s}' % ( a.ID, a.ID, captionsDict['scaffolds'] )
         print '\\label{fig:%sScaffolds}' % a.ID
         print '\\end{figure}'
         print '\\clearpage'
      
   if os.path.exists( os.path.join( options.imagesDir , a.ID+ '.contigPaths.eps') ):
      if options.placeHolders:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\fbox{\\begin{minipage}{6in} '
         print '\\hspace{1in}'
         print '\\vspace{6in}'
         print '\\end{minipage}}'
         print '\\caption[%s contig path cumulative length plot.]{%s contig path cumulative length plot. %s}' % ( a.ID, a.ID, captionsDict['contigPaths'] )
         print '\\label{fig:%sContigPaths}' % a.ID
         print '\\end{figure}'
      else:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\epsfig{file=images/%s.contigPaths.eps, height=6in}' % a.ID
         print '\\caption[%s contig path cumulative length plot.]{%s contig path cumulative length plot. %s}' % ( a.ID, a.ID, captionsDict['contigPaths'] )
         print '\\label{fig:%sContigPaths}' % a.ID
         print '\\end{figure}'
         print '\\clearpage'

   if os.path.exists( os.path.join( options.imagesDir , a.ID+ '.blocks.eps') ):
      if options.placeHolders:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\fbox{\\begin{minipage}{6in} '
         print '\\hspace{1in}'
         print '\\vspace{6in}'
         print '\\end{minipage}}'
         print '\\caption[%s block cumulative length plot.]{%s block cumulative length plot. %s}' % ( a.ID, a.ID, captionsDict['blocks'] )
         print '\\label{fig:%sBlocks}' % a.ID
         print '\\end{figure}'
      else:
         print '\\begin{figure}[htc]'
         print '\\centering'
         print '\\epsfig{file=images/%s.blocks.eps, height=6in}' % a.ID
         print '\\caption[%s block cumulative length plot.]{%s block cumulative length plot. %s}' % ( a.ID, a.ID, captionsDict['blocks'] )
         print '\\label{fig:%sBlocks}' % a.ID
         print '\\end{figure}'
         print '\\clearpage'

def extractRanksString( team ):
   s = ''
   i = 0
   orderedByRank = sorted( team.entries, key=lambda x: x.rank, reverse=False )
   for a in orderedByRank:
      i += 1
      if i == len( team.entries ):
         spacer = ''
      else:
         spacer = ', '
      s += '%d%s' % ( a.rank, spacer )
   return s
   
def printTeam( team, assembliesList, captionsDict, options ):
   print '\\subsection{%s, %s}' % ( team.ID, team.name )
   print '\\noindent Affiliation: %s\\par' % team.affiliations
   print '\\noindent Contact: %s\\par' % team.contact
   print '\\noindent Software: \\textbf{%s}\\par' % team.software
   print '\\noindent Number of entries: %d\\par' % team.numEntries
   print '\\vspace{0.25in}'
   print '\\rowcolors{1}{tableShade}{white}'
   print '\\begin{tabular}{ | r | c | c | c | c | }'
   print '\\hline'
   print 'ID & Total & Hap 1 & Hap 2 & Bac \\\\'
   print '\\hline \\hline'
   orderedByRank = sorted( team.entries, key=lambda x: x.rank, reverse=False )
   for a in orderedByRank:
      printRankRow( a )
   print '\\hline'
   print '\\end{tabular}\par'
   
   print '\\vspace{0.25in}'
   print '\\noindent \\textbf{Assemblies:}\par'
   for a in team.entries:
      printAssembly( a, assembliesList, captionsDict, options )
      print '\n'
   print '\\clearpage'
   print '\n'

def printLatex( teamsDict, assembliesList, captionsDict, options ):
   order = teamsDict.keys()
   order.sort()
   for t in order:
      printTeam( teamsDict[ t ], assembliesList, captionsDict, options )

def sortEntriesLists( teamsDict ):
   for t in teamsDict:
      teamsDict[t].entries = sorted( teamsDict[t].entries, key=lambda x: int(x.ID[1:]), reverse=False)

def readSubStatsDir( assembliesList, options ):
   for a in assembliesList:
      if os.path.exists( os.path.join( options.subStatsDir, a.ID+'.subStats.lower.txt')):
         f = open( os.path.join( options.subStatsDir, a.ID+'.subStats.lower.txt'), 'r' )
         for line in f:
            line = line.strip()
            d = line.split('\t')
            a.subStatsLower[ d[0] ] = d[ 1 ]
         f.close()
      if os.path.exists( os.path.join( options.subStatsDir, a.ID+'.subStats.upper.txt')):
         f = open( os.path.join( options.subStatsDir, a.ID+'.subStats.upper.txt'), 'r' )
         for line in f:
            line = line.strip()
            d = line.split('\t')
            a.subStatsUpper[ d[0] ] = d[ 1 ]
         f.close()

def readSubmissionStatsDir( assembliesList, options ):
   for a in assembliesList:
      if os.path.exists( os.path.join( options.submissionStatsDir, a.ID+'.summary.txt')):
         f = open( os.path.join( options.submissionStatsDir, a.ID+'.summary.txt'), 'r' )
         scaffold = False
         for line in f:
            line = line.strip()
            d = line.split()
            if d[0] == 'n':
               scaffold = not scaffold
               continue
            if scaffold:
               a.sizeStatsScaffold = d
            else:
               a.sizeStatsContigs  = d
         f.close()

def readSubmissionLengthsDir( assembliesList, options ):
   for a in assembliesList:
      if os.path.exists( os.path.join( options.submissionLengthsDir, a.ID+'.scaffolds.txt')):
         f = open( os.path.join( options.submissionLengthsDir, a.ID+'.scaffolds.txt'), 'r' )
         cum = 0
         for line in f:
            line = line.strip()
            cum += int( line )
         f.close()
         a.sizeStatsScaffold.append( cum )
      if os.path.exists( os.path.join( options.submissionLengthsDir, a.ID+'.contigs.txt')):
         f = open( os.path.join( options.submissionLengthsDir, a.ID+'.contigs.txt'), 'r' )
         cum = 0
         for line in f:
            line = line.strip()
            cum += int( line )
         f.close()
         a.sizeStatsContigs.append( cum )

def main():
   usage = ( 'usage: %prog --rankFile=rFile --infoFile=iFile --subStatsDir=path/to/dir/ --submissionStatsDir=path/to/dir/ --submissionLengthsDir=path/to/dir/ [options]\n\n'
             '%prog writes the Individual Section of the report in latex format.\n\n'
             '%prog takes in an assembly rank file ( --rankFile ), an assembly info file\n'
             ' ( --infoFile ), a substitution stats directory ( --subStatsDir ), a submission stats\n'
             'directory ( --submissionStatsDir ), a submissionLengthsDir ( --submissionLengthsDir )\n'
             'and an images directory ( --imagesDir ) where eps images are already stored and prints\n'
             'to STDOUT a latex formated section of text.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )

   teamsDict = readInfoFile( options )
   assembliesList = readRankFile( teamsDict, options )
   readSubStatsDir( assembliesList, options )
   readSubmissionStatsDir( assembliesList, options )
   readSubmissionLengthsDir( assembliesList, options )
   
   sortEntriesLists( teamsDict )

   boilerPlate = ''
   captionsDict = { 'contigs':'The three facets, from bottom to top are the cumulative stacked proportion of alignment column category as a function of contig length, the cumulative number of errors relative to the largest number of errors in all assemblies, and a blow-up region of the first facet. Regions are colored according to alignment columns that are (1) present in both haplotype 1 and 2 but absent in the assembly, (2) present in either haplotype 1 or 2 but not both and absent in the assembly, (3) present in either haplotype 1 or 2 but not both and present in the assembly (4) present in both haplotype 1 and haplotype 2 and the assembly, and (5) which refers to the middle facet are columns absent in both haplotype 1, haplotype 2, the bacterial contamination and present in the assembly.',
                    'scaffolds':'The two facets, on the bottom is the cumulative stacked proportion of alignment column category as a function of scaffold path length, and on top is a blow-up region of the first facet. Regions are colored according to alignment columns that are (1) present in both haplotype 1 and 2 but absent in the assembly, (2) present in either haplotype 1 or 2 but not both and absent in the assembly, (3) present in either haplotype 1 or 2 but not both and present in the assembly and (4) present in both haplotype 1 and haplotype 2 and the assembly.',
                    'contigPaths':'The two facets, on the bottom is the cumulative stacked proportion of alignment column category as a function of contig path length, and on top is a blow-up region of the first facet. Regions are colored according to alignment columns that are (1) present in both haplotype 1 and 2 but absent in the assembly, (2) present in either haplotype 1 or 2 but not both and absent in the assembly, (3) present in either haplotype 1 or 2 but not both and present in the assembly and (4) present in both haplotype 1 and haplotype 2 and the assembly.',
                    'blocks':'The three facets, from bottom to top are the cumulative stacked proportion of alignment column category as a function of block length, the cumulative number of errors relative to the largest number of errors in all assemblies, and a blow-up region of the first facet. Regions are colored according to alignment columns that are (1) present in both haplotype 1 and 2 but absent in the assembly, (2) present in either haplotype 1 or 2 but not both and absent in the assembly, (3) present in either haplotype 1 or 2 but not both and present in the assembly (4) present in both haplotype 1 and haplotype 2 and the assembly, and (5) which refers to the middle facet are columns absent in both haplotype 1, haplotype 2, the bacterial contamination and present in the assembly.' }
   printLatex( teamsDict, assembliesList, captionsDict, options )

if __name__ == '__main__':
   main()
