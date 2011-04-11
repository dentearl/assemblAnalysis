#!/usr/bin/env python
"""
latexWrapper.py
6 April 2011
dent earl, dearl (a) soe ucsc edu

latexWrapper takes a latex stream from STDIN,
wraps it in a hard coded header and prints to
STDOUT the result.

It is useful when you want to check the latex
output of a script, say a latex table, without
compiling the entire project.

"""
from optparse import OptionParser
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   pass

def checkOptions( args, options, parser ):
   pass

def main():
   usage = ( 'usage: %prog [options] < latexTable.tex \n\n'
             '%prog takes a latex file from STDIN and wraps it in a hard coded header\n'
             'and writes to STDOUT a latex compilable document.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   ( options, args ) = parser.parse_args()
   checkOptions( args, options, parser )
   
   header='''\\documentclass[11pt]{report}
\\usepackage{fltpage} % full page figures, captions on next page
\\usepackage{epsfig}
\\usepackage{empheq}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{moreverb} % for verbatimtab                                    
\\usepackage{multirow}
\\usepackage{hyperref} % internal links
%\\usepackage{rotating} % snpStats global Table
%\\usepackage{undertilde} % for utilde                                       
\\usepackage[table]{xcolor} % for table row colors                          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NSF margins                                       %
\\setlength{\\oddsidemargin}{0.5in}                    %
\\setlength{\\evensidemargin}{0.5in}                   %
\\setlength{\\topmargin}{0.0in}                        %
\\setlength{\\headheight}{0.1in}                       %
\\setlength{\\headsep}{0.1in}                          %
\\setlength{\\footskip}{0.5in}                         %
\\setlength{\\textheight}{8.75in}                      %
\\setlength{\\textwidth}{5.5in}                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%\\usepackage{color}                                                        
\\newcommand{\\var}{\\ensuremath{\\text{Var}}}
\\newcommand{\\vect}[1]{\\boldsymbol{#1}}
\\newcommand{\\superscript}[1]{\\ensuremath{^{\\textrm{#1}}}}
\\newcommand{\\subscript}[1]{\\ensuremath{_{\\textrm{#1}}}}
\\long\\def\\symbolfootnote[#1]#2{\\begingroup%                                
\\def\\thefootnote{\\fnsymbol{footnote}}\\footnote[#1]{#2}\\endgroup}
\\definecolor{light}{gray}{.75}
\\definecolor{myGray60}{gray}{.4}
\\definecolor{myGray40}{gray}{.6}
\\definecolor{myGray20}{gray}{.8}
\\definecolor{myLBlue}{rgb}{0.6862745, 0.8745098, 0.8941176} %{175,223,228} 
\\definecolor{myDOrange}{rgb}{0.9529412, 0.4313725, 0.1294118} %{243,110,33}
\\definecolor{tableShade}{HTML}{F1F5FA}   %iTunes                           
\\definecolor{tableShade2}{HTML}{ECF3FE}  %Finder                           

\\begin{document}
'''
   footer = '''
\\end{document}
'''
   
   print header
   for line in sys.stdin:
      sys.stdout.write( line )
   print footer
   

if __name__ == '__main__':
   main()
