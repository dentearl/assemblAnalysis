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
import signal # deal with broken pipes
import sys

signal.signal( signal.SIGPIPE, signal.SIG_DFL ) # broken pipes

def initOptions( parser ):
   parser.add_option( '--noMargins', dest='noMargins',
                      default=False, action='store_true',
                      help=('Turns off the default NSF style margins.'))

def checkOptions( args, options, parser ):
   pass

def main():
   usage = ( 'usage: %prog [options] < latexTable.tex \n\n'
             '%prog takes a latex file from STDIN and wraps it in a hard coded header\n'
             'and writes to STDOUT a latex compilable document.')
   parser = OptionParser( usage=usage )
   initOptions( parser )
   options, args = parser.parse_args()
   checkOptions( args, options, parser )

   NSFMargins = '''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NSF margins                                       %
\\setlength{\\oddsidemargin}{-0.5in}                   %
\\setlength{\\evensidemargin}{-0.5in}                  %
\\setlength{\\topmargin}{0.0in}                        %
\\setlength{\\headheight}{0.1in}                       %
\\setlength{\\headsep}{0.1in}                          %
\\setlength{\\footskip}{0.5in}                         %
\\setlength{\\textheight}{8.75in}                      %
\\setlength{\\textwidth}{7.0in}                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
   if options.noMargins:
      margins = ''
   else:
      margins = NSFMargins
      
   header='''\\documentclass[11pt]{report}
\\usepackage{fltpage} %% full page figures, captions on next page
\\usepackage{epsfig}
\\usepackage{empheq}
\\usepackage{amsmath}
\\usepackage{amssymb}
\\usepackage{moreverb} %% for verbatimtab                                    
\\usepackage{multirow}
\\usepackage{hyperref} %% internal links
%%\\usepackage{rotating} %% snpStats global Table
%%\\usepackage{undertilde} %% for utilde                                       
\\usepackage[table]{xcolor} %% for table row colors                          
%s
%%\\usepackage{color}                                                        
\\newcommand{\\var}{\\ensuremath{\\text{Var}}}
\\newcommand{\\vect}[1]{\\boldsymbol{#1}}
\\newcommand{\\superscript}[1]{\\ensuremath{^{\\textrm{#1}}}}
\\newcommand{\\subscript}[1]{\\ensuremath{_{\\textrm{#1}}}}
\\long\\def\\symbolfootnote[#1]#2{\\begingroup%%                                
\\def\\thefootnote{\\fnsymbol{footnote}}\\footnote[#1]{#2}\\endgroup}
\\definecolor{light}{gray}{.75}
\\definecolor{myGray60}{gray}{.4}
\\definecolor{myGray40}{gray}{.6}
\\definecolor{myGray20}{gray}{.8}
\\definecolor{myLBlue}{rgb}{0.6862745, 0.8745098, 0.8941176} %%{175,223,228} 
\\definecolor{myDOrange}{rgb}{0.9529412, 0.4313725, 0.1294118} %%{243,110,33}
\\definecolor{tableShade}{HTML}{F1F5FA}   %%iTunes                           
\\definecolor{tableShade2}{HTML}{ECF3FE}  %%Finder                           

\\begin{document}
''' % margins
   footer = '''
\\end{document}
'''
   
   print header
   for line in sys.stdin:
      sys.stdout.write( line )
   print footer
   

if __name__ == '__main__':
   main()
