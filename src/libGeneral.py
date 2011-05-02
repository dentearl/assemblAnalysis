# libGeneral.py
# a library for repetitive code in the assemblathon analysis
# dent earl, dearl (a) soe ucsc edu
#
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
def prettyNumber( n ):
   """ pretty number takes a float or int and returns a string
   """
   from libGeneral import prettyFloat, prettyInt
   import re
   if isinstance( n, str ):
      pat = re.compile( '^[0-9]+\.[0-9]+$' )
      if re.match( pat, n ):
         return prettyFloat( float(n), 2 )
      else:
         return prettyInt( int(n) )
   elif isinstance( n, float ):
      return prettyFloat( n, 2 )
   elif isinstance( n, int ):
      return prettyInt( n )
   else:
      sys.stderr.write('Error, libGeneral.py: prettyNumber: unexpected type for n: %s' % n.__class__ )
      sys.exit(1)

def prettyInt( i ):
   s = ''
   r = '%d' % i
   for j in xrange(0, len( r )):
      if j > 0 and not j % 3:
         s = '%s,%s' % ( r[ (len(r) - 1) - j], s )
      else:
         s = '%s%s' % ( r[ (len(r) - 1) - j], s )
   return s

def prettyFloat( f, n ):
   s = ''
   r = '%.*f' % ( n, f )
   for j in xrange(0, len( r ) ):
      if not (j - n - 1) % 3 and j > (n + 1):
         s = '%s,%s' % ( r[ (len(r) - 1) - j], s )
      else:
         s = '%s%s' % ( r[ (len(r) - 1) - j], s )
   return s

##############################
# idMap provides a mapping between the assemly ID numbers used
# in the filesystem to the names displayed in figures and tables.
# If you add your own assembly into the analysis, you should also
# add your ID code and six letter name into this dict.
idMap = { 'A':'Astr', 'B':'Wtsi P', 'C':'Ebi', 'D':'Wtsi S',
          'E':'Cracs', 'F':'Bccgsc', 'G':'Doejgi', 'H':'Irisa',
          'I':'Cshl', 'J':'Dscisu', 'K':'Cbslug', 'L':'Ucsf', 
          'M':'Rhul', 'N':'Gacwt', 'O':'Dcsuoc', 'P':'Bgi', 
          'Q':'Broad', 'V':'Nvelv', 'W':'Nclc', 'X':'Nabyss' }
