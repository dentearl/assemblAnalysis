# libGeneral.py
# a library for repetitive code in the assemblathon analysis
# dent earl, dearl (a) soe ucsc edu
#
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
