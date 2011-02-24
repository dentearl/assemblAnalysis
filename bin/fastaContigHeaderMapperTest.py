#!/usr/bin/env python
import unittest
import os
import sys
myBinDir = os.path.normpath( os.path.dirname(sys.argv[0]) )
sys.path.append(myBinDir + "/../../..")
os.environ["PATH"] = myBinDir + "/../../../../bin:" + os.environ["PATH"]

class RoundTripCheck( unittest.TestCase ):
   import os
   knownValues = (('''>name1
ACGTnnnACGT
>name2
ACGttttttttt
ttttttttt
''','''>contig000001
ACGTnnnACGT
>contig000002
ACGttttttttt
ttttttttt
'''), ('''>apple
ACTGT
>apple2
ACTGTACTGT
>Horrible W0rds and a tab	 4@!#@!!!$&*){}
ACGTACGT
>emptyContig

>Some other stuff, odd extra space.
ACGT

>Last one
TGCATGCAacgt bad characters
''', '''>contig000001
ACTGT
>contig000002
ACTGTACTGT
>contig000003
ACGTACGT
>contig000004

>contig000005
ACGT

>contig000006
TGCATGCAacgt bad characters
'''))
   if not os.path.exists( 'tempTestFiles' ):
      os.mkdir( 'tempTestFiles' )
   def test_oneWay( self ):
      """fastaContigHeaderMapper should produce known results."""
      import subprocess
      for pre, post in self.knownValues:
         # generate map
         cmd = [os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--createMap=%s' %
                os.path.join('tempTestFiles','testMap.map')]
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT )
         ( sout ) = p.communicate( pre )[0]
         # go forward
         cmd = [os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--map=%s' %
                os.path.join('tempTestFiles','testMap.map'),
                '--goForward']
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT )
         ( outFor ) = p.communicate( pre )[0]
         self.assertEqual( post, outFor )

   def test_roundTrip( self ):
      """fastaContigHeaderMapper should be invertible."""
      import subprocess
      for pre, post in self.knownValues:
         # generate map
         cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--createMap=%s' %
                os.path.join('tempTestFiles','testMap.map')]
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT )
         ( sout ) = p.communicate( pre )[0]
         # go forward
         cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--map=%s' %
                os.path.join('tempTestFiles','testMap.map'),
                '--goForward']
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT )
         ( outFor ) = p.communicate( pre )[0]
         self.assertEqual( post, outFor )
         # go backward
         cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--map=%s' %
                os.path.join('tempTestFiles','testMap.map'),
                '--goBackward']
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                               stderr=subprocess.STDOUT )
         ( outBack ) = p.communicate( outFor )[0]
         
         self.assertEqual( pre, outBack )
   def test_roundTripPrefix( self ):
      """fastaContigHeaderMapper should be invertible with prefixes."""
      import random
      import string
      import subprocess
      print ' '
      chars = string.letters + string.digits + ' ' + '\t' + string.punctuation
      for i in range(50):
         prefix = ''.join( random.choice( chars ) for x in range(30))
         for pre, post in self.knownValues:
            #add prefix to post
            post2 = ''
            j = 0
            for p in post.split('\n'):
               j += 1
               p = p.strip()
               if p == '':
                  if j != len( post.split('\n') ):
                     post2 += '\n'
                  continue
               if p[0] == '>':
                  post2+= '>%s.%s\n' % ( prefix, p[1:] )
               else:
                  post2+= '%s\n' % p
            post = post2
            # generate map
            cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--createMap=%s' %
                   os.path.join('tempTestFiles','testMap.map'),
                   '--prefix=%s' % prefix]
            p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT )
            ( sout ) = p.communicate( pre )[0]
            # go forward
            cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--map=%s' %
                   os.path.join('tempTestFiles','testMap.map'),
                   '--goForward']
            p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT )
            ( outFor ) = p.communicate( pre )[0]
            self.assertEqual( post, outFor )
            # go backward
            cmd = [ os.path.join( myBinDir, 'fastaContigHeaderMapper.py'), '--map=%s' %
                   os.path.join('tempTestFiles','testMap.map'),
                   '--goBackward']
            p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT )
            ( outBack ) = p.communicate( outFor )[0]
            
            self.assertEqual( pre, outBack )

if __name__ == '__main__':
   unittest.main()