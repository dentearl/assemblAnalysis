import unittest
import os
import sys
myBinDir = os.path.abspath( os.path.dirname( sys.argv[0] ))

class VerifyKnownInputOutput( unittest.TestCase ):   
   knownValues = (('''>sequence1NNNNNnnnnnNNNnNn
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNGTACGTACGTACGTACGTACGTACGTACGTACG
''', '''>sequence1NNNNNnnnnnNNNnNn
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNNNNNNNNNNNNNNNNNNNNNNNGTACGTACGTA
CGTACGTACGTACGTACGTACG
'''),('''>sequence2!@#$
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGN
NNACGTACGTACGT
''','''>sequence2!@#$
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGN
NNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGT
'''))
   
   def test_knownValuesExpansion( self ):
      """standardizeNumNs should expand values >= --expandAt to 25
      """
      import os
      import subprocess
      cmd = [ os.path.join( myBinDir, 'standardizeNumNs.py' ), 
              '--expandAt=%d' % 3, '--lineLength=%d' % 50 ]
      for IN, expectedOut in self.knownValues:
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT )
         ( streamOut ) = p.communicate( IN )[0]
         self.assertEqual( streamOut, expectedOut )
   knownValues = (('''>sequence1NNNNNnnnnnNNNnNn
ACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNN
''', '''>sequence1NNNNNnnnnnNNNnNn
ACGTNNNNNNNNNNNNNNNNNNNNNNNNN
'''),('''>sequence2!@#$
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
NN
''','''>sequence2!@#$
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNN
'''))
   def test_knownValuesShrinkage( self ):
      """standardizeNumNs should shrink values > 25 to 25
      """
      import os
      import subprocess
      cmd = [ os.path.join( myBinDir, 'standardizeNumNs.py' ), 
              '--expandAt=%d' % 3, '--lineLength=%d' % 50 ]
      for IN, expectedOut in self.knownValues:
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT )
         ( streamOut ) = p.communicate( IN )[0]
         self.assertEqual( streamOut, expectedOut )
   def test_lineLengthSetting( self ):
      """ Length of output lines should be <= --lineLength.
      """
      import os
      import subprocess
      for i in range(30, 61):
         cmd = [ os.path.join( myBinDir, 'standardizeNumNs.py' ), 
                 '--expandAt=%d' % 25, '--lineLength=%d' % i ]
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT )
         IN = '>seq1\n' + 'A' * 300
         ( streamOut ) = p.communicate( IN )[0]
         for s in streamOut:
            if s[0] == '>':
               continue
            self.assertTrue( len( s ) <= i )

if __name__ == '__main__':
   unittest.main()
