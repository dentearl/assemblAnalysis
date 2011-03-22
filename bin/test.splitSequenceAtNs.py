import unittest
import os
import sys
myBinDir = os.path.abspath( os.path.dirname( sys.argv[0] ))

class VerifyKnownInputOutput( unittest.TestCase ):   
   knownValues = (( '''>sequence1NNNNNnnnnnNNNnNn
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNGTACGTACGTACGTACGTACGTACGTACGTACG
''', 
                    '''>sequence1NNNNNnnnnnNNNnNn.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNGTACGTACGTACGTACGTACGTACGTACGTACG
''') , 
                  ( '''>sequence2!@#$
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGN
NNNNNNNNNNNNNNNNNNNNNNNNNACGTACGTACGT
''',
                    '''>sequence2!@#$.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
>sequence2!@#$.split002
ACGTACGTACGT
'''), 
                  ( '''>sequence3
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNACGTACGT
NNNNNNNNNNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
ACGTACGT
''',
                    '''>sequence3.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTACGTACGTACGTACGTACGTACGT
>sequence3.split002
ACGTACGTNNNNNNNNNNACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGT
>sequence3.split003
ACGTACGT
'''),
                  ( '''>sequence4
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTNNNNNNNNNNNNNNNNNNNNNNNNN
''',
                    '''>sequence4.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGT
'''),
                  ( '''>sequence5
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTNNNNNNNNNNNNNNNNNNNNNNNN
''',
                    '''>sequence5.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
GTACGTNNNNNNNNNNNNNNNNNNNNNNNN
'''))
   
   def test_knownValues( self ):
      """splitSequenceAtNs should split sequences with N >= --splitAt
      """
      import os
      import subprocess
      cmd = [ os.path.join( myBinDir, 'splitSequenceAtNs.py' ), 
              '--splitAt=%d' % 25, '--lineLength=%d' % 50,
              '--label=%s' % 'split' ]
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
         cmd = [ os.path.join( myBinDir, 'splitSequenceAtNs.py' ), 
                 '--splitAt=%d' % 25, '--lineLength=%d' % i,
                 '--label=%s' % 'split']
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
