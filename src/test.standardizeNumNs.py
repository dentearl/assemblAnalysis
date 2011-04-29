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
