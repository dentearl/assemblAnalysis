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
   knownValues = (( '''>sequence1NNNNNnnnnnNNNnNn
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNGTACGTACGTACGTACGTACGTACGTACGTACG
>sequenceB
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG
>A1.scaffold003847
TAGCTTGTGTCATTATTAATTCTGTGACTGCTGGTACTTTTAGGTTTACC
AGTGAAATCAACTAACTATAATTACCTACAATATTGAGGCCTACCATATT
TAGGCGGTGTTGGTCCAGA
>A1.scaffold003848
TAGTGTGGGAGGTTGATAGTATATTGTTCCAAAATATTGTTCCCAACAGT
AATATGTTCACGATCTCATATTTCCTTAATCGTTACATGGAATTAACTTG
''', 
                    '''>sequence1NNNNNnnnnnNNNnNn.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
ACGTACGTACGTACNNNGTACGTACGTACGTACGTACGTACGTACGTACG
>sequenceB.split001
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
G
>A1.scaffold003847.split001
TAGCTTGTGTCATTATTAATTCTGTGACTGCTGGTACTTTTAGGTTTACC
AGTGAAATCAACTAACTATAATTACCTACAATATTGAGGCCTACCATATT
TAGGCGGTGTTGGTCCAGA
>A1.scaffold003848.split001
TAGTGTGGGAGGTTGATAGTATATTGTTCCAAAATATTGTTCCCAACAGT
AATATGTTCACGATCTCATATTTCCTTAATCGTTACATGGAATTAACTTG
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
      for i in xrange(30, 61):
         cmd = [ os.path.join( myBinDir, 'splitSequenceAtNs.py' ), 
                 '--splitAt=%d' % 25, '--lineLength=%d' % i,
                 '--label=%s' % 'split']
         p = subprocess.Popen( cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT )
         IN = '>seq1\n' + 'A' * 300
         ( streamOut ) = p.communicate( IN )[0]
         for s in streamOut:
            if s.startswith('>'):
               continue
            self.assertTrue( len( s ) <= i )

if __name__ == '__main__':
   unittest.main()
