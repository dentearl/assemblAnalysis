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
class Data:
   """ A single Data object will be used to store
   small pieces of information that need to be used
   by many different functions.
   """
   def __init__( self ):
      self.chroms            = {} # values will be strings
      self.mafBlocksByChrom  = {} # values will be lists of MafBlock objs
      self.gffRecordsByChrom = {} # values will be dicts of key: annot 
                                  # names, value: list of GffRecord objs
class MafBlock:
   """ a MafBlock object is made up of the reference
   and it's pair. The both have their own coordinates 
   that define the nature of their alignment to one
   another. When buliding a MafBlock from file one can 
   give refEnd and pairEnd values equal to refStart and
   pairStart and then walk through the sequences and 
   use the increment method whenever there is not a gap
   '-' character in the sequence.
   Hap path edge (five, three) codes:
    * 0 = contig ends.
    * 1 = correct adjacency.
    * 2 = error adjacency.
    * 3 = scaffold gap
   """
   def __init__( self ):
      self.refGenome   = ''
      self.refChr      = ''
      self.refStart    = -1
      self.refEnd      = -1
      self.refStrand   = 0
      self.refTotalLength = -1
      self.refSeq      = '' # future use
      self.pairGenome  = ''
      self.pairChr     = ''
      self.pairStart   = -1
      self.pairEnd     = -1
      self.pairStrand  = 0
      self.pairTotalLength = -1
      self.pairSeq     = '' # future use
      self.hpl         = -1
      self.hplStart    = -1
      self.hplEnd      = -1
      self.spl         = -1
   def increment( self ):
      self.refEnd  += self.refStrand
      self.pairEnd += self.pairStrand

class MafLine:
   """ a MafLine object stores the raw input from the maf file
   and is subsequently joined with another MafLine object to 
   become one or more MafBlock objects.
   """
   def __init__( self ):
      self.genome = ''
      self.chr    = ''
      self.start  = -1
      self.length = -1
      self.strand = ''
      self.totalLength = -1
      self.sequence = ''
      # order > -1 implies this sequence is the comparativeGenome
      # and this is the order the sequence appears in the block from
      # (0..n-1) for a block with n sequences.
      self.order = -1 


class GffRecord:
   """ a GffRecord object contains the relevant information needed
   from one line of a .gff file to create a plot of annotation 
   density across the genome.
   """
   def __init__( self ):
      self.chr    = ''
      self.source = ''
      self.type   = ''
      self.start  = -1
      self.end    = -1
      self.score  = -1
      self.strand = ''
      self.frame  = ''
      self.group  = ''

def newMafWigDict( numBins ):
    """ Returns a new mafWigDict object. note that the xAxis starts 
    with a dummy valuable. The xAxis must be filled in with objListUtility_xAxis().
    """ 
    import numpy
    import sys
    if numBins < 1:
        sys.stderr.write('libMafGffPlot.py, numBins=%s is less than one\n' % str(numBins))
        sys.exit( 1 )
    return { 'maf'       : numpy.zeros( shape = ( numBins )),
             'maf1e2'    : numpy.zeros( shape = ( numBins )),
             'maf1e3'    : numpy.zeros( shape = ( numBins )),
             'maf1e4'    : numpy.zeros( shape = ( numBins )),
             'maf1e5'    : numpy.zeros( shape = ( numBins )),
             'maf1e6'    : numpy.zeros( shape = ( numBins )),
             'maf1e7'    : numpy.zeros( shape = ( numBins )),
             'xAxis'     : -1,
             'mafCpl1e2' : numpy.zeros( shape = ( numBins )),
             'mafCpl1e3' : numpy.zeros( shape = ( numBins )),
             'mafCpl1e4' : numpy.zeros( shape = ( numBins )),
             'mafCpl1e5' : numpy.zeros( shape = ( numBins )),
             'mafCpl1e6' : numpy.zeros( shape = ( numBins )),
             'mafCpl1e7' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e2' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e3' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e4' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e5' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e6' : numpy.zeros( shape = ( numBins )),
             'mafCtg1e7' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e2' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e3' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e4' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e5' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e6' : numpy.zeros( shape = ( numBins )),
             'mafSpl1e7' : numpy.zeros( shape = ( numBins )),
             'mafCpEdgeCount'    : numpy.zeros( shape = ( numBins ), dtype=numpy.int ),
             'mafCpEdgeMax'      : 0,
             'mafCpErrorCount'   : numpy.zeros( shape = ( numBins ), dtype=numpy.int ),
             'mafCpErrorMax'     : 0,
             'mafCpScafGapCount' : numpy.zeros( shape = ( numBins ), dtype=numpy.int ),
             'mafCpScafGapMax'   : 0,
             'blockEdgeCount'    : numpy.zeros( shape = ( numBins ), dtype=numpy.int ),
             'blockEdgeMax'      : 0,
             'columnsInBlocks'   : 0 }

def objListToBinnedWiggle( objList, featLen, numBins, filename ):
   """ obj can be either a GffRecord object or a MafBlock object.
   featLen is the length of the chromosome.
   returns a numpy vector of length numBins normalized by the maximum
   possible number of bases per bin.
   """
   from libMafGffPlot import GffRecord
   from libMafGffPlot import MafBlock
   from libMafGffPlot import newMafWigDict
   from libMafGffPlot import objListUtility_xAxis
   import numpy
   import sys
   if objList is None or len( objList ) < 1:
      return None
   if isinstance( objList[0], GffRecord ):
      """ the Gff return is a single numpy vector of numBins length
      """
      data = {}
      # populate xAxis
      data['xAxis'] = objListUtility_xAxis( featLen, numBins )
      annotTypes = set([ 'CDS', 'UTR', 'NXE', 'NGE', 
                         'island', 'tandem', 'repeat' ])
      for t in annotTypes:
         data[ t + 'Count' ] = numpy.zeros( shape = ( numBins ))
         data[ t + 'Max' ]   = 0
        
      for a in objList:
         if a.type not in annotTypes:
            continue
         # verify input
         if a.start > featLen or a.end > featLen:
            sys.stderr.write( 'libMafGffPlot.py: file %s has annotation on chr %s '
                              'with bounds [%d - %d] which are beyond featLen (%d)\n' %
                              ( filename, a.chr, a.start, a.end, featLen ))
            sys.exit( 1 )
         # index position in a 'numBins' length array.
         pos = objListUtility_rangeToPos( a.start, a.end, featLen, numBins )
            
         # tough to follow index hack to get around the fact that numpy will not use
         # += 1 for a list of indices that contain repeats.
         plo, phi = pos.min(), pos.max()
         pbins = numpy.bincount( pos - plo )
         data[ a.type + 'Count' ][ plo:phi + 1 ] += pbins

         for p in pos:
            if data[ a.type + 'Max' ] < data[ a.type + 'Count' ][ p ]:
               data[ a.type + 'Max' ] = data[ a.type + 'Count' ][ p ]
      return data
   elif isinstance( objList[0], MafBlock ):
      """ the Maf return is a dictionary with the following keys
      maf               all maf block bases
      maf1e2            maf blocks 100 or greater
      maf1e3            maf blocks 1,000 or greater
      maf1e4            maf blocks 10,000 or greater
      maf1e5            maf blocks 100,000 or greater
      maf1e6            maf blocks 1,000,000 or greater
      maf1e7            maf blocks 10,000,000 or greater
      xAxis             x Values

      mafCpl1eX         maf contig paths of X or greater

      mafCtg1eX         maf contigs of X or greater. taken from totalLength field of maf.
      
      mafSpl1eX         maf scaffold paths of X or greater
      
      mafCpEdgeCounts   each contig path has two edges, a left and a right
      mafCpEdgeMax      max count
      mafCpErrorCounts  contig paths are made up of segments, segments may have errors at junctions.
      mafCpErrorMax     max count
      mafSpEdgeCounts   Same as above, but for scaffold paths
      mafSpEdgeMax      
      mafSpErrorCounts  
      mafSpErrorMax     
      blockEdgeCounts   each block has two edges, a left and a right
      blockEdgeMax      max count
      
      """
      from libMafGffPlot import objListUtility_addContigPathEdgeErrors
      from libMafGffPlot import objListUtility_addBlockEdges
      from libMafGffPlot import objListUtility_normalizeCategories
      data = newMafWigDict( numBins )
      
      # populate xAxis
      data['xAxis'] = objListUtility_xAxis( featLen, numBins )
      for mb in objList:
         # do block edges
         objListUtility_addBlockEdges( data, mb, featLen, numBins )
         
         # do contige path edges and errors
         objListUtility_addContigPathEdgeErrors( data, mb, featLen, numBins )
            
         # do all of the different maf block flavors
         objListUtility_mafBlockCounts( data, mb, featLen, numBins )
            
      # normalize all categories
      objListUtility_normalizeCategories( data, featLen, numBins )
        
      return data
   # closing the elif isinstance() checks
   else:
      return None

def objListUtility_addContigPathEdgeErrors( data, mb, featLen, numBins ):
   """ Utility function for the MafBlock instance version of 
   libMafGffPlot.objListToBinnedWiggle()
   """ 
   from libMafGffPlot import objListUtility_indexToPos
   import math
   posSt  = objListUtility_indexToPos( mb.refStart, featLen, numBins )
   posEnd = objListUtility_indexToPos( mb.refEnd, featLen, numBins )
   # five prime
   if mb.hplStart == 0:
      # edge
      data['mafCpEdgeCount'][ posSt ] += 1.0
      if data['mafCpEdgeCount'][ posSt ] > data[ 'mafCpEdgeMax' ]:
         data[ 'mafCpEdgeMax' ] = data['mafCpEdgeCount'][ posSt ]
   elif mb.hplStart == 2:
      # error
      data['mafCpErrorCount'][ posSt ] += 1.0
      if data['mafCpErrorCount'][ posSt ] > data[ 'mafCpErrorMax' ]:
         data[ 'mafCpErrorMax' ] = data['mafCpErrorCount'][ posSt ]
   elif mb.hplStart == 3:
      # scaffold gap
      data['mafCpScafGapCount'][ posSt ] += 1.0
      if data['mafCpScafGapCount'][ posSt ] > data[ 'mafCpScafGapMax' ]:
         data[ 'mafCpScafGapMax' ] = data['mafCpScafGapCount'][ posSt ]
   # three prime
   if mb.hplEnd == 0:
      # edge
      data['mafCpEdgeCount'][ posEnd ] += 1.0
      if data['mafCpEdgeCount'][ posEnd ] > data[ 'mafCpEdgeMax' ]:
         data[ 'mafCpEdgeMax' ] = data['mafCpEdgeCount'][ posEnd ]
   elif mb.hplEnd == 2:
      # error
      data['mafCpErrorCount'][ posEnd ] += 1.0
      if data['mafCpErrorCount'][ posEnd ] > data[ 'mafCpErrorMax' ]:
         data[ 'mafCpErrorMax' ] = data['mafCpErrorCount'][ posEnd ]
   elif mb.hplEnd == 3:
      # scaffold gap
      data['mafCpScafGapCount'][ posEnd ] += 1.0
      if data['mafCpScafGapCount'][ posEnd ] > data[ 'mafCpScafGapMax' ]:
         data[ 'mafCpScafGapMax' ] = data['mafCpScafGapCount'][ posEnd ]

def objListUtility_normalizeCategories( data, featLen, numBins ):
   """ Utility function for the MafBlock instance version of 
   libMafGffPlot.objListToBinnedWiggle()
   """ 
   import math
   import sys
   maxPossibleCount = math.ceil( float( featLen ) / float( numBins ))
    
   for r in [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
              'maf1e5', 'maf1e6', 'maf1e7', 'mafCpl1e2',
              'mafCpl1e3', 'mafCpl1e4', 'mafCpl1e5', 
              'mafCpl1e6', 'mafCpl1e7',
              'mafCtg1e2', 'mafCtg1e3', 'mafCtg1e4',
              'mafCtg1e5', 'mafCtg1e6', 'mafCtg1e7', 
              'mafSpl1e2', 'mafSpl1e3', 'mafSpl1e4', 
              'mafSpl1e5', 'mafSpl1e6', 'mafSpl1e7' ]:

      # verify data 1
      if sum( data[ r ] > maxPossibleCount ) > 0:
         sys.stderr.write('libMafGffPlot.py: Error in normalization step, category \'%s\' has elements '
                          'greater than max %d (featLen/numBins = %d/%d)\n' 
                          % ( r, maxPossibleCount, featLen, numBins ))
         # verify data 2
         i = -1
         for d in data[ r ]:
            i += 1
            if d > maxPossibleCount:
               start = math.floor( i * ( float( featLen ) / numBins ))
               end   = math.floor( (i + 1) * ( float( featLen ) / numBins ))
               sys.stderr.write('   i=%d [%d,%d] count %d\n' % ( i, start, end, d ))
               sys.exit( 1 )

      # normalize
      data[ r ] /= float( maxPossibleCount )

def objListUtility_addBlockEdges( data, mb, featLen, numBins ):
   """ Utility function for the MafBlock instance version of 
   libMafGffPlot.objListToBinnedWiggle()
   """ 
   from libMafGffPlot import objListUtility_indexToPos
   import math
   import sys
   for r in [ mb.refStart, mb.refEnd ]:
      if r > featLen:
         sys.stderr.write('libMafGffPlot.py: Error in block edge step, a position is '
                          'greater than the feature length, %d > %d.\n' 
                          % (r, featLen))
         sys.exit( 1 )

      p = objListUtility_indexToPos( r, featLen, numBins )
      if p < 0:
         sys.stderr.write('libMafGffPlot.py: Error in block edge step, a position, %d, is less than 0\n' % p )
      elif p >= len( data['blockEdgeCount'] ):
         sys.stderr.write('a position, %d, is greater than or '
                          'equal to len(data[\'blockEdgeCount\']) %d [%d-%d]\n' 
                          % (p, len(data['blockEdgeCount']), mb.refStart, mb.refEnd))
      data['blockEdgeCount'][ p ] += 1.0
      if data['blockEdgeCount'][ p ] > data[ 'blockEdgeMax' ]:
         data[ 'blockEdgeMax' ] = data['blockEdgeCount'][ p ]

def objListUtility_mafBlockCounts( data, mb, featLen, numBins ):
   """ Utility function for the MafBlock instance version of 
   libMafGffPlot.objListToBinnedWiggle()
   This is by far the most costly routine in the objList creation process
   """ 
   from libMafGffPlot import objListUtility_rangeToPos
   import numpy
   length = mb.refEnd - ( mb.refStart + 1 )
   
   # tough to follow index hack to get around the fact that numpy will not use
   # += 1 for a list of indices that contain repeats.
   pos = objListUtility_rangeToPos( mb.refStart, mb.refEnd, featLen, numBins )
   plo, phi = pos.min(), pos.max()
   pbins = numpy.bincount( pos - plo )
   data['maf'][ plo:phi + 1 ] += pbins
   
   for i in xrange( 2, 8 ):
      if length >= 10 ** i:
         data[ 'maf1e%d' % i ][ plo:phi+1 ] += pbins
      if mb.spl >= 10 ** i:
         data[ 'mafSpl1e%d' % i ][ plo:phi+1 ] += pbins
      if mb.pairTotalLength >= 10 ** i:
         data[ 'mafCtg1e%d' % i ][ plo:phi+1 ] += pbins
      if mb.hpl >= 10 ** i:
         data[ 'mafCpl1e%d' % i ][ plo:phi+1 ] += pbins

def objListUtility_xAxis( featLen, numBins ):
    """ Utility function for the MafBlock instance version of 
    libMafGffPlot.objListToBinnedWiggle()
    """ 
    import numpy
    xaxis = numpy.linspace( start=0, stop=featLen, num=numBins )
    return xaxis

def objListUtility_indexToPos( p, featLen, numBins ):
    """ Utility function for the MafBlock instance version of 
    libMafGffPlot.objListToBinnedWiggle()
    Takes in a single index position in a [1, featLen] range and then provides the
    appropriate position in an array of length numBins. Maps an integer to the space
    [0, numBins].
    """ 
    import math
    import numpy
    # map point p from [1, featLen] to [0, featLen - 1]
    z  = float( p - 1.0 )
    # scale z to [0, featLen)
    z /= featLen
    # map  to [0, numBins)
    z *= float( numBins )
    # return the index in a length numBins array
    return math.floor( z )

def objListUtility_rangeToPos( p1, p2, featLen, numBins ):
   """ Utility function for the MafBlock instance version of 
   libMafGffPlot.objListToBinnedWiggle()
   Takes in two positions, p1 and p2, that are in the range [1, featLen] and 
   returns a numpy array containing indices in the range [0, numBins]. Maps
   the integers p1..p2 to the space [0, numBins].
   NOTE: the array we return can contain repeated indices. If it does and you
   store the result in pos to access another array , a, you cannot
   simply use a[ pos ] += 1 and get what you expect. Instead, use:

   plo, phi = pos.min(), pos.max()
   pbins = numpy.bincount( pos - plo )
   a[ plo:phi+1 ] += pbins
   
   This will be slightly faster than walking through each value in pos and 
   incrementing.
   """ 
   import numpy
   # convert [1, featLen ] to [0, featLen - 1]
   z  = numpy.arange( start=p1, stop=p2 + 1.0, step=1, dtype=numpy.int ) - 1.0
   # scale elements in array to [0, 1)
   z /= featLen
   # map elements in array to [0, numBins)
   z *= float( numBins )
   # return indices for the mapping
   return ( numpy.floor( z ) ).astype('int')

def packData( packMe, filename, options, prot='py23Bin' ):
   """ packData takes some data and a filename and
   packs it away. prot refers to the protocol to use.
   """
   import cPickle
   protocols = { 'ASCII' : 0,
                 'pre23Bin' : 1,
                 'py23Bin'  : 2 }
   f = open( filename, 'wb' )
   cPickle.dump( packMe, f, protocol = protocols[ prot ] )
   f.close()

def unpackData( filename, options, data ):
   """ unpackData unpacks a pickle object
   """
   import cPickle
   import os
   if not os.path.exists( filename ):
      sys.stderr.write( 'Error, %s does not exist.\n' % filename )
      sys.exit( 1 )
   f = open( filename, 'rb' )
   d = cPickle.load( f )
   f.close()
   return d

