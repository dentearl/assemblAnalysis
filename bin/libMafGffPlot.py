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
    import numpy
    return { 'maf'    : numpy.zeros( shape = ( numBins )),
             'maf1e2' : numpy.zeros( shape = ( numBins )),
             'maf1e3' : numpy.zeros( shape = ( numBins )),
             'maf1e4' : numpy.zeros( shape = ( numBins )),
             'maf1e5' : numpy.zeros( shape = ( numBins )),
             'maf1e6' : numpy.zeros( shape = ( numBins )),
             'maf1e7' : numpy.zeros( shape = ( numBins )),
             'xAxis'  : numpy.zeros( shape = ( numBins )),
             'mafHpl1e2' : numpy.zeros( shape = ( numBins )),
             'mafHpl1e3' : numpy.zeros( shape = ( numBins )),
             'mafHpl1e4' : numpy.zeros( shape = ( numBins )),
             'mafHpl1e5' : numpy.zeros( shape = ( numBins )),
             'mafHpl1e6' : numpy.zeros( shape = ( numBins )),
             'mafHpl1e7' : numpy.zeros( shape = ( numBins )),
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
             'mafHpEdgeCount'  : numpy.zeros( shape = ( numBins )),
             'mafHpEdgeMax'    : 0,
             'mafHpErrorCount' : numpy.zeros( shape = ( numBins )),
             'mafHpErrorMax'   : 0,
             'mafHpScafGapCount' : numpy.zeros( shape = ( numBins )),
             'mafHpScafGapMax'   : 0,
             'blockEdgeCount'  : numpy.zeros( shape = ( numBins )),
             'blockEdgeMax'   : 0 }

def objListToBinnedWiggle( objList, featLen, numBins, filename ):
    """ obj can be either a GffRecord object or a MafBlock object.
    featLen is the length of the chromosome.
    returns a numpy vector of length numBins normalized by the maximum
    possible number of bases per bin.
    """
    from libMafGffPlot import GffRecord
    from libMafGffPlot import MafBlock
    from libMafGffPlot import newMafWigDict
    import numpy
    import sys
    if objList == None or len( objList ) < 1:
        return None
    if isinstance( objList[0], GffRecord ):
        """ the Gff return is a single numpy vector of numBins length
        """
        data = { 'xAxis'    : numpy.zeros( shape = ( numBins ))}
        for t in [ 'CDS', 'UTR', 'NXE', 'NGE', 'island', 'tandem', 'repeat' ]:
                 data[ t + 'Count' ] = numpy.zeros( shape = ( numBins ))
                 data[ t + 'Max' ]   = 0
        # populate xAxis
        for i in range( 0, numBins ):
            data['xAxis'][ i ] = (float( i ) / ( numBins - 1.0 )) * float( featLen )

        annotTypes = set([ 'CDS', 'UTR', 'NXE', 'NGE', 
                           'island', 'tandem', 'repeat' ] )
        for a in objList:
            if a.type not in annotTypes:
                continue
            # index position in a 'numBins' length array.
            if a.start > featLen or a.end > featLen:
                sys.stderr.write( 'Error, file %s has annotation on chr %s '
                                  'with bounds [%d - %d] which are beyond featLen (%d)\n' %
                                  ( filename, a.chr, a.start, a.end, featLen ))
                sys.exit( 1 )
            for i in range( a.start, a.end + 1 ):
                pos = int(( float( i ) / (( featLen + 1.0 ) / float( numBins ) )))
                data[ a.type + 'Count' ][ pos ] += 1
                if data[ a.type + 'Max' ] < data[ a.type + 'Count' ][ pos ]:
                   data[ a.type + 'Max' ] = data[ a.type + 'Count' ][ pos ]
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

        mafHpl1eX         maf haplotype paths of X or greater

        mafCtg1eX         maf contigs of X or greater. taken from totalLength field of maf.
        
        mafSpl1eX         maf scaffold paths of X or greater
        
        mafHpEdgeCounts   each haplotype path has two edges, a left and a right
        mafHpEdgeMax      max count
        mafHpErrorCounts  haplotype paths are made up of segments, segments may have errors at junctions.
        mafHpErrorMax     max count
        mafSpEdgeCounts   Same as above, but for scaffold paths
        mafSpEdgeMax      
        mafSpErrorCounts  
        mafSpErrorMax     
        blockEdgeCounts   each block has two edges, a left and a right
        blockEdgeMax      max count
        
        """
        data = newMafWigDict( numBins )
        maxPossibleCount = float( featLen ) / float( numBins )
        # populate xAxis
        for i in range( 0, numBins ):
            data['xAxis'][ i ] = (float( i ) / ( numBins - 1.0 )) * float( featLen )
        for mb in objList:
            # do block edges
            for r in [ mb.refStart, mb.refEnd ]:
                pos = int(( float( r ) / (( featLen + 1.0 ) / float( numBins ) ) ))
                data['blockEdgeCount'][ pos ] += 1.0
                if data['blockEdgeCount'][ pos ] > data[ 'blockEdgeMax' ]:
                    data[ 'blockEdgeMax' ] = data['blockEdgeCount'][ pos ]
            # do haplotype path edges and errors
            posSt  = int(( float( mb.refStart  ) / (( featLen + 1.0 ) / float( numBins ) ) ))
            posEnd = int(( float( mb.refEnd    ) / (( featLen + 1.0 ) / float( numBins ) ) ))
            # five prime
            if mb.hplStart == 0:
                # edge
                data['mafHpEdgeCount'][ posSt ] += 1.0
                if data['mafHpEdgeCount'][ posSt ] > data[ 'mafHpEdgeMax' ]:
                    data[ 'mafHpEdgeMax' ] = data['mafHpEdgeCount'][ posSt ]
            elif mb.hplStart == 2:
                # error
                data['mafHpErrorCount'][ posSt ] += 1.0
                if data['mafHpErrorCount'][ posSt ] > data[ 'mafHpErrorMax' ]:
                    data[ 'mafHpErrorMax' ] = data['mafHpErrorCount'][ posSt ]
            elif mb.hplStart == 3:
                # scaffold gap
                data['mafHpScafGapCount'][ posSt ] += 1.0
                if data['mafHpScafGapCount'][ posSt ] > data[ 'mafHpScafGapMax' ]:
                    data[ 'mafHpScafGapMax' ] = data['mafHpScafGapCount'][ posSt ]
            # three prime
            if mb.hplEnd == 0:
                # edge
                data['mafHpEdgeCount'][ posEnd ] += 1.0
                if data['mafHpEdgeCount'][ posEnd ] > data[ 'mafHpEdgeMax' ]:
                    data[ 'mafHpEdgeMax' ] = data['mafHpEdgeCount'][ posEnd ]
            elif mb.hplEnd == 2:
                # error
                data['mafHpErrorCount'][ posEnd ] += 1.0
                if data['mafHpErrorCount'][ posEnd ] > data[ 'mafHpErrorMax' ]:
                    data[ 'mafHpErrorMax' ] = data['mafHpErrorCount'][ posEnd ]
            elif mb.hplEnd == 3:
                # scaffold gap
                data['mafHpScafGapCount'][ posEnd ] += 1.0
                if data['mafHpScafGapCount'][ posEnd ] > data[ 'mafHpScafGapMax' ]:
                    data[ 'mafHpScafGapMax' ] = data['mafHpScafGapCount'][ posEnd ]
            # do all of the different maf block flavors
            for i in range( mb.refStart, mb.refEnd + 1 ):
                pos = int(( float( i ) / (( featLen + 1.0 ) / float( numBins ) ) ))
                data['maf'][ pos ] += 1
                length = mb.refEnd - mb.refStart
                for i in range( 2, 8 ):
                    if length >= 10 ** i:
                        data[ 'maf1e%d' % i ][ pos ] += 1
                    if mb.spl >= 10 ** i:
                        data[ 'mafSpl1e%d' % i ][ pos ] += 1
                    if mb.pairTotalLength >= 10 ** i:
                        data[ 'mafCtg1e%d' % i ][ pos ] += 1
                    if mb.hpl >= 10 ** i:
                        data[ 'mafHpl1e%d' % i ][ pos ] += 1
                    
        for r in [ 'maf', 'maf1e2', 'maf1e3', 'maf1e4', 
                   'maf1e5', 'maf1e6', 'maf1e7', 'mafHpl1e2',
                   'mafHpl1e3', 'mafHpl1e4', 'mafHpl1e5', 
                   'mafHpl1e6', 'mafHpl1e7',
                   'mafCtg1e2', 'mafCtg1e3', 'mafCtg1e4',
                   'mafCtg1e5', 'mafCtg1e6', 'mafCtg1e7', 
                   'mafSpl1e2',
                   'mafSpl1e3', 'mafSpl1e4', 'mafSpl1e5', 
                   'mafSpl1e6', 'mafSpl1e7' ]:
            for i in range( 0, numBins ):
                data[ r ][ i ] /= float( maxPossibleCount )
        return data
    else:
        return None
