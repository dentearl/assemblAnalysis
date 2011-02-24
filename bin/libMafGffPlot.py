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
   """
   def __init__( self ):
      self.refGenome  = ''
      self.refChr     = ''
      self.refStart   = -1
      self.refEnd     = -1
      self.refStrand  = 0
      self.refSeq     = '' # future use
      self.pairGenome = ''
      self.pairChr    = ''
      self.pairStart  = -1
      self.pairEnd    = -1
      self.pairStrand = 0
      self.pairSeq    = '' # future use
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

def objListToBinnedWiggle( objList, featLen, numBins, filename ):
    """ obj can be either a GffRecord object or a MafBlock object.
    featLen is the length of the chromosome.
    returns a numpy vector of length numBins normalized by the maximum
    possible number of bases per bin.
    """
    from libMafGffPlot import GffRecord
    from libMafGffPlot import MafBlock
    import numpy
    import sys
    if objList == None or len( objList ) < 1:
            return None
    if isinstance( objList[0], GffRecord ):
        """ the Gff return is a single numpy vector of numBins length
        """
        vec = numpy.zeros( shape = ( numBins ))
        if objList == None or len( objList ) < 1:
            return vec
        maxCount = 0
        for a in objList:
            # index position in a 'numBins' length array.
            if a.start > featLen or a.end > featLen:
                sys.stderr.write( 'Error, file %s has annotation on chr %s '
                                  'with bounds [%d - %d] which are beyond featLen (%d)\n' %
                                  ( filename, a.chr, a.start, a.end, featLen ))
                sys.exit( 1 )
            for i in range( a.start, a.end + 1 ):
                pos = int(( float( i ) / (( featLen + 1.0 ) / float( numBins ) )))
                vec[ pos ] += 1
                if vec[ pos ] > maxCount:
                    maxCount = vec[ pos ]
        for i in range( 0, numBins):
            vec[ i ] /= float( maxCount )
        return vec
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
        blockEdgeDensity  each block has two edges, a left and a right
        """
        data = { 'maf'   : numpy.zeros( shape = ( numBins )),
                 'maf1e2': numpy.zeros( shape = ( numBins )),
                 'maf1e3': numpy.zeros( shape = ( numBins )),
                 'maf1e4': numpy.zeros( shape = ( numBins )),
                 'maf1e5': numpy.zeros( shape = ( numBins )),
                 'maf1e6': numpy.zeros( shape = ( numBins )),
                 'maf1e7': numpy.zeros( shape = ( numBins )),
                 'xAxis' : numpy.zeros( shape = ( numBins )),
                 'blockEdgeDensity': numpy.zeros( shape = ( numBins )) }
        maxPossibleCount = float( featLen ) / float( numBins )
        # populate xAxis
        for i in range( 0, numBins ):
            data['xAxis'][ i ] = (float( i ) / ( numBins - 1.0 )) * float( featLen )
        for m in objList:
            # do block edges
            for r in [ m.refStart, m.refEnd ]:
                pos = int(( float( m.refStart ) / (( featLen + 1.0 ) / float( numBins ) ) ))
                data['blockEdgeDensity'][ pos ] += 1
            # do all of the different maf block flavors
            for i in range( m.refStart, m.refEnd + 1 ):
                pos = int(( float( i ) / (( featLen + 1.0 ) / float( numBins ) ) ))
                data['maf'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 100:
                    data['maf1e2'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 1000:
                    data['maf1e3'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 10000:
                    data['maf1e4'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 100000:
                    data['maf1e5'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 1000000:
                    data['maf1e6'][ pos ] += 1
                if ( m.refEnd - m.refStart ) >= 10000000:
                    data['maf1e7'][ pos ] += 1
        for r in ['maf', 'maf1e2', 'maf1e3', 'maf1e4', 'maf1e5',
                  'maf1e6', 'maf1e7', 'blockEdgeDensity' ]:
            for i in range( 0, numBins):
                data[ r ][ i ] /= float( maxPossibleCount )
        return data
    else:
        return None
