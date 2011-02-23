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
    vec = numpy.zeros( shape = ( numBins ))
    if objList == None or len( objList ) < 1:
        return vec
    if isinstance( objList[0], GffRecord ):
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
        maxCount = float( featLen ) / numBins
        for m in objList:
            if m.refStart > featLen or m.refEnd > featLen:
                sys.stderr.write( 'Error, file %s has maf block on chr %s with '
                                  'bounds [%d - %d] which are beyond featLen (%d)\n' %
                                  ( filename, m.refChr, m.refStart, m.refEnd, featLen ))
                sys.exit( 1 )
            # index position in a 'numBins' length array.
            for i in range( m.refStart, m.refEnd + 1 ):
                pos = int(( float( i ) / (( featLen + 1.0 ) / float( numBins ) ) ))
                vec[ pos ] += 1
        for i in range( 0, numBins):
            vec[ i ] /= float( maxCount )
        return vec
    else:
        return None
