# libAssemblySubset.py
# dent earl, dearl (a) soe ucsc edu
# 8 April 2011
#
# Module to handle the creation of a set of assemblies
# that are to be used by the calling script.
# It is up to each script to make use of the set 
# options.assemblySubset in it's own way.
#
# One suggestion is to place code like the following in a readData()
# function:
#
# for input in inputAssemblys:
#    if options.subsetFile:
#       if input.name not in options.assemblySubset:
#          continue
#   ...
#
def initOptions( parser ):
   parser.add_option( '--subsetFile', dest='subsetFile',
                      type='string', help=('Subset file contains one assembly ID per line.'))

def checkOptions( options, parser ):
   import os
   options.assemblySubset = set()
   if options.subsetFile == None:
      return
   if not os.path.exists( options.subsetFile ):
      parser.error('%s does not exist!\n' % ( options.subsetFile ))
   f = open( options.subsetFile, 'r' )
   for line in f:
      line = line.strip()
      if line.startswith('#'):
         continue
      d = line.split()
      assembly = d[0]
      options.assemblySubset.add( assembly )
