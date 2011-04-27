# libPlotting.py
# a library for repetitive plotting code in the assemblathon analysis
# dent earl, dearl (a) soe ucsc edu
#

def initImage( width, height, options, data ):
   """
   initImage takes a width and height and returns
   both a fig and pdf object. options must contain outFormat,
   and dpi
   """
   import matplotlib.backends.backend_pdf as pltBack
   import matplotlib.pyplot as plt
   pdf = None
   if options.outFormat == 'pdf' or options.outFormat == 'all':
      pdf = pltBack.PdfPages( options.out + '.pdf' )
   fig = plt.figure( figsize=( width, height ), dpi=options.dpi, facecolor='w' )
   data.fig = fig
   return ( fig, pdf )

def writeImage( fig, pdf, options ):
   """
   writeImage assumes options contains outFormat and dpi.
   """
   if options.outFormat == 'pdf':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
   elif options.outFormat == 'png':
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
   elif options.outFormat == 'all':
      fig.savefig( pdf, format='pdf' )
      pdf.close()
      fig.savefig( options.out + '.png', format='png', dpi=options.dpi )
      fig.savefig( options.out + '.eps', format='eps' )
   elif options.outFormat == 'eps':
      fig.savefig( options.out + '.eps', format='eps' )
