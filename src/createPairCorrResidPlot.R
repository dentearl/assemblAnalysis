#!/usr/bin/env Rscript
# runMyPairs.R
# dent earl, dearl (a) soe ucsc edu
# 4 May 2011
#
# an R script to produce two plots, correlations.pdf and residuals.pdf
# for all pairs of metrics related to the assemblathon project.
#
# argument 1 should be the data.tab file,
# argument 2 should be my version of pairs()
# argument 3 should be the directory to write the plots to
#
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
#########
args = commandArgs(TRUE)

d.raw = read.table(args[1], header=TRUE)
source( args[2] )

d.raw$S_N50 = log(d.raw$S_N50)
d.raw$C_N50 = log(d.raw$C_N50)
d.raw$N50_CPNG50 = log(d.raw$N50_CPNG50)
d.raw$N50_SPNG50 = log(d.raw$N50_SPNG50)
d.raw$contiguousRanks = log(d.raw$contiguousRanks)
d.raw$substitutionErrors = log(d.raw$substitutionErrors)
d.raw$substitutionErrors[d.raw$substitutionErrors == -Inf] = NA
d.raw$structuralContigPathErrors = log(d.raw$structuralContigPathErrors)
d.raw$copyNumberErrors = log(d.raw$copyNumberErrors)

panel.hist = function(x, ...){ 
   usr = par("usr"); on.exit(par(usr))
   par(usr = c(usr[1:2], 0, 1.5) )
   h = hist(x, plot = FALSE)
   breaks = h$breaks; nB = length(breaks)
   y = h$counts; y = y/max(y)
   rect(breaks[-nB], 0, breaks[-1], y, col="gray85", lty=0, ...)
}

panel.lm = function( x, y, ...){
   usr = par("usr"); on.exit(par(usr))
   points(x, y, pch=21)
   abline(lm(y ~ x), col='red', lwd=1)
}

panel.residuals = function( x, y, ...){
   usr = par("usr");
   fit = lm( y ~ x )
   qqVals = qqnorm(residuals(fit), plot.it=FALSE)
   usr = par(usr = c(min(qqVals$x), max(qqVals$x), min(qqVals$y), max(qqVals$y)))
   points(qqVals$x, qqVals$y)
   qqline(residuals(fit), col='red')
   on.exit(par(usr))
}

panel.cor = function(x, y, digits=2, cex.cor){
   usr = par("usr"); on.exit(par(usr))
   par(usr = c(0, 1, 0, 1))
   r = abs(cor(x, y, use='pairwise', method='pearson'))
   test = cor.test(x,y) 
   # bonferroni correction
   if (test$p.value * choose(length(x), 2) > 1.0) {
      test.p = 1.0
   }else{
      test.p = test$p.value * choose(length(x), 2)
   }
   sigStr = symnum(test.p, corr = FALSE, na = FALSE, 
              cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
              symbols = c("***", "**", "*", ".", " ")) 
   txt = format(c(r, 0.123456789), digits=digits)[1]
   if(missing(cex.cor)) {
      if (3.0 * r > 0.95){
         cex = 3.0*r
      }else{
         cex = 0.95
      }
   }
   text(0.5, 0.5, txt, cex = cex)
   text(.8, .8, sigStr, cex=cex, col=2) 
}

pdf(paste( args[3], 'correlations.pdf', sep=''), height=10, width=12)
pairs(d.raw, lower.panel=panel.lm, upper.panel=panel.cor, diag.panel=panel.hist)
out = dev.off()
pdf(paste( args[3], 'residuals.pdf', sep=''), height=10, width=12)
pairs(d.raw, lower.panel=panel.residuals, upper.panel=panel.cor, diag.panel=panel.hist, printLabels=FALSE)
out = dev.off()
