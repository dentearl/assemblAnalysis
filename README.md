#Assemblathon Analysis Scripts
[Dent Earl](https://github.com/dentearl/) and [Benedict Paten](https://github.com/benedictpaten/)

Scripts to automate the creation of Tables and Figures for the Assemblathon 1 analysis.

##Installation and use
1. Run <code>assemblaScripts</code> on some data.
2. Create a new directory where you'd like to create perform an analysis.
3. Run <code>importCactusResult.py</code> to bring in the necessary files for both the **Scaffold** and **Contig** alignments into your new directory.
4. <code>cp</code> the <code>analysisMakefile</code> into the new directory.
5. Edit <code>analysisMakefile</code> to set the variables for <code>binDir</code>, <code>fltpage</code> and <code>simulationStatsTab</code>. The <code>binDir</code> should be the relative path to this repo's <code>bin/</code>. The <code>simulationStatsTab</code> should point to the file in <code>extra/</code>.
6. Run <code>make -f analysisMakefile</code>. **Pro tip**: The makefile was written to use <code>-j</code> for a speedup.
7. Once the make finishes, all of the results are stored in <code>publication/</code>
8. _There is no step eight!_

##Dependencies
* assemblaScripts: https://github.com/benedictpaten/assemblaScripts which has its own dependencies
* matplotlib: http://matplotlib.sourceforge.net/ for plots
* latex: http://www.latex-project.org/ for latex tables
* fltpage: http://www.ctan.org/tex-archive/macros/latex/contrib/fltpage for formatting the latex tables
* Data: http://hgwdev.cse.ucsc.edu/~benedict/code/Assemblathon_1.html
* Annotations: I'm finding a permanent home for these, will update this soon.
