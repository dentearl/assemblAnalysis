#Assemblathon Analysis Scripts
[Dent Earl](https://github.com/dentearl/) and [Benedict Paten](https://github.com/benedictpaten/)

Scripts to automate the creation of Tables and Figures for the Assemblathon 1 project analysis.

##Dependencies
* assemblaScripts: https://github.com/benedictpaten/assemblaScripts which has its own dependencies
* matplotlib: http://matplotlib.sourceforge.net/ for plots
* latex: http://www.latex-project.org/ for latex tables
* fltpage: http://www.ctan.org/tex-archive/macros/latex/contrib/fltpage for formatting the latex tables
* Data: http://compbio.soe.ucsc.edu/assemblathon1/
* Annotations: http://compbio.soe.ucsc.edu/assemblathon1/

##Installation
1. Download the package. 
2. <code>cd</code> into the directory.
3. Type <code>make</code>.

##Use
1. Run your assembler on the assemblathon dataset and produce a fasta file. Pick a letter to be your ID from {R, S, T, U, Y, Z}. I'll assume you use R for the rest of this list. Edit the <code>idMap</code> dict in <code>src/libGeneral.py</code> to add your ID and your name code. This is used throughout the code to place names on figures and in tables.
2. Create a sequence name map using </code>$ fastaContigHeaderMapper.py --prefix R1 --createMap R1.map < R1.rawAssembly.fa</code> .
3. Use the name map to transform all of the sequence names using <code>$ fastaContigHeaderMapper.py --map R1.map --goForward < R1.rawAssembly.fa > R1.fa</code>
4. Standardize the number of Ns in your sequences to be no greater than 25. For Assemblathon 1 25 Ns marked the difference between contigs within a scaffold. <code>$ standardizeNumNs.py --expandAt 25 < R1.fa > R1_scaffolds.fa</code> . So if your assembler used 4 Ns as a scaffold gap then you would use <code>--expandAt 4</code>, and then if there were 4 or more Ns in a row they were made to have 25 Ns. If you looked at the distribution of Ns in the resulting fasta you'd see runs of Ns of length 1, 2, 3 and 25.
5. Run RepeatMasker on the sequence using the simulated mobile element library available at http://compbio.soe.ucsc.edu/assemblathon1/ . We used <code>$ RepeatMasker -lib MELib.fa -parallel 10 -qq -xsmall -dir tempRepMaskDir/ seq.fa</code> . We used RepeatMasker v1.25.
6. Run trf on the sequence: <code>$ trf fasta.fa 2 7 7 80 10 50 2000 -m -d -h</code>. We used trf v4.00.
7. Soft mask the assembly with the union of the RepeatMasker and trf outputs.
8. Split the assembly into contigs using <code>$ splitSequenceAtNs.py --splitAt 25 < R1_scaffolds.trf.repmask.fa  > R1_contigs.trf.repmask.fa</code> . You now have two fasta files for your assembly, a scaffold and a contig file. You're now ready to run the Cactus aligner.
9. Run <code>assemblaScripts</code> on your data.
10. Create a new directory where you'd like to create perform an analysis.
11. Run <code>importCactusResult.py</code> to bring in the necessary files for both the **Scaffold** and **Contig** alignments into your new directory.
12. <code>cp</code> the <code>analysisMakefile</code> into the new directory.
13. Edit <code>analysisMakefile</code> to set the variables for <code>binDir</code>, <code>fltpage</code> and <code>simulationStatsTab</code>. The <code>binDir</code> should be the relative path to this repo's <code>bin/</code>. The <code>simulationStatsTab</code> should point to the file in <code>extra/</code>.
14. Run <code>make -f analysisMakefile</code>. **Pro tip**: The makefile was written to use <code>-j</code> for a speedup.
15. Once the make finishes, all of the results are stored in <code>publication/</code>
16. _There is no step sixteen!_
