SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
# dent earl dearl (a) soe ucsc edu
# 11 April 2011
##############################
binDir:=../../code/trunk/bin
projectDir:=.
resultsDir:=${projectDir}/publication
fltpage:=/hive/users/dearl/assemblathon/latexSupplemental/trunk/fltpage.sty
chrNames:= chr0,chr1,chr2
hap1ChrLengths:= 76252953,18509915,17699484
hap2ChrLengths:= 76285600,18539192,17710169
bacChrLength:= 4747819
####################
tablesDir:=${resultsDir}/tables
copyNumberPlotsDir:=${resultsDir}/copyNumberPlots
linkagePlotsDir:=${resultsDir}/linkagePlots
coverageSFDir:=${resultsDir}/coverageStackFillPlots
hap1WPSDir:=${projectDir}/picklesMafsScaffoldHap1
hap1WPCDir:=${projectDir}/picklesMafsContigHap1
hap1AWPDir:=${projectDir}/picklesAnnotsHap1
hap2WPSDir:=${projectDir}/picklesMafsScaffoldHap2
hap2WPCDir:=${projectDir}/picklesMafsContigHap2
hap2AWPDir:=${projectDir}/picklesAnnotsHap2
bacWPSDir:=${projectDir}/picklesMafsScaffoldBac
bacWPCDir:=${projectDir}/picklesMafsContigBac
bacAWPDir:=${projectDir}/picklesAnnotsBac
rankingsDir:=${resultsDir}/rankings
rankingsFiles:=$(foreach r,N50_CPN50.tab N50_NA50.tab N50_SPN50.tab copyNumberErrors.tab coverageTotal.tab linkageRanks.tab structuralContigPathErrors.tab substitutionErrors.tab,${rankingsDir}/$r)
hap1CommonCoverageArgs:= --annotPickleDir ${hap1AWPDir} --mafPickleDir ${hap1WPSDir} --chrNames ${chrNames} --chrLengths ${hap1ChrLengths} --referenceGenome hap1 --outFormat all --annotationOrder CDS,UTR,NXE,NGE,island,repeat --annotationCeilings 2000,2000,4000,4200,2000,8500 --subsetFile ${rankingsDir}/subsetTeamTopEntry
hap1CommonCoverageArgsNoSubset:= --annotPickleDir ${hap1AWPDir} --mafPickleDir ${hap1WPSDir} --chrNames ${chrNames} --chrLengths ${hap1ChrLengths} --referenceGenome hap1 --outFormat all --annotationOrder CDS,UTR,NXE,NGE,island,repeat --annotationCeilings 2000,2000,4000,4200,2000,8500
hap2CommonCoverageArgs:= --annotPickleDir ${hap2AWPDir} --mafPickleDir ${hap2WPSDir} --chrNames ${chrNames} --chrLengths ${hap2ChrLengths} --referenceGenome hap2 --outFormat all --annotationOrder CDS,UTR,NXE,NGE,island,repeat --annotationCeilings 2000,2000,4000,4200,2000,8500 --subsetFile ${rankingsDir}/subsetTeamTopEntry
hap2CommonCoverageArgsNoSubset:= --annotPickleDir ${hap2AWPDir} --mafPickleDir ${hap2WPSDir} --chrNames ${chrNames} --chrLengths ${hap2ChrLengths} --referenceGenome hap2 --outFormat all --annotationOrder CDS,UTR,NXE,NGE,island,repeat --annotationCeilings 2000,2000,4000,4200,2000,8500
bacCommonCoverageArgs:= --annotPickleDir ${bacAWPDir} --mafPickleDir ${bacWPSDir} --chrNames bac --chrLengths ${bacChrLength} --referenceGenome bac --outFormat all --annotationOrder repeat --annotationCeilings 50 --subsetFile ${rankingsDir}/subsetTeamTopEntry
bacCommonCoverageArgsNoSubset:= --annotPickleDir ${bacAWPDir} --mafPickleDir ${bacWPSDir} --chrNames bac --chrLengths ${bacChrLength} --referenceGenome bac --outFormat all --annotationOrder repeat --annotationCeilings 50
##############################

.SECONDARY: # leave this blank to force make to keep intermediate files
.PHONY: all clean

all: mafWigPicks annotationWigPicks copyNumberPlots linkagePlots covStackFillPlots ${resultsDir}/facets/n50statsFacets.pdf ${resultsDir}/facets/linkageStatsFacets.pdf ${resultsDir}/facets/copyNumberFacetLinear.pdf ${resultsDir}/facets/subStatsFacets.pdf ${tablesDir}/rankingSummaryTable.pdf ${tablesDir}/rankingSummaryTableTeams.pdf ${tablesDir}/contigPathStatsTable.pdf ${tablesDir}/coverageRankingTableTeams.pdf ${tablesDir}/coverageRankingTable.pdf ${tablesDir}/fltpage.sty ${rankingsDir}/subsetTeamTopEntry
	@echo Work Complete.

clean:
	rm -rf ${resultsDir}

${tablesDir}/fltpage.sty: ${fltpage}
	mkdir -p $(dir $@)
	cp $< $@.tmp
	mv $@.tmp $@

# maf wigPickle files (for coverage stack fill plots)
mafWigPicks: contigMafWigPicks scaffoldMafWigPicks
contigMafWigPicks: hap1ContigMafWigPicks hap2ContigMafWigPicks bacContigMafWigPicks
hap1ContigMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsContigs/*))),${hap1WPCDir}/hap1.$a.maf.pickle)
hap2ContigMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsContigs/*))),${hap2WPCDir}/hap2.$a.maf.pickle)
bacContigMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsContigs/*))),${bacWPCDir}/bac.$a.maf.pickle)
scaffoldMafWigPicks: hap1ScaffoldMafWigPicks hap2ScaffoldMafWigPicks bacScaffoldMafWigPicks
hap1ScaffoldMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsScaffolds/*))),${hap1WPSDir}/hap1.$a.maf.pickle)
hap2ScaffoldMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsScaffolds/*))),${hap2WPSDir}/hap2.$a.maf.pickle)
bacScaffoldMafWigPicks: $(foreach a,$(subst .maf,,$(notdir $(wildcard ${projectDir}/cactusMafsScaffolds/*))),${bacWPSDir}/bac.$a.maf.pickle)
${hap1WPSDir}/hap1.%.maf.pickle: ${projectDir}/cactusMafsScaffolds/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome hap1 --comparisonGenome $* \
		--outDir $(dir $@) --name hap1Tmp.$* --chrNames ${chrNames} --chrLengths ${hap1ChrLengths} --verify
	mv ${hap1WPSDir}/hap1Tmp.$*.maf.pickle $@
${hap1WPCDir}/hap1.%.maf.pickle: ${projectDir}/cactusMafsContigs/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome hap1 --comparisonGenome $* \
		--outDir $(dir $@) --name hap1Tmp.$* --chrNames ${chrNames} --chrLengths ${hap1ChrLengths} --verify
	mv ${hap1WPCDir}/hap1Tmp.$*.maf.pickle $@
${hap2WPSDir}/hap2.%.maf.pickle: ${projectDir}/cactusMafsScaffolds/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome hap2 --comparisonGenome $* \
		--outDir $(dir $@) --name hap2Tmp.$* --chrNames ${chrNames} --chrLengths ${hap2ChrLengths} --verify
	mv ${hap2WPSDir}/hap2Tmp.$*.maf.pickle $@
${hap2WPCDir}/hap2.%.maf.pickle: ${projectDir}/cactusMafsContigs/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome hap2 --comparisonGenome $* \
		--outDir $(dir $@) --name hap2Tmp.$* --chrNames ${chrNames} --chrLengths ${hap2ChrLengths} --verify
	mv ${hap2WPCDir}/hap2Tmp.$*.maf.pickle $@
${bacWPSDir}/bac.%.maf.pickle: ${projectDir}/cactusMafsScaffolds/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome bac --comparisonGenome $* \
		--outDir $(dir $@) --name bacTmp.$* --chrNames bac --chrLengths ${bacChrLength} --verify
	mv ${bacWPSDir}/bacTmp.$*.maf.pickle $@
${bacWPCDir}/bac.%.maf.pickle: ${projectDir}/cactusMafsContigs/%.maf ${binDir}/mafToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/mafToPlotPickles.py --maf $< --referenceGenome bac --comparisonGenome $* \
		--outDir $(dir $@) --name bacTmp.$* --chrNames bac --chrLengths ${bacChrLength} --verify
	mv ${bacWPCDir}/bacTmp.$*.maf.pickle $@

# annotation wigPickle files (for coverage stack fill plots)
annotationWigPicks:  ${hap1AWPDir}/hap1.annots.pickle ${hap2AWPDir}/hap2.annots.pickle ${bacAWPDir}/bac.annots.pickle
${hap1AWPDir}/hap1.annots.pickle: ${projectDir}/annotations/hap1.annots.trf.repmask.gff ${binDir}/gffToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/gffToPlotPickles.py --gff $< \
		--outDir $(dir $@) --prefix hap1Tmp --chrNames ${chrNames} \
		--chrLengths ${hap1ChrLengths}
	mv ${hap1AWPDir}/hap1Tmp.annots.pickle $@
${hap2AWPDir}/hap2.annots.pickle: ${projectDir}/annotations/hap2.annots.trf.repmask.gff ${binDir}/gffToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/gffToPlotPickles.py --gff $< \
		--outDir $(dir $@) --prefix hap2Tmp --chrNames ${chrNames} \
		--chrLengths ${hap2ChrLengths}
	mv ${hap2AWPDir}/hap2Tmp.annots.pickle $@
${bacAWPDir}/bac.annots.pickle: ${projectDir}/annotations/bac.annots.trf.repmask.gff ${binDir}/gffToPlotPickles.py
	mkdir -p $(dir $@)
	${binDir}/gffToPlotPickles.py --gff $< \
		--outDir $(dir $@) --prefix bacTmp --chrNames bac --chrLengths ${bacChrLength}
	mv ${bacAWPDir}/bacTmp.annots.pickle $@

# coverage stack fill plots:
covStackFillPlots: covSFPHap1 covSFPHap2 covSFPBac
covSFPHap1: $(foreach p,cactusHap1 cactusHap1Fill cactusHap1BED cactusHap1StackFillBlocks cactusHap1CP cactusHap1StackFillContigPaths cactusHap1StackFillContigs cactusHap1StackFillScaffPaths, ${coverageSFDir}/$p.pdf)
covSFPHap2: $(foreach p,cactusHap2 cactusHap2Fill cactusHap2BED cactusHap2StackFillBlocks cactusHap2CP cactusHap2StackFillContigPaths cactusHap2StackFillContigs cactusHap2StackFillScaffPaths, ${coverageSFDir}/$p.pdf)
covSFPBac: $(foreach p,cactusBac cactusBacFill cactusBacBED cactusBacStackFillBlocks cactusBacCP cactusBacStackFillContigPaths, ${coverageSFDir}/$p.pdf)
# haplotype 1
${coverageSFDir}/cactusHap1.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@)
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1Tmp.$$p ${coverageSFDir}/cactusHap1.$$p ; \
	done
${coverageSFDir}/all_cactusHap1.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@)
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1Tmp.$$p ${coverageSFDir}/all_cactusHap1.$$p ; \
	done
${coverageSFDir}/cactusHap1Fill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1FillTmp.$$p ${coverageSFDir}/cactusHap1Fill.$$p ; \
	done
${coverageSFDir}/all_cactusHap1Fill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1FillTmp.$$p ${coverageSFDir}/all_cactusHap1Fill.$$p ; \
	done
${coverageSFDir}/cactusHap1BED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1BEDTmp.$$p ${coverageSFDir}/cactusHap1BED.$$p ; \
	done
${coverageSFDir}/all_cactusHap1BED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1BEDTmp.$$p ${coverageSFDir}/all_cactusHap1BED.$$p ; \
	done
${coverageSFDir}/cactusHap1StackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1StackFillBlocksTmp.$$p ${coverageSFDir}/cactusHap1StackFillBlocks.$$p ; \
	done
${coverageSFDir}/all_cactusHap1StackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*)  ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1StackFillBlocksTmp.$$p ${coverageSFDir}/all_cactusHap1StackFillBlocks.$$p ; \
	done
${coverageSFDir}/cactusHap1CP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1CPTmp.$$p ${coverageSFDir}/cactusHap1CP.$$p ; \
	done
${coverageSFDir}/all_cactusHap1CP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*)  ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1CPTmp.$$p ${coverageSFDir}/all_cactusHap1CP.$$p ; \
	done
${coverageSFDir}/cactusHap1StackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1StackFillContigPathsTmp.$$p ${coverageSFDir}/cactusHap1StackFillContigPaths.$$p ; \
	done
${coverageSFDir}/all_cactusHap1StackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*)  ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1StackFillContigPathsTmp.$$p ${coverageSFDir}/all_cactusHap1StackFillContigPaths.$$p ; \
	done
${coverageSFDir}/cactusHap1StackFillContigs.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ContigMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py --annotPickleDir ${hap1AWPDir}/ \
		--mafPickleDir ${hap1WPCDir} --chrNames ${chrNames} \
		--chrLengths ${hap1ChrLengths} --referenceGenome hap1 --outFormat all \
		--annotationOrder CDS,UTR,NXE,NGE,island,repeat \
		--annotationCeilings 2000,2000,4000,4200,2000,8500 --subsetFile ${rankingsDir}/subsetTeamTopEntry \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigs
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1StackFillContigsTmp.$$p ${coverageSFDir}/cactusHap1StackFillContigs.$$p ; \
	done
${coverageSFDir}/all_cactusHap1StackFillContigs.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${hap1AWPDir}/hap1.annots.pickle hap1ContigMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py --annotPickleDir ${hap1AWPDir}/ \
		--mafPickleDir ${hap1WPCDir} --chrNames ${chrNames} \
		--chrLengths ${hap1ChrLengths} --referenceGenome hap1 --outFormat all \
		--annotationOrder CDS,UTR,NXE,NGE,island,repeat \
		--annotationCeilings 2000,2000,4000,4200,2000,8500 \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigs
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1StackFillContigsTmp.$$p ${coverageSFDir}/all_cactusHap1StackFillContigs.$$p ; \
	done
${coverageSFDir}/cactusHap1StackFillScaffPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillScaffPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap1StackFillScaffPathsTmp.$$p ${coverageSFDir}/cactusHap1StackFillScaffPaths.$$p ; \
	done
${coverageSFDir}/all_cactusHap1StackFillScaffPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap1WPSDir}/*) ${hap1AWPDir}/hap1.annots.pickle hap1ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap1CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillScaffPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap1StackFillScaffPathsTmp.$$p ${coverageSFDir}/all_cactusHap1StackFillScaffPaths.$$p ; \
	done
# haplotype 2
${coverageSFDir}/cactusHap2.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) 
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2Tmp.$$p ${coverageSFDir}/cactusHap2.$$p ; \
	done
${coverageSFDir}/all_cactusHap2.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*)  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) 
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2Tmp.$$p ${coverageSFDir}/all_cactusHap2.$$p ; \
	done
${coverageSFDir}/cactusHap2Fill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2FillTmp.$$p ${coverageSFDir}/cactusHap2Fill.$$p ; \
	done
${coverageSFDir}/all_cactusHap2Fill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2FillTmp.$$p ${coverageSFDir}/all_cactusHap2Fill.$$p ; \
	done
${coverageSFDir}/cactusHap2BED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2BEDTmp.$$p ${coverageSFDir}/cactusHap2BED.$$p ; \
	done
${coverageSFDir}/all_cactusHap2BED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2BEDTmp.$$p ${coverageSFDir}/all_cactusHap2BED.$$p ; \
	done
${coverageSFDir}/cactusHap2StackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2StackFillBlocksTmp.$$p ${coverageSFDir}/cactusHap2StackFillBlocks.$$p ; \
	done
${coverageSFDir}/all_cactusHap2StackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2StackFillBlocksTmp.$$p ${coverageSFDir}/all_cactusHap2StackFillBlocks.$$p ; \
	done
${coverageSFDir}/cactusHap2CP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2CPTmp.$$p ${coverageSFDir}/cactusHap2CP.$$p ; \
	done
${coverageSFDir}/all_cactusHap2CP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2CPTmp.$$p ${coverageSFDir}/all_cactusHap2CP.$$p ; \
	done
${coverageSFDir}/cactusHap2StackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2StackFillContigPathsTmp.$$p ${coverageSFDir}/cactusHap2StackFillContigPaths.$$p ; \
	done
${coverageSFDir}/all_cactusHap2StackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2StackFillContigPathsTmp.$$p ${coverageSFDir}/all_cactusHap2StackFillContigPaths.$$p ; \
	done
${coverageSFDir}/cactusHap2StackFillContigs.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ContigMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py --annotPickleDir ${hap2AWPDir}/ \
		--mafPickleDir ${hap2WPCDir} --chrNames ${chrNames} \
		--chrLengths ${hap2ChrLengths} --referenceGenome hap2 --outFormat all \
		--annotationOrder CDS,UTR,NXE,NGE,island,repeat \
		--annotationCeilings 2000,2000,4000,4200,2000,8500 --subsetFile ${rankingsDir}/subsetTeamTopEntry \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigs
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2StackFillContigsTmp.$$p ${coverageSFDir}/cactusHap2StackFillContigs.$$p ; \
	done
${coverageSFDir}/all_cactusHap2StackFillContigs.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ContigMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py --annotPickleDir ${hap2AWPDir}/ \
		--mafPickleDir ${hap2WPCDir} --chrNames ${chrNames} \
		--chrLengths ${hap2ChrLengths} --referenceGenome hap2 --outFormat all \
		--annotationOrder CDS,UTR,NXE,NGE,island,repeat \
		--annotationCeilings 2000,2000,4000,4200,2000,8500 \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigs
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2StackFillContigsTmp.$$p ${coverageSFDir}/all_cactusHap2StackFillContigs.$$p ; \
	done
${coverageSFDir}/cactusHap2StackFillScaffPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${rankingsDir}/subsetTeamTopEntry  ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillScaffPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusHap2StackFillScaffPathsTmp.$$p ${coverageSFDir}/cactusHap2StackFillScaffPaths.$$p ; \
	done
${coverageSFDir}/all_cactusHap2StackFillScaffPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${hap2WPSDir}/*) ${hap2AWPDir}/hap2.annots.pickle hap2ScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${hap2CommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillScaffPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusHap2StackFillScaffPathsTmp.$$p ${coverageSFDir}/all_cactusHap2StackFillScaffPaths.$$p ; \
	done
# bacterial
${coverageSFDir}/cactusBac.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) 
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacTmp.$$p ${coverageSFDir}/cactusBac.$$p ; \
	done
${coverageSFDir}/all_cactusBac.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) 
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacTmp.$$p ${coverageSFDir}/all_cactusBac.$$p ; \
	done
${coverageSFDir}/cactusBacFill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacFillTmp.$$p ${coverageSFDir}/cactusBacFill.$$p ; \
	done
${coverageSFDir}/all_cactusBacFill.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --fill
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacFillTmp.$$p ${coverageSFDir}/all_cactusBacFill.$$p ; \
	done
${coverageSFDir}/cactusBacBED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacBEDTmp.$$p ${coverageSFDir}/cactusBacBED.$$p ; \
	done
${coverageSFDir}/all_cactusBacBED.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --edgeErrorCeiling 300
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacBEDTmp.$$p ${coverageSFDir}/all_cactusBacBED.$$p ; \
	done
${coverageSFDir}/cactusBacStackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacStackFillBlocksTmp.$$p ${coverageSFDir}/cactusBacStackFillBlocks.$$p ; \
	done
${coverageSFDir}/all_cactusBacStackFillBlocks.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillBlocks
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacStackFillBlocksTmp.$$p ${coverageSFDir}/all_cactusBacStackFillBlocks.$$p ; \
	done
${coverageSFDir}/cactusBacCP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacCPTmp.$$p ${coverageSFDir}/cactusBacCP.$$p ; \
	done
${coverageSFDir}/all_cactusBacCP.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --contigPathErrorDensity --contigPathEdgeDensity \
		--edgeErrorCeiling 30
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacCPTmp.$$p ${coverageSFDir}/all_cactusBacCP.$$p ; \
	done
${coverageSFDir}/cactusBacStackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${rankingsDir}/subsetTeamTopEntry ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgs} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/cactusBacStackFillContigPathsTmp.$$p ${coverageSFDir}/cactusBacStackFillContigPaths.$$p ; \
	done
${coverageSFDir}/all_cactusBacStackFillContigPaths.pdf: ${binDir}/plotPicklesToPlot.py $(wildcard ${bacWPSDir}/*) ${bacAWPDir}/bac.annots.pickle bacScaffoldMafWigPicks
	mkdir -p $(dir $@)
	${binDir}/plotPicklesToPlot.py ${bacCommonCoverageArgsNoSubset} \
		--out $(subst .pdf,Tmp.pdf,$@) --stackFillContigPaths
	for p in png eps pdf; do \
		mv ${coverageSFDir}/all_cactusBacStackFillContigPathsTmp.$$p ${coverageSFDir}/all_cactusBacStackFillContigPaths.$$p ; \
	done

# n50 facetPlot
${resultsDir}/facets/n50statsFacets.pdf: ${binDir}/createN50StatsPlot.py $(wildcard ${projectDir}/contigPathStatsScaffolds/*) ${rankingsDir}/subsetTeamTopEntry
	mkdir -p $(dir $@)
	${binDir}/createN50StatsPlot.py --contigPathStatsDir ${projectDir}/contigPathStatsScaffolds/ \
		--out ${resultsDir}/facets/n50statsFacetsTmp.pdf --outFormat pdf \
	--subsetFile ${rankingsDir}/subsetTeamTopEntry
	mv ${resultsDir}/facets/n50statsFacetsTmp.pdf $@

# linkage individual plots
${linkagePlotsDir}/%.linkageStats.pdf: ${projectDir}/linkageStatsScaffolds/%.linkageStats.xml ${binDir}/createLinkageStatsPlot.py
	mkdir -p $(dir $@)
	${binDir}/createLinkageStatsPlot.py --title "$* Copy Number Statistics" --out ${linkagePlotsDir}/$*.linkageStats.tmp.pdf --outFormat pdf $<
	mv ${linkagePlotsDir}/$*.linkageStats.tmp.pdf $@
# linkage individual plots
linkagePlots: $(foreach f,$(wildcard ${projectDir}/linkageStatsScaffolds/*.xml), $(join ${linkagePlotsDir}/, $(subst xml,pdf, $(notdir $f)) ))

# linkage facet plot
${resultsDir}/facets/linkageStatsFacets.pdf: ${binDir}/createLinkageStatsPlot.py $(wildcard ${projectDir}/linkageStatsScaffolds/*) ${rankingsDir}/subsetTeamTopEntry
	mkdir -p $(dir $@)
	${binDir}/createLinkageStatsPlot.py --out ${resultsDir}/facets/linkageStatsFacetsTmp --outFormat pdf --title 'Linkage Statistics' ${projectDir}/linkageStatsScaffolds/* --subsetFile ${rankingsDir}/subsetTeamTopEntry
	mv ${resultsDir}/facets/linkageStatsFacetsTmp.pdf $@

# copy number individual plots
${copyNumberPlotsDir}/%.pdf: ${projectDir}/copyNumberPathStatsScaffolds/%.xml ${binDir}/createCopyNumberStatsPlot.py
	mkdir -p $(dir $@)
	${binDir}/createCopyNumberStatsPlot.py --title "$* Copy Number Statistics" --out ${copyNumberPlotsDir}/$*.tmp.pdf --outFormat pdf $<
	mv ${copyNumberPlotsDir}/$*.tmp.pdf $@
# copy number individual plots
copyNumberPlots: $(foreach f,$(wildcard ${projectDir}/copyNumberPathStatsScaffolds/*.xml), $(join ${copyNumberPlotsDir}/, $(subst xml,pdf, $(notdir $f)) ))

# copy number facet plot
${resultsDir}/facets/copyNumberFacetLinear.pdf: ${binDir}/createCopyNumberStatsFacetedPlot.py $(wildcard ${projectDir}/copyNumberPathStatsScaffolds/*) ${rankingsDir}/subsetTeamTopEntry
	mkdir -p $(dir $@)
	${binDir}/createCopyNumberStatsFacetedPlot.py --dir ${projectDir}/copyNumberPathStatsScaffolds/ \
		--subsetFile ${rankingsDir}/subsetTeamTopEntry --out ${resultsDir}/facets/copyNumberFacetLinearTmp --outFormat pdf
	mv ${resultsDir}/facets/copyNumberFacetLinearTmp.pdf $@

# substitution stats facet plot
${resultsDir}/facets/subStatsFacets.pdf: ${binDir}/createSubStatsMegaPlot.py $(wildcard ${projectDir}/subStatsScaffolds/*) ${rankingsDir}/subsetTeamTopEntry
	mkdir -p $(dir $@)
	${binDir}/createSubStatsMegaPlot.py --subStatsDir ${projectDir}/subStatsScaffolds/ --subsetFile \
		${rankingsDir}/subsetTeamTopEntry --out ${resultsDir}/facets/subStatsFacetsTmp --outFormat pdf
	mv ${resultsDir}/facets/subStatsFacetsTmp.pdf $@

# hap1 coverage rankings
${rankingsDir}/coverageHap1.txt: ${binDir}/extractPlotPickleCounts.py $(wildcard ${hap1WPSDir}/*)
	mkdir -p $(dir $@)
	for f in ${hap1WPSDir}/*.pickle ; do \
		a=$$(basename $$f | perl -ple 's/^hap\d\.(.{2,3})\..*/$$1/;') ; \
		echo $$a $$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${hap1WPSDir}/hap1.$$a.* | awk '{sum += $$1} END{ print sum }' ) \
		$$( echo $$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${hap1WPSDir}/hap1.$$a.* | awk '{sum += $$1} END{ print sum }' )/112462352 | bc -l ) >> $@.tmp ; \
		done
	mv $@.tmp $@
# hap2 coverage rankings
${rankingsDir}/coverageHap2.txt: ${binDir}/extractPlotPickleCounts.py $(wildcard ${hap2WPSDir}/*)
	mkdir -p $(dir $@)
	for f in ${hap2WPSDir}/*.pickle ; do \
		a=$$(basename $$f | perl -ple 's/^hap\d\.(.{2,3})\..*/$$1/;') ; \
		echo $$a $$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${hap2WPSDir}/hap2.$$a.* | awk '{sum += $$1} END{ print sum }' ) \
		$$( echo $$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${hap2WPSDir}/hap2.$$a.* | awk '{sum += $$1} END{ print sum }' )/112462352 | bc -l ) >> $@.tmp ; \
		done
	mv $@.tmp $@
# bac coverage rankings
${rankingsDir}/coverageBac.txt: ${binDir}/extractPlotPickleCounts.py $(wildcard ${bacWPSDir}/*)
	mkdir -p $(dir $@)
	for f in ${bacWPSDir}/*.pickle ; do \
		a=$$(basename $$f | perl -ple 's/^bac\.(.{2,3})\..*/$$1/;') ; \
		echo $$a $$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${bacWPSDir}/bac.$$a.* | awk '{sum += $$1} END{ print sum }' ) $$( echo \
		$$( ${binDir}/extractPlotPickleCounts.py --key columnsInBlocks \
		${bacWPSDir}/bac.$$a.* | awk '{sum += $$1} END{ print sum }' )/112462352 | bc -l ) >> $@.tmp ; \
		done
	mv $@.tmp $@
# total coverage ranking
${rankingsDir}/coverageTotal.tab: ${rankingsDir}/coverageHap1.txt ${rankingsDir}/coverageHap2.txt ${rankingsDir}/coverageBac.txt
	mkdir -p $(dir $@)
	paste $^ > ${rankingsDir}/coverageTotal.txt
	echo -e '#Assembly\tHap Total\tHap 1\tHap 2\tBac' > ${rankingsDir}/coverageTotalSummary.tab
	cat ${rankingsDir}/coverageTotal.txt | awk '{printf "%s\t%f\t%f\t%f\t%f\n", $$1, (($$2+$$5)/224997313), $$2/112462352, $$5/112534961, $$8/4747819}' >> ${rankingsDir}/coverageTotalSummary.tab
	awk 'NR==1; NR>1 {print $$0 | "sort -rn -k2"}' ${rankingsDir}/coverageTotalSummary.tab | awk 'NR==1; NR>1 {printf "%s\t%.4f\t%.4f\t%.4f\t%.4f\n", $$1, $$2, $$3, $$4, $$5 }'  > $@.tmp
	rm ${rankingsDir}/coverageTotalSummary.tab
	mv $@.tmp $@

# coverage rankings table, all
${tablesDir}/coverageRankingTable.pdf: ${rankingsDir}/coverageTotal.tab ${binDir}/tab2tex.py ${binDir}/latexWrapper.py ${tablesDir}/fltpage.sty
	mkdir -p $(dir $@)
	cat ${rankingsDir}/coverageTotal.tab | ${binDir}/tab2tex.py | \
		${binDir}/latexWrapper.py > ${tablesDir}/coverageRankingTableTmp.tex
	pushd ${tablesDir} && latex coverageRankingTableTmp.tex && latex coverageRankingTableTmp.tex \
		latex coverageRankingTableTmp.tex && dvips -Ppdf coverageRankingTableTmp.dvi && \
		ps2pdf coverageRankingTableTmp.ps && popd
	mv ${tablesDir}/coverageRankingTableTmp.pdf $@
# coverage rankings table, subset
${tablesDir}/coverageRankingTableTeams.pdf: ${rankingsDir}/coverageTotal.tab ${binDir}/tab2tex.py ${binDir}/latexWrapper.py ${binDir}/subsetTabs.py ${rankingsDir}/subsetTeamTopEntry ${tablesDir}/fltpage.sty
	mkdir -p $(dir $@)
	cat ${rankingsDir}/coverageTotal.tab | ${binDir}/subsetTabs.py --subsetFile \
		${rankingsDir}/subsetTeamTopEntry | ${binDir}/tab2tex.py | \
		${binDir}/latexWrapper.py > ${tablesDir}/coverageRankingTableTeamsTmp.tex
	pushd ${tablesDir} && latex coverageRankingTableTeamsTmp.tex && latex coverageRankingTableTeamsTmp.tex \
		latex coverageRankingTableTeamsTmp.tex && dvips -Ppdf coverageRankingTableTeamsTmp.dvi && \
		ps2pdf coverageRankingTableTeamsTmp.ps && popd
	mv ${tablesDir}/coverageRankingTableTeamsTmp.pdf $@

# Contig Path (haplotype path) / structural error table
${tablesDir}/contigPathStatsTable.pdf: ${binDir}/createContigPathStatsTable.py ${binDir}/latexWrapper.py ${rankingsDir}/subsetTeamTopEntry $(wildcard ${projectDir}/contigPathStatsScaffolds/*) ${tablesDir}/fltpage.sty
	mkdir -p $(dir $@)
	${binDir}/createContigPathStatsTable.py --subsetFile ${rankingsDir}/subsetTeamTopEntry \
		 --contigPathStatsDir ${projectDir}/contigPathStatsScaffolds/ | ${binDir}/latexWrapper.py > \
		${tablesDir}/contigPathStatsTableTmp.tex
	pushd ${tablesDir} && latex contigPathStatsTableTmp.tex && latex contigPathStatsTableTmp.tex \
		latex contigPathStatsTableTmp.tex && dvips -Ppdf contigPathStatsTableTmp.dvi && \
		ps2pdf contigPathStatsTableTmp.ps && popd
	mv ${tablesDir}/contigPathStatsTableTmp.pdf $@

# linkage rankings
${rankingsDir}/linkageRanks.tab: ${binDir}/createLinkageStatsPlot.py $(wildcard ${projectDir}/linkageStatsScaffolds/*)
	mkdir -p $(dir $@)
	${binDir}/createLinkageStatsPlot.py --outputRanks ${projectDir}/linkageStatsScaffolds/* > $@.tmp
	mv $@.tmp $@

# structural (contig path / haplotype path) errors rankings
${rankingsDir}/structuralContigPathErrors.tab: ${binDir}/createContigPathStatsTable.py $(wildcard ${projectDir}/contigPathStatsCaffolds/*)
	mkdir -p $(dir $@)
	${binDir}/createContigPathStatsTable.py --contigPathStatsDir ${projectDir}/contigPathStatsScaffolds/ --outputRanks > $@.tmp
	mv $@.tmp $@

# base substitution errors rankings
${rankingsDir}/substitutionErrors.tab: ${binDir}/createSubStatsMegaPlot.py $(wildcard ${projectDir}/subStatsScaffolds/*)
	mkdir -p $(dir $@)
	${binDir}/createSubStatsMegaPlot.py --subStatsDir ${projectDir}/subStatsScaffolds/ --outputRanks | \
		awk 'NR==1; NR>1 {printf "%s\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\t%.3e\n",$$1, $$2, $$3, $$4, $$5, $$6, $$7, $$8, $$9}' > $@.tmp
	mv $@.tmp $@

# copy number errors rankings
${rankingsDir}/copyNumberErrors.tab: ${binDir}/createCopyNumberStatsFacetedPlot.py $(wildcard ${projectDir}/copyNumberStatsScaffolds/*)
	mkdir -p $(dir $@)
	${binDir}/createCopyNumberStatsFacetedPlot.py --dir ${projectDir}/copyNumberPathStatsScaffolds --outputRanks > $@.tmp
	mv $@.tmp $@

# all n50 stats
${rankingsDir}/N50All.txt: ${binDir}/createN50StatsPlot.py $(wildcard ${projectDir}/contigPathStatsScaffolds/*)
	mkdir -p $(dir $@)
	${binDir}/createN50StatsPlot.py --contigPathStatsDir ${projectDir}/contigPathStatsScaffolds/ --outputRanks --sortOn contigNA50 > $@.tmp
	mv $@.tmp $@

# NA50 rankings
${rankingsDir}/N50_NA50.tab: ${rankingsDir}/N50All.txt
	mkdir -p $(dir $@)
	awk -F='\t' '{print $$1,"\t",$$2}' $< | awk 'NR==1; NR>1 {print $$0 | "sort -rn -k2"}' > $@.tmp
	mv $@.tmp $@

# contig path N50 rankings
${rankingsDir}/N50_CPN50.tab: ${rankingsDir}/N50All.txt
	mkdir -p $(dir $@)
	awk -F='\t' '{print $$1,"\t",$$3}' $< | awk 'NR==1; NR>1 {print $$0 | "sort -rn -k2"}' > $@.tmp
	mv $@.tmp $@

# scaffold path N50 rankings
${rankingsDir}/N50_SPN50.tab: ${rankingsDir}/N50All.txt
	mkdir -p $(dir $@)
	awk -F='\t' '{print $$1,"\t",$$4}' $< | awk 'NR==1; NR>1 {print $$0 | "sort -rn -k2"}' > $@.tmp
	mv $@.tmp $@

# The top assembly from each team
${rankingsDir}/subsetTeamTopEntry: ${binDir}/summarizeRankings.py ${rankingsFiles} ${binDir}/summarizedRankingToTeamTop.py
	mkdir -p $(dir $@)
	${binDir}/summarizeRankings.py ${rankingsDir}/*tab | ${binDir}/summarizedRankingToTeamTop.py | \
		 sed '/^#/d' | sed 's/.linkageStats.xml//' | awk '{print $$1}' > $@.tmp
	mv $@.tmp $@

# rankings table, all
${tablesDir}/rankingSummaryTable.pdf: ${binDir}/summarizeRankings.py ${rankingsFiles} ${binDir}/tab2tex.py ${tablesDir}/fltpage.sty
	mkdir -p $(dir $@)
	${binDir}/summarizeRankings.py ${rankingsDir}/*tab --retainValues | awk 'NR==1; NR>1 {print $$0 | "sort -n -k2"}' | \
	 	sed '/^$$/d' | ${binDir}/tab2tex.py > ${tablesDir}/rankingSummaryTable.tex
	pushd ${tablesDir} && ../../${binDir}/latexWrapper.py < rankingSummaryTable.tex > \
		rankingSummaryTableWrapped.tex && latex rankingSummaryTableWrapped.tex && latex rankingSummaryTableWrapped.tex \
		latex rankingSummaryTableWrapped.tex && dvips -Ppdf rankingSummaryTableWrapped.dvi && \
		ps2pdf rankingSummaryTableWrapped.ps && popd
	mv ${tablesDir}/rankingSummaryTableWrapped.pdf $@
# rankings table, subset
${tablesDir}/rankingSummaryTableTeams.pdf: ${binDir}/summarizeRankings.py ${rankingsFiles} ${binDir}/summarizedRankingToTeamTop.py ${binDir}/tab2tex.py ${rankingsDir}/subsetTeamTopEntry ${tablesDir}/fltpage.sty
	mkdir -p $(dir $@)
	${binDir}/summarizeRankings.py ${rankingsDir}/*tab --retainValues | ${binDir}/summarizedRankingToTeamTop.py | \
		awk 'NR==1; NR>1 {print $$0 | "sort -n -k2"}' | sed '/^$$/d' | ${binDir}/tab2tex.py > \
		${tablesDir}/rankingSummaryTableTeams.tex
	pushd ${tablesDir} && ../../${binDir}/latexWrapper.py < rankingSummaryTableTeams.tex > \
		rankingSummaryTableTeamsWrapped.tex && latex rankingSummaryTableTeamsWrapped.tex && latex rankingSummaryTableTeamsWrapped.tex \
		latex rankingSummaryTableTeamsWrapped.tex && dvips -Ppdf rankingSummaryTableTeamsWrapped.dvi && \
		ps2pdf rankingSummaryTableTeamsWrapped.ps && popd
	mv ${tablesDir}/rankingSummaryTableTeamsWrapped.pdf $@

