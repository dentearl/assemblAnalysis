SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
# dent earl dearl (a) soe ucsc edu
# 26 Jan 2011
##############################
# EDIT THESE:
PROJECT_DIR:=/hive/users/dearl/assemblathon/prod
BIN_DIR:=/hive/users/dearl/assemblathon/code/trunk/bin
##############################
# It is assumed that all assemblies will be in fasta format, gzipped
# and named according to:
#            mySweetAssembly.fa.gz
# All assemblies go in the ${PROJECT_DIR}/archives directory.
# Assumed initial directory structure:
# ${PROJECT_DIR}/
#     |------ /haplotypes/hap1.trf.repmask.2bit  # hap1, trf'd and repeatMasked
#     |------ /haplotypes/hap2.trf.repmask.2bit  # hap2. 
#     |------ /archives/*    # all *.fa.gz assembly entrants
#     |------ /MElibrary/MELib.fa    # The mobile element library to be used in repeat masking
# The rule structure below explicity contains hap1 and hap2 rules to make
# the Makefile simpler, though less general.
#
# DEPENDENCIES
#   RepeatMasker
#   trf - tandem repeat finder
#   partitionSequence.pl
#   constructLiftFile.pl
#   runLastzChain.sh
#   runPara.sh
#   rsync
#   removeEmptyContigs.py
#
##############################
# DO NOT EDIT BELOW THIS LINE
.SECONDARY: # leave this blank to force make to keep intermediate files
HAPS_DIR:=${PROJECT_DIR}/haplotypes
SUBMISSION_DIR:=${PROJECT_DIR}/submissions
RAW_DIR:=${PROJECT_DIR}/archives
CHAINS_DIR:=${PROJECT_DIR}/chains
CHAINSCRIPTS_DIR:=${PROJECT_DIR}/chainScripts
MAFS_DIR:=${PROJECT_DIR}/mafs
ASSEMBLIES_DIR:=${PROJECT_DIR}/assemblies
REPMASK_DIR:=${PROJECT_DIR}/repeatMasking
TRF_DIR:=${PROJECT_DIR}/tandemRepeatFinder
MELIB:=${PROJECT_DIR}/MELibrary/MELib.fa
ASSEMBLIES:=$(patsubst %.fa.gz,%,$(notdir $(wildcard ${RAW_DIR}/*contigs.fa.gz)))
HAPLOTYPES:=$(patsubst %.trf.repmask.2bit,%,$(notdir $(wildcard ${HAPS_DIR}/*.trf.repmask.2bit)))
UNIQ:=$(strip $(shell date '+%s' | perl -ple 's/^\d{6}//;')) # last 4 digits of time.
REPEATMASKER:=/scratch/data/RepeatMasker/RepeatMasker
HOST=$(shell hostname)
PPID=$(shell echo $$PPID)
tmpExt=${HOST}.${PPID}.tmp
REPMASKTMP_DIR:=${TMPDIR}/dearlRepMask${tmpExt}
TRFTMP_DIR:=${TMPDIR}/dearlTRF${tmpExt}

all: repMask chains
	@echo "All work complete."

downloadSubmissions:
	mkdir -p ${SUBMISSION_DIR}
	${BIN_DIR}/grabAssemblies.sh ${SUBMISSION_DIR}

updateSubmissions: downloadSubmissions 
	mkdir -p ${RAW_DIR}
	rsync -av ${SUBMISSION_DIR}/*contigs.fa.gz ${RAW_DIR}

# verify the fasta has unique ids
# use only the first 40 characters on the header line as
# RepeatMasker freaks out if you have 50 characters or more.
# We reserve 9 characters for future use.
${RAW_DIR}/%.map: ${RAW_DIR}/%.fa.gz
	zcat $< | ${BIN_DIR}/fastaContigHeaderMapper.py --createMap $@.tmp
	mv $@.tmp $@

# extract files, remove everything in the header line after the unique int id
${ASSEMBLIES_DIR}/%.fa: ${RAW_DIR}/%.fa.gz ${RAW_DIR}/%.map
	mkdir -p $(dir $@)
	zcat $< | ${BIN_DIR}/fastaContigHeaderMapper.py --map ${RAW_DIR}/${*F}.map --goForward > $@.${tmpExt}3 # change headers to contain only unique IDs
	${BIN_DIR}/removeEmptyContigs.py < $@.${tmpExt}3 > $@.${tmpExt}2  # remove the empty contigs from the sequence file
	faFilter -minSize=100 $@.${tmpExt}2 $@.${tmpExt}1 # throw away contigs < 100
	perl -ple 'if(! m/^>/){ s/[^ACGTacgt]/N/g;};' < $@.${tmpExt}1 > $@.${tmpExt} # mask weird IUPACs
	rm $@.${tmpExt}1 $@.${tmpExt}2 $@.${tmpExt}3
	mv $@.${tmpExt} $@

# # run trf on fasta, get the .bed
${TRF_DIR}/%.trf.bed: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p $(dir $@)
	mkdir -p ${TRFTMP_DIR}${*F}
	cd ${TRFTMP_DIR}${*F} && trfBig -bedAt=${TRFTMP_DIR}${*F}/$(notdir $@).tmp $< ${TRFTMP_DIR}${*F}/${*F}.trf.fa -tempDir=${TRFTMP_DIR}${*F}
	mv ${TRFTMP_DIR}${*F}/${*F}.trf.fa ${TRF_DIR}/${*F}.trf.fa
	mv ${TRFTMP_DIR}${*F}/$(notdir $@).tmp $@
	rm -rf ${TRFTMP_DIR}${*F}/

# stub used when you just want trf
trf: $(foreach a, ${ASSEMBLIES}, $(join ${TRF_DIR}/${a},.trf.bed))
	echo "TRF complete."

# run repeatMasker on trf'd assemblies
${REPMASK_DIR}/%/seq.repmask.fa: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p $(dir $@)
	mkdir -p ${REPMASKTMP_DIR}${*F}
	cd ${REPMASKTMP_DIR}${*F} && ${REPEATMASKER} -lib ${MELIB} -parallel 10 -qq -xsmall -dir ${REPMASK_DIR}/${*F} $<
	mv ${REPMASK_DIR}/${*F}/${*F}.fa.masked $@	
	rm -rf ${REPMASKTMP_DIR}${*F}/

# soft add the trf masking via .bed to the already repeat masked fasta file
${ASSEMBLIES_DIR}/%.trf.repmask.fa: ${TRF_DIR}/%.trf.bed ${REPMASK_DIR}/%/seq.repmask.fa
	maskOutFa -softAdd ${REPMASK_DIR}/${*F}/seq.repmask.fa $< $@.${tmpExt}
	mv $@.${tmpExt} $@

# convert to .2bit to facilitate lastz/chain pipeline
${ASSEMBLIES_DIR}/%.trf.repmask.2bit: ${ASSEMBLIES_DIR}/%.trf.repmask.fa
	faToTwoBit $< $@.${tmpExt}
	mv $@.${tmpExt} $@

# stub used when you just want trf and repmask
repMask: $(foreach a, ${ASSEMBLIES}, $(join ${ASSEMBLIES_DIR}/${a},.trf.repmask.2bit))
	@echo "RepeatMasking complete."

# create chainJobs script
${CHAINSCRIPTS_DIR}/hap1/%/chainJobs.csh: ${ASSEMBLIES_DIR}/%.trf.repmask.2bit
	mkdir -p $(dir $@)
	cd $(dir $@) && ${BIN_DIR}/runLastzChain.sh $(dir $@) ${HAPS_DIR}/hap1.trf.repmask.2bit $<
${CHAINSCRIPTS_DIR}/hap2/%/chainJobs.csh: ${ASSEMBLIES_DIR}/%.trf.repmask.2bit
	mkdir -p $(dir $@)
	cd $(dir $@) && ${BIN_DIR}/runLastzChain.sh $(dir $@) ${HAPS_DIR}/hap2.trf.repmask.2bit $<
${CHAINSCRIPTS_DIR}/bac/%/chainJobs.csh: ${ASSEMBLIES_DIR}/%.trf.repmask.2bit
	mkdir -p $(dir $@)
	cd $(dir $@) && ${BIN_DIR}/runLastzChain.sh $(dir $@) ${HAPS_DIR}/bac.trf.repmask.2bit $<

# create the parasol jobList
${CHAINSCRIPTS_DIR}/hap1/%/jobList: ${CHAINSCRIPTS_DIR}/hap1/%/chainJobs.csh
	gensub2 $(dir $@)target.list $(dir $@)query.list $(dir $@)template $@.${tmpExt}
	mv $@.${tmpExt} $@
${CHAINSCRIPTS_DIR}/hap2/%/jobList: ${CHAINSCRIPTS_DIR}/hap2/%/chainJobs.csh
	gensub2 $(dir $@)target.list $(dir $@)query.list $(dir $@)template $@.${tmpExt}
	mv $@.${tmpExt} $@
${CHAINSCRIPTS_DIR}/bac/%/jobList: ${CHAINSCRIPTS_DIR}/bac/%/chainJobs.csh
	gensub2 $(dir $@)target.list $(dir $@)query.list $(dir $@)template $@.${tmpExt}
	mv $@.${tmpExt} $@

# run parasol on swarm
${CHAINSCRIPTS_DIR}/hap1/%/para-complete: ${CHAINSCRIPTS_DIR}/hap1/%/jobList
	ssh swarm ${BIN_DIR}/runPara.sh $(dir $<)
	touch $@
${CHAINSCRIPTS_DIR}/hap2/%/para-complete: ${CHAINSCRIPTS_DIR}/hap2/%/jobList
	ssh swarm ${BIN_DIR}/runPara.sh $(dir $<)
	touch $@
${CHAINSCRIPTS_DIR}/bac/%/para-complete: ${CHAINSCRIPTS_DIR}/bac/%/jobList
	ssh swarm ${BIN_DIR}/runPara.sh $(dir $<)
	touch $@

# create chains
${CHAINS_DIR}/hap1.%.all.chain.gz: ${CHAINSCRIPTS_DIR}/hap1/%/para-complete
	mkdir -p $(dir $@)
	cd $(dir $<) && ./chainJobs.csh
	cp $(dir $<)hap1.${*F}.all.chain.gz $@
${CHAINS_DIR}/hap2.%.all.chain.gz: ${CHAINSCRIPTS_DIR}/hap2/%/para-complete
	mkdir -p $(dir $@)
	cd $(dir $<) && ./chainJobs.csh
	cp $(dir $<)hap2.${*F}.all.chain.gz $@
${CHAINS_DIR}/bac.%.all.chain.gz: ${CHAINSCRIPTS_DIR}/bac/%/para-complete
	mkdir -p $(dir $@)
	cd $(dir $<) && ./chainJobs.csh
	cp $(dir $<)bac.${*F}.all.chain.gz $@

chains: $(foreach a, ${ASSEMBLIES}, $(foreach b, ${HAPLOTYPES}, $(join ${CHAINS_DIR}/${b}.${a},.all.chain.gz)))
	@echo "Chain generation complete."

# chainToAxt will not operate on stdin
${CHAINS_DIR}/hap1.%.all.chain: ${CHAINS_DIR}/hap1.%.all.chain.gz
	zcat $< > $@.tmp
	mv $@.tmp $@
${CHAINS_DIR}/hap2.%.all.chain: ${CHAINS_DIR}/hap2.%.all.chain.gz
	zcat $< > $@.tmp
	mv $@.tmp $@
${CHAINS_DIR}/bac.%.all.chain: ${CHAINS_DIR}/bac.%.all.chain.gz
	zcat $< > $@.tmp
	mv $@.tmp $@

# create axts
${CHAINS_DIR}/hap1.%.axt: ${CHAINS_DIR}/hap1.%.all.chain
	chainToAxt $< ${HAPS_DIR}/hap1.trf.repmask.2bit ${ASSEMBLIES_DIR}/${*F}.trf.repmask.2bit $@.tmp
	mv $@.tmp $@
${CHAINS_DIR}/hap2.%.axt: ${CHAINS_DIR}/hap2.%.all.chain
	chainToAxt $< ${HAPS_DIR}/hap2.trf.repmask.2bit ${ASSEMBLIES_DIR}/${*F}.trf.repmask.2bit $@.tmp
	mv $@.tmp $@
${CHAINS_DIR}/bac.%.axt: ${CHAINS_DIR}/bac.%.all.chain
	chainToAxt $< ${HAPS_DIR}/bac.trf.repmask.2bit ${ASSEMBLIES_DIR}/${*F}.trf.repmask.2bit $@.tmp
	mv $@.tmp $@

# create mafs
${MAFS_DIR}/hap1.%.maf: ${CHAINS_DIR}/hap1.%.axt
	mkdir -p $(dir $@)
	axtToMaf $< ${CHAINSCRIPTS_DIR}/hap1/${*F}/hap1.chrom.sizes ${CHAINSCRIPTS_DIR}/hap1/${*F}/${*F}.chrom.sizes $@.tmp
	mv $@.tmp $@
${MAFS_DIR}/hap2.%.maf: ${CHAINS_DIR}/hap2.%.axt
	mkdir -p $(dir $@)
	axtToMaf $< ${CHAINSCRIPTS_DIR}/hap2/${*F}/hap2.chrom.sizes ${CHAINSCRIPTS_DIR}/hap2/${*F}/${*F}.chrom.sizes $@.tmp
	mv $@.tmp $@
${MAFS_DIR}/bac.%.maf: ${CHAINS_DIR}/bac.%.axt
	mkdir -p $(dir $@)
	axtToMaf $< ${CHAINSCRIPTS_DIR}/bac/${*F}/bac.chrom.sizes ${CHAINSCRIPTS_DIR}/bac/${*F}/${*F}.chrom.sizes $@.tmp
	mv $@.tmp $@

mafs: $(foreach a, ${ASSEMBLIES}, $(foreach b, ${HAPLOTYPES}, $(join ${MAFS_DIR}/${b}.${a},.maf)))
	@echo "MAF generation complete."

##############################

chainClean:
	mkdir -p trash${tmpExt}
	if [ -e ${CHAINS_DIR} ]; then mv ${CHAINS_DIR} trash${tmpExt}/chains; fi
	if [ -e ${CHAINSCRIPTS_DIR} ]; then mv ${CHAINSCRIPTS_DIR} trash${tmpExt}/chainScripts; fi
	rm -rf trash${tmpExt}/chains trash${tmpExt}/chainScripts

trfRepClean:
	mkdir -p trash${tmpExt}
	if [ -e ${REPMASK_DIR} ]; then mv ${REPMASK_DIR} trash${tmpExt}/repmask; fi
	if [ -e ${TRF_DIR} ]; then mv ${TRF_DIR} trash${tmpExt}/trf; fi
	rm -rf trash${tmpExt}/repmask trash${tmpExt}/trf

fullClean: chainClean trfRepClean
	mkdir -p trash${tmpExt}
	if [ -e ${RAW_DIR} ]; then mv ${RAW_DIR} trash${tmpExt}/raw; fi
	if [ -e ${ASSEMBLIES_DIR} ]; then mv ${ASSEMBLIES_DIR} trash${tmpExt}/assemblies; fi
	if [ -e ${MAFS_DIR} ]; then mv ${MAFS_DIR} trash${tmpExt}/mafs; fi
	rm -rf trash${tmpExt}/ &
