SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
# dent earl dearl (a) soe ucsc edu
# 26 Jan 2011
##############################
# EDIT ONLY PROJECT_DIR 
PROJECT_DIR:=/hive/users/dearl/assemblathon
BIN_DIR:=/hive/users/dearl/assemblathon/assemblathon/trunk/bin
# It is assumed that all assemblies will be in fasta format, gzipped
# and named according to:
#            mySweetAssembly.fa.gz
# All assemblies go in the ${PROJECT_DIR}/assemblyArchives directory.
# Assumed initial directory structure:
# ${PROJECT_DIR}/
#     |------ /haplotypes/hap1.trf.repmask.2bit  # hap1, trf'd and repeatMasked
#     |------ /haplotypes/hap2.trf.repmask.2bit  # hap2. 
#     |------ /assemblyArchives/*    # all *.fa.gz assembly entrants
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
#
##############################
# DO NOT EDIT BELOW THIS LINE
.SECONDARY: # leave this blank to force make to keep intermediate files
HAPS_DIR:=${PROJECT_DIR}/haplotypes
RAW_DIR:=${PROJECT_DIR}/assemblyArchives
CHAINS_DIR:=${PROJECT_DIR}/chains
CHAINSCRIPTS_DIR:=${PROJECT_DIR}/chainScripts
ASSEMBLIES_DIR:=${PROJECT_DIR}/assemblies
REPMASK_DIR:=${PROJECT_DIR}/repeatMasking
TRF_DIR:=${PROJECT_DIR}/tandemRepeatFinder
MELIB:=${PROJECT_DIR}/MELibrary/MELib.fa
ASSEMBLIES:=$(patsubst %.fa.gz,%,$(notdir $(wildcard ${RAW_DIR}/*.gz)))
HAPLOTYPES:=$(patsubst %.trf.repmask.2bit,%,$(notdir $(wildcard ${HAPS_DIR}/*.trf.repmask.2bit)))
UNIQ:=$(strip $(shell date '+%s' | perl -ple 's/^\d{6}//;')) # last 4 digits of time.
REPEATMASKER:=/scratch/data/RepeatMasker/RepeatMasker
HOST=$(shell hostname)
PPID=$(shell echo $$PPID)
tmpExt=${HOST}.${PPID}.tmp
REPMASKTMP_DIR:=${TMPDIR}/dearlRepMask${tmpExt}
TRFTMP_DIR:=${TMPDIR}/dearlTRF${tmpExt}

all: chains
	@echo "Work complete"

updateSubmissions:
	mkdir -p submissions
	${BIN_DIR}/grabAssemblies.sh ${PROJECT_DIR}/submissions

${RAW_DIR}/%.fa-verified: ${RAW_DIR}/%.fa.gz
	if [[ ! -z $$(zcat $< | grep '>' - | perl -ple 's/>(\S+).*/$$1/;' | head | sort -n | uniq -d ) ]]; then \
		echo "File $< contains duplicate ids, exiting" >&2 ; \
		exit 1; \
	fi
	touch $@

# extract files, remove everything in the header line after the unique int id
${ASSEMBLIES_DIR}/%.fa: ${RAW_DIR}/%.fa.gz ${RAW_DIR}/%.fa-verified
	mkdir -p $(dir $@)
	zcat $< | perl -ple 's/>(\S+).*/>$$1/g;' > $@.${tmpExt}1
	faFilter -minSize=100 $@.${tmpExt}1 $@.${tmpExt}
	mv $@.${tmpExt} $@

# # run trf on fasta, get the .bed
${TRF_DIR}/%.trf.bed: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p $(dir $@)
	mkdir -p ${TRFTMP_DIR}${*F}
	cd ${TRFTMP_DIR}${*F} && trfBig -bedAt=$@.tmp $< ${TRF_DIR}/${*F}.trf.fa
	mv $@.tmp $@
	rm -rf ${TRFTMP_DIR}${*F}

# run repeatMasker on trf'd assemblies
${REPMASK_DIR}/%/seq.repmask.fa: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p $(dir $@)
	mkdir -p ${REPMASKTMP_DIR}${*F}
	cd ${REPMASKTMP_DIR}${*F} && ${REPEATMASKER} -lib ${MELIB} -parallel 10 -qq -xsmall -dir ${REPMASK_DIR}/${*F} $<
	mv ${REPMASK_DIR}/${*F}/${*F}.fa.masked $@
	rm -rf ${RMKSITMP_DIR}${*F}

# soft add the trf masking via .bed to the already repeat masked fasta file
${ASSEMBLIES_DIR}/%.trf.repmask.fa: ${TRF_DIR}/%.trf.bed ${REPMASK_DIR}/%/seq.repmask.fa
	maskOutFa -softAdd ${REPMASK_DIR}/${*F}/seq.repmask.fa $< $@.${tmpExt}
	mv $@.${tmpExt} $@

# convert to .2bit to facilitate lastz/chain pipeline
${ASSEMBLIES_DIR}/%.trf.repmask.2bit: ${ASSEMBLIES_DIR}/%.trf.repmask.fa
	faToTwoBit $< $@.${tmpExt}
	mv $@.${tmpExt} $@

# create chainJobs script
${CHAINSCRIPTS_DIR}/hap1/%/chainJobs.csh: ${ASSEMBLIES_DIR}/%.trf.repmask.2bit
	mkdir -p $(dir $@)
	cd $(dir $@) && ${BIN_DIR}/runLastzChain.sh $(dir $@) ${HAPS_DIR}/hap1.trf.repmask.2bit $<

${CHAINSCRIPTS_DIR}/hap2/%/chainJobs.csh: ${ASSEMBLIES_DIR}/%.trf.repmask.2bit
	mkdir -p $(dir $@)
	cd $(dir $@) && ${BIN_DIR}/runLastzChain.sh $(dir $@) ${HAPS_DIR}/hap2.trf.repmask.2bit $<

# create the parasol jobList
${CHAINSCRIPTS_DIR}/hap1/%/jobList: ${CHAINSCRIPTS_DIR}/hap1/%/chainJobs.csh
	gensub2 $(dir $@)target.list $(dir $@)query.list $(dir $@)template $@.${tmpExt}
	mv $@.${tmpExt} $@

${CHAINSCRIPTS_DIR}/hap2/%/jobList: ${CHAINSCRIPTS_DIR}/hap2/%/chainJobs.csh
	gensub2 $(dir $@)target.list $(dir $@)query.list $(dir $@)template $@.${tmpExt}
	mv $@.${tmpExt} $@

# run parasol on swarm
${CHAINSCRIPTS_DIR}/hap1/%/para-complete: ${CHAINSCRIPTS_DIR}/hap1/%/jobList
	aoeu ssh swarm ${BIN_DIR}/runPara.sh $(dir $<)
	touch $@

${CHAINSCRIPTS_DIR}/hap2/%/para-complete: ${CHAINSCRIPTS_DIR}/hap2/%/jobList
	aoeu ssh swarm ${BIN_DIR}/runPara.sh $(dir $<)
	touch $@

${CHAINS_DIR}/hap1.%.all.chain.gz: ${CHAINSCRIPTS_DIR}/hap1/%/para-complete
	mkdir -p $(dir $@)
	cd $(dir $<) && ./chainJobs.csh
	cp $(dir $<)hap1.${*F}.all.chain.gz $@

${CHAINS_DIR}/hap2.%.all.chain.gz: ${CHAINSCRIPTS_DIR}/hap2/%/para-complete
	mkdir -p $(dir $@)
	cd $(dir $<) && ./chainJobs.csh
	cp $(dir $<)hap2.${*F}.all.chain.gz $@

chains: $(foreach a, ${ASSEMBLIES}, $(foreach b, ${HAPLOTYPES}, $(join ${CHAINS_DIR}/${b}.${a},.all.chain.gz)))
	@echo $^

chainClean:
	mkdir -p trash
	mv ${CHAINS_DIR} trash/${CHAINS_DIR}
	mv ${CHAINSCRIPTS_DIR} trash/${CHAINSCRIPTS_DIR}
	rm -rf trash/${CHAINS_DIR} trash/${CHAINSCRIPTS_DIR}
	rm -rf trash/

fullClean:
	mkdir -p trash
	mv ${ASSEMBLIES_DIR} trash/${ASSEMBLIES_DIR}
	mv ${REPMASK_DIR} trash/${REPMASK_DIR}
	mv ${TRF_DIR} trash/${TRF_DIR}
	mv ${CHAINS_DIR} trash/${CHAINS_DIR}
	mv ${CHAINSCRIPTS_DIR} trash/${CHAINSCRIPTS_DIR}
	rm -rf trash/${ASSEMBLIES_DIR} trash/${REPMASK_DIR} trash/${TRF_DIR} $(wildcard ${RAW_DIR}/*.fa-verified) trash/${CHAINS_DIR} trash/${CHAINSCRIPTS_DIR}
	rm -rf trash/
