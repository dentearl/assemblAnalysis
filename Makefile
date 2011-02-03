SHELL:=/bin/bash -e
export SHELLOPTS=pipefail
# dent earl dearl (a) soe ucsc edu
# 26 Jan 2011
##############################
# EDIT ONLY PROJECT_DIR 
PROJECT_DIR:=/hive/users/dearl/assemblathon
# It is assumed that all assemblies will be in fasta format, gzipped
# and named according to:
#            mySweetAssembly.fa.gz
# All assemblies go in the ${PROJECT_DIR}/assemblyArchives directory.
# Assumed initial directory structure:
# ${PROJECT_DIR}/
#     |------ /haplotypes/hap*.fa  # the two haplotypes, repeat masked
#     |------ /assemblyArchives/*  # all *.fa.gz assembly entrants
#     |------ /MElibrary/MELib.fa  # The mobile element library to be used in repeat masking
#
# DEPENDENCIES
#   repeatMasking_doCluster.py
#   repeatMasker
#   trf - tandem repeat finder
#
##############################
# DO NOT EDIT BELOW THIS LINE
.SECONDARY: # leave this blank to force make to keep intermediate files
HAPSDIR:=${PROJECT_DIR}/haplotypes
MELIB:=${PROJECT_DIR}/MELibrary/MELib.fa
RAW_DIR:=${PROJECT_DIR}/assemblyArchives
ASSEMBLIES_DIR:=${PROJECT_DIR}/assemblies
REPMASK_DIR:=${PROJECT_DIR}/repeatMasking
TRF_DIR:=${PROJECT_DIR}/tandemRepeatFinder
ASSEMBLIES:=$(patsubst %.fa.gz,%,$(notdir $(wildcard ${RAW_DIR}/*.gz)))
UNIQ:=$(strip $(shell date '+%s' | perl -ple 's/^\d{6}//;')) # last 4 digits of time.
REPEATMASKER:=/scratch/data/RepeatMasker/RepeatMasker
HOST=$(shell hostname)
PPID=$(shell echo $$PPID)
TMPEXT=${HOST}.${PPID}.tmp
REPMASKTMP_DIR:=${TMPDIR}/dearlRepMask${TMPEXT}
TRFTMP_DIR:=${TMPDIR}/dearlTRF${TMPEXT}

all: repeatMask
	@echo "Work complete"

${RAW_DIR}/%.fa-verified: ${RAW_DIR}/%.fa.gz
	@if [[ ! -z $$(zcat $< | grep '>' - | perl -ple 's/>(\S+).*/$$1/;' | head | sort -n | uniq -d ) ]]; then \
		echo "File $< contains duplicate ids, exiting" >&2 ; \
		exit 1; \
	fi
	touch $@

# extract files, remove everything in the header line after the unique int id
${ASSEMBLIES_DIR}/%.fa: ${RAW_DIR}/%.fa.gz ${RAW_DIR}/%.fa-verified
	mkdir -p $(dir $@)
	zcat $< | perl -ple 's/>(\S+).*/>$$1/g;' > $@.${TMPEXT}
	mv $@.${TMPEXT} $@

# # run trf on fasta, get the .bed
${TRF_DIR}/%.trf.bed: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p ${TRF_DIR}
	mkdir -p ${TRFTMP_DIR}${*F}
	cd ${TRFTMP_DIR}${*F} && trfBig -bedAt=$@.tmp $< ${TRF_DIR}/${*F}.trf.fa
	mv $@.tmp $@
	rm -rf ${TRFTMP_DIR}${*F}

# run repeatMasker on trf'd assemblies
${REPMASK_DIR}/%/seq.repmask.fa: ${ASSEMBLIES_DIR}/%.fa
	mkdir -p ${REPMASK_DIR}/${*F}
	mkdir -p ${REPMASKTMP_DIR}${*F}
	cd ${REPMASKTMP_DIR}${*F} && ${REPEATMASKER} -lib ${MELIB} -parallel 10 -qq -xsmall -dir ${REPMASK_DIR}/${*F} $<
	mv ${REPMASK_DIR}/${*F}/${*F}.fa.masked $@
	rm -rf ${RMKSITMP_DIR}${*F}

# soft add the trf masking via .bed to the already repeat masked fasta file
${ASSEMBLIES_DIR}/%.trf.repmask.fa: ${TRF_DIR}/%.trf.bed ${REPMASK_DIR}/%/seq.repmask.fa
	maskOutFa -softAdd ${REPMASK_DIR}/${*F}/seq.repmask.fa $< $@.${TMPEXT}
	mv $@.${TMPEXT} $@

# convert to .2bit to facilitate lastz/chain pipeline
${ASSEMBLIES_DIR}/%.trf.repmask.2bit: ${ASSEMBLIES_DIR}/%.trf.repmask.fa
	faToTwoBit $< $@.${TMPEXT}
	mv $@.${TMPEXT} $@

repeatMask: $(foreach a, ${ASSEMBLIES}, $(join ${ASSEMBLIES_DIR}/${a}, .trf.repmask.2bit ))
	@echo $^

clean:
	rm -rf ${ASSEMBLIES_DIR} ${REPMASK_DIR} ${TRF_DIR} $(wildcard ${RAW_DIR}/*.fa-verified)
