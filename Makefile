SHELL:=/bin/bash
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

# extract files
${ASSEMBLIES_DIR}/%.fa: ${RAW_DIR}/%.fa.gz
	mkdir -p ${ASSEMBLIES_DIR}
	gunzip < $< > $@.tmp
ifneq ($(strip $(shell grep '>' $@.tmp | perl -ple 's/>(\S+).*/$$1/;' | head | sort -n - | uniq -d )),)
	@echo "File ${*F}.fa.tmp contains duplicate ids, exiting"
	@exit 1
endif
	mv $@.tmp $@

# remove everything in the header line after the unique int id
${ASSEMBLIES_DIR}/%.preTRF.fa: ${ASSEMBLIES_DIR}/%.fa
	perl -ple 's/>(\S+).*/>$$1/g;' < $< > $@.tmp
	mv $@.tmp $@

# run trf on fasta
${TRF_DIR}/%.trf.fa: ${ASSEMBLIES_DIR}/%.preTRF.fa
	mkdir -p ${TRF_DIR}
	trfBig $< $@.tmp
	mv $@.tmp $@

# run repeatMasker on trf output
${REPMASK_DIR}/%/rmsk.2bit: ${TRF_DIR}/%.trf.fa
	mkdir -p ${REPMASK_DIR}
	repeatMasking_doCluster.py --genome $< --lib ${MELIB} --screenOverride --maxJob 400 --chunkSize 75000 --workDir ${REPMASK_DIR}/${*F}
	mv ${REPMASK_DIR}/${*F}/${*F}.rmsk.2bit $@

# convert 2bit back to fasta
${ASSEMBLIES_DIR}/%.repmask.fa: ${REPMASK_DIR}/%/rmsk.2bit
	twoBitToFa $< $@.tmp
	mv $@.tmp $@

repeatMask: $(foreach a, ${ASSEMBLIES}, $(join ${ASSEMBLIES_DIR}/${a}, .repmask.fa ))
	echo $^

clean:
	rm -rf ${ASSEMBLIES_DIR} ${REPMASK_DIR} ${TRF_DIR}
