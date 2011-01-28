SHELL:=/bin/bash
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
##############################
# DO NOT EDIT BELOW THIS LINE
.SECONDARY: # leave this blank to force make to keep intermediate files
HAPSDIR:=${PROJECT_DIR}/haplotypes
MELIB:=${PROJECT_DIR}/MELibrary/MELib.fa
RAW_DIR:=${PROJECT_DIR}/assemblyArchives
ASSEMBLIES_DIR:=${PROJECT_DIR}/assemblies
REPMASK_DIR:=${PROJECT_DIR}/repeatMasking
ASSEMBLIES:=$(patsubst %.fa.gz,%,$(notdir $(wildcard ${RAW_DIR}/*.gz)))

${ASSEMBLIES_DIR}/%.fa: ${RAW_DIR}/%.fa.gz
	mkdir -p ${ASSEMBLIES_DIR}
	gunzip < $< > $@.tmp
	mv $@.tmp $@

${REPMASK_DIR}/%/seq.rmsk.2bit: ${ASSEMBLIES_DIR}/%.fa
	repeatMasking_doCluster.py --genome $< --lib ${MELIB} --screenOverride --maxJob 500 --chunkSize 50000 --workDir ${REPMASK_DIR}/${*F}

${ASSEMBLIES_DIR}/%.repmask.fa: ${REPMASK_DIR}/%/seq.rmsk.2bit
	twoBitToFa $< $@.tmp
	mv $@.tmp $@

repeatMask: $(foreach a, ${ASSEMBLIES}, $(join ${ASSEMBLIES_DIR}/${a}, .repmask.fa ))
	@echo $^

clean:
	rm -rf ${ASSEMBLIES_DIR} ${REPMASK_DIR}
