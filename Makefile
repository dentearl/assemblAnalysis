SHELL:=/bin/bash -e
export SHELLOPTS=pipefail

binDir:=bin
srcDir:=src

.PHONY: all clean
all: $(foreach f,$(notdir $(wildcard ${srcDir}/*py)),${binDir}/$f) ${binDir}/createPairCorrResidPlot.R ${binDir}/myPairs.R

${binDir}/%.py: ${srcDir}/%.py
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	if [[ $* != lib* ]] && [[ $* != test* ]] ; then chmod 755 $@.tmp ; fi
	mv $@.tmp $@

${binDir}/%.R: ${srcDir}/%.R
	@mkdir -p $(dir $@)
	cp $< $@.tmp
	if [[ $* != lib* ]] && [[ $* != test* ]] ; then chmod 755 $@.tmp ; fi
	mv $@.tmp $@

clean:
	rm -rf ${binDir}
