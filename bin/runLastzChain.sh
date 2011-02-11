#!/bin/bash -e
set -beEu -o pipefail
# runLastzChain.sh
# dent earl dearl (a) soe ucsc edu
# 9 Feb 2011
# adapted from a shell script of the same name by Hiram Clawson hiram (a) soe
#
# shell script to prepare scripts to run lastz/chain pipeline
# this script produces, among others, the following files in the
# supplied workDir:
#  * runLastz       wrapper script to be run on the cluster
#  * template       used by gensub2 to create jobList
#  * jobList        jobList used by parasol
#  * chainJobs.csh  script to be run following the cluster job
#
# usage: runLastzChain.sh workDir targetFile.2bit queryFile.2bit [near|medium|far]
#
# select one of three different parameter sets
# near == genomes close to each other
# medium == genomes at middle distance from each other
# far == genomes distant from each other
##############################

usage()
{
    echo -e "$0\nusage: runLastzChain.sh workDir targetFile.2bit queryFile.2bit [near|medium|far]"
    exit 2
}
if [[ $# -lt 3 ]] || [[ $# -gt 4 ]]; then
    usage
fi
export chainNear="-minScore=5000 -linearGap=medium"
export chainMedium="-minScore=3000 -linearGap=medium"
export chainFar="-minScore=5000 -linearGap=loose"
export lastzNear="B=0 C=0 E=150 H=0 K=4500 L=3000 M=254 O=600 Q=/scratch/data/blastz/human_chimp.v2.q T=2 Y=15000"
export lastzMedium="B=0 C=0 E=30 H=0 K=3000 L=3000 M=50 O=400 T=1 Y=9400"
export lastzFar="B=0 C=0 E=30 H=2000 K=2200 L=6000 M=50 O=400 Q=/scratch/data/blastz/HoxD55.q T=2 Y=3400"

if [[ $# -eq 4 ]]; then
    dist=$(echo $4 | tr '[:upper:]' '[:lower:]')
    if [[ $dist != 'near' && $dist != 'medium' && $dist != 'far' ]]; then
        usage
    fi
else
    dist='near'
fi
case $dist in
    near ) export chainParams="$chainNear"
        export lastzParams="$lastzNear";;
    medium ) export chainParams="$chainMedium"
        export lastzParams="$lastzMedium";;
    far ) export chainParams="$chainFar"
        export lastzParams="$lastzFar";;
esac
if [ ! -e $1 ] || [ ! -d $1 ]; then
    usage
else
    export WRKDIR=$1
    export WRKDIR=${WRKDIR%/} # remove trailing slash
fi
if [ ! -e $2 ] || [ ! -e $3 ]; then
    usage
else
    export TNAME=$( basename $2 .trf.repmask.2bit )
    export TARGET=$2
    export TARGET_DIR=$( dirname $TARGET )
    export QNAME=$( basename $3 .trf.repmask.2bit )
    export QUERY=$3
    export QUERY_DIR=$( dirname $QUERY )
fi

ls -ld $TARGET $QUERY

if [ ! -s ${TNAME}.chrom.sizes ]; then
    twoBitInfo ${TARGET} stdout | sort -k2nr > ${TNAME}.chrom.sizes
    rm -rf ${TNAME}PartList ${TNAME}.part.list
    mkdir ${TNAME}PartList
fi
if [ ! -s ${QNAME}.chrom.sizes ]; then
    twoBitInfo ${QUERY} stdout | sort -k2nr > ${QNAME}.chrom.sizes
    rm -rf ${QNAME}PartList ${QNAME}.part.list
    mkdir ${QNAME}PartList
fi

if [ ! -s ${TNAME}.part.list ]; then
    #partitionSequence.pl 10000000 10000 ${TARGET} ${TNAME}.chrom.sizes 1 \
    # Hiram's call: 10 Mb chunks with 10,000 bp overlap with 1 contig per partition
    partitionSequence.pl 10000000 10000 ${TARGET} ${TNAME}.chrom.sizes 1 \
    -lstDir ${TNAME}PartList > ${TNAME}.part.list
fi
if [ ! -s ${QNAME}.part.list ]; then
    #partitionSequence.pl 20000000 0 ${QUERY} ${QNAME}.chrom.sizes 1 \
    partitionSequence.pl 20000000 0 ${QUERY} ${QNAME}.chrom.sizes 2000 \
    -lstDir ${QNAME}PartList > ${QNAME}.part.list
fi

grep --invert-match PartList ${TNAME}.part.list > target.list
if [ $(ls -1A ${TNAME}PartList/ | wc -l) -gt 0 ]; then
    for F in ${TNAME}PartList/*.lst
      do
      cat ${F}
    done >> target.list
fi

####
# we may potentially do away with this section
# the following conditional is to prevent the script from exiting if the
# grep comes back with no matches.
if grep --invert-match PartList ${QNAME}.part.list > query.list ; then
    echo 'matches found' >> /dev/null
    else 
    echo 'no matches found' >> /dev/null
fi
if [ $(ls -1A ${QNAME}PartList/ | wc -l) -gt 0 ]; then
    for F in ${QNAME}PartList/*.lst
      do
      # cat ${F}
      echo ${F}
    done >> query.list
fi
# 
####
echo "constructLiftFile.pl ${TNAME}.chrom.sizes target.list > target.lift"
constructLiftFile.pl ${TNAME}.chrom.sizes target.list > target.lift
echo "constructLiftFile.pl ${QNAME}.chrom.sizes query.list > query.lift"
constructLiftFile.pl ${QNAME}.chrom.sizes query.list > query.lift

echo "#LOOP" > template
echo './runLastz $(file1) $(file2) { check out exists+ psl/$(file1).$(file2).psl.gz }' >> template
echo "#ENDLOOP" >> template

# This script will be run on the cluster.
cat <<_EOF_ > runLastz
#!/bin/bash -e
set -beEu -o pipefail

export WDIR=$WRKDIR
export TDIR=$TARGET_DIR
export QDIR=$QUERY_DIR
export FT=\$1
export FQ=\$2
export tmpDir=/scratch/tmp/\${FT}
mkdir -p raw psl \${tmpDir}
twoBitToFa \${TDIR}/\${FT} \${tmpDir}/\${FT}.fa
if [ \$(echo \$FQ | perl -ple 's/.*(\..*?)\$/\$1/;') == '.lst' ]; then
   cat \${WDIR}/*PartList/\${FQ} \
       | while read filename; do 
           twoBitToFa \$filename \${tmpDir}/\${FQ}.fa.tmp
           cat \${tmpDir}/\${FQ}.fa.tmp >> \${tmpDir}/\${FQ}.fa
        done
else
   twoBitToFa \${QDIR}/\${FQ} \${tmpDir}/\${FQ}.fa
fi
/cluster/bin/penn/lastz-distrib-1.02.00/bin/lastz \${tmpDir}/\${FT}.fa \
    \${tmpDir}/\${FQ}.fa \
    ${lastzParams} \
    > raw/\${FT}.\${FQ}.lav
lavToPsl raw/\${FT}.\${FQ}.lav stdout \
    | liftUp -type=.psl stdout target.lift error stdin \
    | liftUp -nohead -pslQ -type=.psl stdout query.lift carry stdin \
    | gzip -c > psl/\${FT}.\${FQ}.psl.gz
rm -f \${tmpDir}/\${FT}.fa \${tmpDir}/\${FQ}.fa
rmdir --ignore-fail-on-non-empty \${tmpDir}
_EOF_

chmod 755 $WRKDIR/runLastz

mkdir -p chain
echo "#!/bin/csh -fe" > chainJobs.csh
for T in `cat target.list | sed -e "s#${TARGET_DIR}/##"` ; do
  echo    "zcat psl/${T}.*.psl.gz \\"
  echo    "    | axtChain -psl -verbose=0 ${chainParams} \\"
  echo    "      stdin ${TARGET} ${QUERY} stdout \\"
  echo    "    | chainAntiRepeat ${TARGET} ${QUERY} stdin chain/${T}.chain"
done >> chainJobs.csh

chmod 755 $WRKDIR/chainJobs.csh

echo "find ./chain -name \"*.chain\" | chainMergeSort -inputList=stdin | gzip -c > ${TNAME}.${QNAME}.all.chain.gz" >> chainJobs.csh
