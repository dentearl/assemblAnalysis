#!/bin/bash -e
set -beEu -o pipefail
usage()
{
    echo -e "usage: $(basename $0) workDir"
    echo "workDir should contain a viable parasol jobList called jobList."
    exit 2
}
if [[ $# -ne 1 ]]; then
    usage
fi
if [[ ! -d $1 ]]; then
    usage
fi
if [[ ! -e $1/jobList ]]; then
    usage
fi

workDir=$1
cd $workDir
para create jobList
para try jobList
para check jobList
para shove jobList -retries=4
exit 0
