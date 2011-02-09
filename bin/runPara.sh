#!/bin/bash -e
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
para shove jobList -retries=10
exit 0
