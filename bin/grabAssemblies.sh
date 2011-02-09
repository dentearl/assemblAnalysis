#!/bin/bash -e
DEST_DIR=${1%/}
if [ ! -e $DEST_DIR ]; then
    echo "Local Destination directory $DEST_DIR does not exist!"
    exit 2
fi
if [ ! -d $DEST_DIR ]; then
    echo "Local Destination directory $DEST_DIR is not a directory!"
    exit 2
fi

REMOTE=http://bioshare.bioinformatics.ucdavis.edu/Data/h4dys0s28v

files=$(curl -s $REMOTE/ | grep 'gz' | perl -ple 's/^.*?\.fa\.gz">(\S+?\.fa\.gz)<\/a><\/td>.*$/$1/;')

i=1
tot=${#files[@]}
for f in $files; do
    if [ ! -e $DEST_DIR/$f ] ; then
        echo "$f $i of $tot"
        curl $REMOTE/$f -o $DEST_DIR/$f
        i=((i++))
    fi
done
