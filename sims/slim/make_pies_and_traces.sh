#!/bin/bash

USAGE="\
Usage:
    $0 (xval directory name)\
"

if [ $# -lt 1 ]
then
    echo "$USAGE"
    exit 1
fi

for DIR in "$@"
do
    for K in 2 3 4 5 6 7
    do
        FILES=$(ls -t ${DIR}/xval_sp_rep*K${K}_{pie,trace}*.pdf)
        OUTFILE="${DIR}/K${K}_pies_and_traces.pdf"
        pdfjoin --outfile "$OUTFILE" --rotateoversize false --nup 2x1 $(for x in $FILES; do echo "$x 1"; done)
    done
done
