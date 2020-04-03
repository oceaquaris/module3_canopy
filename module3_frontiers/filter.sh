#!/bin/bash

outdir='../filtered'

for i in *.gz
do
    bname=${i%.csv.gz}
    zcat $i | cut -f 3,7,8 -d ',' | grep -v ',NA' | grep -v 'fill' > "${outdir}/${bname}.filter.csv"
    gzip "${outdir}/${bname}.filter.csv"
done
