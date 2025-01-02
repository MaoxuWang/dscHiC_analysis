#!/usr/bin/bash
set -e
set -o pipefail

sampleid=$1
outdir=$2

hickit_addReadName=$(grep "hickit_addReadName" $config | sed 's/.*=//')
hickit_addReadName=$(grep "hickit_js" $config | sed 's/.*=//')

$hickit_js sam2seg \
    $outdir/align/${sampleid}.sam.gz |
    $hickit_js chronly -y - | bgzip \
    >$outdir/contacts/${sampleid}.contacts.seg.gz

# dup-dist=100, remove duplicates by cell barcode, start, end
$hickit_addReadName --dup-dist=100 \
    -i $outdir/contacts/${sampleid}.contacts.seg.gz -o $outdir/contacts/${sampleid}.contacts.pairs