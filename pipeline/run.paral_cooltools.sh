#!/usr/bin/bash

list=$1
indir=$2
outdir=$3

if [[ ! -d $outdir ]];then
    mkdir -p $outdir
fi 

function run_cooltools(){
    sampleid=$1
    indir=$2
    outdir=$3
    species=$4
    cooltools="/share/home/wangmx/anaconda3/envs/cooler/bin/cooltools"

    path_to_cpg="/share/home/wangmx/software/dip-c/color"
    ## step 1. calculate 500Kb compartment

    res="100000","250000","500000","1000000"

    $cooltools insulation \
        --threshold Li \
        -o $outdir/$sampleid.insulation.10K.tsv \
        $indir/$sampleid.mcool::resolutions/10000 10000 50000


    $cooltools eigs-cis \
        -o $outdir/$sampleid.100k \
        --n-eigs 1 \
        --phasing-track $path_to_cpg/${species}.cpg.cooltools.100k.txt \
        $indir/$sampleid.mcool::resolutions/100000 \
        --view /share/home/wangmx/software/dip-c/color/mm10.view.df.txt


    $cooltools eigs-cis \
        -o $outdir/$sampleid.500k \
        --n-eigs 1 \
        --phasing-track $path_to_cpg/${species}.cpg.cooltools.500k.txt \
        $indir/$sampleid.mcool::resolutions/500000 \
        --view /share/home/wangmx/software/dip-c/color/mm10.view.df.txt


    $cooltools expected-cis \
        $indir/$sampleid.mcool::resolutions/500000 \
        -o $outdir/$sampleid.expected.cis.500k.tsv

    $cooltools insulation \
        --threshold Li \
        -o $outdir/$sampleid.insulation.500K.tsv \
        $indir/$sampleid.mcool::resolutions/250000 250000 500000

    $cooltools saddle \
        --qrange 0.02 0.98 \
        --fig png \
        -o $outdir/$sampleid.saddle.cis.500K \
        $indir/$sampleid.mcool::resolutions/500000 \
        $outdir/$sampleid.500k.cis.vecs.tsv $outdir/$sampleid.expected.cis.500k.tsv

    # $cooltools dots \
    #     --nproc $thread \
    #     --max-loci-separation 10000000 \
    #     --fdr 0.01 \
    #     -o $outdir/${sampleid}.dots_10k.tsv \
    #     --view /share/home/wangmx/software/dip-c/color/mm10.view.df.txt \
    #     $indir/${sampleid}.mcool::resolutions/10000 \
    #     $outdir/${sampleid}.expected.cis.10k.tsv

}

run_cooltools="/share/home/wangmx/software/pipelines/dscHiC-process/scripts/pip.cooltoos.sh"
species="mm10"

export -f run_cooltools
parallel --progress --verbose run_cooltools ::: $(cat $list | xargs) ::: $indir ::: $outdir ::: $species
