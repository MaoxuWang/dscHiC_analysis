#!/usr/bin/bash

set -e
set -o pipefail

helpdoc() {
    cat <<EOF
Description:
    [main::] Demultiplexed single cell Hi-C or dscHiC-multiome data processing pipeline created by Maoxu Wang, Xie Lab, Peking University.
    Last update: 2023.06.07
Usage:
    $run [options]
        -h print user manual information
Options:
    --indir             -- input directory
    --outdir            -- output directory
    --sampleid          -- sample name specified 
    --reference         -- reference genome folder with bwa indexeds
    --run_mode          -- [pipe | mapping | merge | impute | mixing_mapping | mixing_merge]
    --config            -- software configuration (path) 
    --ref_cpg           -- CpG frequency file of certain resolution 
    --calLibrarySize    -- whether to calculate library size (True or False, default: False)
    --rmsk              -- repeats file used in AMULET pipe
    --autosome          -- chromosomes file 
    --chromSize         -- chromosome size file (TAB sperated)
    --debug             -- debug mode
    --library           -- dscHiC or dscHiC-multiome or dscHiC-deep
    --help              -- print manual
EOF
}
ARGS=$(getopt -o o:i:h -l indir:,outdir:,sampleid:,reference:,calLibrarySize:,run_mode:,debug:,rmsk:,chromSize:,calib_file:,cytoBand:,config:,autosome:,ref_cpg:,library:,min_contactsN:,help -- "$@")

if [ $? != 0 ]; then
    echo "terminating..." >&2
    exit 1
fi

eval set -- "$ARGS"

while true; do
    case "$1" in
    -i | --indir)
        indir=$2
        shift 2
        ;;
    -o | --outdir)
        outdir=$2
        shift 2
        ;;
    --sampleid)
        sampleid=$2
        shift 2
        ;;
    --reference)
        reference=$2
        shift 2
        ;;
    --debug)
        debug=$2
        shift 2
        ;;
    --run_mode)
        run_mode=$2
        shift 2
        ;;
    --cytoBand)
        cytoBand=$2
        shift 2
        ;;
    --calLibrarySize)
        calLibrarySize=$2
        shift 2
        ;;
    --calib_file)
        calib_file=$2
        shift 2
        ;;
    --ref_cpg)
        ref_cpg=$2
        shift 2
        ;;
    --rmsk)
        rmsk=$2
        shift 2
        ;;
    --chromSize)
        chromSize=$2
        shift 2
        ;;
    --autosome)
        autosome=$2
        shift 2
        ;;
    --config)
        config=$2
        shift 2
        ;;
    --library)
        library=$2
        shift 2
        ;;
    --min_contactsN)
        min_contactsN=$2
        shift 2
        ;;
    -h | --help)
        helpdoc
        exit 1
        ;;
    --)
        shift
        echo "Parameters recevied."
        break
        ;;
    *)
        echo "unknown parameter:" $1
        helpdoc
        exit 1
        ;;
    esac
done

## inspect arguments

if [[ ! -n $outdir ]]; then outdir=output; fi
if [[ ! -n $library ]]; then library="dscHiC"; fi
echo "Library: $library"
if [[ ! -n $calLibrarySize ]]; then calLibrarySize=False; fi
if [[ ! -n $debug ]]; then debug=False; fi

if [[ ! -n $ref_cpg ]]; then
    if [[ ! $run_mode == "mapping" ]]; then
        echo "Please specify CpG frequency file"
        exit
    fi
fi

if [[ ! -n $sampleid ]]; then
    echo "Please specify sample name"
    exit
fi

if [[ ! -n $config ]]; then
    echo "Please specify software configuration file "
    exit
fi

if [[ -n $min_contactsN ]];then
    echo "Filter contacts > $min_contactsN ..."    
    downstream_cb="$outdir/QC/$sampleid.downstream.cellbarcode"
else
    min_contactsN=1000
    downstream_cb="$outdir/QC/$sampleid.CellCalled.txt"
fi

species=$(basename -s ".fa" $reference)

src_dir=$(echo $config | sed 's/software.config/scripts/' | xargs)

## check software
bwa=$(grep bwa $config | sed 's/.*=//')
samtools=$(grep samtools $config | sed 's/.*=//')
dscHiCtools=$(grep "dscHiCtools" $config | sed 's/.*=//')
Rscript=$(grep "Rscript" $config | sed 's/.*=//')
hickit_addReadName=$(grep "hickit_addReadName" $config | sed 's/.*=//')
hickit_js=$(grep "hickit_js" $config | sed 's/.*=//')
AMULET_dir=$(grep "AMULET_dir" $config | sed 's/.*=//')
bedtools=$(grep "bedtools" $config | sed 's/.*=//')
python=$(grep "python" $config | sed 's/.*=//')
hicluster=$(grep "hicluster" $config | sed 's/.*=//')

## package scripts
EstimateLibrarySize="$src_dir/calculateLibrarySize.R"
CellCalling="$src_dir/CellCalling.R"
py_scAB="$src_dir/AB_compartments.py"
AB_embedding="$src_dir/AB_embedding.R"
contacts_distribution="$src_dir/contacts_distribution.R"
mixingCount="$src_dir/mixingCount.py"
mixingPlot="$src_dir/mixingPlot.R"
metrics="$src_dir/metrics.R"
scVI_input="$src_dir/prepareInput.R"

# ------- pipeline ----------
echo "-----Demultiplexed single cell Hi-C analysis-----"
echo -e "Start processing :\n" \
    "\tindex: $reference\n" \
    "\tCpG frequency: $ref_cpg"
echo "Outdir is " $outdir
echo "CPUS used: $SLURM_CPUS_PER_TASK"
echo "Sample is  $sampleid"

if [[ $SLURM_CPUS_PER_TASK > 20 ]];then 
    sort_thread=20
else
    sort_thread=$SLURM_CPUS_PER_TASK
fi 

 if [[ ! -d $outdir/countCB ]]; then
        mkdir -p $outdir/countCB
fi

## check input files
Read1=$(ls $indir/${sampleid}_R1_001.fastq.gz)
Read2=$(ls $indir/${sampleid}_I2_001.fastq.gz)
Read3=$(ls $indir/${sampleid}_R2_001.fastq.gz)

if [ ! -s $Read1 ]; then
    echo $Read1 "Read not exist"
    exit
fi

if [ ! -s $Read2 ]; then
    echo $Read2 "Read not exist"
    exit
fi

if [ ! -s $Read3 ]; then
    echo $Read3 "Read not exist"
    exit
fi

if [[ $library == "dscHiC-multiome" ]];then 
    start_idx=8
    end_idx=24
else
    start_idx=0
    end_idx=16
fi

if [[ ! -s $outdir/countCB/${sampleid}.R1.barcoded.fastq.gz ]];then 
    $dscHiCtools mapBarcode \
        --read1 $Read1 \
        --index $Read2 \
        --read2 $Read3 \
        --start_idx $start_idx \
        --end_idx $end_idx \
        --outdir $outdir/countCB \
        --sampleName $sampleid \
        --threads $SLURM_CPUS_PER_TASK \
        --n_lines -1 \
        --chemistry $chemistry
    echo "[mapBarcode::] finished!"
    time_n=`date`
    echo -e "$time_n\nfinished" > $outdir/countCB/mapBarcode.done
else 
    echo "[mapBarcode::] finished!"
fi 

 if [[ ! -d $outdir/align ]]; then
        mkdir -p $outdir/align
    fi

## prepare input file
if [[ ! -s $outdir/align/${sampleid}.sam.gz ]]; then
    barcodedRead1="$outdir/countCB/${sampleid}.R1.barcoded.fastq.gz"
    barcodedRead2="$outdir/countCB/${sampleid}.R2.barcoded.fastq.gz"

            ## mapping using bwa-mem
    $bwa mem -5SPM -t $SLURM_CPUS_PER_TASK \
        $reference $barcodedRead1 $barcodedRead2 \
        | pigz  >$outdir/align/${sampleid}.sam.gz

    echo "[mapping::] finished!"
else
    echo "[mapping::] finished!"
fi