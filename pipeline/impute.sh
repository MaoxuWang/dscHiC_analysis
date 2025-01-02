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


function scAB() {
    if [[ ! -s $outdir/embedding/AB.done ]];then 
        if [[ ! -d $outdir/embedding/AB ]]; then
            mkdir -p $outdir/embedding/AB
        fi

        if [[ ! -s $outdir/contacts/${sampleid}.contacts.pairs.txt.gz ]]; then
            zcat $outdir/contacts/${sampleid}.contacts.pairs.gz |
                grep -v "#" |
                awk '{print $1,$2,$3,$4,$5}' OFS="\t" \
                    | pigz >$outdir/contacts/${sampleid}.contacts.pairs.txt.gz
        fi

        ## default res: 1M        
        ### using all cells
        $python $py_scAB \
            --input_contact_pairs $outdir/contacts/${sampleid}.contacts.pairs.txt.gz \
            --output $outdir/embedding/AB/$sampleid.AB_value.txt.gz \
            --weighted True \
            --ref_CpG $ref_cpg \
            --barcode_kept $downstream_cb \
            --cores $SLURM_CPUS_PER_TASK \
            --input_format "raw" 


        $Rscript $AB_embedding \
            $outdir/embedding/AB/$sampleid.AB_value.txt.gz \
            $outdir/QC/$sampleid.metadata.txt \
            "AB"

        echo "[scAB::] finished!"
        echo "Finshed" >$outdir/embedding/AB.done
    else 
        echo "[scAB::] finished!"
    fi 
}

function impute_scHiCluster(){
    dir=$1
    resolution=$2 # 100000
    cell_per_process=2000
    res_symbol=$(echo "scale=0; $resolution / 1000" | bc | xargs)"K"

    if [[ $resolution == 10000 ]];then 
        $hicluster prepare-impute \
            --cell_table $dir/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 2 \
            --cpu_per_job $SLURM_CPUS_PER_TASK \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 5050000 \
            --window_size 30000000 \
            --step_size 10000000 \
            --output_dir $dir/impute/$res_symbol/ \
            --chrom_size_path $chromSize \
            --resolution $resolution
    fi 

    if [[ $resolution == 100000 ]];then
        $hicluster prepare-impute \
            --cell_table $dir/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 1 \
            --cpu_per_job $SLURM_CPUS_PER_TASK \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 500000000 \
            --window_size 500000000 \
            --step_size 500000000 \
            --output_dir $dir/impute/$res_symbol/ \
            --chrom_size_path $chromSize \
            --resolution $resolution
    fi 
    
    if [[ $resolution == 25000 ]];then
        $hicluster prepare-impute \
            --cell_table $dir/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 1 \
            --cpu_per_job $SLURM_CPUS_PER_TASK \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 5050000 \
            --window_size 500000000 \
            --step_size 500000000 \
            --output_dir $dir/impute/$res_symbol/ \
            --chrom_size_path $chromSize \
            --resolution $resolution
    fi 
}


function scHiCluster_prepare(){
    hg38_blacklist_scHiCluster="/share/home/wangmx/software/scHiCluster/files/blacklist/hg38-blacklist.v2.bed.gz"
    mm10_blacklist_scHiCluster="/share/home/wangmx/software/scHiCluster/files/blacklist/mm10-blacklist.v2.bed.gz"
    species=$(basename -s ".fa" $reference)
    declare -A blacklists
    blacklists["hg38"]=$hg38_blacklist_scHiCluster
    blacklists["mm10"]=$mm10_blacklist_scHiCluster
    blacklist=${blacklists[$species]}

    if [[ ! -d $outdir/embedding/scHiCluster/$sampleid ]];then 
        mkdir -p $outdir/embedding/scHiCluster/$sampleid
    fi 
    
    scHiCluster_wd=$outdir/embedding/scHiCluster/$sampleid

    if [[ ! -d $scHiCluster_wd/sc_contacts ]];then 
        mkdir -p $scHiCluster_wd/sc_contacts
    fi 
    ## split into single cell contact pairs file
    $dscHiCtools splitContacts \
        --input_contacts $outdir/contacts/${sampleid}.contacts.pairs.txt.gz \
        --outdir $scHiCluster_wd/sc_contacts \
        --threads $SLURM_CPUS_PER_TASK \
        --barcode_file $downstream_cb

    ## filter contacts

    cat $downstream_cb \
        | awk -v p=$scHiCluster_wd/sc_contacts '{printf("%s/%s\n", p, $0".contacts.pairs.txt.gz")}' \
        > $scHiCluster_wd/all.contact_table.txt

    paste <(awk -F'/' '{print $NF}' $scHiCluster_wd/all.contact_table.txt | cut -d. -f1) $scHiCluster_wd/all.contact_table.txt | sort -k1,1 > $scHiCluster_wd/all.contact_table.tsv
    rm $scHiCluster_wd/all.contact_table.txt

    if [[ ! -d $scHiCluster_wd/rmbkl ]];then 
        mkdir -p $scHiCluster_wd/rmbkl
    fi 

    $hicluster filter-contact \
        --output_dir $scHiCluster_wd/rmbkl/ \
        --blacklist_1d_path $blacklist \
        --chr1 1 \
        --pos1 2 \
        --chr2 3 \
        --pos2 4 \
        --chrom_size_path $chromSize \
        --contact_table $scHiCluster_wd/all.contact_table.tsv
    
    ## impute in different resolution

    cat $downstream_cb \
        | awk -v p=$scHiCluster_wd/rmbkl '{printf("%s/%s\n", p, $0".contact.rmbkl.tsv.gz")}' > $scHiCluster_wd/contact_table_rmbkl.txt
    paste <(awk -F'/' '{print $NF}' $scHiCluster_wd/contact_table_rmbkl.txt | cut -d. -f1) $scHiCluster_wd/contact_table_rmbkl.txt | sort -k1,1 > $scHiCluster_wd/contact_table_rmbkl.tsv
    rm $scHiCluster_wd/contact_table_rmbkl.txt
    
    # for resolution in "100000" "25000" "10000";do
    for resolution in "10000";do

        impute_scHiCluster $scHiCluster_wd $resolution
    done 


}

function scVI-3d() {
    contacts_filter=$min_contactsN
    if [[ ! -n $contacts_filter ]]; then contacts_filter=7000; fi

    if [[ ! -d $outdir/embedding/scVI_3D/${sampleid} ]]; then
        mkdir -p $outdir/embedding/scVI_3D/${sampleid}
    fi

    python="/share/home/wangmx/anaconda3/envs/scvi-3d/bin/python"
    scVI_3D_maoxu="/share/home/wangmx/software/scVI-3D_maoxu/scripts/scVI-3D.py"

    # resolution=1000000
    resolution=500000

    bandMax=10

    if [[ -s $outdir/QC/$sampleid.metadata.txt ]]; then
        cat $outdir/QC/$sampleid.metadata.txt |
            awk '$4 < 0.45' |
            awk -v contacts_filter=$contacts_filter \
                '{if($3 > contacts_filter){print $1}}' \
                >$outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.txt
    else
        echo "Please run main process first!"
        exit
    fi

    ## prepare Input
    if [[ ! -s $outdir/contacts/${sampleid}.contacts.pairs.txt.gz ]]; then
        zcat $outdir/contacts/${sampleid}.contacts.pairs.gz |
            grep -v "#" |
            awk '{print $1,$2,$3,$4,$5}' OFS="\t" \
                | pigz >$outdir/contacts/${sampleid}.contacts.pairs.txt.gz
    fi

    if [[ ! -s $outdir/embedding/scVI_3D/${sampleid}/${sampleid}.scVI-3d.input.txt ]]; then
        $Rscript $scVI_input \
            $outdir/contacts/${sampleid}.contacts.pairs.txt.gz \
            $resolution \
            $outdir/embedding/scVI_3D/${sampleid} \
            $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.txt
    fi

    if [[ -s $outdir/splitBygenotype/clusters.cellline.tsv ]]; then
        grep -v "celline" $outdir/splitBygenotype/clusters.cellline.tsv |
            awk '{print $1, $2}' OFS="\t" | sort -k1,1 >$outdir/embedding/scVI_3D/${sampleid}/CellData_summary.tmp

        sort -k1,1 $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.txt >$outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.tmp

        join -j 1 $outdir/embedding/scVI_3D/${sampleid}/CellData_summary.tmp $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.tmp |
            awk 'BEGIN{print "name\tcell_type"}{print $1, $2}' OFS="\t" \
                >$outdir/embedding/scVI_3D/${sampleid}/CellData_summary.txt
        rm $outdir/embedding/scVI_3D/${sampleid}/CellData_summary.tmp
        rm $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.tmp
    else
        awk 'BEGIN{print "name\tcell_type"}{print $1, "Unknown"}' OFS="\t" $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.txt \
            >$outdir/embedding/scVI_3D/${sampleid}/CellData_summary.txt
    fi

    cpuN=1

    chromosomes="$(cat $autosome | sort -k1,1V | tr "\n" "," | xargs)"

    if [[ ! -s $outdir/embedding/scVI_3D/${sampleid}/scVI_3D.done ]]; then
        $python $scVI_3D_maoxu \
            -c $chromosomes \
            -r $resolution \
            -i $outdir/embedding/scVI_3D/${sampleid}/$sampleid.scVI-3d.input.txt \
            -o $outdir/embedding/scVI_3D/${sampleid} \
            -cs $outdir/embedding/scVI_3D/${sampleid}/CellData_summary.txt \
            -g $chromSize \
            -n 100 -p $cpuN -pca 50 -s -v -gpu
        echo "Finshed" >$outdir/embedding/scVI_3D/${sampleid}/scVI_3D.done
    fi

    ## downstream analysis
    $python $py_scAB \
        --input_contact_pairs $outdir/embedding/scVI_3D/${sampleid}/$sampleid.scVI_imputed.contacts.paris.txt.gz \
        --output $outdir/embedding/scVI_3D/${sampleid}/$sampleid.AB_value.txt.gz \
        --weighted True \
        --ref_CpG $ref_cpg \
        --barcode_kept $outdir/embedding/scVI_3D/${sampleid}/cellbarcodes.scVI.txt \
        --input_format "scVI-3d" \
        --cores $SLURM_CPUS_PER_TASK 

    $Rscript $AB_embedding \
        $outdir/embedding/scVI_3D/${sampleid}/$sampleid.AB_value.txt.gz \
        $outdir/QC/$sampleid.metadata.txt \
        "scVI"
}

function higashi() {
    res=$1
    contacts_filter=$min_contactsN
    if [[ ! -n $contacts_filter ]]; then contacts_filter=7000; fi
    # res="10K" #  100K, 200K, 1Mb
    # res="500K" #  100K, 250K, 1Mb
    calib_file="/share/home/wangmx/software/Higashi/calib/cpg_density.${res}.${species}.txt"
    if [[ ! -d $outdir/embedding/higashi_${res}/ ]]; then
        mkdir -p $outdir/embedding/higashi_${res}/
    fi
    python="/share/home/wangmx/anaconda3/envs/higashi/bin/python"
    fast_higashi="/share/home/wangmx/software/pipelines/dscHiC-process/scripts/run_FastHigashi.py"
    Process="/share/home/wangmx/software/Higashi/higashi/Process.py"
    main_cell="/share/home/wangmx/software/Higashi/higashi/main_cell.py"
    scCompartment="/share/home/wangmx/software/Higashi/higashi/scCompartment.py"
    scTAD="/share/home/wangmx/software/Higashi/higashi/scTAD.py"
    if [[ -s $outdir/QC/$sampleid.metadata.txt ]]; then
        cat $outdir/QC/$sampleid.metadata.txt \
            | awk '$4 < 0.45' \
            | awk -v contacts_filter=$contacts_filter \
                '{if($3 > contacts_filter){print $1}}' \
            | awk '{print $1, NR-1}' OFS="\t" > $outdir/embedding/higashi_${res}/${sampleid}.cellbarcodes.higashi.txt
    else
        echo "Please run main process first!"
        exit
    fi

    if [[ ! -s $outdir/contacts/${sampleid}.contacts.pairs.txt.gz ]]; then
        zcat $outdir/contacts/${sampleid}.contacts.pairs.gz |
            grep -v "#" |
            awk '{print $1,$2,$3,$4,$5}' OFS="\t" \
                | pigz >$outdir/contacts/${sampleid}.contacts.pairs.txt.gz
    fi

    if [[ ! -d $outdir/embedding/higashi_${res}/${sampleid} ]];then 
        mkdir -p $outdir/embedding/higashi_${res}/${sampleid}
    fi 

    if [[ ! -s $outdir/embedding/higashi_${res}/${sampleid}/data.txt ]]; then
        zcat $outdir/contacts/${sampleid}.contacts.pairs.txt.gz | sort -k1,1 >$outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.contacts.pairs.tmp

        join -j 1 $outdir/embedding/higashi_${res}/${sampleid}.cellbarcodes.higashi.txt \
            $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.contacts.pairs.tmp |
            tr " " "\t" |
            awk '$3 == $5' |
            $bedtools groupby -g 1,2,3,4,5,6 -c 6 -o count -i - |
            awk 'BEGIN{print "cell_id\tchrom1\tpos1\tchrom2\tpos2\tcount"}{print $2,$3,$4,$5,$6,$7}' OFS="\t" >$outdir/embedding/higashi_${res}/${sampleid}/data.txt

            rm $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.contacts.pairs.tmp
    fi

    ## prepare config
    chromosomes=$(cat $autosome | tr "\n" "," | sed 's/,/\\",\\"/g' | sed 's/,""//g' | xargs)

    # cat /share/home/wangmx/software/Higashi/template.JSON |
    cat /share/home/wangmx/software/Higashi/template_${res}.fh.JSON |
        sed "s|NAME|${sampleid}|" |
        sed "s|CHRLIST|${chromosomes}|" |
        sed "s|chromSIZE|${chromSize}|" |
        sed "s|DIRECTORY|$outdir/embedding/higashi_${res}/${sampleid}|" |
        sed "s|CYTOBAND|${cytoBand}|" |
        sed 's/,""//g' >$outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config

    ## runing this on gpu node
    if [[ ! -s $outdir/embedding/higashi_${res}/higashi.done ]]; then
        $python $fast_higashi \
            -c $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config \
            --cellbarcode $outdir/embedding/higashi_${res}/${sampleid}.cellbarcodes.higashi.txt \
            --filter \
            --do_rwr \
            --preprocess
        # $python $Process -c $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config
        $python $main_cell -c $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config
        echo "Finshed" >$outdir/embedding/higashi_${res}/higashi.done
        echo "Higashi finished!"
    else
        echo "Higashi finished!"
    fi 

    ## cpu is also OK for this
    if [[ ! -s $outdir/embedding/higashi_${res}/${sampleid}.AB.hdf5 ]];then 
        if [[ -n $calib_file ]]; then
            $python $scCompartment -c $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config \
                --calib_file $calib_file \
                --neighbor \
                -o ${sampleid}.AB.hdf5
        else
            echo "calib file not specified!"
            exit
        fi
    fi 

    ## downstream analysis
    $Rscript $AB_embedding \
        $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.AB.hdf5 \
        $outdir/QC/$sampleid.metadata.txt \
        "higashi"

    if [[ $res == "250K" ]];then
        $python $scTAD -c $outdir/embedding/higashi_${res}/${sampleid}/${sampleid}.config \
            --neighbor \
            -o ${sampleid}.TAD.hdf5
    fi 
}


scHiCluster_prepare

scAB

scVI-3d

higashi "500K"
higashi "250K"
