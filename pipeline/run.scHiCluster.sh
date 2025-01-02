#!/usr/bin/bash
set -e
set -o pipefail


hicluster="/share/home/wuhg/miniforge3/envs/schicluster/bin/hicluster"
mm10_blacklist="/share/home/wangmx/software/scHiCluster/files/blacklist/mm10-blacklist.v2.bed.gz"
chrom_size="/share/home/wangmx/database/chromosomeSize/mm10.autosome.sizes"
dscHiCtools="/share/home/wangmx/anaconda3/bin/dscHiCtools"
ddate=$(date -R | awk '{print $2$3}')

cell_per_process=2000
contacts_filter=10000
wd=`pwd`


function filter_contacts(){
    barcode_file=$1
    outdir=$2

    cat $barcode_file | awk -v p=$outdir/sc_contacts '{printf("%s/%s\n", p, $0".contacts.pairs.txt.gz")}' > $outdir/chunk_id/all.contact_table.txt
    paste <(awk -F'/' '{print $NF}' $outdir/chunk_id/all.contact_table.txt | cut -d. -f1) $outdir/chunk_id/all.contact_table.txt | sort -k1,1 > $outdir/chunk_id/all.contact_table.tsv
    rm $outdir/chunk_id/all.contact_table.txt


    if [[ ! -d $outdir/rmbkl ]];then 
        mkdir -p $outdir/rmbkl
    fi 
    thread=40
    job="rmbkl_all"
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=160g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00

    $hicluster filter-contact \
        --output_dir $outdir/rmbkl/ \
        --blacklist_1d_path $mm10_blacklist \
        --chr1 1 \
        --pos1 2 \
        --chr2 3 \
        --pos2 4 \
        --chrom_size_path $chrom_size \
        --contact_table $outdir/chunk_id/all.contact_table.tsv
RUN
}


function prepare_impute(){
    outdir=$1
    resolution=$2 # 100000
    barcode_file=$outdir/chunk_id/cellbarcodes.scHiCluster.txt

    res_symbol=$(echo "scale=0; $resolution / 1000" | bc | xargs)"K"
    cpu_per_job=40
    cat $barcode_file | awk -v p=$outdir/rmbkl '{printf("%s/%s\n", p, $0".contact.rmbkl.tsv.gz")}' > $outdir/chunk_id/contact_table_rmbkl.txt
    paste <(awk -F'/' '{print $NF}' $outdir/chunk_id/contact_table_rmbkl.txt | cut -d. -f1) $outdir/chunk_id/contact_table_rmbkl.txt | sort -k1,1 > $outdir/chunk_id/contact_table_rmbkl.tsv

    if [[ $resolution == 10000 ]];then 
        $hicluster prepare-impute \
            --cell_table $outdir/chunk_id/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 2 \
            --cpu_per_job $cpu_per_job \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 5050000 \
            --window_size 30000000 \
            --step_size 10000000 \
            --output_dir $outdir/impute/$res_symbol/ \
            --chrom_size_path $chrom_size \
            --resolution $resolution
    fi 

    if [[ $resolution == 100000 ]];then
        $hicluster prepare-impute \
            --cell_table $outdir/chunk_id/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 1 \
            --cpu_per_job $cpu_per_job \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 500000000 \
            --window_size 500000000 \
            --step_size 500000000 \
            --output_dir $outdir/impute/$res_symbol/ \
            --chrom_size_path $chrom_size \
            --resolution $resolution
    fi 
    
    if [[ $resolution == 25000 ]];then
        $hicluster prepare-impute \
            --cell_table $outdir/chunk_id/contact_table_rmbkl.tsv \
            --batch_size $cell_per_process \
            --pad 1 \
            --cpu_per_job $cpu_per_job \
            --chr1 1 \
            --pos1 2 \
            --chr2 3 \
            --pos2 4 \
            --output_dist 10050000 \
            --window_size 500000000 \
            --step_size 500000000 \
            --output_dir $outdir/impute/$res_symbol/ \
            --chrom_size_path $chrom_size \
            --resolution $resolution
    fi 

}


function embedding(){
    indir=$1 #$outdir/impute/$res_symbol
    resolution=$2 # 100000

    # Compute embedding
    find $indir/ | grep cool > $indir/cell_table.txt
    paste <(awk -F'/' '{print $NF}' $indir/cell_table.txt | cut -d. -f1) $indir/cell_table.txt | sort -k1,1 > $indir/cell_table.tsv
    rm $indir/cell_table.txt

    if [[ ! -d $indir/embedding ]];then 
        mkdir -p $indir/embedding
    fi 
    thread=20
    job="embedding"
    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=110g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00


    $hicluster embedding \
        --cell_table_path $indir/cell_table.tsv \
        --output_dir $indir/embedding \
        --dim 50 \
        --dist 1000000 \
        --resolution $resolution \
        --scale_factor 100000 \
        --norm_sig \
        --save_raw \
        --cpu 20
RUN

}


function gene_score(){
    indir=$1 # 
    sampleid=$2

    # Compute gene contact score
    find $indir/ | grep cool > $indir/cell_table.txt 
    paste <(awk -F'/' '{print $NF}' $indir/cell_table.txt  | cut -d. -f1) $indir/cell_table.txt  | sort -k1,1 > $indir/cell_table.tsv 
    rm $indir/cell_table.txt 
    job="impute_gene_score"
    thread=80
    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $job
#SBATCH -o logs/${job}.${ddate}.log
#SBATCH -e logs/${job}.${ddate}.log
#SBATCH --mem=400g
#SBATCH --cpus-per-task=$thread
#SBATCH --partition=compute_fat
#SBATCH -t 7-00:00:00


    $hicluster gene-score \
        --cell_table_path $indir/cell_table.tsv  \
        --gene_meta_path /share/home/wangmx/database/GTF/mm10/gencode.vM32.annotation.gene.sorted.bed.gz \
        --resolution 10000 \
        --output_hdf_path $indir/${sampleid}.geneimputescore.hdf \
        --chrom_size_path $chrom_size \
        --cpu $thread \
        --mode impute
RUN

}

ls output/QC/*.metadata.txt | while read metadata;do
    sampleid=$(basename -s ".metadata.txt" $metadata)
    echo "$sampleid is processing..."
    outdir=$wd/output/embedding/scHiCluster/$sampleid
    # prepare $metadata $outdir
    # filter_contacts $outdir/chunk_id/cellbarcodes.scHiCluster.txt $outdir
    # resolution=100000
    for resolution in "100000" "25000" "10000";do
        res_symbol=$(echo "scale=0; $resolution / 1000" | bc | xargs)"K" 
        prepare_impute $outdir $resolution
        cat $outdir/impute/$res_symbol/snakemake_cmd.txt | while read CMD; do
            echo $CMD
            sbatch --wrap="$CMD" \
                --mem=120g \
                --cpus-per-task=40 \
                -t 7-00:00:00 \
                -o logs/snakemake_${res_symbol}.${ddate}.log \
                -e logs/snakemake_${res_symbol}.${ddate}.log \
                -J snakemake_impute_${res_symbol} \
                --partition=compute_fat
        done 
    done 
    # embedding $outdir/impute/$res_symbol $resolution
    # gene_score $outdir/impute/10K $sampleid
done 
