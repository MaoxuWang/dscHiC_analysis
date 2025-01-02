library(pbapply)
library(parallel)
library(tidyverse)
library(data.table)
library(scales)
library("unixtools")
library('bedtoolsr')
library(pheatmap)

args = commandArgs(T)
wd = args[1]
repli_chip_file = args[2]
expand_ = 10000

set.tempdir("/dshare/xielab/analysisdata/IVF/dashA/temp")
cl = 8

split_mtx <- function(matrix){
    ## split & merge contacts
    matrix.p1 <- matrix[, c("cellbarcode", "chrA", "startA")]
    colnames(matrix.p1) <- c("cellbarcode", "chrom", "site")
    matrix.p2 <- matrix[, c("cellbarcode", "chrB", "startB")]
    colnames(matrix.p2) <- c("cellbarcode", "chrom", "site")
    merged.mtx <- bind_rows(matrix.p1, matrix.p2)
    return(merged.mtx)
}

if (! "merged.contacts.Rdata" %in% list.files(wd)){
    merged.mtx = data.frame()
    for (input in list.files(path = wd, pattern = "*.contacts.pairs.txt.gz")){
        print(input)
        sample = gsub(".contacts.pairs.txt.gz", "",  input)
        suffix = gsub("20230427.mESC", "", sample) 
        cb = paste0(sample, ".singlecells.tsv")
        matrix <- read_delim(input, delim = "\t", col_names = F)
        
        # matrix <- read_delim(input, delim = "\t", col_names=F) 

        colnames(matrix) <- c("cellbarcode", "chrA", "startA","chrB", "startB")
        print(head(matrix))

        cellbarcode_f <- read_delim(cb, delim = "\t", col_names = F)
        colnames(cellbarcode_f) <- "cellbarcode"

        matrix %>%
            inner_join(cellbarcode_f) -> matrix
        
        matrix$cellbarcode <- paste0(matrix$cellbarcode, suffix)
        merged.mtx <- bind_rows(merged.mtx, matrix)
    }

    save(merged.mtx, file = "merged.contacts.Rdata")
}else{
    print("Data loading...")

    load("merged.contacts.Rdata")
}

### --contacts decay--
break_log = round(1000 * 2 ** (37:142 * 0.125), 0)

merged.mtx %>%
    filter(chrA == chrB) %>%
    mutate(distance = (startB - startA + 1)) %>%
    mutate(bin = cut(distance, breaks = break_log)) %>%
    filter(bin != "NA") %>%
    group_by(cellbarcode, bin) %>%
    summarise(count=n()) %>%
    ungroup() %>%
    group_by(cellbarcode) %>%
    summarise(bin=bin, count=count, valid=sum(count)) %>%
    ungroup() %>%
    mutate(freq = round(count / valid * 100, 2)) %>%
    select(cellbarcode, bin, freq) %>%
    pivot_wider(names_from = cellbarcode, values_from = freq, values_fill = 0) %>%
    separate(col = bin, sep = ",", into=c("start", "end")) %>%
    mutate(start = as.numeric(gsub("\\(", "", start))) %>%
    arrange(desc(start)) %>%
    select(-end) -> freq.table

merged.mtx %>%
    filter(chrA == chrB) %>%
    mutate(distance = (startB - startA + 1)) %>%
    filter(distance >= (1000 * 2 ** (37 * 0.125)),
            distance <= (1000 * 2 ** (142 * 0.125))) %>%
    group_by(cellbarcode) %>%
    summarise(near=mean((
            distance >= (1000 * 2 ** (37 * 0.125))
        ) &
        (
            distance <= (1000 * 2 ** (88 * 0.125))
        )) * 100,
            mitotic=mean((
            distance >= (1000 * 2 ** (89 * 0.125))
        ) &
        (
            distance <= (1000 * 2 ** (108 * 0.125))
        )) * 100) -> near_mitotic.percent


## farAvg
break_log_farAvgDist = round(1000 * 2 ** (97:142 * 0.125), 0)

merged.mtx %>%
    filter(chrA == chrB) %>%
    mutate(distance = (startB - startA + 1)) %>%
    mutate(bin = cut(distance, breaks = break_log_farAvgDist)) %>%
    filter(bin != "NA") %>%
    group_by(cellbarcode) %>%
    summarise(farAvgDist=mean(distance)) %>%
    select(cellbarcode, farAvgDist) -> farAvgDist.mtx

## repli-score
split.mtx <- split_mtx(merged.mtx)
split.mtx %>%
    mutate(start=site, end = site + 1) %>%
    select(chrom, start,end, cellbarcode) -> split.mtx


repli_chip <- read_delim(repli_chip_file, delim = "\t", col_names = F)
colnames(repli_chip) <- c("chr", "start", "end","signal")
repli_chip$start <- repli_chip$start - expand_
repli_chip$end <- repli_chip$end + expand_
result <- bedtoolsr::bt.intersect(split.mtx, repli_chip, wa = T, wb = T) 

result %>%
    group_by(V4) %>%
    summarise(raw_repli_score = mean(V8 > 0) ) -> repli_score.mtx
colnames(repli_score.mtx) <- c("cellbarcode", "raw_repli_score")

near_mitotic.percent %>%
    inner_join(farAvgDist.mtx) %>%
    inner_join(repli_score.mtx) -> cellbarcode.repli

cellbarcode.repli %>%
    mutate(group = case_when(
    ( mitotic >= 30 & near <= 50 ) ~ "Post-M",
    ( near > 50 & (near + 1.8 * mitotic) > 100 ) ~ "Pre-M",
    ( near <= 63 ) ~ "G1",
    ( near > 63 & near <= 78.5 ) ~ "Early-S",
    ( near > 78.5 ) ~ "Late S-G2",
    TRUE ~ "blank")) -> cellbarcode.repli

## ordering
cellbarcode.repli %>%
    filter(group == "Post-M") %>%
    arrange(desc(mitotic)) -> Post_M.repli

cellbarcode.repli %>%
    filter(group == "G1") %>%
    mutate(order_index = scale(near) + scale(farAvgDist)) %>%
    arrange(order_index) %>%
    select(-order_index) -> G1.repli

cellbarcode.repli %>%
    filter(group == "Early-S") %>%
    mutate(order_index = (near / var(near) + raw_repli_score / var(raw_repli_score))) %>%
    arrange(order_index) %>%
    select(-order_index) -> early_mid_S.repli

cellbarcode.repli %>%
    filter(group == "Late S-G2") %>%
    mutate(order_index = (near / var(near) - raw_repli_score / var(raw_repli_score))) %>%
    arrange(order_index) %>%
    select(-order_index) -> mid_S_G2.repli

cellbarcode.repli %>%
    filter(group == "Pre-M") %>%
    arrange(mitotic) -> Pre_M.repli

ordered.mtx <- bind_rows(
    Post_M.repli,
    G1.repli,
    early_mid_S.repli,
    mid_S_G2.repli,
    Pre_M.repli
)

freq.table.ordered <- freq.table[,ordered.mtx$cellbarcode]
freq.table.ordered$start <- freq.table$start

write.table(cellbarcode.repli, file="mESC.repli_score.txt", quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)

write.table(freq.table.ordered, file="mESC.contacts_frequency.ordered.mtx", quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)

# ## plot
df@meta.data[,c("seurat_clusters",source)] %>%
    arrange(seurat_clusters)  -> temp

colnames(temp) <- c("seurat_clusters",source)
my_mat_num <- as.matrix.noquote(temp[,source])

my_mat_num.scale <- t(my_mat_num)

cellbarcode.repli %>%
    dplyr::select(cellbarcode, group) -> annotation_col_data

rownames(annotation_col_data) <- annotation_col_data$cellbarcode


break1 = seq(-0.025, 0, length.out=25)

break2 = seq(0, 0.025, length.out=26)

the_break <- c(break1, break2[-1])
pheatmap(freq.table.ordered, 
        border_color=NA,
        cluster_cols = F,
        cluster_row = F,
        show_colnames = F,  
        annotation_col = annotation_col_data,
        cellwidth = 0.2, 
        cellheight = 100,
        color = parula(n = 50),
        # breaks = the_break,
        filename = out_name
        )