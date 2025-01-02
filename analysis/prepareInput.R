library(tidyverse)
args = commandArgs(T)
input = args[1]
res = args[2]
outdir = args[3] 
valid_cb = args[4]

sampleid = gsub("\\.contacts.pairs.txt.gz", "", gsub(".*\\/", "", input))



res = as.numeric(res)

Input2scVI_3D <- function(input, df_cb, res=1000000){
    ## Input: .paris.txt
    df <- read_delim(input, col_names=F,delim="\t")
    colnames(df) <- c("cell","chromA","binA","chromB", "binB")
    df$binA <- floor(df$binA / res)*res
    df$binB <- floor(df$binB / res)*res

    df %>%
        inner_join(df_cb) %>%
        filter(str_detect(chromA, "chr")) %>%
        filter(str_detect(chromB, "chr")) %>%
        filter(chromA == chromB) %>%
        group_by(cell, chromA, binA, chromB, binB) %>%    
        summarise(count = n()) %>%
        ungroup() -> count_df
    return(count_df)

}


df_cb <- read_delim(valid_cb, col_names=F,delim="\t")
colnames(df_cb) <- "cell"
print("cell barcode whitelist is loaded.")
df.total <- Input2scVI_3D(input, df_cb, res)
print("Total counts matrix is ready!")


write.table( df.total, file=paste0(outdir, "/", sampleid, ".scVI-3d.input.txt"), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)

