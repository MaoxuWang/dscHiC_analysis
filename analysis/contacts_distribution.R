library(pbapply)
library(parallel)
library(tidyverse)
library(data.table)
library(scales)
library("unixtools")


args = commandArgs(T)
input = args[1] # contact.pairs.gz
wd = args[2] # working directory
cb = args[3] # called cell 

sampleid = gsub(".contacts.pairs.gz", "", gsub(".*\\/", "", input))
set.tempdir("/dshare/xielab/analysisdata/IVF/dashA/temp")

cl = 8

sample = paste0(wd, "/", sampleid)

assignContacts <- function(array_n){
    chrA = array_n[2]
    chrB = array_n[4]
    startA = as.double(array_n[3])
    startB = as.double(array_n[5])

    if(chrA == chrB){
        distance = abs(startA - startB)
        if (distance <= 20000){
            return("Intra (<20 kb)")
        }else{
            return("Intra (>20 kb)")
        }
    }else{
        return("Interchromosomal")
    }
}

mitoticContacts <- function(array_n){
    chrA = array_n[2]
    chrB = array_n[4]
    startA = as.double(array_n[3])
    startB = as.double(array_n[5])

    if(chrA == chrB){
        distance = abs(startA - startB)
        if (distance >= 2000000 & distance <= 12000000){
            return("Mitotic_band")
        }else{
            return("Non_mitotic_band")
        }
    }else{
        return("Non_mitotic_band")
    }
}

matrix.raw <- read_delim(input, delim = "\t", skip=27, col_names = F)
# matrix <- read_delim(input, delim = "\t", col_names=F) 

colnames(matrix.raw) <- c("cellbarcode", "chrA", "startA","chrB", "startB", "minus", "plus")
print(head(matrix.raw))
kept_barcode <- fread(cb, sep = "\t", header=F)
colnames(kept_barcode) <- c("cellbarcode")
print("cell barcode loaded!")

matrix.raw %>%
    inner_join(kept_barcode) -> matrix


matrix$assign <- apply(matrix, 1, FUN = assignContacts)

print(head(matrix))
matrix %>%
    group_by(cellbarcode) %>%
    count(assign) %>%
    ungroup() %>%
    group_by(cellbarcode) %>%
    pivot_wider(names_from = "assign", values_from = n) -> df.mtx 

df.mtx <- as.data.frame(df.mtx)
rownames(df.mtx) <- df.mtx$cellbarcode


df.mtx <- df.mtx[-1]

df.mtx.distribution <-  df.mtx / rowSums(df.mtx)

# write.table(df.mtx.distribution, file=paste0(sample, ".contacts_distribution.txt"), quote=FALSE, sep = "\t",row.names = TRUE,col.names = TRUE)

df.mtx.distribution %>%
    pivot_longer(cols = colnames(df.mtx.distribution), names_to = "assign", values_to = "n") -> df.mtx.distribution.plot

## draw Plot

plot.distribution <- ggplot(df.mtx.distribution.plot, aes(x=factor(assign, levels = c("Intra (<20 kb)", "Intra (>20 kb)", "Interchromosomal")), y=n)) +
geom_violin(position = position_dodge(width=0.8), scale = "area") + 
geom_boxplot(position = position_dodge(width=0.8),
        outlier.size = 0, width = 0.1, show.legend = F)+
scale_y_continuous(limits = c(0,0.8),labels = scales::percent) + 
# geom_text(aes(label=paste0(round(mean,2),"%"), y=mean-0.15), position=position_dodge(0.9), vjust=-0.4) +
# scale_fill_manual(values=c("euploidy"='#f6ad49',"aneuploidy"='#84a2d4')) +
ylab("% Contacts")+ xlab("") + labs("")+
theme_bw()+
theme(panel.grid = element_line(colour = 'white'),
    # legend.position = c(0.82,0.82),
    # legend.title = element_blank(),
        axis.title.y=element_text(size=16),
        axis.text.y=element_text(size=10,colour = "black"),
        axis.text.x=element_text(size=16,colour = "black"))

ggsave(plot.distribution, filename = paste0(sample, ".contacts_distribution.violinBox_plot.png"), width = 210, height = 200,units = "mm", dpi = 600)

ggsave(plot.distribution, filename = paste0(sample, ".contacts_distribution.violinBox_plot.pdf"), width = 210, height = 200,units = "mm", dpi = 600)


## mitotic band

matrix.raw %>%
    inner_join(kept_barcode) -> matrix.mito

matrix.mito$assign <- apply(matrix.mito, 1, FUN = mitoticContacts)

 ## mito
        
matrix.mito %>%
    group_by(cellbarcode) %>%
    count(assign) %>%
    ungroup() %>%
    group_by(cellbarcode) %>%
    pivot_wider(names_from = "assign", values_from = n) -> df.matrix.mito 

df.matrix.mito  <- as.data.frame(df.matrix.mito )
rownames(df.matrix.mito ) <- df.matrix.mito $cellbarcode


df.matrix.mito  <- df.matrix.mito [-1]

df.matrix.mito.distribution <-  df.matrix.mito  / rowSums(df.matrix.mito )
df.mtx.distribution <- as.data.frame(df.mtx.distribution)

df.matrix.mito.distribution <- as.data.frame(df.matrix.mito.distribution)

# write.table(df.matrix.mito.distribution, file=paste0(sample, ".mitotic_band.txt"), quote=FALSE, sep = "\t",row.names = TRUE,col.names = TRUE)
df.mtx.distribution$cellbarcode <- rownames(df.mtx.distribution)
df.matrix.mito.distribution$cellbarcode <- rownames(df.matrix.mito.distribution)

df.mtx.distribution %>%
    inner_join(df.matrix.mito.distribution) %>%
    arrange(cellbarcode) %>%
    select(cellbarcode, Interchromosomal, "Intra (<20 kb)", "Intra (>20 kb)", Mitotic_band, Non_mitotic_band) -> df.result

write.table(df.result, file=paste0(sample, ".contacts_distribution.txt"), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
