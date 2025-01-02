library(tidyverse)
library(umap)
library(tsne)
library("bnstruct")
library(Seurat)
library(paletteer)
library(ggsci)
library(BuenColors)
library(RColorBrewer)
library(rhdf5)


args = commandArgs(T)
AB_matrix = args[1]

metadata = args[2]
mode = args[3] ## whether from higashi

# 1. prepare metadata
metadata <- read_delim(metadata, 
                col_names=T,
                delim="\t")
metadata %>%
    mutate(Intra_ge20kb = Intra_ge20kb * 100,
            Interchromosomal = Interchromosomal * 100,
            Intra_le20kb = Intra_le20kb * 100,
            Mitotic_band = Mitotic_band * 100,
            Non_mitotic_band = Non_mitotic_band * 100
    ) -> metadata
metadata$cis_over_trans <- metadata$Intra_ge20kb / metadata$Interchromosomal

calMissing <- function(line){
    return(sum(is.na(line)) / (length(line)-1))
}

scaleOne <- function(x){
    return((x - 1) / (sum(!is.na(x)) - 1) )
}


# 2. read inpute matrix
## higashi: z-score AB comparet value
## AB and scVI-3d_imputed AB: normalized count averaged AB value
print(paste0("Processing AB matrix in ", mode, "..."))
if (mode == "higashi" ){
    sampleid = gsub(".AB.hdf5", "", AB_matrix) # out/sampleid
    sampleName = gsub(".*\\/", "", sampleid)
    wd = gsub(paste0("/",sampleName), "", sampleid)

    h5f = H5Fopen(AB_matrix)

    temp = as.data.frame(h5f$'compartment_zscore')
    print("input matrix loaded!")
    temp$bin_index <- paste(temp$bin.chrom, temp$bin.start, sep = "_")
    rownames(temp) <- temp$bin_index

    temp %>%
        mutate(bin.chrom=NULL, bin.start=NULL, bin.end=NULL, bin_index=NULL) -> df.fillNA

    cb = paste0(wd, "/", sampleName, ".cellbarcodes.higashi.txt")

    if(file.exists(cb)){
        dict_cb <- read_delim(cb, col_names = F,delim = "\t")
        colnames(dict_cb) <- c("cellbarcode", "ID")
    }else{
        stop(paste0("cellbarcode file not found at ", cb))
    }
    
    data.frame(ID=as.double(gsub('cell_', "", colnames(df.fillNA)))) %>%
        inner_join(dict_cb) %>%
        mutate(ID=paste0("cell_", ID)) -> new

    df.fillNA = df.fillNA[,new$ID]
    colnames(df.fillNA) <- new$cellbarcode
}else{
    sampleid = gsub(".AB_value.txt.gz", "", AB_matrix) # out/sampleid
    sampleName = gsub(".*\\/", "", sampleid)
    wd = gsub(paste0("/",sampleName), "", sampleid)

    AB <- read_delim(AB_matrix, 
                    col_names=F,
                    delim="\t")
    print("input matrix loaded!")
    colnames(AB) = c("cell", "bin", "value")

    AB %>%
        filter(!str_detect(bin, "chrM")) %>%
        filter(!str_detect(bin, "chrY")) %>%
        filter(!str_detect(bin, "chrX")) %>%
        pivot_wider(names_from = "bin", values_from = "value", values_fill =NA) -> full_AB

    ## retain autosmal
    auto_AB <- full_AB[,c(1,grep("chr",colnames(full_AB)))]

    # save(auto_AB, 
    #     file=paste0(wd, "/", sampleName, ".auto_AB.Rdata")
    # )


    cellNames <- auto_AB$cell
    auto_AB <- auto_AB[,-1]

    df.rank <- apply(auto_AB, 1, rank, na.last = "keep")

    df.rank.scaled <- apply(df.rank, 2, FUN = scaleOne)

    df <- as.data.frame(df.rank.scaled)

    df %>%
        replace(is.na(.), 0) -> df.fillNA

    colnames(df.fillNA) <- cellNames

}


# 3. Seurat analysis 

### ---- Only retain bins that is missing at most 1% of cells
# genes.filter <- rownames(df.fillNA[which(rowSums(df.fillNA == 0) / ncol(df.fillNA) <= 0.01 ),])

# df.fillNA <- df.fillNA[genes.filter,]
##### 

dscHiC <- CreateSeuratObject(counts = df.fillNA, paste0(sampleName, "_AB"))

dscHiC@meta.data$cellbarcode <- rownames(dscHiC@meta.data)
dscHiC@meta.data %>%
    left_join(metadata) -> dscHiC@meta.data

rownames(dscHiC@meta.data) <- dscHiC@meta.data$cellbarcode

dscHiC <- subset(dscHiC, Interchromosomal != "NA")


saveRDS(dscHiC, file = paste0(sampleid, ".", mode, ".AB.raw.rds"))

qc_plot <- VlnPlot(dscHiC, features = c("nFeature_RNA", "nCount_RNA","Interchromosomal", "Mitotic_band"), ncol = 4)

ggsave(qc_plot, filename = paste0(sampleid, ".", mode, ".AB.QC.png"), dpi = 600, width = 40, height = 10)
ggsave(qc_plot, filename = paste0(sampleid, ".", mode, ".AB.QC.pdf"), dpi = 600, width = 40, height = 10)

### filter inter > 45%
dscHiC.filtered <- subset(dscHiC, Interchromosomal <= 45 )

dscHiC.filtered <- FindVariableFeatures(dscHiC.filtered, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
dscHiC.filtered <- ScaleData(dscHiC.filtered, verbose = FALSE)
dscHiC.filtered <- RunPCA(dscHiC.filtered, npcs = 50, verbose = FALSE)

elbow.plot <- ElbowPlot(dscHiC.filtered, ndims = 50)
ggsave(elbow.plot, filename = paste0(sampleid, ".", mode, ".AB.Elowbow.pdf"), dpi = 600, width = 12, height = 10)
ggsave(elbow.plot, filename = paste0(sampleid, ".", mode, ".AB.Elowbow.png"), dpi = 600, width = 12, height = 10)

### UMAP adjusted parameters
dscHiC.filtered <- RunUMAP(dscHiC.filtered, dims = 1:10, verbose = F, min.dist = 0.2, n.neighbors = 50)
dscHiC.filtered <- FindNeighbors(dscHiC.filtered, dims = 1:10, verbose = F, k.param = 50)
dscHiC.filtered <- FindClusters(dscHiC.filtered, resolution = 0.2)

pal <- jdb_palette("corona")
pal <- pal[1:length(unique(dscHiC.filtered@meta.data$seurat_clusters))]

p_cluster <- DimPlot(dscHiC.filtered, reduction = "umap",
                group.by = "seurat_clusters", 
                pt.size = 2.0, 
                label = TRUE, 
                label.size = 8, 
                repel = TRUE, 
                shuffle = TRUE, cols = pal) +
            NoLegend() + 
            theme(axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_text(size=16),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.title.x = element_text(size=16)) +
            theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                        
p_cluster <- p_cluster + labs(subtitle = paste0("Cell number: ", nrow(dscHiC.filtered@meta.data)))

contactsN.plot <- FeaturePlot(object = dscHiC.filtered,
                        features = c("contactsN"),
                        pt.size=2.0,
                        reduction="umap",
                        ncol=1,
                        label=F, 
                        cols = c("lightgrey", "firebrick3")) +
                        theme(axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_text(size=16),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.title.x = element_text(size=16)) +
                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

Mitotic_band.plot <- FeaturePlot(object = dscHiC.filtered,
                        features = c("Mitotic_band"),
                        pt.size=2.0,
                        reduction="umap",
                        ncol=1,
                        label=F, 
                        cols = c("lightgrey", "firebrick3")) +
                        theme(axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_text(size=16),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.title.x = element_text(size=16)) +
                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
Interchromosomal.plot <- FeaturePlot(object = dscHiC.filtered,
                        features = c("Interchromosomal"),
                        pt.size=2.0,
                        reduction="umap",
                        ncol=1,
                        label=F, 
                        cols = c("lightgrey", "firebrick3")) +
                        theme(axis.text.y = element_blank(), 
                              axis.ticks.y = element_blank(),
                              axis.title.y = element_text(size=16),
                              axis.ticks.x = element_blank(),
                              axis.text.x = element_blank(),
                              axis.title.x = element_text(size=16)) +
                        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

singlets.plot <-  p_cluster | contactsN.plot | Mitotic_band.plot | Interchromosomal.plot
ggsave(singlets.plot, filename = paste0(sampleid, ".", mode, ".AB.seurat_clusters.pdf"), dpi = 600, width = 40, height = 12)
ggsave(singlets.plot, filename = paste0(sampleid, ".", mode, ".AB.seurat_clusters.png"), dpi = 600, width = 40, height = 12)

saveRDS(dscHiC.filtered, file = paste0(sampleid, ".", mode, ".AB.filtered.rds"))