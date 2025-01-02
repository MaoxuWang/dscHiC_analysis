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
library(pheatmap)
library(pals)
options(scipen = 200)

args = commandArgs(T)
rds_path = args[1]
marker_file = args[2] # brain.marker_gene.Top100.txt
harmony = as.logical(args[3]) ## TRUE or FALSE
species = args[4]

sampleid = gsub(".rds", "", rds_path) # out/sampleid


### ----- function -----

#### draw plot
drawPlot <- function(Sample.harmony.filtered, harmony){
    pal <- jdb_palette("corona")
    pal <- pal[1:length(unique(Sample.harmony.filtered@meta.data$seurat_clusters))]

    p_cluster <- DimPlot(Sample.harmony.filtered, reduction = "umap",
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

    p_cluster <- p_cluster + labs(subtitle = paste0("Cell number: ", nrow(Sample.harmony.filtered@meta.data)))

    pal <- pal[1:length(unique(Sample.harmony.filtered@meta.data$celline))]

    pal <- pal[1:length(unique(Sample.harmony.filtered@meta.data$orig.ident))]

    if (harmony){

        atch_cluster <- DimPlot(Sample.harmony.filtered, reduction = "umap", group.by = "orig.ident", pt.size = 2.0, label = TRUE, label.size = 8, repel = TRUE, shuffle = TRUE, cols = pal) + 
                NoLegend() +
                FontSize(x.title = 16, y.title = 16) + 
                theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) +
                theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    }
    contactsN.plot <- FeaturePlot(object = Sample.harmony.filtered,
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


    Mitotic_band.plot <- FeaturePlot(object = Sample.harmony.filtered,
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
    Interchromosomal.plot <- FeaturePlot(object = Sample.harmony.filtered,
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
    
    if (harmony){
        final_plot.Harmony <- (p_cluster | batch_cluster | contactsN.plot | Mitotic_band.plot | Interchromosomal.plot )
    }else{
        final_plot.Harmony <- (p_cluster | contactsN.plot | Mitotic_band.plot | Interchromosomal.plot )
    }
    return(final_plot.Harmony)
}


#### calculate marker AB value
calAB_by_geneSets <- function(df, source_string, marker, species, binSize=1000000){
    feature_gene_map <- read_delim(paste0("/share/home/wangmx/database/GTF/", species, "/genomic_elements/genes.", species, ".region.bed"), col_names = F)
    colnames(feature_gene_map) <- c("chr", "start", "end", "gene")
    feature_gene_map$bin <- floor(((feature_gene_map$start + feature_gene_map$end ) / 2) / binSize) * binSize
    feature_gene_map$feature <- paste(feature_gene_map$chr, as.character(feature_gene_map$bin), sep = "-")
    feature_gene_map %>%
        filter(chr != "chrX") %>%
        filter(chr != "chrY") %>%
        filter(chr != "chrM") %>%
        dplyr::select(feature, gene)-> feature_gene_map
    
    feature_gene_map %>%
        inner_join(marker, by = "gene") %>%
        filter(source == source_string) -> marker.AB_feature

    print(head(marker.AB_feature))
    assay.data <- GetAssayData(object = df, slot = 'count')
    scAB.df <- as.data.frame(assay.data)
    scAB.df$feature <- rownames(scAB.df)
    
    scAB.df %>%
        inner_join(marker.AB_feature, by = "feature") %>%
        mutate(feature = NULL, gene=NULL, source = NULL) -> df.source

    df <- AddMetaData(
                object = df,
                metadata = data.frame("value"=apply(df.source, MARGIN = 2, FUN = mean) - mean(apply(df.source, MARGIN = 2, FUN = mean))),
                col.name = source_string
                )
    return(df)
}

#### plot AB heatmap
plotHeatmap_AB <- function(df, source, out_name){
    df@meta.data[,c("seurat_clusters",source)] %>%
        arrange(seurat_clusters)  -> temp

    colnames(temp) <- c("seurat_clusters",source)
    my_mat_num <- as.matrix.noquote(temp[,source])

    my_mat_num.scale <- t(my_mat_num)

    temp %>%
        dplyr::select(seurat_clusters) -> annotation_col_data

    break1 = seq(-0.025, 0, length.out=25)

    break2 = seq(0, 0.025, length.out=26)

    the_break <- c(break1, break2[-1])
    pheatmap(my_mat_num.scale, 
             border_color=NA,
             cluster_cols = F,
             cluster_row = F,
             show_colnames = F,  
             annotation_col = annotation_col_data,
             cellwidth = 0.2, 
             cellheight = 100,
             color = parula(n = 50),
             breaks = the_break,
             filename = out_name
            )
}


df <- readRDS(rds_path)

if (harmony){
    df %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>%
        ScaleData(verbose = FALSE) %>%
        RunPCA(npcs = 50, verbose = FALSE) -> df
        df <- suppressMessages({ harmony::RunHarmony(df, "orig.ident") })

        df %>%
            RunUMAP(dims = 1:10, verbose = F, min.dist = 0.2, n.neighbors = 50, reduction = 'harmony') %>%
            FindNeighbors(dims = 1:10, verbose = F, k.param = 50, reduction = 'harmony') %>%
            FindClusters(resolution = 0.2) -> df

        final.plot.Harmony <- drawPlot(df.P0)
        
        ggsave(final.plot.Harmony, filename = paste0(sampleid, ".Harmony.seurat_clusters.pdf"), dpi = 600, width = 48, height = 12)
        ggsave(final.plot.Harmony, filename = paste0(sampleid, ".Harmony.seurat_clusters.png"), dpi = 600, width = 48, height = 12)
        saveRDS(df, file = paste0(sampleid, ".Harmony.rds"))
}

marker <- read_delim(marker_file, col_names = T, delim = "\t")
for (source in unique(marker$source)){
    print(paste0(source, " is processing..."))
    df <- calAB_by_geneSets(df, source, marker, species)
    print("Done.")
}

heatmap_file_out = paste0(sampleid, ".AB.heatmap.pdf")


plotHeatmap_AB(df, unique(marker$source), heatmap_file_out)


