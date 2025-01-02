library(ComplexHeatmap)
library(circlize)
library(pals)
library(rhdf5)
library(foreach)
library(ggrastr)
library(bedtoolsr)
source("/dshare/xielab/analysisdata/IVF/wangmx/dscHi-C/scripts/plot_function.R")

## calculate gene sets scAB value


queryGenesByRegion <- function(region_df, species, by="gene_body", extend=NA){
    options(scipen = 200)
    feature_gene_map <- read_delim(paste0("/share/home/wangmx/database/GTF/", species, "/gencode.vM25.annotation.gene.sorted.bed.gz"), col_names = F)
    colnames(feature_gene_map) <- c("chr", "start", "end", "gene", "strand")
    feature_gene_map %>% filter(!duplicated(gene)) -> feature_gene_map
    if (by == "gene_body"){
        feature_gene_map %>%
            mutate(gene_length = end - start) %>%
            dplyr::select(chr, start, end, gene, gene_length) -> feature_gene_map
    }else if(by == "TSS"){
        feature_gene_map %>%
            mutate(gene_length = end - start) %>%
            mutate(TSS = if_else(strand == "+", start, end)) %>%
            mutate(start = TSS - 1, end = TSS + 1) %>%
            dplyr::select(chr, start,end,gene, gene_length) -> feature_gene_map
    }else if(by == "TES"){
        feature_gene_map %>%
            mutate(gene_length = end - start) %>%
            mutate(TES = if_else(strand == "+", end, start)) %>%
            mutate(start = TES - 1, end = TES + 1) %>%
            dplyr::select(chr ,start, end, gene, gene_length) -> feature_gene_map
    }else{
        stop("Wrong method specified ! {gene_body, TSS, TES}")
    }
    if(! is.na(extend)){
        region_df %>%
            mutate(start = start - extend, end = end + extend) -> region_df
    }
    feature_gene_map %>%
        bedtoolsr::bt.intersect(region_df, wa = TRUE, wb = TRUE) -> gene.list

    return(gene.list)
}



calAB_by_geneSets <- function(df, source_string, marker, species, method="heatmap", binSize=1000000, assay = "RNA"){
    options(scipen = 200)
    feature_gene_map <- read_delim(paste0("/share/home/wangmx/database/GTF/", species, "/genomic_elements/genes.", species, ".region.bed"), col_names = F)
    colnames(feature_gene_map) <- c("chr", "start", "end", "gene")
    feature_gene_map$bin <- floor(((feature_gene_map$start + feature_gene_map$end ) / 2) / binSize) * binSize
    feature_gene_map$feature <- paste(feature_gene_map$chr, as.character(feature_gene_map$bin), sep = "-")
    feature_gene_map %>%
        filter(chr != "chrX") %>%
        filter(chr != "chrY") %>%
        filter(chr != "chrM") %>%
        dplyr::select(feature, gene)-> feature_gene_map
    
    colnames(marker) <- c("gene", "source")
    feature_gene_map %>%
        inner_join(marker, by = "gene") %>%
        filter(source == source_string) -> marker.AB_feature

    print(head(marker.AB_feature))
    assay.data <- GetAssayData(object = df, slot = 'count', assay = assay)
    scAB.df <- as.data.frame(assay.data)
    scAB.df$feature <- rownames(scAB.df)
    
    scAB.df %>%
        inner_join(marker.AB_feature, by = "feature") %>%
        mutate(feature = NULL, gene=NULL, source = NULL) -> df.source

    if(method == "heatmap"){
        df <- AddMetaData(
                    object = df,
                    metadata = data.frame("value"=apply(df.source, MARGIN = 2, FUN = mean) - mean(apply(df.source, MARGIN = 2, FUN = mean))),
                    col.name = source_string
                    )

    }else{
        df <- AddMetaData(
                object = df,
                metadata = data.frame("value"=apply(df.source, MARGIN = 2, FUN = mean)),
                col.name = source_string
                )
    }

    return(df)
}



scaleOne <- function(x){
    return((x - 1) / (sum(!is.na(x)) - 1) )
}

cal_na <- function(x){
    return(sum(is.na(x)) / length(x))
}


createscABObj <- function(scAB_value_file, min_missing_feature = 0.6, mode='layer', assay='RNA'){

    scAB.mtx <- read_delim(scAB_value_file, delim = "\t")
    colnames(scAB.mtx) = c("cell", "bin", "value")

    scAB.mtx %>%
        filter(!str_detect(bin, "chrM")) %>%
        filter(!str_detect(bin, "chrY")) %>%
        filter(!str_detect(bin, "chrX")) %>%
        pivot_wider(names_from = "bin", values_from = "value", values_fill =NA) -> full_AB

    cellNames <- full_AB$cell
    full_AB <- full_AB[,-1]

    df.rank <- apply(full_AB, 1, rank, na.last = "keep")

    df.rank.scaled <- apply(df.rank, 2, FUN = scaleOne)

    df <- as.data.frame(df.rank.scaled)
    df <- df[apply(df, 1, cal_na) < min_missing_feature, ]
    df %>%
        replace(is.na(.), 0) -> df.fillNA
    colnames(df.fillNA) <- cellNames

    if(mode == "layer"){
        obj <- CreateAssayObject(counts = df.fillNA)
    }else{
        obj <- CreateSeuratObject(counts = df.fillNA, assay = assay)
    }
    return(obj)
}


convertFormat <- function(temp){
    temp$bin_index <- paste(temp$bin.chrom, temp$bin.start, sep = "_")
    rownames(temp) <- temp$bin_index

    temp %>%
        mutate(bin.chrom=NULL, bin.start=NULL, bin.end=NULL, bin_index=NULL) -> out.df
    return(out.df) 
}


convertID <- function(obj.df, dict_cb){

    data.frame(ID=as.double(gsub('cell_', "", colnames(obj.df)))) %>%
            inner_join(dict_cb) %>%
            mutate(ID=paste0("cell_", ID)) -> new

    obj.df = obj.df[,new$ID]
    colnames(obj.df) <- new$cellbarcode 
    return(obj.df)
}


## calculate bin-level scAB or insulation score by higashi imputation
createHigashiObj <- function(
    HDF5_file,
    cb, #".cellbarcodes.higashi.txt"
    mode="compartment",
    method="raw",
    assay=NULL
){
    print(paste0(mode, "as input."))

    h5f = H5Fopen(HDF5_file)  
    print("input matrix loaded!")
    if(file.exists(cb)){
        dict_cb <- read_delim(cb, col_names = F,delim = "\t")
        colnames(dict_cb) <- c("cellbarcode", "ID")
    }else{
        stop(paste0("cellbarcode file not found at ", cb))
    }

    if (mode == "compartment"){
        if (method == "z_score"){
            print("Using Z-score")
            compartment.df = as.data.frame(h5f$'compartment_zscore') %>% convertFormat()
            print("input matrix loaded!")

        }else if(method == "raw"){
            print("Using raw score")
            compartment.df = as.data.frame(h5f$'compartment_raw') %>% convertFormat()
        }else{
            stop("Wrong method specified ! {z_score, raw} ")
        }

        compartment.df <- convertID(compartment.df, dict_cb)
        if (! is.null(assay)){
            layer_compartment <- CreateSeuratObject(counts = compartment.df, assay = assay)

        }else{
            layer_compartment <- CreateAssayObject(counts = compartment.df)
            layer_compartment@scale.data <- as.matrix(layer_compartment@data)

        }

        return(layer_compartment)
    }else{
        insu.df <- as.data.frame(h5f$'insulation') %>% dplyr::select(-bulk) %>% convertFormat() %>% convertID(dict_cb=dict_cb) 
        tads.df <- as.data.frame(h5f$'calib_tads') %>% convertFormat() %>% convertID(dict_cb=dict_cb)

        if (! is.null(assay)){
            layer_insu <- CreateSeuratObject(counts = insu.df, assay = paste0(assay, "_score"))
            layer_tads <- CreateSeuratObject(counts = tads.df, assay = paste0(assay, "_boundary"))

        }else{
            layer_insu <- CreateAssayObject(counts = insu.df)
            layer_tads <- CreateAssayObject(counts = tads.df)
        }
        layer_tads@scale.data <- as.matrix(layer_tads@data)
        layer_insu@scale.data <- as.matrix(layer_insu@data)

        return(c('insulation'=layer_insu, 'tads'=layer_tads))

    }

}



createEmbedLayer <- function(
    seurat_obj,
    npy_file, # row cell; col data;
    layer_name,
    rownames_, #cb_higashi$cellbarcode,
    n_dim = -1,
    process = "raw",
    regress_out = FALSE
){
    
    n_last <- 3
    format_str <- substr(npy_file, nchar(npy_file) - n_last + 1, nchar(npy_file))
    if(format_str == "npy" || format_str == "npz"){
        library(reticulate)
        use_python("/share/home/wangmx/anaconda3/envs/bioinfo/bin/python")
        np <- import("numpy")
        if (format_str == "npz"){
            mtx = np$load(npy_file)['arr_0']
        }else{
            mtx = np$load(npy_file)
        }

        colnames(mtx) <- paste("Dim", 1:ncol(mtx))
        rownames(mtx ) <- rownames_
        mtx <- mtx[Cells(seurat_obj),]
        ## select significant features
        mtx <- mtx[, names(sort(apply(mtx, 2, sd), decreasing = TRUE))]
        colnames(mtx) <- paste("Dim", 1:ncol(mtx))
        mtx <- t(mtx)

    }else if(format_str == "hdf"){
        temp = H5Fopen(npy_file)
        mtx <- as.data.frame(temp$data$block0_values)
        colnames(mtx) <- temp$data$axis1
        rownames(mtx) <- temp$data$axis0

        mtx <- mtx[,Cells(seurat_obj)]
    }


    if (n_dim == "variable"){
        print("Selecting variable dims...")
        n_dim = feature_test(mtx, p_cutoff = 0.01)
    }else if(n_dim == -1){
        print("Extracting all dims...")
        n_dim = ncol(mtx)
    }
    print(paste0("Using ", n_dim, " for downstream analysis..."))

    fh_l2_norm <- CreateAssayObject(counts = mtx)
    print(paste0("Assay Created: ", layer_name))
    seurat_obj <- seurat_obj
    seurat_obj[[layer_name]] = fh_l2_norm
    DefaultAssay(seurat_obj) <- layer_name
    print("Run Clustering...")
    if(process == "pca"){
        print("Performing on pca matrix.")

        seurat_obj[[layer_name]]@scale.data <- as.matrix(seurat_obj[[layer_name]]@data)
        seurat_obj %>%
            RunPCA(features = rownames(seurat_obj), npcs = 50, verbose = FALSE) %>%
            RunUMAP(dims = 1:n_dim, verbose = F) %>%
            FindNeighbors(dims = 1:n_dim, verbose = F) %>%
            FindClusters(resolution = seq(0.2,0.6,0.2), verbose = F) -> seurat_obj

    }else if(process == "raw"){
        print("Performing on raw matrix.")

        seurat_obj %>%
            RunUMAP(dims = NULL, slot = "data", features = rownames(seurat_obj)[1:n_dim], verbose = F) %>%
            FindNeighbors(dims = NULL, features = rownames(seurat_obj)[1:n_dim], verbose = F) %>%
            FindClusters(resolution = seq(0.2,0.6,0.2), verbose = F) -> seurat_obj
    }else if(process == "gene_score"){
        print("Performing gene_score workflow.")

        sd_feature <- data.frame(apply(mtx, 1, sd))
        colnames(sd_feature) <- c("sd")

        target_sum = round(median(colSums(mtx)), 0)
        sd_feature %>%
            mutate(gene = rownames(sd_feature)) %>%
            arrange(desc(sd)) %>%
            top_n(wt = sd, n = 2000) %>%
            pull(gene) -> var.features
        
        VariableFeatures(seurat_obj) <- var.features

        if(regress_out){
            seurat_obj %>%
                NormalizeData(normalization.method = "RC", scale.factor = target_sum) %>%
                ScaleData(scale.max = 10, vars.to.regress = "contactsN") %>% 
                RunPCA(npcs = 50, verbose = FALSE) %>%
                RunUMAP(dims = 1:n_dim, verbose = F) %>%
                FindNeighbors(dims = 1:n_dim, verbose = F) %>%
                FindClusters(resolution = seq(0.2,0.6,0.2), verbose = F) -> seurat_obj
        }else{
            seurat_obj %>%
                NormalizeData(normalization.method = "RC", scale.factor = target_sum) %>%
                ScaleData(scale.max = 10) %>%
                RunPCA(npcs = 50, verbose = FALSE) %>%
                RunUMAP(dims = 1:n_dim, verbose = F) %>%
                FindNeighbors(dims = 1:n_dim, verbose = F) %>%
                FindClusters(resolution = seq(0.2,0.6,0.2), verbose = F) -> seurat_obj
        }
      

    }

    print("Done.")

    return(seurat_obj)

}
#### GIVEN two cell type information of same cells, calculate the overlap ratio
calOverLapScore <- function(label_df, rownames_, colnames_){
    colnames(label_df) <- c("target", "query")

    label_df.mtx <- foreach(target_celltype = levels(label_df$target),
                    .combine = "rbind", .errorhandling = "stop") %dopar% {
    label_df %>%
        filter(target == target_celltype) %>%
        group_by(target) %>%
        count(query) %>%
        mutate(overlap= n / sum(n)) %>%
        right_join(data.frame(query = as.factor(levels(label_df$query)))) %>%
        arrange(query) %>%
        mutate(overlap = replace_na(overlap, 0))
    }

    label_df.mtx %>%
        pull(overlap) %>%
        matrix(nrow = length(unique(label_df$target)),
            ncol = length(unique(label_df$query))) -> mtx
    rownames(mtx) <- rownames_
    colnames(mtx) <-colnames_
    return(mtx)

}

#### select differential insulation feature; only consider bondary
test_diff <- function(df, target_cell_type){
    data.frame(t(df)) %>%
        dplyr::select(contains(target_cell_type)) -> target_df
    target_value <- as.vector(unlist(target_df[1,]))

    data.frame(t(df)) %>%
        dplyr::select(!contains(target_cell_type)) -> background_df
    background_value <- as.vector(unlist(background_df[1,]))
    fitted = fitdistr(background_value, "normal")

    pvalue = pnorm(target_value, mean = fitted$estimate[1], sd = fitted$estimate[2])
    diff = target_value - median(background_value)
    return(paste0(pvalue, ":", diff, ":", target_value, ":", median(background_value)))
}


## cooltools output as input
getDiffBin <- function(insulation_table, boundary_table, cell_type){
    insulation_table %>%
        mutate(bin_id=rownames(insulation_table)) %>%
        dplyr::select(bin_id) %>%
        separate(bin_id, sep = "_", into = c("chr", "start", "end")) %>%
        mutate(bin_id=rownames(insulation_table)) -> base
    
    boundary_table %>%
        mutate(bin_id = rownames(boundary_table)) %>%
        pivot_longer(cols = colnames(boundary_table), names_to = "celltype", values_to = "is_boundary") %>%
        filter(str_detect(celltype, cell_type)) %>%
        filter(is_boundary == TRUE) %>%
        inner_join(base) %>%
        dplyr::select(-bin_id) -> boundary.df

    base$test = apply(insulation_table, 1, FUN = test_diff, target_cell_type = cell_type)
    base %>%
        separate(test, sep=":", into = c("p_val", "delta", "target_val", "median_background_val"), convert = TRUE) %>%
        mutate(target = cell_type) -> base
    return(base)

}



## call true cells based on multiome data
joint_cell_calling <- function(df,
    min_contacstN,
    min_RNAN,
    save_path=NULL
    ){
    df %>%
        mutate(cell_indicator = 
            if_else(contactsN >= min_contacstN & nCount_RNA >= min_RNAN, "Cells", "Non-cells")) -> df

    plot <- ggplot(df, aes(x = contactsN, y =nCount_RNA, colour = cell_indicator)) +
        geom_point_rast(alpha = 1) +
        scale_y_continuous(breaks = c(0, 10^(0:5)),
                    labels = c("0", "1", "10", "100", "1000", "10K", "100K"),
                    trans="log10") +
        scale_x_continuous(breaks = c(0, 10^(0:5)),
                    labels = c("0", "1", "10", "100", "1000", "10K", "100K"),
                    trans="log10") + 
        scale_colour_manual(values = c("Cells" = "orange",
                                    "Non-cells" = "grey"))+
        # geom_bin2d(binwidth = c(0.08, 0.08)) + 
        # scale_fill_distiller(palette="")
                                            # ,
                                    # labels=c(paste0('Cells(', sum(df$cell_indicator == "Cells"), ")") ,
                                    #         paste0('Non-cells(', sum(df$cell_indicator == "Non-cells"), ")"))) + 
        theme_bw()+
        coord_fixed() + fixed_coordinates +
        geom_hline(yintercept= min_RNAN, linetype = 'dashed', color = 'firebrick3', size = 1.2) +
        geom_vline(xintercept= min_contacstN, linetype = 'dashed', color = 'firebrick3', size = 1.2) +
        theme(panel.grid = element_line(colour = 'white'), 
            legend.title=element_blank(),
            # legend.position = 'middle',
            legend.text = element_text(size=12),
        axis.title.x =element_text(size=18), 
            axis.title.y =element_text(size=18),
            axis.text.x = element_text(size=10,colour = "black"),
            axis.text.y = element_text(size=10,colour = "black"))+
        labs(x="Number of Unique Contacts per barcode", y="Number of UMIs per barcode")

    if( ! is.null(save_path )){
       ggsave(filename = save_path,
        plot,
        width=180,
        height = 150, units = "mm",dpi = 600) 
    }
    print(paste0('Cells: ', sum(df$cell_indicator == "Cells")))
    print(paste0('Non-cells: ', sum(df$cell_indicator == "Non-cells")))

    return(plot)

}


calRNAChange <- function(
    seurat_obj,
    ident.1,
    ident.2,
    cluster_name,
    logfc.threshold = 0.25,
    min.pct = 0
){ 
    cells.1 <- WhichCells(object = seurat_obj, idents = ident.1)
    cells.2 <- WhichCells(object = seurat_obj, idents = ident.2)


    seurat_obj@meta.data %>%
        mutate(cellbarcode = rownames(seurat_obj@meta.data)) %>%
        mutate(new_group = case_when(
            cellbarcode %in% cells.1 ~ "target",
            cellbarcode %in% cells.2 ~ "background",
            TRUE ~ "Not_associated"
        )) -> seurat_obj@meta.data
    seurat_obj@meta.data$new_group <- factor(seurat_obj@meta.data$new_group, levels = unique(seurat_obj@meta.data$new_group))
    AverageExpression(
        seurat_obj,
        assays = "RNA",
        group.by = "new_group",
        layer = "data"
        )[[1]] %>% as.data.frame() %>%
    dplyr::select(target, background) -> RNA.result
    RNA.result$gene <- rownames(RNA.result)
    return(RNA.result)

}

## saddle_plot
saddle_plot <- function(mtx, title, col_fun){
    ht2 <- Heatmap(mtx,cluster_columns = FALSE,
        cluster_rows=FALSE,
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill){

            if ( i == 50 & j == 50){
            grid.text(AA, x = x, y = y)
            }else if(i == 3 & j == 3){
            grid.text(BB, x = x, y = y)
            }else if (i == 3 & j == 50){
            grid.text(BA, x = x, y = y)
            }else if (i == 50 & j == 3){
            grid.text(AB, x = x, y = y)
            }
        },
        column_title = title)
}