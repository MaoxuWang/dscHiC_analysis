library(tidyverse)
library(kneedle)
library(scales)


args = commandArgs(T)

input = args[1]
wd = args[2]


selectBestSensitivity <- function(qc.matrix.sorted){
    ratios = c()
    for (i in 1:10){
        out <- kneedle(qc.matrix.sorted$rank, 
            qc.matrix.sorted$contactsN,
            sensitivity = i)
        ratio <- sum(qc.matrix.sorted[qc.matrix.sorted$rank <= out[1],]$contactsN) / sum(qc.matrix.sorted$contactsN)
        ratios <- c(ratios, ratio)
    }
    idx = which(abs(ratios - 0.85) == min(abs(ratios-0.85)))
    if(length(idx) != 1){
        idx = idx[1]
    }
    return(idx)
}


sampleid = gsub(".contacts_count.txt", "", gsub(".*\\/", "", input))
sample = paste0(wd, "/", sampleid,".")

qc.matrix <- read_delim(input, delim = "\t", col_names = F)
colnames(qc.matrix) <- c("cellbarcode", "contactsN")
qc.matrix %>%
    arrange(desc(contactsN)) -> qc.matrix.sorted

qc.matrix.sorted$rank <- 1:nrow(qc.matrix.sorted)

idx = selectBestSensitivity(qc.matrix.sorted)

out <- kneedle(qc.matrix.sorted$rank, 
                qc.matrix.sorted$contactsN,
                sensitivity = idx)

qc.matrix.sorted %>%
     mutate(source=case_when(
    rank <= out[1]  ~ "Cells",
    rank > out[1] ~ "Background")) -> qc.matrix.filtered


cell_calling <- ggplot(qc.matrix.filtered, aes(x=rank, y=contactsN, colour = source)) + 
    geom_line(linewidth = 1.5) + 
    scale_colour_manual(values =c("Cells"="#007b43", "Background"="gray")) + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x))) +
    theme_bw()+
    theme(panel.grid = element_line(colour = 'white'),
          legend.title=element_blank(),
          legend.position = c(0.85, 0.85),
          axis.title.x =element_text(size=16), 
        axis.title.y =element_text(size=16),
        axis.text.x = element_text(size=10,colour = "black"),
        axis.text.y = element_text(size=10,colour = "black"))+
    labs(x="Barcodes", y="Contact Counts", subtitle = paste0("Estimated Cells: ", sum(qc.matrix.filtered$source == "Cells")))


ggsave(cell_calling, filename = paste0(sample, "CellCalling.png"), width = 140, height = 120,units = "mm", dpi = 600)

## output cellbarcodes.txt

qc.matrix.filtered %>%
    filter(source == "Cells") %>%
    select(cellbarcode) -> cellbarcode

write.table(cellbarcode, file=paste0(sample, "CellCalled.txt"), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
