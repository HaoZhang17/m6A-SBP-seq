# 1. Import library
```{r}
library(BiocParallel)
library(Biostrings)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggsci)
library(GenomicAlignments)
library(GenomicRanges)
library(Guitar)
library(openxlsx)
library(paletteer)
library(reshape2)
library(rtracklayer)
library(Rsamtools)
library(tidyverse)
```


# 2. Import peak files
```{r}
dir_path <- '/results/peaks/'

# for bed files
c293_peak_bed_list <- list()
for (i in c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG')){
  filepath_forward <- paste0(dir_path, 'HEK293T_', i, '_forward_strand/peaks.bed')
  filepath_reverse <- paste0(dir_path, 'HEK293T_', i, '_reverse_strand/peaks.bed')
  bed_file_forward <- read.table(filepath_forward, sep = '\t')
  bed_file_reverse <- read.table(filepath_reverse, sep = '\t')
  bed_file_forward_filtered <- bed_file_forward %>% dplyr::filter(V6 == '+')
  bed_file_reverse_filtered <- bed_file_reverse %>% dplyr::filter(V6 == '-')
  c293_peak_bed_list[[i]] <- rbind(bed_file_forward_filtered, bed_file_reverse_filtered)
}

# for csv files
c293_peak_csv_list <- list()
for (i in c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG')){
  filepath_forward <- paste0(dir_path, 'HEK293T_', i, '_forward_strand/peaks.csv')
  filepath_reverse <- paste0(dir_path, 'HEK293T_', i, '_reverse_strand/peaks.csv')
  csv_file_forward <- read.csv(filepath_forward, header = T)
  csv_file_reverse <- read.csv(filepath_reverse, header = T)
  csv_file_forward_filtered <- csv_file_forward %>% dplyr::filter(strand == '+')
  csv_file_reverse_filtered <- csv_file_reverse %>% dplyr::filter(strand == '-')
  c293_peak_csv_list[[i]] <- rbind(csv_file_forward_filtered, csv_file_reverse_filtered)
}
```

# 3. Filter peaks
```{r}
# filter peaks with conditions: IP/INPUT >= 2, RPM.IP >= 1
peak_filter <- function(peak_csv, peak_bed) {
  peak_csv$Fold <- peak_csv$RPM.IP / peak_csv$RPM.input
  peak_csv$peakID <- paste(peak_csv$chr, peak_csv$chromStart, peak_csv$chromEnd, peak_csv$strand, sep = '_')
  peak_bed$peakID <- paste(peak_bed$V1, peak_bed$V2, peak_bed$V3, peak_bed$V6, sep = '_')
  peak_csv_filtered <- peak_csv %>% dplyr::filter(Fold >= 2 & RPM.IP >= 1) 
  peak_bed_filtered <- peak_bed[which(peak_bed$peakID %in% peak_csv_filtered$peakID), ]
  peak_bed_filtered_join <- merge(peak_bed_filtered, peak_csv_filtered, by = 'peakID', all.x = T)
  peak_bed_filtered_join <- peak_bed_filtered_join[, c(2, 3, 4, 25, 6, 7, 8, 9, 10, 11, 12, 13)]
  return(list(peak_csv_filtered, peak_bed_filtered_join))
}


Em6ABP1_peak_list <- peak_filter(c293_peak_csv_list[['Em6ABP1']], c293_peak_bed_list[['Em6ABP1']])
Em6ABP1_peak_csv_filtered <- Em6ABP1_peak_list[[1]]
Em6ABP1_peak_bed_filtered <- Em6ABP1_peak_list[[2]]

SYSY_peak_list <- peak_filter(c293_peak_csv_list[['SYSY']], c293_peak_bed_list[['SYSY']])
SYSY_peak_csv_filtered <- SYSY_peak_list[[1]]
SYSY_peak_bed_filtered <- SYSY_peak_list[[2]]

Millipore_peak_list <- peak_filter(c293_peak_csv_list[['Millipore']], c293_peak_bed_list[['Millipore']])
Millipore_peak_csv_filtered <- Millipore_peak_list[[1]]
Millipore_peak_bed_filtered <- Millipore_peak_list[[2]]

Active_Motif_peak_list <- peak_filter(c293_peak_csv_list[['Active_Motif']], c293_peak_bed_list[['Active_Motif']])
Active_Motif_peak_csv_filtered <- Active_Motif_peak_list[[1]]
Active_Motif_peak_bed_filtered <- Active_Motif_peak_list[[2]]

IgG_peak_list <- peak_filter(c293_peak_csv_list[['IgG']], c293_peak_bed_list[['IgG']])
IgG_peak_csv_filtered <- IgG_peak_list[[1]]
IgG_peak_bed_filtered <- IgG_peak_list[[2]]


write.table(Em6ABP1_peak_bed_filtered, 'results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(SYSY_peak_bed_filtered, 'results/peaks/HEK293T_SYSYSY/peaks_filtered.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(Millipore_peak_bed_filtered, 'results/peaks/HEK293T_Millipore/peaks_filtered.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(Active_Motif_peak_bed_filtered, 'results/peaks/HEK293T_Active_Motif/peaks_filtered.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(IgG_peak_bed_filtered, 'results/peaks/HEK293T_IgG/peaks_filtered.bed', sep = '\t', col.names = F, row.names = F, quote = F)
```


# 4. Obtain sample-specific peaks
## 4.1 bedtools intersect
## >>> bash

## 4.2 Load and filter peaks
```{r}
dir_path <- 'results/peaks/intersected/'
samples <- c("Em6ABP1", "SYSY", "Active_Motif", "Millipore", "IgG")
lapply(samples, function(sample) {
  input_file <- file.path(dir_path, paste0(sample, "_peaks_intersected.bed"))
  peak_data <- read.table(input_file, sep = "\t", header = FALSE)
  peak_specific <- peak_data %>% filter(V13 == '.')
  peak_specific_bed12 <- peak_specific[, 1:12]
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  write.table(peak_specific, 
              file.path(dir_path, paste0(sample, "_peak_specific.bed")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(peak_specific_bed12, 
              file.path(dir_path, paste0(sample, "_peak_specific.bed12")),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
})
```


# 5. Metagene analysis [Fig. 4b, Fig. 4g, Supp. Fig. 8c]
## Definition of the plotting function
```{r}
create_metagene_plot <- function(data, title, y_max, text_y) {
  color_value <- c("#529e3fff", "#C53A32", "#3B76AF", "#EF8636")
  
  ggplot(data, aes(x = x, y = density, fill = group, color = group)) +
    geom_line(linewidth = 1.2, alpha = 0.8) +
    geom_ribbon(aes(ymin = 0, ymax = density), 
                alpha = 0.1, size = 0, show.legend = FALSE) + 
    scale_color_manual(name = '', values = color_value) +
    scale_fill_manual(name = '', values = color_value) +
    scale_x_continuous(
      name = NULL,
      breaks = c(0, 0.142, 0.588, 1),
      labels = c("", "Start", "Stop", ""),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      name = "Frequency (%)",
      breaks = seq(0, y_max, by = 0.5),
      limits = c(0, y_max),
      expand = c(0, 0)
    ) +
    geom_vline(xintercept = c(0.142, 0.588), 
               linetype = "dashed", color = "gray50", 
               alpha = 0.7, size = 0.5) +
    annotate("text", x = 0.07, y = text_y, 
             label = "5'UTR", size = 4, fontface = "bold") +
    annotate("text", x = 0.35, y = text_y, 
             label = "CDS", size = 4, fontface = "bold") +
    annotate("text", x = 0.70, y = text_y, 
             label = "3'UTR", size = 4, fontface = "bold") +
    theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold",
                                margin = margin(b = 15)),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 10),
      axis.title.y = element_text(margin = margin(r = 10)),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    ggtitle(title)
}
```{r}


## 5.1 all peaks [Fig. 4b, Supp. Fig. 8c]
```{r}
options(stringsAsFactors = F)
txdb <- txdbmaker::makeTxDbFromGFF(file = "resources/reference/Homo_sapiens.GRCh38.95.gtf", 
                        format="gtf", 
                        dataSource="Ensembl", 
                        organism="Homo sapiens")

# Em6ABP1
metagene_all_Em6ABP1 <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = "results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed", 
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("Em6ABP1"))

# Antibodies
stBedFiles_all_antibody <- list("results/peaks/HEK293T_SYSYSY/peaks_filtered.bed",
                               "results/peaks/HEK293T_Active_Motif/peaks_filtered.bed",
                               "results/peaks/HEK293T_Millipore/peaks_filtered.bed")

metagene_all_antibody <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = stBedFiles_antibody, 
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("SY", "AM", "Milli"))


metagene_all_Em6ABP1_data <- metagene_all_Em6ABP1$data
metagene_all_antibody_data <- metagene_all_antibody$data

metagene_all_Em6ABP1_plot <- create_metagene_plot(
  data = metagene_all_Em6ABP1_data,
  title = "All peaks_Em6ABP1",
  y_max = 3.5,
  text_y = 3.3
)

metagene_all_antibody_plot <- create_metagene_plot(
  data = metagene_all_antibody_data,
  title = "All peaks_Antibodies",
  y_max = 4,
  text_y = 3.5
)

plot_grid(plotlist = list(metagene_all_Em6ABP1_plot, metagene_all_antibody_plot), ncol = 2, align = 'h')
```


## 5.2 specific peaks [Fig. 4g]
```{r}
dir_path <- 'results/peaks/intersected/'
options(stringsAsFactors = F)

txdb <- txdbmaker::makeTxDbFromGFF(file = "resources/reference/Homo_sapiens.GRCh38.95.gtf", 
                        format="gtf", 
                        dataSource="Ensembl", 
                        organism="Homo sapiens")

metagene_sp_Em6ABP1 <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = paste(dir_path, "Em6ABP1_peak_specific.bed", sep = ''),
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("Em6ABP1"))
metagene_sp_SY <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = paste(dir_path, "SYSY_peak_specific.bed", sep = ''),
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("SY"))
metagene_sp_AM <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = paste(dir_path, "Active_Motif_peak_specific.bed", sep = ''),
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("AM"))
metagene_sp_Milli <- GuitarPlot(txTxdb = txdb, 
                        stBedFiles = paste(dir_path, "Millipore_peak_specific.bed", sep = ''),
                        headOrtail = FALSE,
                        enableCI = FALSE, 
                        mapFilterTranscript = TRUE, 
                        pltTxType = c("mrna"), 
                        stGroupName = c("Milli"))

metagene_sp_Em6ABP1_data <- metagene_sp_Em6ABP1$data
metagene_sp_SY_data <- metagene_sp_SY$data
metagene_sp_AM_data <- metagene_sp_AM$data
metagene_sp_Milli_data <- metagene_sp_Milli$data

metagene_sp_Em6ABP1_plot <- create_metagene_plot(
  data = metagene_sp_Em6ABP1_data,
  title = "Specific peaks_Em6ABP1",
  y_max = 4,
  text_y = 3.5
)

metagene_sp_SY_plot <- create_metagene_plot(
  data = metagene_sp_SY_data,
  title = "Specific peaks_SYSY",
  y_max = 4,
  text_y = 3.5
)

metagene_sp_AM_plot <- create_metagene_plot(
  data = metagene_sp_AM_data,
  title = "Specific peaks_Active_Motif",
  y_max = 4,
  text_y = 3.5
)

metagene_sp_Milli_plot <- create_metagene_plot(
  data = metagene_sp_Milli_data,
  title = "Specific peaks_Millipore",
  y_max = 4,
  text_y = 3.5
)

plot_grid(plotlist = list(metagene_sp_Em6ABP1_plot, metagene_sp_SY_plot, metagene_sp_AM_plot, metagene_sp_Milli_plot), ncol = 2, align = 'h')
```


# 6. Distribution of peaks across transcript regions [Fig. 4c, Supp. Fig. 8e]
```{r}
process_anno <- function(file_path) {
  stats <- fread(file_path, sep = "\t", header = TRUE)
  
  stats_1 <- unique(stats[seq(1, min(13, nrow(stats))), .(Annotation, Number.of.peaks)])
  stats_1 <- stats_1[Number.of.peaks != '0.0']
  stats_1$Number.of.peaks <- as.numeric(stats_1$Number.of.peaks)

  stats_1$Annotation[stats_1$Annotation == 'Promoter'] <- 'TSS'
  stats_1$Annotation[stats_1$Annotation == 'Exon'] <- 'CDS'
  
  ncRNA_categories <- c('pseudo', 'miRNA', 'ncRNA', 'rRNA')
  existing_ncRNA <- ncRNA_categories[ncRNA_categories %in% stats_1$Annotation]
  
  if (length(existing_ncRNA) > 0) {
    ncRNA_sum <- sum(stats_1$Number.of.peaks[stats_1$Annotation %in% existing_ncRNA])
    stats_1 <- stats_1[!Annotation %in% existing_ncRNA]
        stats_1 <- rbind(stats_1, data.table(Annotation = 'ncRNA', Number.of.peaks = ncRNA_sum))
  }
  
  stats_1$Number.of.peaks <- as.numeric(stats_1$Number.of.peaks)
  stats_1$per <- round(stats_1$Number.of.peaks / sum(stats_1$Number.of.peaks) * 100, digits = 2)
  stats_1$Anno_new <- paste0(stats_1$Annotation, ' (', stats_1$per, '%)')
  
  return(stats_1)
}

create_pie <- function(data, title) {
  new_order <- c('TSS', '5UTR', 'CDS', '3UTR', 'Intron', 'TTS', 'Intergenic', 'ncRNA')
  data$Annotation <- factor(data$Annotation, levels = new_order)
  
  ggplot(data, aes(x = "", y = Number.of.peaks, fill = Annotation)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(aes(label = Anno_new), position = position_stack(vjust = 0.5), size = 3) +
    scale_fill_paletteer_d("ggsci::nrc_npg") +
    ggtitle(title)
}



dir_path_all <- 'results/peaks/'
dir_path_sp <- 'results/peaks/intersected/'

samples <- c("Em6ABP1", "SYSY", "Active_Motif", "Millipore", "IgG")

all_plots <- lapply(samples, function(s) {
  file <- paste0(dir_path_all, "HEK293T_", s, "/peak_filtered_anno.stats.txt")
  data <- process_anno(file)
  create_pie(data, paste(s, 'all peaks'))
})

sp_plots <- lapply(samples, function(s) {
  file <- paste0(dir_path_sp, s, "_specific_peak_anno.stats.txt")
  data <- process_anno(file)
  create_pie(data, paste(s, 'sp peaks'))
})

plot_grid(plotlist = c(all_plots, sp_plots), ncol = 2, align = 'v')
```


# 7. Peak overlap analysis [Fig. 4d]
```{r}
dir_path <- 'results/peaks/intersected/'
Em6ABP1_peak <- read.table(paste0(dir_path, 'Em6ABP1_peaks_intersected.bed'), sep = '\t', header = F)
SYSY_peak <- read.table(paste0(dir_path, 'SYSY_peaks_intersected.bed'), sep = '\t', header = F)
Active_Motif_peak <- read.table(paste0(dir_path, 'Active_Motif_peaks_intersected.bed'), sep = '\t', header = F)
Millipore_peak <- read.table(paste0(dir_path, 'Millipore_peaks_intersected.bed'), sep = '\t', header = F)
IgG_peak <- read.table(paste0(dir_path, 'IgG_peaks_intersected.bed'), sep = '\t', header = F)


Em6ABP1_peak$Em6ABP1_peakID <- paste(Em6ABP1_peak$V1, Em6ABP1_peak$V2, Em6ABP1_peak$V3, Em6ABP1_peak$V6, sep = '_')
SYSY_peak$SYSY_peakID <- paste(SYSY_peak$V1, SYSY_peak$V2, SYSY_peak$V3, SYSY_peak$V6, sep = '_')
Active_Motif_peak$Active_Motif_peakID <- paste(Active_Motif_peak$V1, Active_Motif_peak$V2, Active_Motif_peak$V3, Active_Motif_peak$V6, sep = '_')
Millipore_peak$Millipore_peakID <- paste(Millipore_peak$V1, Millipore_peak$V2, Millipore_peak$V3, Millipore_peak$V6, sep = '_')
IgG_peak$IgG_peakID <- paste(IgG_peak$V1, IgG_peak$V2, IgG_peak$V3, IgG_peak$V6, sep = '_')


# Em6ABP1 vs SY
# from Em6ABP1's side
Em6ABP1_total_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID)) # 22203
Em6ABP1_SYSY_overlap_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID[Em6ABP1_peak$V13 == 'SYSY'])) # 14910
Em6ABP1_sp_peak_num <- Em6ABP1_total_peak_num - Em6ABP1_SYSY_overlap_peak_num # 7293
# from SY's side
SYSY_total_peak_num <- length(unique(SYSY_peak$SYSY_peakID)) # 18452
SB_Em6ABP1_overlap_peak_num <- length(unique(SYSY_peak$SYSY_peakID[SYSY_peak$V13 == 'Em6ABP1'])) # 15538
SYSY_sp_peak_num <- SYSY_total_peak_num - SB_Em6ABP1_overlap_peak_num # 2914


# Em6ABP1 vs AM
# from Em6ABP1's side
Em6ABP1_total_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID)) # 22203
Em6ABP1_Active_Motif_overlap_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID[Em6ABP1_peak$V13 == 'Active_Motif'])) # 3516
Em6ABP1_sp_peak_num <- Em6ABP1_total_peak_num - Em6ABP1_Active_Motif_overlap_peak_num # 18687
# from AM's side
Active_Motif_total_peak_num <- length(unique(Active_Motif_peak$Active_Motif_peakID)) # 4428
SB_Em6ABP1_overlap_peak_num <- length(unique(Active_Motif_peak$Active_Motif_peakID[Active_Motif_peak$V13 == 'Em6ABP1'])) # 4187
Active_Motif_sp_peak_num <- Active_Motif_total_peak_num - SB_Em6ABP1_overlap_peak_num # 241


# Em6ABP1 vs Milli
# from Em6ABP1's side
Em6ABP1_total_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID)) # 22203
Em6ABP1_Millipore_overlap_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID[Em6ABP1_peak$V13 == 'Millipore'])) # 3817
Em6ABP1_sp_peak_num <- Em6ABP1_total_peak_num - Em6ABP1_Millipore_overlap_peak_num # 18386
# from Milli's side
Millipore_total_peak_num <- length(unique(Millipore_peak$Millipore_peakID)) # 5055
SB_Em6ABP1_overlap_peak_num <- length(unique(Millipore_peak$Millipore_peakID[Millipore_peak$V13 == 'Em6ABP1'])) # 4611
Millipore_sp_peak_num <- Millipore_total_peak_num - SB_Em6ABP1_overlap_peak_num # 444


# Em6ABP1 vs IgG
# from Em6ABP1's side
Em6ABP1_total_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID)) # 22203
Em6ABP1_IgG_overlap_peak_num <- length(unique(Em6ABP1_peak$Em6ABP1_peakID[Em6ABP1_peak$V13 == 'IgG'])) # 208
Em6ABP1_sp_peak_num <- Em6ABP1_total_peak_num - Em6ABP1_IgG_overlap_peak_num # 21995
# from IgG's side
IgG_total_peak_num <- length(unique(IgG_peak$IgG_peakID)) # 655
SB_Em6ABP1_overlap_peak_num <- length(unique(IgG_peak$IgG_peakID[IgG_peak$V13 == 'Em6ABP1'])) # 225
IgG_sp_peak_num <- IgG_total_peak_num - SB_Em6ABP1_overlap_peak_num # 430
```


# 8. DRACH motif in peaks [Fig. 4e]
```{r}
motif_scan_in_peak <- function(peak, genome_fasta){
  # extract peak sequence
  peak_seqs <- scanFa(genome_fasta, param = peak)
  names(peak_seqs) <- peak$name
  
  # search DRACH motif
  drach_pattern <- DNAString("DRACH")
  
  # calculate the number of peaks containing DRACH
  contains_drach <- vapply(peak_seqs, function(seq) {
    matches <- matchPattern(drach_pattern, seq, fixed = FALSE)
    length(matches) > 0
  }, logical(1))
  
  
  results <- data.frame(
    peak_id = names(peak_seqs),
    chr = as.character(seqnames(peak)),
    start = start(peak),
    end = end(peak),
    width = width(peak),
    has_rrach = contains_drach,
    drach_count = drach_counts,
    sequence = as.character(peak_seqs),
    stringsAsFactors = FALSE
  )
  
  # output
  cat("Total peak number:", length(peak_seqs), "\n")
  cat("Number of peaks with DRACH:", sum(contains_drach), "\n")
  cat("Number of peaks without DRACH:", sum(!contains_drach), "\n")
  cat("Porportion of peaks containing DRACH:", 
      round(mean(contains_drach) * 100, 2), "%\n")
  
  return(mean(contains_drach))
}

fasta_file <- "/gold1/zhanghao/Refseq/Homo_sapiens/Sequence_Ensembl/GRCh38.fa" 
fa <- FaFile(fasta_file)



Em6ABP1_all_peaks <- import("results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed")
SYSY_all_peaks <- import("rresults/peaks/HEK293T_SYSY/peaks_filtered.bed")
Active_Motif_all_peaks <- import("results/peaks/HEK293T_Active_Motif/peaks_filtered.bed")
Millipore_all_peaks <- import("results/peaks/HEK293T_Millipore/peaks_filtered.bed")
IgG_all_peaks <- import("results/peaks/HEK293T_IgG/peaks_filtered.bed")

Em6ABP1_sp_peaks <- import("results/peaks/intersected/Em6ABP1_peak_specific.bed")
SYSY_sp_peaks <- import("results/peaks/intersected/SYSY_peak_specific.bed")
Active_Motif_sp_peaks <- import("results/peaks/intersected/Active_Motif_peak_specific.bed")
Millipore_sp_peaks <- import("results/peaks/intersected/Millipore_peak_specific.bed")
IgG_sp_peaks <- import("results/peaks/intersected/IgG_peak_specific.bed")



Em6ABP1_all_peak_DRACH_ratio <- motif_scan_in_peak(Em6ABP1_all_peaks, fa)
SYSY_all_peak_DRACH_ratio <- motif_scan_in_peak(SYSY_all_peaks, fa)
Active_Motif_all_peak_DRACH_ratio <- motif_scan_in_peak(Active_Motif_all_peaks, fa)
Millipore_all_peak_DRACH_ratio <- motif_scan_in_peak(Millipore_all_peaks, fa)
IgG_all_peak_DRACH_ratio <- motif_scan_in_peak(IgG_all_peaks, fa)


Em6ABP1_sp_peak_DRACH_ratio <- motif_scan_in_peak(Em6ABP1_sp_peaks, fa)
SYSY_sp_peak_DRACH_ratio <- motif_scan_in_peak(SYSY_sp_peaks, fa)
Active_Motif_sp_peak_DRACH_ratio <- motif_scan_in_peak(Active_Motif_sp_peaks, fa)
Millipore_sp_peak_DRACH_ratio <- motif_scan_in_peak(Millipore_sp_peaks, fa)
IgG_sp_peak_DRACH_ratio <- motif_scan_in_peak(IgG_sp_peaks, fa)


# Ratio of peaks with DRACH
all_peak_DRACH_ratio <- data.frame(sample = rep(c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore'), each = 2),
                                   type = rep(c('With DRACH', 'Without DRACH'), time = 4), 
                                   ratio = c(Em6ABP1_all_peak_DRACH_ratio, 1 - Em6ABP1_all_peak_DRACH_ratio,
                                             SYSY_all_peak_DRACH_ratio, 1 - SYSY_all_peak_DRACH_ratio,
                                             Active_Motif_all_peak_DRACH_ratio, 1 - Active_Motif_all_peak_DRACH_ratio,
                                             Millipore_all_peak_DRACH_ratio, 1 - Millipore_all_peak_DRACH_ratio))

sp_peak_DRACH_ratio <- data.frame(sample = rep(c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore'), each = 2),
                                   type = rep(c('With DRACH', 'Without DRACH'), time = 4), 
                                   ratio = c(Em6ABP1_sp_peak_DRACH_ratio, 1 - Em6ABP1_sp_peak_DRACH_ratio,
                                             SYSY_sp_peak_DRACH_ratio, 1 - SYSY_sp_peak_DRACH_ratio,
                                             Active_Motif_sp_peak_DRACH_ratio, 1 - Active_Motif_sp_peak_DRACH_ratio,
                                             Millipore_sp_peak_DRACH_ratio, 1 - Millipore_sp_peak_DRACH_ratio))


all_peak_DRACH_ratio$sample <- factor(all_peak_DRACH_ratio$sample, levels = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore'))
sp_peak_DRACH_ratio$sample <- factor(sp_peak_DRACH_ratio$sample, levels = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore'))
all_peak_DRACH_ratio$ratio <- round(all_peak_DRACH_ratio$ratio, 2)
sp_peak_DRACH_ratio$ratio <- round(sp_peak_DRACH_ratio$ratio, 2)


p1_motif_peaks <- ggplot(data = all_peak_DRACH_ratio, aes(x = sample,y = ratio, fill = type))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,   
           color = "black", 
           linewidth = 0.5) +
  geom_text(aes(label = ratio), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Red_Blue_Brown") +
  ggtitle("All peaks")


p2_motif_peaks <- ggplot(data = sp_peak_DRACH_ratio, aes(x = sample,y = ratio, fill = type))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,   
           color = "black", 
           linewidth = 0.5) +
  geom_text(aes(label = ratio), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Red_Blue_Brown") +
  ggtitle("Sample-specific peaks")

pp <- list(p1, p2)
plot_grid(plotlist = pp, ncol = 2, align = 'h')


# Number of peaks with DRACH
sp_peak_DRACH_num <- data.frame(Sample = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG'),
                       Sp_total = c(nrow(Em6ABP1_sp_peaks),
                                    nrow(SYSY_sp_peaks),
                                    nrow(Active_Motif_sp_peaks),
                                    nrow(Millipore_sp_peaks),
                                    nrow(IgG_peak_specific)),
                       Sp_DRACH = c(nrow(Em6ABP1_sp_peaks) * Em6ABP1_sp_peak_DRACH_ratio,
                                    nrow(SYSY_sp_peaks) * SYSY_sp_peak_DRACH_ratio,
                                    nrow(Active_Motif_sp_peaks) * Active_Motif_sp_peak_DRACH_ratio,
                                    nrow(Millipore_sp_peaks) * Millipore_sp_peak_DRACH_ratio,
                                    nrow(IgG_peak_specific) * IgG_sp_peak_DRACH_ratio),
                       Sp_no_DRACH = c(nrow(Em6ABP1_sp_peaks) * (1 - Em6ABP1_sp_peak_DRACH_ratio),
                                    nrow(SYSY_sp_peaks) * (1 - SYSY_sp_peak_DRACH_ratio),
                                    nrow(Active_Motif_sp_peaks) * (1 - Active_Motif_sp_peak_DRACH_ratio),
                                    nrow(Millipore_sp_peaks) * (1 - Millipore_sp_peak_DRACH_ratio),
                                    nrow(IgG_peak_specific) * (1 - IgG_sp_peak_DRACH_ratio)))
sp_peak_DRACH_num_melt <- melt(sp_peak_DRACH_num, id = c('Sample', 'Sp_total'))


sp_peak_DRACH_num_melt$Sample <- factor(sp_peak_DRACH_num_melt$Sample, levels = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG'))
sp_peak_DRACH_num_melt$variable <- factor(sp_peak_DRACH_num_melt$variable, levels = c('Sp_no_DRACH', 'Sp_DRACH'))


ggplot(data = sp_peak_DRACH_num_melt, aes(x = Sample, y = value, fill = variable)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_fill_manual(name = '', values = c('#F4A582', '#B30000')) + 
  labs(x = '', y = 'Peak number')
```


# 9. Peak number [Supp. Fig. 8a]
```{r}
peak_num <- data.frame(Sample = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG'),
                       Total = c(nrow(Em6ABP1_peak_bed_filtered),
                                 nrow(SYSY_peak_bed_filtered),
                                 nrow(Active_Motif_peak_bed_filtered),
                                 nrow(Millipore_peak_bed_filtered),
                                 nrow(IgG_peak_bed_filtered)),
                       Specific = c(nrow(Em6ABP1_sp_peaks),
                                    nrow(SYSY_sp_peaks),
                                    nrow(Active_Motif_sp_peaks),
                                    nrow(Millipore_sp_peaks),
                                    nrow(IgG_peak_specific)))

peak_num$Sample <- factor(peak_num$Sample, levels = c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore', 'IgG'))

color_value = c("#C53A32", "#EF8636", "#529E3F", "#3B76AF", "#979998")
total_peak_num_plot <- ggplot(data = peak_num, aes(x = Sample, y = Total, fill = Sample)) + 
  geom_bar(stat = "identity", size = 1.5, width = 0.7) +
  geom_text(aes(label = Total), size = 4.5,vjust = -0.5, position = position_dodge(0.3)) +
  theme_classic(base_size = 15) + 
  theme(legend.position = "none") +
  scale_fill_manual(name = '', values = color_value) + 
  labs(x = '', y = 'Total peak number')

sp_peak_num_plot <- ggplot(data = peak_num, aes(x = Sample, y = Specific, fill = Sample)) + 
  geom_bar(stat = "identity", size = 1.5, width = 0.7) +
  geom_text(aes(label = Specific), size = 4.5, vjust = -0.5, position = position_dodge(0.3)) +
  theme_classic(base_size = 15) +
  theme(legend.position = "none") +
  scale_fill_manual(name = '', values = color_value) + 
  labs(x = '', y = 'Sample-specific peak number')

plot_grid(plotlist = list(total_peak_num_plot, sp_peak_num_plot), ncol = 1, align = 'h') # 4 * 5
```


# 10. DRACH motif location in read [Supp. Fig. 8d]
```{r}
Em6ABP1_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_Em6ABP1.sorted.cp_ds.bam'
SYSY_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_SYSYSY.sorted.cp_ds.bam'
Active_Motif_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_Active_Motif.sorted.cp_ds.bam'
Millipore_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_Millipore.sorted.cp_ds.bam'
input_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_Input.sorted.cp_ds.bam'
IgG_bam_file <- 'results/mapping/mapped/rep_merged/HEK293T_IgG.sorted.cp_ds.bam'


motif_read <- function(bam_file, peaks){
  param <- ScanBamParam(what = c("qname", "seq", "qwidth", "pos", "rname", "strand"))
  aln <- scanBam(bam_file, param = param)[[1]]
  reads_gr <- GRanges(
    seqnames = aln$rname,
    ranges   = IRanges(start = aln$pos, width = aln$qwidth),
    strand   = aln$strand
  )
  
  # overlap with provided peaks
  ov <- findOverlaps(reads_gr, peaks, ignore.strand = FALSE)
  reads_idx <- unique(queryHits(ov))
  
  
  # prepare sequence for scanning
  seqs <- DNAStringSet(aln$seq[reads_idx])
  lens <- aln$qwidth[reads_idx]
  strs <- aln$strand[reads_idx]
  
  # extract reversely complementary sequence for reads mapped to the minus strand
  neg <- strs == "-"
  seqs[neg] <- reverseComplement(seqs[neg])
  
  # scan for DRACH
  hits <- vmatchPattern("DRACH", seqs, fixed = FALSE)
  
  # obtain reads containing DRACH
  nr <- elementNROWS(hits)
  has_rr <- nr > 0
  ratio <- sum(nr > 0) / length(nr > 0)
  
  return(ratio)
}


Em6ABP1_sp_motif_contained_ratio <- motif_read(Em6ABP1_bam_file, Em6ABP1_sp_peaks)
SYSY_sp_motif_contained_ratio <- motif_read(SYSY_bam_file, SYSY_sp_peaks)
Active_Motif_sp_motif_contained_ratio <- motif_read(Active_Motif_bam_file, Active_Motif_sp_peaks)
Millipore_sp_motif_contained_ratio <- motif_read(Millipore_bam_file, Millipore_sp_peaks)

Em6ABP1_sp_na_ratio <- 1 - Em6ABP1_sp_motif_contained_ratio
SYSY_sp_na_ratio <- 1 - SYSY_sp_motif_contained_ratio
Active_Motif_sp_na_ratio <- 1 - Active_Motif_sp_motif_contained_ratio
Millipore_sp_na_ratio <- 1 - Millipore_sp_motif_contained_ratio


sp_motif_mid_ratio <- data.frame(sample = rep(c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore'), each = 2),
                                   type = rep(c('With DRACH', 'Without DRACH'), times = 4),
                                   value = c(Em6ABP1_sp_motif_contained_ratio, 
                                             Em6ABP1_sp_na_ratio, 
                                             SYSY_sp_motif_contained_ratio, 
                                             SYSY_sp_na_ratio, 
                                             Active_Motif_sp_motif_contained_ratio, 
                                             Active_Motif_sp_na_ratio, 
                                             Millipore_sp_motif_contained_ratio, 
                                             Millipore_sp_na_ratio))

sp_motif_mid_ratio$value <- round(sp_motif_mid_ratio$value, 2)
sample_order <- c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore')
sp_motif_mid_ratio$sample <- factor(sp_motif_mid_ratio$sample, levels = sample_order)


ggplot(data = sp_motif_mid_ratio, aes(x = sample,y = value, fill = type))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,   
           color = "black", 
           linewidth = 0.5) +
  geom_text(aes(label = value), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Red_Blue_Brown") +
  ggtitle("Sample-specific peaks")
```


# 11. Peak abundance distribution analysis [Supp. Fig. 8f]
```{r}
Em6ABP1_sp_peaks$peakID <- paste(Em6ABP1_sp_peaks$V1, Em6ABP1_sp_peaks$V2, Em6ABP1_sp_peaks$V3, Em6ABP1_sp_peaks$V6, sep = '_')
SYSY_sp_peaks$peakID <- paste(SYSY_sp_peaks$V1, SYSY_sp_peaks$V2, SYSY_sp_peaks$V3, SYSY_sp_peaks$V6, sep = '_')
Active_Motif_sp_peaks$peakID <- paste(Active_Motif_sp_peaks$V1, Active_Motif_sp_peaks$V2, Active_Motif_sp_peaks$V3, Active_Motif_sp_peaks$V6, sep = '_')
Millipore_sp_peaks$peakID <- paste(Millipore_sp_peaks$V1, Millipore_sp_peaks$V2, Millipore_sp_peaks$V3, Millipore_sp_peaks$V6, sep = '_')

Em6ABP1_sp_peaks_csv_filtered <- Em6ABP1_peak_csv_filtered[which(Em6ABP1_peak_csv_filtered$peakID %in% Em6ABP1_sp_peaks$peakID), ]
SYSY_sp_peaks_csv_filtered <- SYSY_peak_csv_filtered[which(SYSY_peak_csv_filtered$peakID %in% SYSY_sp_peaks$peakID), ]
Active_Motif_sp_peaks_csv_filtered <- Active_Motif_peak_csv_filtered[which(Active_Motif_peak_csv_filtered$peakID %in% Active_Motif_sp_peaks$peakID), ]
Millipore_sp_peaks_csv_filtered <- Millipore_peak_csv_filtered[which(Millipore_peak_csv_filtered$peakID %in% Millipore_sp_peaks$peakID), ]


Em6ABP1_all_abun <- Em6ABP1_peak_csv_filtered$RPM.IP
SYSY_all_abun <- SYSY_peak_csv_filtered$RPM.IP
Active_Motif_all_abun <- Active_Motif_peak_csv_filtered$RPM.IP
Millipore_all_abun <- Millipore_peak_csv_filtered$RPM.IP

Em6ABP1_sp_abun <- Em6ABP1_sp_peaks_csv_filtered$RPM.IP
SYSY_sp_abun <- SYSY_sp_peaks_csv_filtered$RPM.IP
Active_Motif_sp_abun <- Active_Motif_sp_peaks_csv_filtered$RPM.IP
Millipore_sp_abun <- Millipore_sp_peaks_csv_filtered$RPM.IP


max_length_all <- max(
  length(Em6ABP1_all_abun),
  length(SYSY_all_abun),
  length(Active_Motif_all_abun),
  length(Millipore_all_abun)
  )

max_length_sp <- max(
  length(Em6ABP1_sp_abun),
  length(SYSY_sp_abun),
  length(Active_Motif_sp_abun),
  length(Millipore_sp_abun)
  )


abun_all_df <- data.frame(
  No = seq(max_length_all),
  Em6ABP1_all = c(Em6ABP1_all_abun, rep(NA, max_length_all - length(Em6ABP1_all_abun))),
  SYSY_all = c(SYSY_all_abun, rep(NA, max_length_all - length(SYSY_all_abun))),
  Active_Motif_all = c(Active_Motif_all_abun, rep(NA, max_length_all - length(Active_Motif_all_abun))),
  Millipore_all = c(Millipore_all_abun, rep(NA, max_length_all - length(Millipore_all_abun)))
  )

abun_sp_df <- data.frame(
  No = seq(max_length_sp),
  Em6ABP1_sp = c(Em6ABP1_sp_abun, rep(NA, max_length_sp - length(Em6ABP1_sp_abun))),
  SYSY_sp = c(SYSY_sp_abun, rep(NA, max_length_sp - length(SYSY_sp_abun))),
  Active_Motif_sp = c(Active_Motif_sp_abun, rep(NA, max_length_sp - length(Active_Motif_sp_abun))),
  Millipore_sp = c(Millipore_sp_abun, rep(NA, max_length_sp - length(Millipore_sp_abun)))
  )


abun_all_df <- melt(abun_all_df, id = 'No') %>% drop_na() %>% mutate(log2abun = log2(value))
abun_sp_df <- melt(abun_sp_df, id = 'No') %>% drop_na() %>% mutate(log2abun = log2(value))


abun_all_df <- abun_all_df %>%
  mutate(
    level = case_when(
      value >= 1 & value < 5 ~ "1<=x<5",
      value >= 5 & value < 10 ~ "5<=x<10",
      value >= 10 & value < 50 ~ "10<=x<50",
      value >= 50 & value < 100 ~ "50<=x<100",
      value >= 100 ~ "x>=100"
    )
  )

abun_sp_df <- abun_sp_df %>%
  mutate(
    level = case_when(
      value >= 1 & value < 5 ~ "1<=x<5",
      value >= 5 & value < 10 ~ "5<=x<10",
      value >= 10 & value < 50 ~ "10<=x<50",
      value >= 50 & value < 100 ~ "50<=x<100",
      value >= 100 ~ "x>=100"
    )
  )


abun_all_df2 <- abun_all_df %>% group_by(variable, level) %>% tally()
abun_analysis_all <- abun_all_df2 %>%
  group_by(variable) %>%
  mutate(variable_total = sum(n),
    prop_in_variable = round(n / variable_total, 2))


abun_sp_df2 <- abun_sp_df %>% group_by(variable, level) %>% tally()
abun_analysis_sp <- abun_sp_df2 %>%
  group_by(variable) %>%
  mutate(variable_total = sum(n),
    prop_in_variable = round(n / variable_total, 2))


abun_analysis_all$level <- factor(abun_analysis_all$level, levels = c('1<=x<5', '5<=x<10', '10<=x<50', '50<=x<100', 'x>=100'))
abun_analysis_sp$level <- factor(abun_analysis_sp$level, levels = c('1<=x<5', '5<=x<10', '10<=x<50', '50<=x<100', 'x>=100'))


p1_abun <- ggplot(data = abun_analysis_all, aes(x = variable,y = prop_in_variable, fill = level))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,   
           color = "black", 
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "Peak abundance (RPM)", "ggthemes::Tableau_10") # 4 * 5.5

p2_abun <- ggplot(data = abun_analysis_sp, aes(x = variable,y = prop_in_variable, fill = level))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,   
           color = "black", 
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "Peak abundance (RPM)", "ggthemes::Tableau_10") # 4 * 5.5


plot_grid(plotlist = list(p1_abun, p2_abun), ncol = 1, align = 'h')
```


# 12. Overlap with GLORI m6A sites [Supp. Fig. 9a, Supp. Fig. 9b]
## 12.1 Proportion of GLORI m6A sites overlapped with all peaks [Supp. Fig. 9a]
### 12.1.1 GLORI RPM>=1
```{r}
# import bedtools result
GLORI_others_intersect <- read.table('results/peaks/intersect_with_GLORI/GLORI_others_intersected_1.bed', sep = '\t')
GLORI_others_intersect$GLORI_siteID <- paste(GLORI_others_intersect$V1, GLORI_others_intersect$V2, GLORI_others_intersect$V3, GLORI_others_intersect$V6, sep = '_')

# Em6ABP1 vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Em6ABP1_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Em6ABP1']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Em6ABP1_overlap_site_num
GLORI_Em6ABP1_df <- data.frame(sample = c('Em6ABP1', 'Em6ABP1'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Em6ABP1_overlap_site_num))

# SY vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_SYSY_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'SYSY'])) 
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_SYSY_overlap_site_num 
GLORI_SYSY_df <- data.frame(sample = c('SYSY', 'SYSY'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_SYSY_overlap_site_num))

# AM vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_Active_Motif_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Active_Motif']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Active_Motif_overlap_site_num 
GLORI_Active_Motif_df <- data.frame(sample = c('Active_Motif', 'Active_Motif'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Active_Motif_overlap_site_num))

# Milli vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Millipore_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Millipore']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Millipore_overlap_site_num
GLORI_Millipore_df <- data.frame(sample = c('Millipore', 'Millipore'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Millipore_overlap_site_num))

# integrate data
GLORI_df <- rbind(GLORI_Em6ABP1_df, GLORI_SYSY_df, GLORI_Active_Motif_df, GLORI_Millipore_df)


GLORI_df_more1 <- GLORI_df %>%
  group_by(sample) %>%
  mutate(variable_total = sum(value),
    prop_in_variable = round(value / variable_total, 2))

order <- c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore')
order2 <- c('Overlapped sites', 'GLORI-specific sites')
GLORI_df_more1$sample <-  factor(GLORI_df_more1$sample, levels = order)
GLORI_df_more1$type <- factor(GLORI_df_more1$type, levels = order2)


p1_GLORI <- ggplot(data = GLORI_df_more1, aes(x = sample, y = prop_in_variable, fill = type))+
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Tableau_20") +
  ggtitle('RPM>1')
```


### 12.1.2 GLORI 0.5<RPM<=1
```{r}
# import bedtools result
GLORI_others_intersect <- read.table('results/peaks/intersect_with_GLORI/GLORI_others_intersected_0.5_1.bed', sep = '\t')
GLORI_others_intersect$GLORI_siteID <- paste(GLORI_others_intersect$V1, GLORI_others_intersect$V2, GLORI_others_intersect$V3, GLORI_others_intersect$V6, sep = '_')

# Em6ABP1 vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Em6ABP1_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Em6ABP1'])) 
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Em6ABP1_overlap_site_num
GLORI_Em6ABP1_df <- data.frame(sample = c('Em6ABP1', 'Em6ABP1'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Em6ABP1_overlap_site_num))


# SY vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_SYSY_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'SYSY'])) 
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_SYSY_overlap_site_num 
GLORI_SYSY_df <- data.frame(sample = c('SYSY', 'SYSY'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_SYSY_overlap_site_num))


# AM vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_Active_Motif_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Active_Motif']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Active_Motif_overlap_site_num 
GLORI_Active_Motif_df <- data.frame(sample = c('Active_Motif', 'Active_Motif'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Active_Motif_overlap_site_num))


# Milli vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Millipore_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Millipore']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Millipore_overlap_site_num
GLORI_Millipore_df <- data.frame(sample = c('Millipore', 'Millipore'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Millipore_overlap_site_num))


# integrate data
GLORI_df <- rbind(GLORI_Em6ABP1_df, GLORI_SYSY_df, GLORI_Active_Motif_df, GLORI_Millipore_df)


GLORI_df_2 <- GLORI_df %>%
  group_by(sample) %>%
  mutate(variable_total = sum(value),
    prop_in_variable = round(value / variable_total, 2))

order <- c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore')
order2 <- c('Overlapped sites', 'GLORI-specific sites')
GLORI_df_2$sample <-  factor(GLORI_df_2$sample, levels = order)
GLORI_df_2$type <- factor(GLORI_df_2$type, levels = order2)


p2_GLORI <- ggplot(data = GLORI_df_2, aes(x = sample, y = prop_in_variable, fill = type))+
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Tableau_20") + 
  ggtitle('0.5<RPM<=1')
```


### 12.1.3 GLORI 0.25<RPM<=0.5
```{r}
# import bedtools result
GLORI_others_intersect <- read.table('results/peaks/intersect_with_GLORI/GLORI_others_intersected_0.25_0.5.bed', sep = '\t')
GLORI_others_intersect$GLORI_siteID <- paste(GLORI_others_intersect$V1, GLORI_others_intersect$V2, GLORI_others_intersect$V3, GLORI_others_intersect$V6, sep = '_')

# Em6ABP1 vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Em6ABP1_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Em6ABP1']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Em6ABP1_overlap_site_num
GLORI_Em6ABP1_df <- data.frame(sample = c('Em6ABP1', 'Em6ABP1'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Em6ABP1_overlap_site_num))

# SY vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_SYSY_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'SYSY'])) 
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_SYSY_overlap_site_num 
GLORI_SYSY_df <- data.frame(sample = c('SYSY', 'SYSY'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_SYSY_overlap_site_num))

# AM vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID)) 
GLORI_Active_Motif_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Active_Motif']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Active_Motif_overlap_site_num 
GLORI_Active_Motif_df <- data.frame(sample = c('Active_Motif', 'Active_Motif'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Active_Motif_overlap_site_num))

# Milli vs GLORI
GLORI_total_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID))
GLORI_Millipore_overlap_site_num <- length(unique(GLORI_others_intersect$GLORI_siteID[GLORI_others_intersect$V8 == 'Millipore']))
GLORI_sp_site_num <- GLORI_total_site_num - GLORI_Millipore_overlap_site_num
GLORI_Millipore_df <- data.frame(sample = c('Millipore', 'Millipore'), type = c('GLORI-specific sites', 'Overlapped sites'), value = c(GLORI_sp_site_num, GLORI_Millipore_overlap_site_num))

# integrate data
GLORI_df <- rbind(GLORI_Em6ABP1_df, GLORI_SYSY_df, GLORI_Active_Motif_df, GLORI_Millipore_df)


GLORI_df_2 <- GLORI_df %>%
  group_by(sample) %>%
  mutate(variable_total = sum(value),
    prop_in_variable = round(value / variable_total, 2))

order <- c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore')
order2 <- c('Overlapped sites', 'GLORI-specific sites')
GLORI_df_2$sample <-  factor(GLORI_df_2$sample, levels = order)
GLORI_df_2$type <- factor(GLORI_df_2$type, levels = order2)


p3_GLORI <- ggplot(data = GLORI_df_2, aes(x = sample, y = prop_in_variable, fill = type))+
  geom_bar(stat = "identity",
           position = "fill",
           width = 0.8,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "", "ggthemes::Tableau_20") + 
  ggtitle('0.25<RPM<=0.5')
```

### 12.1.4 GLORI result integration
```{r}
p_aligned <- align_plots(p3_GLORI, p2_GLORI, p1_GLORI, align = "h", axis = "l")
plot_grid(p_aligned[[1]], p_aligned[[2]], p_aligned[[3]], ncol = 3) # 14 * 9.5
```


## 12.2 Peaks overlapped with GLORI m6A sites [Supp. Fig. 9b]
```{r}
GLORI_data <- read.xlsx('/gold1/zhanghao/ZR/other_study/GLORI_HEK293T_single_m6A_site.xlsx', colNames = T)

GLORI_data$Chr <- gsub('chr', '', GLORI_data$Chr)
GLORI_data$start <- as.character(GLORI_data$Sites)
GLORI_data$end <- as.character(GLORI_data$Sites + 1)
GLORI_data$ACov_rep1 <- round(GLORI_data$AGCov_rep1 * GLORI_data$m6A_level_rep1, 0)
GLORI_data$ACov_rep2 <- round(GLORI_data$AGCov_rep2 * GLORI_data$m6A_level_rep2, 0)
GLORI_data$ACPM_rep1 <- round(GLORI_data$ACov_rep1 / 360, 2)
GLORI_data$ACPM_rep2 <- round(GLORI_data$ACov_rep2 / 360, 2)
GLORI_data$ACPM_avg <- rowMeans(dplyr::select(GLORI_data, c(ACPM_rep1, ACPM_rep2)))

GLORI_data_slim <- GLORI_data %>% dplyr::select(Chr, start, end, Gene, Cluster_info, Strand, ACPM_avg)


# separate GLORI m6A sites by abundance
GLORI_data_slim_0_0.1 <- GLORI_data_slim %>% filter(ACPM_avg > 0 & ACPM_avg <= 0.1)
GLORI_data_slim_0.1_0.25 <- GLORI_data_slim %>% filter(ACPM_avg > 0.1 & ACPM_avg <= 0.25)
GLORI_data_slim_0.25_0.5 <- GLORI_data_slim %>% filter(ACPM_avg > 0.25 & ACPM_avg <= 0.5)
GLORI_data_slim_0.5_1 <- GLORI_data_slim %>% filter(ACPM_avg > 0.5 & ACPM_avg <= 1)
GLORI_data_slim_1_5 <- GLORI_data_slim %>% filter(ACPM_avg > 1 & ACPM_avg <= 5)
GLORI_data_slim_5 <- GLORI_data_slim %>% filter(ACPM_avg > 5)


write.table(GLORI_data_slim, 'results/peaks/intersect_with_GLORI/GLORI_slim.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_0_0.1, 'results/peaks/intersect_with_GLORI/GLORI_slim_0_0.1.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_0.1_0.25, 'results/peaks/intersect_with_GLORI/GLORI_slim_0.1_0.25.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_0.25_0.5, 'results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_0.5_1, 'results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_1_5, 'results/peaks/intersect_with_GLORI/GLORI_slim_1_5.bed', sep = '\t', col.names = F, row.names = F, quote = F)
write.table(GLORI_data_slim_5, 'results/peaks/intersect_with_GLORI/GLORI_slim_5.bed', sep = '\t', col.names = F, row.names = F, quote = F)


# >>> to bash
# bedtools intersect


# import bedtools result
Em6ABP1_others_intersect <- read.table('results/peaks/intersect_with_GLORI/Em6ABP1_others_intersected.bed', sep = '\t')
SYSY_others_intersect <- read.table('results/peaks/intersect_with_GLORI/SYSY_others_intersected.bed', sep = '\t')
Active_Motif_others_intersect <- read.table('results/peaks/intersect_with_GLORI/Active_Motif_others_intersected.bed', sep = '\t')
Millipore_others_intersect <- read.table('results/peaks/intersect_with_GLORI/Millipore_others_intersected.bed', sep = '\t')


# stats by GLORI abundance
Em6ABP1_others_intersect$Em6ABP1_peakID <- paste(Em6ABP1_others_intersect$V1, Em6ABP1_others_intersect$V2, Em6ABP1_others_intersect$V3, Em6ABP1_others_intersect$V6, sep = '_')
SYSY_others_intersect$SYSY_peakID <- paste(SYSY_others_intersect$V1, SYSY_others_intersect$V2, SYSY_others_intersect$V3, SYSY_others_intersect$V6, sep = '_')
Active_Motif_others_intersect$Active_Motif_peakID <- paste(Active_Motif_others_intersect$V1, Active_Motif_others_intersect$V2, Active_Motif_others_intersect$V3, Active_Motif_others_intersect$V6, sep = '_')
Millipore_others_intersect$Millipore_peakID <- paste(Millipore_others_intersect$V1, Millipore_others_intersect$V2, Millipore_others_intersect$V3, Millipore_others_intersect$V6, sep = '_')


# Em6ABP1 vs GLORI
Em6ABP1_sp_others_intersect <- Em6ABP1_others_intersect[Em6ABP1_others_intersect$Em6ABP1_peakID %in% Em6ABP1_sp_peaks$peakID, ]
Em6ABP1_sp_GLORI_intersect <- Em6ABP1_sp_others_intersect %>% filter(V13 == 'GLORI')
Em6ABP1_sp_GLORI_overlap_peak_num <- length(unique(Em6ABP1_sp_GLORI_intersect$Em6ABP1_peakID))
GLORI_Em6ABP1_sp_df <- data.frame(sample = c('Em6ABP1', 'Em6ABP1'), type = c('Other peaks', 'Intersected with GLORI data'), value = c(length(unique(Em6ABP1_sp_others_intersect$Em6ABP1_peakID)) - Em6ABP1_sp_GLORI_overlap_peak_num, Em6ABP1_sp_GLORI_overlap_peak_num))


# SY vs GLORI
SYSY_sp_others_intersect <- SYSY_others_intersect[SYSY_others_intersect$SYSY_peakID %in% SYSY_sp_peaks$peakID, ]
SYSY_sp_GLORI_intersect <- SYSY_sp_others_intersect %>% filter(V13 == 'GLORI')
SYSY_sp_GLORI_overlap_peak_num <- length(unique(SYSY_sp_GLORI_intersect$SYSY_peakID))
SYSY_sp_GLORI_df <- data.frame(sample = c('SYSY', 'SYSY'), type = c('Other peaks', 'Intersected with GLORI data'), value = c(length(unique(SYSY_sp_others_intersect$SYSY_peakID)) - SYSY_sp_GLORI_overlap_peak_num, SYSY_sp_GLORI_overlap_peak_num))


# AM vs GLORI
Active_Motif_sp_others_intersect <- Active_Motif_others_intersect[Active_Motif_others_intersect$Active_Motif_peakID %in% Active_Motif_sp_peaks$peakID, ]
Active_Motif_sp_GLORI_intersect <- Active_Motif_sp_others_intersect %>% filter(V13 == 'GLORI')
Active_Motif_sp_GLORI_overlap_peak_num <- length(unique(Active_Motif_sp_GLORI_intersect$Active_Motif_peakID))
Active_Motif_sp_GLORI_df <- data.frame(sample = c('Active_Motif', 'Active_Motif'), type = c('Other peaks', 'Intersected with GLORI data'), value = c(length(unique(Active_Motif_sp_others_intersect$Active_Motif_peakID)) - Active_Motif_sp_GLORI_overlap_peak_num, Active_Motif_sp_GLORI_overlap_peak_num))


# Milli vs GLORI
Millipore_sp_others_intersect <- Millipore_others_intersect[Millipore_others_intersect$Millipore_peakID %in% Millipore_sp_peaks$peakID, ]
Millipore_sp_GLORI_intersect <- Millipore_sp_others_intersect %>% filter(V13 == 'GLORI')
Millipore_sp_GLORI_overlap_peak_num <- length(unique(Millipore_sp_GLORI_intersect$Millipore_peakID))
Millipore_sp_GLORI_df <- data.frame(sample = c('Millipore', 'Millipore'), type = c('Other peaks', 'Intersected with GLORI data'), value = c(length(unique(Millipore_sp_others_intersect$Millipore_peakID)) - Millipore_sp_GLORI_overlap_peak_num, Millipore_sp_GLORI_overlap_peak_num))


# integrate data and plot
sp_GLORI_df <- rbind(Em6ABP1_sp_GLORI_df, SYSY_sp_GLORI_df, Active_Motif_sp_GLORI_df, Millipore_sp_GLORI_df)

order <- c('Em6ABP1', 'SYSY', 'Active_Motif', 'Millipore')
order2 <- c('Peaks overlapped with GLORI m6A sites', 'Other peaks')

sp_GLORI_df$sample <-  factor(sp_GLORI_df$sample, levels = order)
sp_GLORI_df$type <- factor(sp_GLORI_df$type, levels = order2)

ggplot(data = sp_GLORI_df, aes(x = sample, y = value, fill = type)) +
      geom_bar(stat = "identity", 
               position = "fill",
               width = 0.8,   
               color = "black", 
               linewidth = 0.5) +
      geom_text(aes(label = value), position=position_fill(vjust=0.5)) + 
      theme_bw(base_size = 15) + labs(x = '', y = 'Peak number') +
      scale_fill_paletteer_d(name = "", "ggthemes::Tableau_20") +
      ggtitle('Sample-specific peaks')
```

