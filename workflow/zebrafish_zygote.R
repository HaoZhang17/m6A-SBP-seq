# 1. Import library
```{r}
library(AnnotationDbi)
library(cowplot)
library(dplyr)
library(ggplot2)
library(ggsci)
library(Guitar)
library(openxlsx)
library(org.Dr.eg.db)
library(paletteer)
library(reshape2)
library(tidyverse)
library(topGO)
```


# 2. Import peak files
```{r}
dir_path = 'results/peaks'

ZF_10zygotes_num <- 1:5
ZF_single_zygote_num <- 1:10

# for bed file
ZF_peak_F_bed_list <- list()
ZF_peak_R_bed_list <- list()

## forward strand
for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_forward_strand')
  filepath <- paste0(dir_path, 'ZF_single_zygote_', i, '_forward_strand/peaks.bed')
  bed_file <- read.table(filepath, sep = '\t')
  ZF_peak_F_bed_list[[sample_name]] <- bed_file %>% filter(V6 == '+')
}

for (i in ZF_10zygotes_num){
  sample_name <- paste0('ZF_10zygotes_', i, '_forward_strand')
  filepath <- paste0(dir_path, 'ZF_10zygotes_', i, '_forward_strand/peaks.bed')
  bed_file <- read.table(filepath, sep = '\t')
  ZF_peak_F_bed_list[[sample_name]] <- bed_file %>% filter(V6 == '+')
}

## reverse strand
for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_reverse_strand')
  filepath <- paste0(dir_path, 'ZF_single_zygote_', i, '_reverse_strand/peaks.bed')
  bed_file <- read.table(filepath, sep = '\t')
  ZF_peak_R_bed_list[[sample_name]] <- bed_file %>% filter(V6 == '-')
}

for (i in ZF_10zygotes_num){
  sample_name <- paste0('ZF_10zygotes_', i, '_reverse_strand')
  filepath <- paste0(dir_path, 'ZF_10zygotes_', i, '_reverse_strand/peaks.bed')
  bed_file <- read.table(filepath, sep = '\t')
  ZF_peak_R_bed_list[[sample_name]] <- bed_file %>% filter(V6 == '-')
}


## combine forward and reverse strands
ZF_peak_bed_list <- list()
for (i in 1:15){
  ZF_peak_bed_list[[i]] <- rbind(ZF_peak_F_bed_list[[i]], ZF_peak_R_bed_list[[i]])
}
names(ZF_peak_bed_list)[1:10] <- paste0('ZF_single_zygote_', ZF_single_zygote_num, '_peak_bed')
names(ZF_peak_bed_list)[11:15] <- paste0('ZF_10zygotes_', ZF_10zygotes_num, '_peak_bed')


# for csv files
ZF_peak_F_csv_list <- list()
ZF_peak_R_csv_list <- list()

## forward strand
for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_forward_strand')
  filepath <- paste0(dir_path, 'ZF_single_zygote_', i, '_forward_strand/peaks.csv')
  csv_file <- read.csv(filepath, header = T)
  ZF_peak_F_csv_list[[sample_name]] <- csv_file %>% filter(strand == '+')
}

for (i in ZF_10zygotes_num){
  sample_name <- paste0('ZF_10zygotes_', i, '_forward_strand')
  filepath <- paste0(dir_path, 'ZF_10zygotes_', i, '_forward_strand/peaks.csv')
  csv_file <- read.csv(filepath, header = T)
  ZF_peak_F_csv_list[[sample_name]] <- csv_file %>% filter(strand == '+')
}

## reverse strand
for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_reverse_strand')
  filepath <- paste0(dir_path, 'ZF_single_zygote_', i, '_reverse_strand/peaks.csv')
  csv_file <- read.csv(filepath, header = T)
  ZF_peak_R_csv_list[[sample_name]] <- csv_file %>% filter(strand == '-')
}

for (i in ZF_10zygotes_num){
  sample_name <- paste0('ZF_10zygotes_', i, '_reverse_strand')
  filepath <- paste0(dir_path, 'ZF_10zygotes_', i, '_reverse_strand/peaks.csv')
  csv_file <- read.csv(filepath, header = T)
  ZF_peak_R_csv_list[[sample_name]] <- csv_file %>% filter(strand == '-')
}


## combine forward and reverse strands
ZF_peak_csv_list <- list()
for (i in 1:15){
  sample_name <- names(ZF_peak_F_csv_list)[i]
  ZF_peak_csv_list[[sample_name]] <- rbind(ZF_peak_F_csv_list[[i]], ZF_peak_R_csv_list[[i]])
}
```


# 3. Filter peaks
```{r}
# 3.1. Filter peaks with condition: read count > 10
# ratio = number of uniquely mapped reads/1,000,000
ratio <- c(5.194970, 5.193945, 5.187488, 5.197082, 5.197815, 5.195337, 5.194416, 5.193312, 5.194450, 5.196878, 8.880170, 8.885278,  8.881920, 8.875920, 8.883994)

for (i in 1:15){
  ZF_peak_csv_list[[i]] <- ZF_peak_csv_list[[i]] %>% filter(RPM.IP >= 10 / ratio[i])
}


# 3.2. Filter peaks with conditions: IP/INPUT >= 2, RPM.IP >= 1
peak_filter <- function(peak_csv, peak_bed) {
  peak_csv$Fold <- peak_csv$RPM.IP / peak_csv$RPM.input
  peak_csv$peakID <- paste(peak_csv$chr, peak_csv$chromStart, peak_csv$chromEnd, peak_csv$strand, sep = '_')
  peak_bed$peakID <- paste(peak_bed$V1, peak_bed$V2, peak_bed$V3, peak_bed$V6, sep = '_')
  peak_csv_filtered <- peak_csv %>% filter(Fold >= 2 & RPM.IP >= 1)
  peak_bed_filtered <- peak_bed[which(peak_bed$peakID %in% peak_csv_filtered$peakID), ]
  peak_bed_filtered_join <- merge(peak_bed_filtered, peak_csv_filtered, by = 'peakID', all.x = T)
  peak_bed_filtered_join <- peak_bed_filtered_join[, c(2, 3, 4, 1, 6, 7, 8, 9, 10, 11, 12, 13)]
  return(list(peak_csv_filtered, peak_bed_filtered_join))
}


ZF_filtered_peak_bed_list <- list()
ZF_filtered_peak_csv_list <- list()

sample_name1 <- paste0('ZF_single_zygote_', ZF_single_zygote_num)
sample_name2 <- paste0('ZF_10zygotes_', ZF_10zygotes_num)
sample_name <- c(sample_name1, sample_name2)
dir_path <- paste0(dir_path, 'ZF_combined/')
for (i in 1:15){
  ZF_filtered_peak_list <- peak_filter(ZF_peak_csv_list[[i]], ZF_peak_bed_list[[i]])
  
  ZF_filtered_peak_csv_list[i] <- ZF_filtered_peak_list[1]
  names(ZF_filtered_peak_csv_list)[i] <- names(ZF_peak_csv_list)[i]
  
  ZF_filtered_peak_bed_list[i] <- ZF_filtered_peak_list[2]
  names(ZF_filtered_peak_bed_list)[i] <- names(ZF_peak_bed_list)[i]
  
  write.table(ZF_filtered_peak_list[2], paste(dir_path, sample_name[i], '_peak_filtered.bed', sep = ''), sep = '\t', col.names = F, row.names = F, quote = F)
}
```


# 4. Peak number plot [Fig. 5a]
```{r}
#  peaks
peak_num_df <- data.frame(sample = c(paste0("ZF_10zygotes_", ZF_10zygotes_num), paste0("ZF_single_zygote_", ZF_single_zygote_num)),
                       peak_num = c(13058, 12459, 12646, 12416, 12460, 14315, 13353, 14142, 15024, 12581, 12332, 12232, 14419, 11696, 13932),
                       type = c(rep('10IP', 5), rep('scIP', 10)))


peak_num_df$sample <- factor(peak_num_df$sample, levels = c(paste0("ZF_single_zygote_", ZF_single_zygote_num), paste0("ZF_10zygotes_", ZF_10zygotes_num)))

ggplot(data = peak_num_df, aes(x = sample, y = peak_num, fill = type)) + 
  geom_bar(stat = "identity", size = 1.5, width = 0.7) +
  geom_text(aes(label = peak_num), size = 4.5, vjust = -0.5, position = position_dodge(0.3)) +
  theme_classic(base_size = 15) + 
  theme(legend.position = "none", axis.text.x = element_text(angle=45, hjust=1)) +
  labs(x = '', y = 'Peak number') + 
  scale_fill_paletteer_d(name = '', "ggthemes::Classic_Cyclic") +
  scale_y_continuous(
    limits = c(0, 17000),
    expand = c(0, 0)
  )
```


# 5. Metagene analysis [Fig. 5b]
```{r}
options(stringsAsFactors = F)

txdb_zebrafish <- txdbmaker::makeTxDbFromGFF(file = "resources/reference/Danio_rerio.GRCz11.115.gtf", 
                        format="gtf", 
                        dataSource="Ensembl", 
                        organism="Danio rerio")


stBedFiles <- list(paste0(dir_path, 'ZF_10zygotes_1_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_10zygotes_2_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_10zygotes_3_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_10zygotes_4_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_10zygotes_5_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_1_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_2_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_3_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_4_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_5_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_6_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_7_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_8_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_9_peak_filtered.bed'),
                   paste0(dir_path, 'ZF_single_zygote_10_peak_filtered.bed'))


p_mRNA_all <- GuitarPlot(txTxdb = txdb_zebrafish, 
                stBedFiles = stBedFiles,
                headOrtail = FALSE,
                enableCI = FALSE, 
                mapFilterTranscript = TRUE, 
                pltTxType = c("mrna"), 
                stGroupName = c(paste0("ZF_10zygotes_", ZF_10zygotes_num), 
                                paste0("ZF_single_zygote_", ZF_single_zygote_num)))


p_mRNA_all_data <- p_mRNA_all$data

p_mRNA_all_data$group <- factor(p_mRNA_all_data$group, 
                                levels = c(paste0("ZF_single_zygote_", ZF_single_zygote_num), 
                                          paste0("ZF_10zygotes_", ZF_10zygotes_num)))

ggplot(p_mRNA_all_data, aes(x = x, y = density, fill = group, color = group)) +
  geom_line(linewidth = 1.2, alpha = 1) +
  geom_ribbon(aes(ymin = 0, ymax = density), 
              alpha = 0.1,
              size = 0,
              show.legend = FALSE) + 
  scale_x_continuous(
    name = NULL,
    breaks = c(0, 0.120, 0.652, 1),
    labels = c("", "Start", "Stop", ""),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Frequency (%)",
    breaks = seq(0, 5, by = 1),
    limits = c(0, 5),
    expand = c(0, 0)
  ) +
  geom_vline(xintercept = c(0.120, 0.652), 
             linetype = "dashed", 
             color = "gray50", 
             alpha = 0.7, 
             size = 0.5) +
  annotate("text", x = 0.07, y = 4.5, 
           label = "5'UTR", size = 4, fontface = "bold") +
  annotate("text", x = 0.35, y = 4.5, 
           label = "CDS", size = 4, fontface = "bold") +
  annotate("text", x = 0.75, y = 4.5, 
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
  )
```


# 6. Distribution of peaks across transcript regions [Fig. 5c]
```{r}
process_annotation <- function(sample_name, dir_path) {
  stats_file <- file.path(dir_path, paste0(sample_name, "_peak_anno.stats.txt"))
  anno_stats <- read.table(stats_file, sep = "\t", header = TRUE, quote = "")
  
  anno_stats_1 <- anno_stats[1:11, c(1, 2)] %>%
    filter(ZF_annotation != '0.0') %>%
    mutate(ZF_annotation = as.numeric(ZF_annotation))
  
  anno_stats_1$Annotation[anno_stats_1$Annotation == 'Promoter'] <- 'TSS'
  anno_stats_1$Annotation[anno_stats_1$Annotation == 'Exon'] <- 'CDS'
  
  ncRNA_categories <- c('pseudo', 'miRNA', 'ncRNA', 'rRNA')
  ncRNA_sum <- sum(anno_stats_1$ZF_annotation[anno_stats_1$Annotation %in% ncRNA_categories])
  anno_stats_1 <- anno_stats_1[!anno_stats_1$Annotation %in% ncRNA_categories, ]
  anno_stats_1 <- rbind(anno_stats_1, data.frame(Annotation = 'ncRNA', ZF_annotation = ncRNA_sum))
  
  total_peaks <- sum(anno_stats_1$ZF_annotation)
  anno_stats_1 <- anno_stats_1 %>%
    mutate(per = round(ZF_annotation / total_peaks * 100, digits = 2),
           Anno_new = paste(Annotation, ' (', per, '%)', sep = ''))
  
  anno_stats_melt <- anno_stats_1 %>%
    select(Annotation, ZF_annotation) %>%
    rename(value = ZF_annotation)
  
  anno_file <- file.path(dir_path, paste0(sample_name, "_peak_anno.txt"))
  anno_data <- read.table(anno_file, sep = "\t", header = TRUE, quote = "")

  return(anno_stats_melt)
}


all_samples <- c(paste0("ZF_10zygotes_", 1:5), paste0("ZF_single_zygote_", 1:10))

# change the column name "value" to corresponding sample name
results <- lapply(all_samples, function(sample) {
  res <- process_annotation(sample, dir_path)
  res$stats <- res$stats %>%
    rename(!!sym(sample) := value)
  res
})

# merge stat tables for all samples
ZF_annotation <- results %>%
  lapply(function(x) x$stats) %>%
  reduce(full_join, by = "Annotation")

ZF_annotation_melt <- ZF_annotation %>%
  pivot_longer(cols = -Annotation, names_to = "variable", values_to = "value")

new_order <- c('TSS', '5UTR', 'CDS', '3UTR', 'Intron', 'TTS', 'Intergenic', 'ncRNA')
ZF_annotation_melt$Annotation <- factor(ZF_annotation_melt$Annotation, levels = new_order)

ZF_plot_data <- ZF_annotation_melt %>%
  group_by(variable) %>%
  arrange(desc(Annotation)) %>%
  mutate(cum_value = cumsum(value)) %>%
  ungroup()

ggplot(data = ZF_annotation_melt, aes(x = variable, y = value, group = Annotation)) +
  geom_area(aes(fill = Annotation), position = "stack", alpha = 0.4) +
  geom_line(aes(color = Annotation), position = position_stack(), linewidth = 0.8) +
  geom_point(aes(color = Annotation), position = position_stack(), size = 1.5) +
  theme_bw(base_size = 15) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_paletteer_d("ggthemes::Tableau_10") +
  scale_fill_paletteer_d("ggthemes::Tableau_10") +
  labs(y = "Cumulative Proportion", x = "Sample")
```



# 7. Peak overlap analysis [Fig. 5e]
## Intervene
```{r}
dir_path <- 'results/peaks/'
peak_overlap_data <- as.matrix(read.table(paste(dir_path, "Intervene_results/Intervene_pairwise_frac_matrix.txt", sep = ''), header = TRUE, row.names = 1, sep = "\t"))

# to percentage
peak_overlap_data <- peak_overlap_data * 100

colnames(peak_overlap_data) <- c(paste0('ZF_10zygotes_', ZF_10zygotes_num), paste0('ZF_single_zygote_', ZF_single_zygote_num))
rownames(peak_overlap_data) <- c(paste0('ZF_10zygotes_', ZF_10zygotes_num), paste0('ZF_single_zygote_', ZF_single_zygote_num))

sample_name <- c(paste0('ZF_10zygotes_', ZF_10zygotes_num), paste0('ZF_single_zygote_', ZF_single_zygote_num))
peak_overlap_data_2 <- peak_overlap_data[match(sample_name, colnames(peak_overlap_data)), match(sample_name, rownames(peak_overlap_data))]


peak_overlap_data_melt <- reshape2::melt(peak_overlap_data_2)
colnames(peak_overlap_data_melt) <- c("Row", "Column", "Value")
peak_overlap_data_melt$Row <- factor(peak_overlap_data_melt$Row, levels = colnames(peak_overlap_data_2))
peak_overlap_data_melt$Column <- factor(peak_overlap_data_melt$Column, levels = colnames(peak_overlap_data_2))

ggplot(peak_overlap_data_melt, aes(x = Column, y = Row, fill = Value)) +
  geom_tile(color = "white", size = 0.5) +
  geom_tile(color="white", size=0.5) +
  geom_text(aes(label=paste0(round(Value,0))), size=3) +
  scale_fill_gradientn(
    colours = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    limits = c(0, 100)) +
  theme_minimal(base_size = 15) +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  coord_fixed() +
  labs(x = '', y = '', title = 'Peak overlap')
```



# 8. Comparison with other studies [Fig. 6a]
## 8.1. m6A-SBP-seq
```{r}
# obtain all m6A-modified RNAs in all single zygotes 
combined_gene_SBPdata <- vector()
for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_forward_strand')
  ZF_single_zygote_i_gene <- unique(ZF_filtered_peak_csv_list[[sample_name]]$geneID)
  
  ZF_single_zygote_i_gene_num <- length(ZF_single_zygote_i_gene)
  combined_gene_SBPdata <- c(combined_gene_SBPdata, ZF_single_zygote_i_gene)
}
combined_gene_SBPdata <- unique(combined_gene_SBPdata)
```


## 8.2. Zhao et al
```{r}
# The bed file was downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2088167
ZF_Zhao_peak <- read.table('other_study/Zhao/GSM2088167_m6A_0h.peaks.bed', sep = '\t')
ZF_Zhao_peak_geneID <- ZF_Zhao_peak$V1

GeneInfoTable <- select(org.Dr.eg.db, keys = ZF_Zhao_peak_geneID, columns = c("SYMBOL", "ENTREZID", "CHR", "ENSEMBL"), keytype = "ACCNUM" )

ZF_Zhao_peak_merged <- merge(ZF_Zhao_peak, GeneInfoTable, by.x = 'V1', by.y = 'ACCNUM')
ZF_Zhao_peak_merged <- unique(ZF_Zhao_peak_merged)
ZF_Zhao_peak_gene <- unique(ZF_Zhao_peak_merged$ENSEMBL)
```


## 8.3. Li et al
```{r}
ZF_Li_single_zygote_num <- c(1, 2, 3, 4, 5)

ZF_Li_peak_bed_list <- list()
ZF_Li_peak_csv_list <- list()

dir_path = 'other_study/Li/'

for (i in ZF_Li_single_zygote_num){
  sample_name <- paste0('Li_zIP1_', i)
  filepath_bed <- paste0(dir_path, 'Li_zIP1_', i, '/peaks.bed')
  ZF_Li_peak_bed_list[[sample_name]] <- read.table(filepath_bed, sep = '\t')

  filepath_csv <- paste0(dir_path, 'Li_zIP1_', i, '/peaks.csv')
  ZF_Li_peak_csv_list[[sample_name]] <- read.csv(filepath_csv, header = T)
}


# Filter peaks with condition: read count > 10
# ratio = number of uniquely mapped reads/1,000,000
ratio <- c(1.08, 1.08, 1.08, 1.08, 1.08)

for (i in 1:5){
  ZF_Li_peak_csv_list[[i]] <- ZF_Li_peak_csv_list[[i]] %>% filter(RPM.IP >= 10 / ratio[i])
}


# Filter peaks with conditions: IP/INPUT >= 2, RPM.IP >= 1
peak_filter <- function(peak_csv, peak_bed) {
  peak_csv$Fold <- peak_csv$RPM.IP / peak_csv$RPM.input
  peak_csv$peakID <- paste(peak_csv$chr, peak_csv$chromStart, peak_csv$chromEnd, peak_csv$strand, sep = '_')
  peak_bed$peakID <- paste(peak_bed$V1, peak_bed$V2, peak_bed$V3, peak_bed$V6, sep = '_')
  peak_csv_filtered <- peak_csv %>% filter(Fold >= 2 & RPM.IP >= 1)
  peak_bed_filtered <- peak_bed[which(peak_bed$peakID %in% peak_csv_filtered$peakID), ]
  peak_bed_filtered_join <- merge(peak_bed_filtered, peak_csv_filtered, by = 'peakID', all.x = T)
  peak_bed_filtered_join <- peak_bed_filtered_join[, c(2, 3, 4, 1, 6, 7, 8, 9, 10, 11, 12, 13)]
  return(list(peak_csv_filtered, peak_bed_filtered_join))
}


ZF_Li_filtered_peak_bed_list <- list()
ZF_Li_filtered_peak_csv_list <- list()
ZF_Li_peak_gene <- ""

sample_name <- paste('Li_zIP1_', ZF_Li_single_zygote_num, sep = '')
dir_path <- paste0(dir_path, 'ZF_combined/')
for (i in 1:5){
  ZF_Li_filtered_peak_list <- peak_filter(ZF_Li_peak_csv_list[[i]], ZF_Li_peak_bed_list[[i]])
  
  ZF_Li_filtered_peak_csv_list[i] <- ZF_Li_filtered_peak_list[1]
  names(ZF_Li_filtered_peak_csv_list)[i] <- names(ZF_Li_peak_csv_list)[i]
  
  ZF_Li_filtered_peak_bed_list[i] <- ZF_Li_filtered_peak_list[2]
  names(ZF_Li_filtered_peak_bed_list)[i] <- names(ZF_Li_peak_bed_list)[i]
  
  ZF_Li_peak_gene <-  c(ZF_Li_peak_gene, ZF_Li_filtered_peak_list[[1]]$geneID)
  
  write.table(ZF_Li_filtered_peak_list[2], paste(dir_path, sample_name[i], '_peak_filtered.bed', sep = ''), sep = '\t', col.names = F, row.names = F, quote = F)
}

ZF_Li_peak_gene <- unique(ZF_Li_peak_gene)
```


## 8.4. Venn plot of m6A-modified RNAs [Fig. 6a]
```{R}
library(eulerr)
ZF_Li_peak_gene_clean <- ZF_Li_peak_gene[!grepl(",", ZF_Li_peak_gene)]
ZF_Zhao_peak_gene_clean <- ZF_Zhao_peak_gene[!grepl(",", ZF_Zhao_peak_gene)]
ZF_Zhao_peak_gene_clean <- ZF_Zhao_peak_gene_clean[!is.na(ZF_Zhao_peak_gene_clean)]
combined_gene_SBPdata_clean <- combined_gene_SBPdata[!grepl(",", combined_gene_SBPdata)]

ZF_gene_list <- list()
ZF_gene_list[['Li']] <- ZF_Li_peak_gene_clean
ZF_gene_list[['Zhao']] <- ZF_Zhao_peak_gene_clean
ZF_gene_list[['SBP-seq']] <- combined_gene_SBPdata_clean

plot(euler(
     ZF_gene_list,
     shape = "circle"),                    
     quantities = list(type = c("percent","counts"), cex= 1),
     labels = list(cex= 1),
     edges = list(col = c("#659757", "#8790B1", "#C9706A"), lex = 2), 
     fills = list(fill = c("#659757", "#8790B1", "#C9706A"), alpha=0.4)
)
```


## 9. Peak abundance analysis [Fig. 6b]
```{r}
ZF_peak_abun_specific <- list()
ZF_peak_abun_common <- list()

ZF_peak_abun_specific_inte <- vector()
ZF_peak_geneID_specific_inte <- vector()
ZF_peak_abun_common_inte <- vector()
ZF_peak_geneID_common_inte <- vector()


for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_forward_strand')
  ZF_single_zygote_peak_sample_i <- ZF_filtered_peak_csv_list[[sample_name]]
  ZF_peak_sample_i_slim <- ZF_single_zygote_peak_sample_i[ZF_single_zygote_peak_sample_i$geneID %in% ZF_SBP_specific_gene, c('RPM.IP', 'geneID')]
  ZF_peak_abun_specific[[sample_name]] <- ZF_peak_sample_i_slim
  ZF_peak_abun_specific_inte <- c(ZF_peak_abun_specific_inte, ZF_peak_sample_i_slim$RPM.IP)
  ZF_peak_geneID_specific_inte <- c(ZF_peak_geneID_specific_inte, ZF_peak_sample_i_slim$geneID)
}


for (i in ZF_single_zygote_num){
  sample_name <- paste0('ZF_single_zygote_', i, '_forward_strand')
  ZF_single_zygote_peak_sample_i <- ZF_filtered_peak_csv_list[[sample_name]]
  ZF_peak_sample_i_slim <- ZF_single_zygote_peak_sample_i[ZF_single_zygote_peak_sample_i$geneID %in% ZF_SBP_common_gene, c('RPM.IP', 'geneID')]
  ZF_peak_abun_common[[sample_name]] <- ZF_peak_sample_i_slim
  ZF_peak_abun_common_inte <- c(ZF_peak_abun_common_inte, ZF_peak_sample_i_slim$RPM.IP)
  ZF_peak_geneID_common_inte <- c(ZF_peak_geneID_common_inte, ZF_peak_sample_i_slim$geneID)
}


ZF_peak_abun_df <- data.frame(type = c(rep('specific', length(ZF_peak_abun_specific_inte)), 
                                       rep('common', length(ZF_peak_abun_common_inte))), 
                                 value = c(ZF_peak_abun_specific_inte, ZF_peak_abun_common_inte),
                                 geneID = c(ZF_peak_geneID_specific_inte, ZF_peak_geneID_common_inte))

# count peak number
ZF_peak_abun_df2 <- ZF_peak_abun_df %>%
  mutate(
    level = case_when(
      value >= 1 & value < 5 ~ "1<=x<5",
      value >= 5 & value < 10 ~ "5<=x<10",
      value >= 10 & value < 50 ~ "10<=x<50",
      value >= 50 & value < 100 ~ "50<=x<100",
      value >= 100 ~ "x>=100"
    )
  )

# calculate percentage
ZF_peak_abun_df3 <- ZF_peak_abun_df2 %>% group_by(type, level) %>% tally()
ZF_peak_abun_df3 <- ZF_peak_abun_df3 %>%
  group_by(type) %>%
  mutate(variable_total = sum(n),
    prop_in_variable = round(n / variable_total,2 ))


ZF_peak_abun_df3$level <- factor(ZF_peak_abun_df3$level, levels = c('1<=x<5', '5<=x<10', '10<=x<50', '50<=x<100', 'x>=100'))
ZF_peak_abun_df3$type <- factor(ZF_peak_abun_df3$type, levels = c('specific', 'common'))


ggplot(data = ZF_peak_abun_df3, aes(x = type,y = prop_in_variable, fill = level))+
  geom_bar(stat = "identity", 
           position = "fill",
           width = 0.8,
           color = "black",
           linewidth = 0.5) +
  geom_text(aes(label = prop_in_variable), position=position_fill(vjust=0.5)) + 
  theme_bw(base_size = 15) + labs(x = '', y = 'Proportion') +
  scale_fill_paletteer_d(name = "Peak abundance (RPM)", "ggthemes::Tableau_10") # 6 * 5
```


## 10. Function of RNAs with high-abundance peaks (RPM > 50) [Fig. 6c]
```{r}
# obtain ID for RNAs with high-abundance peaks (RPM > 50)
ZF_peak_geneID_sp_50_100 <- unique(ZF_peak_abun_df2$geneID[ZF_peak_abun_df2$value < 100 & ZF_peak_abun_df2$value >= 50 & ZF_peak_abun_df2$type == 'specific']) # 749
ZF_peak_geneID_sp_100 <- unique(ZF_peak_abun_df2$geneID[ZF_peak_abun_df2$value > 100 & ZF_peak_abun_df2$type == 'specific']) # 317
ZF_peak_geneID_sp_high_abun <- c(ZF_peak_geneID_sp_50_100, ZF_peak_geneID_sp_100)

# Functional enrichment analysis for RNAs with high-abundance peaks (RPM > 50) was conducted using Metascape
# GO_high_abundance_peak.xlsx is the report file from Metascape
ZF_sp_peak_50_func_GO <- read.xlsx('other_study/GO_high_abundance_peak.xlsx', 'Enrichment')

ZF_sp_peak_50_func_GO <- ZF_sp_peak_50_func_GO %>% filter(grepl('Summary', GroupID))
ZF_sp_peak_50_func_GO$Gene_count <- as.numeric(sub("/.*", "", ZF_sp_peak_50_func_GO$InTerm_InList))

ggplot(ZF_sp_peak_50_func_GO) + 
    geom_col(aes(x = reorder(Description, -LogP), y = -LogP), fill = '#7AA6DC') +
    geom_text(aes(x = reorder(Description, -LogP), y = 0.1, label = Description), size = 3.3, hjust = 0) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.ticks.length = unit(1.5, "mm"), axis.text = element_text(size = 10, color = 'black')) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 9)) +
    labs(y = '-log10(P value)', x = 'GO category (biological process)') # 5 * 6
```
