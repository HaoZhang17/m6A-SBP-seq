library(exomePeak2)
library(BiocParallel)
register(SerialParam())

set.seed(1)

args <- commandArgs(trailingOnly = TRUE)

GENE_ANNO_GTF <- args[1]
INPUT_BAM <- c(args[2])
IP_BAM <- c(args[3])
EXP_NAME <- args[4]

print(IP_BAM)

exomePeak2(bam_input = INPUT_BAM,
           bam_ip = IP_BAM,
           gff = GENE_ANNO_GTF,
           strandness = "unstrand",
           fragment_length = 100,
           bin_size = 25,
           step_size = 25,
           p_cutoff = 1e-5,
           experiment_name = EXP_NAME,
           mode = "exon")