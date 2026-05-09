configfile: "config.yaml"


# 1. Remove adapter
rule qc_cutadapt:
    params:
        min_length = 20
    input:
        read1 = "resources/fastq/{sample}_1.fq.gz",
        read2 = "resources/fastq/{sample}_2.fq.gz"
    output:
        read1 = temp("resources/fastq/{sample}_qc_1.fq.gz"),
        read2 = temp("resources/fastq/{sample}_qc_2.fq.gz")
    threads: 10
    log:
        "logs/cutadapt/{sample}_cutadapt.log"
    shell:
        """
        cutadapt -j {threads} \
            -q 20 \
            -m {params.min_length} \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
            -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGT \
            --nextseq-trim=20 \
            --poly-a \
            -o {output.read1} \
            -p {output.read2} \
            {input.read1} \
            {input.read2}
            2> {log}
        """


# 2. Remove polyG and low-quality bases
rule qc_fastp:
    input:
        read1 = "resources/fastq/{sample}_qc_1.fq.gz",
        read2 = "resources/fastq/{sample}_qc_2.fq.gz"
    output:
        read1 = temp("resources/fastq/{sample}_qc_2nd_1.fq.gz"),
        read2 = temp("resources/fastq/{sample}_qc_2nd_2.fq.gz")
    threads: 10
    log:
        "logs/fastp/{sample}_fastp.log"
    shell:
        """
        fastp -i {input.read1} \
            -I {input.read2} \
            -o {output.read1} \
            -O {output.read2} \
            --trim_poly_g \
            --trim_poly_x \
            --detect_adapter_for_pe \
            --length_required 15 \
            -q 20 \
            --thread {threads} \
            > {log}
        """


# 3. Remove reads belonging to rRNA, tRNA sequences
rule remove_contaminants:
    input:
        read_1 = "resources/fastq/{sample}_qc_2nd_1.fq.gz",
        read_2 = "resources/fastq/{sample}_qc_2nd_2.fq.gz",
    output:
        read_1 = "resources/fastq/{sample}_qc_3rd_1.fq.gz",
        read_2 = "resources/fastq/{sample}_qc_3rd_2.fq.gz",
    log:
        "logs/QC/{sample}_remove_contaminants.log"
    threads: 10
    params:
        ref = config['contamination_fa']
    shell:
        """
        bbduk.sh in={input.read_1} \
            in2={input.read_2} \
            out={output.read_1} \
            out2={output.read_2} \
            threads={threads} \
            ref={params.ref} \
            2> {log}
        """


# 4. Check QC result - fastqc
rule fastqc:
    input:
        "resources/fastq/{sample}_qc_3rd_1.fq.gz",
        "resources/fastq/{sample}_qc_3rd_2.fq.gz"
    output:
        "results/fastqc/{sample}_qc_3rd_1_fastqc.zip",
        "results/fastqc/{sample}_qc_3rd_2_fastqc.zip"
    log:
        "logs/fastqc/{sample}_fastqc.log"
    threads: 10
    shell:
        "fastqc -o results/fastqc/ -t {threads} {input} 2> {log}"


# 5. Check QC result - multiqc
rule multiqc:
    input:
        "results/fastqc/"
    output:
        "results/fastqc/multiqc_report_filtered.html"
    shell:
        "multiqc {input} -n {output}"


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


# 6. Map reads to the reference genome

## genome version: GRCz11
## strandness: fr-firststrand

rule mapping:
    input:
        hisat2_index = config["ZF_genome_index"],
        read_1 = "resources/fastq/{sample}_qc_3rd_1.fq.gz", 
        read_2 = "resources/fastq/{sample}_qc_3rd_2.fq.gz"
    output:
        sam = temp("results/mapping/{sample}.sam")
    threads: 10
    log:
        "logs/hisat2/{sample}.log"
    shell:
        """
        hisat2 -x {input.hisat2_index} \
                -1 {input.read_1} \
                -2 {input.read_2} \
                --rna-strandness RF \
                --dta \
                -p {threads} \
                -S {output.sam} \
                2> {log}
        """


# 7. Sort reads by coordinates
rule sort:
    input:
        "results/mapping/{sample}.sam"
    output:
        temp("results/mapping/{sample}.sorted.bam")
    threads: 10
    shell:
        "samtools sort -O BAM -o {output} -@ {threads} {input}"


# 8. Add read groups
rule add_rg:
    input:
        "results/mapping/{sample}.sorted.bam"
    output:
        temp("results/mapping/{sample}.rg.bam")
    shell:
        "picard AddOrReplaceReadGroups \
                I={input} \
                O={output} \
                RGID={wildcards.sample} \
                RGLB=lib1 \
                RGPL=ILLUMINA \
                RGPU=unit1 \
                RGSM={wildcards.sample}"


# 9. Mark duplicates
rule mark_duplicate:
    input:
        "results/mapping/{sample}.rg.bam"
    output:
        bam = "results/mapping/{sample}.mkdup.bam",
        metrics = "logs/mark_duplicate/{sample}.mkdup.metrics"
    log:
        "logs/mark_duplicate/{sample}.mkdup.log"
    shell:
        "picard MarkDuplicates \
                ASSUME_SORTED=true \
                I={input} \
                O={output.bam} \
                M={output.metrics}"


# 10. Remove duplicates and secondary alignments
rule remove_dup:
    input:
        "results/mapping/{sample}.mkdup.bam"
    output:
        "results/mapping/{sample}.uniq_sorted.bam"
    threads: 10
    shell:
        "(samtools view -H {input}; samtools view -@ {threads} -F 3332 {input} | "
        "grep -w NH:i:1) | samtools view -@ {threads} -bS - > {output}"


# 11.1 Calculate metrics for all reads
rule stat_reads_before_revdup:
    input:
        "results/mapping/{sample}.mkdup.bam"
    output:
        "logs/reads_stats_before/{sample}.mkdup_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 11.2 Calculate metrics for clean reads
rule stat_reads_after_revdup:
    input:
        "results/mapping/{sample}.uniq_sorted.bam"
    output:
        "logs/reads_stats_after/{sample}.uniq_sorted_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 12. Build index for bam
rule index_bam:
    input:
        "results/mapping/{sample}.uniq_sorted.bam"
    output:
        "results/mapping/{sample}.uniq_sorted.bam.bai"
    threads: 10
    shell:
        "samtools index -@ {threads} {input} {output}"


# 13. Convert bam to bigwig
rule bam_to_bw:
    input:
        bam = "results/mapping/{sample}.uniq_sorted.bam",
        index = "results/mapping/{sample}.uniq_sorted.bam.bai"
    output:
        "results/mapping/{sample}_CPM.bw"
    threads: 20
    log:
        "logs/mapping/{sample}.convertbw.txt"
    shell:
        """
        bamCoverage --bam {input.bam} \
        -o {output} \
        --binSize 5 \
        --outFileFormat bigwig \
        --effectiveGenomeSize 1368780147 \
        --normalizeUsing CPM \
        -p {threads} \
        2> {log}
        """

# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# 14. Sample scaling
# 10zygotes: ZF_10zygotes_1 as the standard，single zygote: ZF_single_zygote_2 as the standard

# 14.1 Increase library size by merging
# directly run in shell, not by snakemake
'''
samtools merge -@ 10 -o ZF_10zygotes_2.cp.bam \
    ZF_10zygotes_2.uniq_sorted.bam \
    ZF_10zygotes_2.uniq_sorted.bam

samtools merge -@ 10 -o ZF_10zygotes_3.cp.bam \
    ZF_10zygotes_3.uniq_sorted.bam \
    ZF_10zygotes_3.uniq_sorted.bam

samtools merge -@ 10 -o ZF_10zygotes_4.cp.bam \
    ZF_10zygotes_4.uniq_sorted.bam \
    ZF_10zygotes_4.uniq_sorted.bam

samtools merge -@ 10 -o ZF_10zygotes_5.cp.bam \
    ZF_10zygotes_5.uniq_sorted.bam \
    ZF_10zygotes_5.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_1.cp.bam \
    ZF_single_zygote_1.uniq_sorted.bam \
    ZF_single_zygote_1.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_3.cp.bam \
    ZF_single_zygote_3.uniq_sorted.bam \
    ZF_single_zygote_3.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_4.cp.bam \
    ZF_single_zygote_4.uniq_sorted.bam \
    ZF_single_zygote_4.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_5.cp.bam \
    ZF_single_zygote_5.uniq_sorted.bam \
    ZF_single_zygote_5.uniq_sorted.bam \
    ZF_single_zygote_5.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_6.cp.bam \
    ZF_single_zygote_6.uniq_sorted.bam \
    ZF_single_zygote_6.uniq_sorted.bam \
    ZF_single_zygote_6.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_7.cp.bam \
    ZF_single_zygote_7.uniq_sorted.bam \
    ZF_single_zygote_7.uniq_sorted.bam \
    ZF_single_zygote_7.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_8.cp.bam \
    ZF_single_zygote_8.uniq_sorted.bam \
    ZF_single_zygote_8.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_9.cp.bam \
    ZF_single_zygote_9.uniq_sorted.bam \
    ZF_single_zygote_9.uniq_sorted.bam

samtools merge -@ 10 -o ZF_single_zygote_10.cp.bam \
    ZF_single_zygote_10.uniq_sorted.bam \
    ZF_single_zygote_10.uniq_sorted.bam


# 14.2 Downsample
# ZF_10zygotes_1 -> No downsampling
cp ZF_10zygotes_1.uniq_sorted.bam ZF_10zygotes_1.cp_ds.bam
# ZF_10zygotes_2  8,880,170/(5,161,874*2)=0.8606
samtools view -b --subsample-seed 1 --subsample 0.8606 ZF_10zygotes_2.cp.bam > ZF_10zygotes_2.cp_ds.bam
# ZF_10zygotes_3  8,880,170/(6,456,755*2)=0.6877
samtools view -b --subsample-seed 1 --subsample 0.6878 ZF_10zygotes_3.cp.bam > ZF_10zygotes_3.cp_ds.bam
# ZF_10zygotes_4  8,880,170/(7,921,503*2)=0.5605
samtools view -b --subsample-seed 1 --subsample 0.5605 ZF_10zygotes_4.cp.bam > ZF_10zygotes_4.cp_ds.bam
# ZF_10zygotes_5  8,880,170/(6,587,662*2)=0.6740
samtools view -b --subsample-seed 1 --subsample 0.6740 ZF_10zygotes_5.cp.bam > ZF_10zygotes_5.cp_ds.bam

# ZF_single_zygote_1  5,193,945/(4,404,088*2)=0.5897
samtools view -b --subsample-seed 1 --subsample 0.5897 ZF_single_zygote_1.cp.bam > ZF_single_zygote_1.cp_ds.bam
# ZF_single_zygote_2 -> No downsampling 
cp ZF_single_zygote_2.uniq_sorted.bam ZF_single_zygote_2.cp_ds.bam
# ZF_single_zygote_3  5,193,945/(3,978,336*2)=0.6526
samtools view -b --subsample-seed 1 --subsample 0.6526 ZF_single_zygote_3.cp.bam > ZF_single_zygote_3.cp_ds.bam
# ZF_single_zygote_4  5,193,945/(4,912,775*2)=0.5287
samtools view -b --subsample-seed 1 --subsample 0.5287 ZF_single_zygote_4.cp.bam > ZF_single_zygote_4.cp_ds.bam
# ZF_single_zygote_5 5,193,945/(2,210,208*3)=0.7837
samtools view -b --subsample-seed 1 --subsample 0.7837 ZF_single_zygote_5.cp.bam > ZF_single_zygote_5.cp_ds.bam
# ZF_single_zygote_6 5,193,945/(2,271,317*3)=0.7623
samtools view -b --subsample-seed 1 --subsample 0.7623 ZF_single_zygote_6.cp.bam > ZF_single_zygote_6.cp_ds.bam
# ZF_single_zygote_7 5,193,945/(2,244,127*3)=0.7715
samtools view -b --subsample-seed 1 --subsample 0.7715 ZF_single_zygote_7.cp.bam > ZF_single_zygote_7.cp_ds.bam
# ZF_single_zygote_8 5,193,945/(4,082,700*2)=0.6361
samtools view -b --subsample-seed 1 --subsample 0.6361 ZF_single_zygote_8.cp.bam > ZF_single_zygote_8.cp_ds.bam
# ZF_single_zygote_9 5,193,945/(2,610,619*2)=0.9948
samtools view -b --subsample-seed 1 --subsample 0.9948 ZF_single_zygote_9.cp.bam > ZF_single_zygote_9.cp_ds.bam
# ZF_single_zygote_10 5,193,945/(4,270,948*2)=0.6081
samtools view -b --subsample-seed 1 --subsample 0.6081 ZF_single_zygote_10.cp.bam > ZF_single_zygote_10.cp_ds.bam
'''


# 15. Separate reads mapped to different strands
rule separate_strand_reads:
    input:
        "results/mapping/{sample}.cp_ds.bam"
    output:
        forward_strand_1 = temp("results/mapping/{sample}.forward_strand.1.bam"),
        forward_strand_2 = temp("results/mapping/{sample}.forward_strand.2.bam"),
        reverse_strand_1 = temp("results/mapping/{sample}.reverse_strand.1.bam"),
        reverse_strand_2 = temp("results/mapping/{sample}.reverse_strand.2.bam")
    threads: 10
    shell:
        """
        samtools view -@ {threads} -b -f 83 {input} > {output.forward_strand_1}
        samtools view -@ {threads} -b -f 163 {input} > {output.forward_strand_2}
        samtools view -@ {threads} -b -f 99 {input} > {output.reverse_strand_1}
        samtools view -@ {threads} -b -f 147 {input} > {output.reverse_strand_2}
        """

rule merge_strand_reads:
    input:
        forward_strand_1 = "results/mapping/{sample}.forward_strand.1.bam",
        forward_strand_2 = "results/mapping/{sample}.forward_strand.2.bam",
        reverse_strand_1 = "results/mapping/{sample}.reverse_strand.1.bam",
        reverse_strand_2 = "results/mapping/{sample}.reverse_strand.2.bam"
    output:
        forward_strand = temp("results/mapping/{sample}.forward_strand.bam"),
        reverse_strand = temp("results/mapping/{sample}.reverse_strand.bam")
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output.forward_strand} {input.forward_strand_1} {input.forward_strand_2}
        samtools merge -@ {threads} {output.reverse_strand} {input.reverse_strand_1} {input.reverse_strand_2}
        """

rule sort_merged_strand_reads:
    input:
        forward_strand = "results/mapping/{sample}.forward_strand.bam",
        reverse_strand = "results/mapping/{sample}.reverse_strand.bam"
    output:
        forward_strand_sorted = "results/mapping/{sample}.forward_strand.sorted.bam",
        reverse_strand_sorted = "results/mapping/{sample}.reverse_strand.sorted.bam"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.forward_strand_sorted} {input.forward_strand}
        samtools sort -@ {threads} -o {output.reverse_strand_sorted} {input.reverse_strand}
        samtools index {output.forward_strand_sorted}
        samtools index {output.reverse_strand_sorted}
        """


# 16. Call peaks on separated strands
samtools merge -@ 10 ZF_Input.forward_strand.sorted.bam ZF_Input_1.forward_strand.sorted.bam ZF_Input_2.forward_strand.sorted.bam
samtools merge -@ 10 ZF_Input.reverse_strand.sorted.bam ZF_Input_1.reverse_strand.sorted.bam ZF_Input_2.reverse_strand.sorted.bam


## forward strand
rule call_peaks_exomePeak_forward_strand:
    input:
        gtf = config["ZF_genome_anno"],
        INPUT = "results/mapping/ZF_Input.forward_strand.sorted.bam",
        IP = "results/mapping/{sample}.forward_strand.sorted.bam"
    output:
        "results/peaks/{sample}_forward_strand/peaks.bed"
    log:
        "logs/exomePeak/{sample}_forward_strand.log"
    shell:
        """
        Rscript module/exomePeak2.R \
        {input.gtf} \
        {input.INPUT} \
        {input.IP} \
        {output} 2> {log}
        """

## reverse strand
rule call_peaks_exomePeak_reverse_strand:
    input:
        gtf = config["ZF_genome_anno"],
        INPUT = "results/mapping/ZF_Input.reverse_strand.sorted.bam",
        IP = "results/mapping/{sample}.reverse_strand.sorted.bam"
    output:
        "results/peaks/{sample}_reverse_strand/peaks.bed"
    log:
        "logs/exomePeak/{sample}_reverse_strand.log"
    shell:
        """
        Rscript module/exomePeak2.R \
        {input.gtf} \
        {input.INPUT} \
        {input.IP} \
        {output} 2> {log}
        """


# 17. Filter peaks using R
# >>> to R

# 18. Peak overlap analysis
# direct run in shell, not by snakemake
'''
intervene pairwise -i results/peaks/*_peak_filtered.bed --filenames --compute frac --htype tribar
# >>> to R
'''


# 19. Motif and distribution analyses for all peaks [Fig. 5c and 5d]
# direct run in shell, not by snakemake
'''
cut -f 1,2 resources/reference/Danio_rerio.GRCz11.dna_sm.primary_assembly.fa.fai > resources/reference/Danio_rerio.GRCz11.dna_sm.primary_assembly.chrom.size
cat resources/reference/Danio_rerio.GRCz11.115.gtf | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' > resources/reference/Danio_rerio.GRCz11.115_TranscriptOnly.bed
'''

# 19.1 Shuffle transcripts for the background of motif enrichment analysis
rule motif_bgshuffle:
    input:
        "results/peaks/{sample}_peak_filtered.bed"
    output:
        "results/peaks/{sample}_peak_filtered_shuffled.bed"
    params:
        genome_size = "resources/reference/Danio_rerio.GRCz11.dna_sm.primary_assembly.chrom.size",
        transcriptome_bed = "resources/reference/Danio_rerio.GRCz11.115_TranscriptOnly.bed"
    shell:
        """
        bedtools shuffle \
            -i {input} \
            -g {params.genome_size} \
            -incl {params.transcriptome_bed} \
            -excl {input} \
            -chrom \
            -seed 1 \
            > {output}
        """

# 19.2 Add "chr" to the chromosome name for HOMER
rule add_chr_to_bed_for_shuffled_peaks:
    input:
        "results/peaks/{sample}_peak_filtered_shuffled.bed"
    output:
        "results/peaks/{sample}_peak_filtered_shuffled_homer.bed"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """

rule add_chr_to_bed_for_peaks:
    input:
        "results/peaks/{sample}_peak_filtered.bed"
    output:
        "results/peaks/{sample}_peak_filtered_homer.bed"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """

# 19.3 Motif analysis [Fig. 5d]
rule motif_analysis:
    input:
        peak_file = "results/peaks/{sample}_peak_filtered_homer.bed",
        bg_file = "results/peaks/{sample}_peak_filtered_shuffled_homer.bed"
    output:
        directory("results/peaks/motif/{sample}_homer_bgshuffled")
    log:
        "logs/homer/{sample}_peaks_motif_bgshuffled.log"
    threads: 10
    shell:
        """
        findMotifsGenome.pl {input.peak_file} danRer11 {output} -bg {input.bg_file} -len 7,8 -rna -p {threads} 2> {log}
        """

# Motif visualization
# direct run in shell, not by snakemake
'''
python module/transform_pwm_transfac.py -i motif1.pwm -o motif1.transfac

weblogo < motif1.transfac > motif1.pdf \
        --format pdf \
        --datatype transfac \
        --color-scheme classic \
        --units bits \
        --first-index 1 \
        --errorbars False \
        --scale-width NO
'''


# 19.4 Peak distribution analysis [Fig. 5c]
rule annotate_peaks_homer:
    input:
        "results/peaks/{sample}_peak_filtered_homer.bed"
    output:
        all_out = "results/peaks/annotation/{sample}_peak_anno.txt",
        stats = "results/peaks/annotation/{sample}_peak_anno.stats.txt"
    log:
        "logs/homer/{sample}_peak_anno.log"
    shell:
        """
        annotatePeaks.pl {input} danRer11 -annStats {output.stats} > {output.all_out} 2> {log}
        """