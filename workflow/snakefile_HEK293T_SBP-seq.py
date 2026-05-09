configfile: "config.yaml"


# 1. Merge replicates
rule merge_fastq:
    input:
        r1 = lambda wildcards: [config["HEK293T_samples"][sample_i]["r1"] for sample_i in config["merged_groups"][wildcards.sample]],
        r2 = lambda wildcards: [config["HEK293T_samples"][sample_i]["r2"] for sample_i in config["merged_groups"][wildcards.sample]]
    output:
        r1 = "results/merged_fastq/{sample}_1.fastq.gz",
        r2 = "results/merged_fastq/{sample}_2.fastq.gz"
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """


# 2. Remove adapter
rule qc_cutadapt:
    params:
        min_length = 20
    input:
        read1 = "results/merged_fastq/{sample}_1.fastq.gz",
        read2 = "results/merged_fastq/{sample}_2.fastq.gz"
    output:
        read1 = temp("results/merged_fastq/{sample}_qc_1.fastq.gz"),
        read2 = temp("results/merged_fastq/{sample}_qc_2.fastq.gz")
    threads: 6
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


# 3. Remove polyG and low-quality bases
rule qc_fastp:
    input:
        read1 = "results/merged_fastq/{sample}_qc_1.fastq.gz",
        read2 = "results/merged_fastq/{sample}_qc_2.fastq.gz"
    output:
        read1 = temp("results/merged_fastq/{sample}_qc_2nd_1.fastq.gz"),
        read2 = temp("results/merged_fastq/{sample}_qc_2nd_2.fastq.gz")
    threads: 6
    log:
        "logs/fastp/{sample}.log"
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
            --thread {threads} \
            > {log}
        """


# 4. Remove reads belonging to rRNA, tRNA sequences
rule remove_contaminants:
    input:
        read_1 = "results/merged_fastq/{sample}_qc_2nd_1.fastq.gz",
        read_2 = "results/merged_fastq/{sample}_qc_2nd_2.fastq.gz",
    output:
        read_1 = "results/merged_fastq/{sample}_qc_3rd_1.fq.gz",
        read_2 = "results/merged_fastq/{sample}_qc_3rd_2.fq.gz"
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


# 5. Check QC result - fastqc
rule fastqc:
    input:
        "results/merged_fastq/{sample}_qc_3rd_1.fq.gz",
        "results/merged_fastq/{sample}_qc_3rd_2.fq.gz"
    output:
        "results/fastqc/{sample}_qc_3rd_1_fastqc.zip",
        "results/fastqc/{sample}_qc_3rd_2_fastqc.zip"
    log:
        "logs/fastqc/{sample}_fastqc.log"
    threads: 4
    shell:
        "fastqc -o results/fastqc -t {threads} {input} 2> {log}"


# 6. Check QC result - multiqc
rule multiqc:
    input:
        "results/fastqc/"
    output:
        "results/fastqc/multiqc_report_filtered.html"
    shell:
        "multiqc {input} -n {output}"


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# 7. Map reads to the reference genome 

## genome version: GRCh38
## strandness: fr-firststrand

rule mapping:
    input:
        hisat2_index = config["human_genome_index"],
        read_1 = "results/merged_fastq/{sample}_qc_3rd_1.fq.gz", 
        read_2 = "results/merged_fastq/{sample}_qc_3rd_2.fq.gz"
    output:
        temp("results/mapping/{sample}.sam")
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
                -S {output} \
                2> {log}
        """


# 8. Sort reads by coordinates
rule sort:
    input:
        "results/mapping/{sample}.sam"
    output:
        temp("results/mapping/{sample}.sorted.bam")
    threads: 10
    shell:
        "samtools sort -O BAM -o {output} -@ {threads} {input}"


# 9. Add read groups
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


# 10. Mark duplicates
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
                M={output.metrics} \
                2> {log}"


# 11. Remove duplicates and secondary alignments
rule remove_dup:
    input:
        "results/mapping/{sample}.mkdup.bam"
    output:
        "results/mapping/{sample}.uniq_sorted.bam"
    threads: 10
    shell:
        "(samtools view -H {input}; samtools view -@ {threads} -F 3332 {input} | "
        "grep -w NH:i:1) | samtools view -@ {threads} -bS - > {output}"


# 12.1 Calculate metrics for all reads
rule stat_reads_before_revdup:
    input:
        "results/mapping/{sample}.mkdup.bam"
    output:
        "logs/reads_stats_before/{sample}.mkdup_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 12.2 calculate metrics for clean reads
rule stat_reads_after_revdup:
    input:
        "results/mapping/{sample}.uniq_sorted.bam"
    output:
        "logs/reads_stats_after/{sample}.uniq_sorted_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 13. Build index for bam
rule index_bam:
    input:
        "results/mapping/{sample}.uniq_sorted.bam"
    output:
        "results/mapping/{sample}.uniq_sorted.bam.bai"
    threads: 10
    shell:
        "samtools index -@ {threads} {input} {output}"


# 14. Convert bam to bigwig
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
        --effectiveGenomeSize 2913022398 \
        --normalizeUsing CPM \
        -p {threads} \
        2> {log}
        """



# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# 15. Sample scaling
# 15.1 Merge bam files of technical replicates
rule merge_technical_replicates:
    input:
        rep1 = "results/mapping/{sample}-1_merged.uniq_sorted.bam",
        rep2 = "results/mapping/{sample}-2_merged.uniq_sorted.bam"
    output:
        "results/mapping/rep_merged/{sample}.uniq_sorted.bam"
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output} {input.rep1} {input.rep2}
        samtools index {output}
        """


# 15.2 Increase library size by merging
# directly run in shell, not by snakemake
'''
samtools merge -@ 10 -o results/mapping/rep_merged/HEK293T_SYSY.cp.bam \
    results/mapping/rep_merged/HEK293T_SYSY.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_SYSY.uniq_sorted.bam

samtools merge -@ 10 -o results/mapping/rep_merged/HEK293T_Active_Motif.cp.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam

samtools merge -@ 10 -o results/mapping/rep_merged/HEK293T_Millipore.cp.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam \
    results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam


samtools stats results/mapping/rep_merged/HEK293T_Em6ABP1.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 32863952
samtools stats results/mapping/rep_merged/HEK293T_SYSY.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 17581379
samtools stats results/mapping/rep_merged/HEK293T_Active_Motif.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 5851884
samtools stats results/mapping/rep_merged/HEK293T_Millipore.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 5359138


# 15.3 Downsample
# Em6ABP1 -> No downsampling
cp HEK293T_Em6ABP1.uniq_sorted.bam HEK293T_Em6ABP1.sorted.cp_ds.bam
cp HEK293T_Em6ABP1.uniq_sorted.bam.bai HEK293T_Em6ABP1.sorted.cp_ds.bam.bai
# IgG -> No downsampling
cp HEK293T_Input.uniq_sorted.bam HEK293T_Input.sorted.cp_ds.bam
cp HEK293T_Input.uniq_sorted.bam.bai HEK293T_Input.sorted.cp_ds.bam.bai
# AM (0.9360)
samtools view -b --subsample-seed 1 --subsample 0.9360 HEK293T_Active_Motif.cp.bam > HEK293T_Active_Motif.cp_ds.bam
# Milli (0.8760)
samtools view -b --subsample-seed 1 --subsample 0.8760 HEK293T_Millipore.cp.bam > 293-mIP-Millli.cp_ds.bam
# SY (0.9346)
samtools view -b --subsample-seed 1 --subsample 0.9346 HEK293T_SYSY.cp.bam > HEK293T_SYSY.cp_ds.bam


samtools sort -O BAM -o HEK293T_Input.sorted.cp_ds.bam -@ 50 HEK293T_Input.uniq_sorted.bam
samtools sort -O BAM -o HEK293T_Active_Motif.sorted.cp_ds.bam -@ 50 HEK293T_Active_Motif.cp_ds.bam
samtools sort -O BAM -o HEK293T_Millipore.sorted.cp_ds.bam -@ 50 293-mIP-Millli.cp_ds.bam
samtools sort -O BAM -o HEK293T_SYSY.sorted.cp_ds.bam -@ 50 HEK293T_SYSY.cp_ds.bam

samtools index HEK293T_Active_Motif.sorted.cp_ds.bam
samtools index HEK293T_Millipore.sorted.cp_ds.bam
samtools index HEK293T_SYSY.sorted.cp_ds.bam

samtools stats results/mapping/rep_merged/HEK293T_Em6ABP1.sorted.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 32863880
samtools stats results/mapping/rep_merged/HEK293T_SYSY.sorted.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 32861224
samtools stats results/mapping/rep_merged/HEK293T_Active_Motif.sorted.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 32872062
samtools stats results/mapping/rep_merged/HEK293T_Millipore.sorted.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 32853401
'''


# 16. Separate reads mapped to different strands
rule separate_strand_reads:
    input:
        "results/mapping/rep_merged/{sample}.sorted.cp_ds.bam"
    output:
        forward_strand_1 = temp("results/mapping/rep_merged/{sample}.forward_strand.1.bam"),
        forward_strand_2 = temp("results/mapping/rep_merged/{sample}.forward_strand.2.bam"),
        reverse_strand_1 = temp("results/mapping/rep_merged/{sample}.reverse_strand.1.bam"),
        reverse_strand_2 = temp("results/mapping/rep_merged/{sample}.reverse_strand.2.bam")
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
        forward_strand_1 = "results/mapping/rep_merged/{sample}.forward_strand.1.bam",
        forward_strand_2 = "results/mapping/rep_merged/{sample}.forward_strand.2.bam",
        reverse_strand_1 = "results/mapping/rep_merged/{sample}.reverse_strand.1.bam",
        reverse_strand_2 = "results/mapping/rep_merged/{sample}.reverse_strand.2.bam"
    output:
        forward_strand = temp("results/mapping/rep_merged/{sample}.forward_strand.bam"),
        reverse_strand = temp("results/mapping/rep_merged/{sample}.reverse_strand.bam")
    threads: 10
    shell:
        """
        samtools merge -@ {threads} {output.forward_strand} {input.forward_strand_1} {input.forward_strand_2}
        samtools merge -@ {threads} {output.reverse_strand} {input.reverse_strand_1} {input.reverse_strand_2}
        """

rule sort_merged_strand_reads:
    input:
        forward_strand = "results/mapping/rep_merged/{sample}.forward_strand.bam",
        reverse_strand = "results/mapping/rep_merged/{sample}.reverse_strand.bam"
    output:
        forward_strand_sorted = "results/mapping/rep_merged/{sample}.forward_strand.sorted.bam",
        reverse_strand_sorted = "results/mapping/rep_merged/{sample}.reverse_strand.sorted.bam"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} -o {output.forward_strand_sorted} {input.forward_strand}
        samtools sort -@ {threads} -o {output.reverse_strand_sorted} {input.reverse_strand}
        samtools index {output.forward_strand_sorted}
        samtools index {output.reverse_strand_sorted}
        """


# 17. Call peaks on separated strands
## forward strand
rule call_peaks_forward:
    input:
        gtf = config["human_genome_anno"],
        INPUT = "results/mapping/rep_merged/HEK293T_Input.forward_strand.sorted.bam",
        IP = "results/mapping/rep_merged/{sample}.forward_strand.sorted.bam"
    output:
        "results/peaks/{sample}_forward_strand/peaks.bed"
    log:
        "logs/exomePeak2/{sample}_forward_strand.log"
    shell:
        """
        Rscript module/exomePeak2.R \
        {input.gtf} \
        {input.INPUT} \
        {input.IP} \
        {output} 2> {log}
        """

## reverse strand
rule call_peaks_reverse:
    input:
        gtf = config["human_genome_anno"],
        INPUT = "results/mapping/rep_merged/HEK293T_Input.reverse_strand.sorted.bam",
        IP = "results/mapping/rep_merged/{sample}.reverse_strand.sorted.bam"
    output:
         "results/peaks/{sample}_reverse_strand/peaks.bed"
    log:
        "logs/exomePeak2/{sample}_reverse_strand.log"
    shell:
        """
        Rscript module/exomePeak2.R \
        {input.gtf} \
        {input.INPUT} \
        {input.IP} \
        {output} 2> {log}
        """


# 18. Filter peaks using R
# >>> to R


# 19. Intersect peaks
# directly run in shell, not by snakemake
'''
# Em6ABP1 peaks
bedtools intersect -a results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    -b results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names SY AM Milli IgG | \
    awk -F"\t" 'BEGIN{OFS="\t"} {for(i=1; i<=25 && i<=NF; i++) {
        printf "%s", $i
        if(i < 25 && i < NF) {printf "%s", OFS}} printf "\n"}'  \
    > results/peaks/intersected/Em6ABP1_peaks_intersected.bed

# SY peaks
bedtools intersect -a results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 AM Milli IgG | \
    awk -F"\t" 'BEGIN{OFS="\t"} {for(i=1; i<=25 && i<=NF; i++) {
        printf "%s", $i
        if(i < 25 && i < NF) {printf "%s", OFS}} printf "\n"}'  \
    > results/peaks/intersected/SY_peaks_intersected.bed

# AM peaks
bedtools intersect -a results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY Milli IgG | \
    awk -F"\t" 'BEGIN{OFS="\t"} {for(i=1; i<=25 && i<=NF; i++) {
        printf "%s", $i
        if(i < 25 && i < NF) {printf "%s", OFS}} printf "\n"}'  \
    > results/peaks/intersected/AM_peaks_intersected.bed

# Milli peaks
bedtools intersect -a results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM IgG | \
    awk -F"\t" 'BEGIN{OFS="\t"} {for(i=1; i<=25 && i<=NF; i++) {
        printf "%s", $i
        if(i < 25 && i < NF) {printf "%s", OFS}} printf "\n"}'  \
    > results/peaks/intersected/Milli_peaks_intersected.bed

# IgG peaks
bedtools intersect -a results/peaks/HEK293T_Input/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM Milli | \
    awk -F"\t" 'BEGIN{OFS="\t"} {for(i=1; i<=25 && i<=NF; i++) {
        printf "%s", $i
        if(i < 25 && i < NF) {printf "%s", OFS}} printf "\n"}'  \
    > results/peaks/intersected/IgG_peaks_intersected.bed
'''


# 20. Motif analysis, annotation for all peaks [Fig. 4a, Supp. Fig. 8b]
# direct run in shell, not by snakemake
'''
cut -f 1,2 resources/reference/GRCh38.fa.fai > resources/reference/GRCh38.chrom.size
cat resources/reference/Homo_sapiens.GRCh38.95.gtf | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d '";' > resources/reference/Homo_sapiens.GRCh38.95_TranscriptOnly.bed
'''

# 20.1 Shuffle transcripts for the background of motif enrichment analysis
rule motif_bgshuffle_all:
    input:
        "results/peaks/{sample}/peaks_filtered.bed"
    output:
        "results/peaks/{sample}/peaks_filtered_shuffled.bed"
    params:
        genome_size = "resources/reference/GRCh38.chrom.size",
        transcriptome_bed = "resources/reference/Homo_sapiens.GRCh38.95_TranscriptOnly.bed"
    shell:
        """
        bedtools shuffle \
            -i {input} \
            -g {params.genome_size} \
            -incl {params.transcriptome_bed} \
            -excl {input} \
            -chrom \
            -seed 100 \
            > {output}
        """


# 20.2 Add "chr" to the chromosome name for HOMER
rule add_chr_to_bed_for_all_peaks:
    input:
        "results/peaks/{sample}/peaks_filtered_shuffled.bed"
    output:
        "results/peaks/{sample}/peaks_filtered_shuffled_homer.bed"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """

rule add_chr_to_bed_for_all_peaks2:
    input:
        "results/peaks/{sample}/peaks_filtered.bed"
    output:
        "results/peaks/{sample}/peaks_filtered_homer.bed"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """


# 20.3 Annotate all peaks
rule annotate_all_peaks:
    input:
        "results/peaks/{sample}/peak_filtered_homer.bed"
    output:
        all_out = "results/peaks/{sample}/peak_filtered_anno.txt",
        stats = "results/peaks/{sample}/peak_filtered_anno.stats.txt"
    log:
        "logs/homer/{sample}_all_peak_anno.log"
    shell:
        """
        annotatePeaks.pl {input} hg38 -annStats {output.stats} > {output.all_out} 2> {log}
        """

# 20.4 Motif analysis
rule motif_analysis_all:
    input:
        peak_file = "results/peaks/{sample}/peaks_filtered_homer.bed",
        bg_file = "results/peaks/{sample}/peaks_filtered_shuffled_homer.bed"
    output:
        directory("results/peaks/motif/{sample}_all_homer_bgshuffled")
    log:
        "logs/homer/{sample}_all_peaks_motif_bgshuffled.log"
    threads: 10
    shell:
        """
        findMotifsGenome.pl {input.peak_file} hg38 {output} \
            -bg {input.bg_file} \
            -len 7,8 -rna -p {threads} 2> {log}
        """


# 21. Motif and distribution analyses for specific peaks [Fig. 4f]
# 21.1 Shuffle transcripts for the background of motif enrichment analysis
rule motif_bgshuffle_sp:
    input:
        "results/peaks/intersected/{sample}_peak_specific.bed12"
    output:
        "results/peaks/intersected/{sample}_peak_specific_shuffled.bed12"
    params:
        genome_size = "resources/reference/GRCh38.chrom.size",
        transcriptome_bed = "resources/reference/Homo_sapiens.GRCh38.95_TranscriptOnly.bed"
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


# 21.2 Add "chr" to the chromosome name for HOMER
rule add_chr_to_bed_for_sp_peaks:
    input:
        "results/peaks/intersected/{sample}_peak_specific_shuffled.bed12"
    output:
        "results/peaks/intersected/{sample}_peak_specific_shuffled_homer.bed12"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """

rule add_chr_to_bed_for_sp_peaks2:
    input:
        "results/peaks/intersected/{sample}_peak_specific.bed12"
    output:
        "results/peaks/intersected/{sample}_peak_specific_homer.bed12"
    shell:
        """
        python3 module/add_chr_for_homer.py {input} {output}
        """


# 21.3 Annotate specific peaks
rule annotate_sp_peaks:
    input:
        "results/peaks/intersected/{sample}_peak_specific_homer.bed12"
    output:
        all_out = "results/peaks/intersected/{sample}_specific_peak_anno.txt",
        stats = "results/peaks/intersected/{sample}_specific_peak_anno.stats.txt"
    log:
        "logs/homer/{sample}_specific_peak_anno.log"
    shell:
        """
        annotatePeaks.pl {input} hg38 -annStats {output.stats} > {output.all_out} 2> {log}
        """


# 21.4 Motif analysis
rule motif_analysis_sp:
    input:
        peak_file = "results/peaks/intersected/{sample}_peak_specific_homer.bed12",
        bg_file = "results/peaks/intersected/{sample}_peak_specific_shuffled_homer.bed12"
    output:
        directory("results/peaks/motif/{sample}_sp_homer_bgshuffled")
    log:
        "logs/homer/{sample}_sp_peaks_motif_bgshuffled.log"
    threads: 10
    shell:
        """
        findMotifsGenome.pl {input.peak_file} hg38 {output} \
            -bg {input.bg_file} \
            -len 7,8 -rna -p {threads} 2> {log}
        """


# 21.5 Motif visualization
# direct run in shell, not by snakemake
'''
python module/transform_pwm_transfac.py -i motif1.motif -o motif1.transfac

weblogo < motif1.transfac > \
        results/peaks/motif/motifname.pdf \
        --format pdf \
        --datatype transfac \
        --color-scheme classic \
        --units bits \
        --first-index 1 \
        --errorbars False \
        --scale-width NO
'''



# 22.1 Intersect with GLORI data - GLORI m6A site level [Supp. Fig. 9a]
# direct run in shell, not by snakemake
'''
## 0.25 < RPM <= 0.5
bedtools intersect -a results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM Milli IgG | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/GLORI_others_intersected_0.25_0.5.bed

bedtools intersect -a results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    -b results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed \
    -s -loj -wao -header -names SY AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Em6ABP1_others_intersected_0.25_0.5.bed


bedtools intersect -a results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed \
    -s -loj -wao -header -names Em6ABP1 AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/SY_others_intersected_0.25_0.5.bed


bedtools intersect -a results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed \
    -s -loj -wao -header -names Em6ABP1 SY Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/AM_others_intersected_0.25_0.5.bed


bedtools intersect -a results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.25_0.5.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Milli_others_intersected_0.25_0.5.bed


## 0.5 < RPM <= 1
bedtools intersect -a results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM Milli IgG | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/GLORI_others_intersected_0.5_1.bed

bedtools intersect -a results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    -b results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed \
    -s -loj -wao -header -names SY AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Em6ABP1_others_intersected_0.5_1.bed


bedtools intersect -a results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed \
    -s -loj -wao -header -names Em6ABP1 AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/SY_others_intersected_0.5_1.bed


bedtools intersect -a results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed \
    -s -loj -wao -header -names Em6ABP1 SY Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/AM_others_intersected_0.5_1.bed


bedtools intersect -a results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_0.5_1.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Milli_others_intersected_0.5_1.bed


## RPM > 1
bedtools intersect -a results/peaks/intersect_with_GLORI/GLORI_slim_1.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM Milli IgG | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/GLORI_others_intersected_1.bed

bedtools intersect -a results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    -b results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_1.bed \
    -s -loj -wao -header -names SY AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Em6ABP1_others_intersected_1.bed

bedtools intersect -a results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_1.bed \
    -s -loj -wao -header -names Em6ABP1 AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/SY_others_intersected_1.bed

bedtools intersect -a results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_1.bed \
    -s -loj -wao -header -names Em6ABP1 SY Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/AM_others_intersected_1.bed

bedtools intersect -a results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim_1.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Milli_others_intersected_1.bed


# 24.2 Intersect with GLORI data - peak level [Supp. Fig. 9b]
bedtools intersect -a results/peaks/intersect_with_GLORI/GLORI_slim.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM Milli IgG | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/GLORI_others_intersected.bed


bedtools intersect -a results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    -b results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim.bed \
    -s -loj -wao -header -names SY AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Em6ABP1_others_intersected.bed


bedtools intersect -a results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim.bed \
    -s -loj -wao -header -names Em6ABP1 AM Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/SYSY_others_intersected.bed


bedtools intersect -a results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim.bed \
    -s -loj -wao -header -names Em6ABP1 SY Milli IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Active_Motif_others_intersected.bed


bedtools intersect -a results/peaks/HEK293T_Millipore/peaks_filtered.bed \
    -b results/peaks/HEK293T_Em6ABP1/peaks_filtered.bed \
    results/peaks/HEK293T_SYSY/peaks_filtered.bed \
    results/peaks/HEK293T_Active_Motif/peaks_filtered.bed \
    results/peaks/HEK293T_Input/peaks_filtered.bed \
    results/peaks/intersect_with_GLORI/GLORI_slim.bed \
    -s -loj -wao -header -names Em6ABP1 SY AM IgG GLORI | \
    cut -f1-19 > results/peaks/intersect_with_GLORI/Millipore_others_intersected.bed
'''
