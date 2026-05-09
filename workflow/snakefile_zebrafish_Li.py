configfile: "config.yaml"

SAMPLES = config['ZF_Li_single_IP']

rule all:
    input:
        expand("{sample}/peaks.bed", sample = SAMPLES)


# Li et al
# 1. Remove adapter
rule qc_cutadapt:
    params:
        min_length = 20
    input:
        read1 = lambda wildcards: config["ZF_Li_samples"][wildcards.sample]["r1"],
        read2 = lambda wildcards: config["ZF_Li_samples"][wildcards.sample]["r2"]
    output:
        read1 = temp("other_study/Li/resources/{sample}_qc_1.fq.gz"),
        read2 = temp("other_study/Li/resources/{sample}_qc_2.fq.gz")
    threads: 10
    log:
        "other_study/Li/logs/cutadapt/{sample}_cutadapt.log"
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
        read1 = "other_study/Li/resources/{sample}_qc_1.fq.gz",
        read2 = "other_study/Li/resources/{sample}_qc_2.fq.gz"
    output:
        read1 = temp("other_study/Li/resources/{sample}_qc_2nd_1.fq.gz"),
        read2 = temp("other_study/Li/resources/{sample}_qc_2nd_2.fq.gz")
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
            > {log} 2>&1
        """


# 3. Remove reads belonging to rRNA, tRNA sequences
rule remove_contaminants:
    input:
        read_1 = "other_study/Li/resources/{sample}_qc_2nd_1.fq.gz",
        read_2 = "other_study/Li/resources/{sample}_qc_2nd_2.fq.gz",
    output:
        read_1 = "other_study/Li/resources/{sample}_qc_3rd_1.fq.gz",
        read_2 = "other_study/Li/resources/{sample}_qc_3rd_2.fq.gz",
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
        "other_study/Li/resources/{sample}_qc_3rd_1.fq.gz",
        "other_study/Li/resources/{sample}_qc_3rd_2.fq.gz"
    output:
        "other_study/Li/results/fastqc/{sample}_qc_3rd_1_fastqc.zip",
        "other_study/Li/results/fastqc/{sample}_qc_3rd_2_fastqc.zip"
    log:
        "other_study/Li/logs/fastqc/{sample}_fastqc.log"
    threads: 10
    shell:
        "fastqc -o other_study/Li/results/fastqc/ -t {threads} {input} 2> {log}"


# 5. Check QC result - multiqc
rule multiqc:
    input:
        "other_study/Li/results/fastqc/"
    output:
        "other_study/Li/results/fastqc/multiqc_report_filtered.html"
    shell:
        "multiqc {input} -n {output}"


# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------


# 6. Map reads to the reference genome

## genome version: GRCz11
## strandness: unstranded

rule mapping:
    input:
        hisat2_index = config["ZF_genome_index"],
        read_1 = "other_study/Li/resources/{sample}_qc_3rd_1.fq.gz", 
        read_2 = "other_study/Li/resources/{sample}_qc_3rd_2.fq.gz"
    output:
        sam = temp("other_study/Li/results/mapping/{sample}.sam")
    threads: 10
    log:
        "logs/hisat2/{sample}.log"
    shell:
        """
        hisat2 -x {input.hisat2_index}/ \
                -1 {input.read_1} \
                -2 {input.read_2} \
                --dta \
                -p {threads} \
                -S {output.sam} \
                2> {log}
        """


# 7. Sort reads by coordinates
rule sort:
    input:
        "other_study/Li/results/mapping/{sample}.sam"
    output:
        temp("other_study/Li/results/mapping/{sample}.sorted.bam")
    threads: 10
    shell:
        "samtools sort -O BAM -o {output} -@ {threads} {input}"


# 8. Add read groups
rule add_rg:
    input:
        "other_study/Li/results/mapping/{sample}.sorted.bam"
    output:
        temp("other_study/Li/results/mapping/{sample}.rg.bam")
    log:
        "logs/add_rg/{sample}.add_rg.log"
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
        "other_study/Li/results/mapping/{sample}.rg.bam"
    output:
        bam = "other_study/Li/results/mapping/{sample}.mkdup.bam",
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
        "other_study/Li/results/mapping/{sample}.mkdup.bam"
    output:
        "other_study/Li/results/mapping/{sample}.uniq_sorted.bam"
    threads: 10
    shell:
        "(samtools view -H {input}; samtools view -@ {threads} -F 3332 {input} | "
        "grep -w NH:i:1) | samtools view -@ {threads} -bS - > {output}"


# 11.1 Calculate metrics for all reads
rule stat_reads_before_revdup:
    input:
        "other_study/Li/results/mapping/{sample}.mkdup.bam"
    output:
        "logs/reads_stats_before/{sample}.mkdup_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 11.2 Calculate metrics for clean reads
rule stat_reads_after_revdup:
    input:
        "other_study/Li/results/mapping/{sample}.uniq_sorted.bam"
    output:
        "logs/reads_stats_after/{sample}.uniq_sorted_metrics.log"
    threads: 10
    shell:
        "samtools flagstat -@ {threads} {input} > {output}"


# 12. Build index for bam
rule index_bam:
    input:
        "other_study/Li/results/mapping/{sample}.uniq_sorted.bam"
    output:
        "other_study/Li/results/mapping/{sample}.uniq_sorted.bam.bai"
    threads: 10
    shell:
        "samtools index -@ {threads} {input} {output}"


# 13. Convert bam to bigwig
rule bam_to_bw:
    input:
        bam = "other_study/Li/results/mapping/{sample}.uniq_sorted.bam",
        index = "other_study/Li/results/mapping/{sample}.uniq_sorted.bam.bai"
    output:
        "other_study/Li/results/mapping/{sample}_CPM.bw"
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


# 14. Sample scaling
# single zygote: Li_zIP1_5 as the standard
'''
samtools stats Li_zIP1_1.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 370306
samtools stats Li_zIP1_2.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 528227
samtools stats Li_zIP1_3.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 637921
samtools stats Li_zIP1_4.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 1016270
samtools stats Li_zIP1_5.uniq_sorted.bam | grep '^SN' | grep 'reads mapped:' # 1077146


# 14.1. Increase library size by merging
# directly run in shell, not by snakemake
samtools merge -@ 10 Li_zInput10.uniq_sorted.bam Li_zInput10_1.uniq_sorted.bam Li_zInput10_2.uniq_sorted.bam

samtools merge -@ 10 -o Li_zIP1_1.cp.bam \
    Li_zIP1_1.uniq_sorted.bam \
    Li_zIP1_1.uniq_sorted.bam \
    Li_zIP1_1.uniq_sorted.bam

samtools merge -@ 10 -o Li_zIP1_2.cp.bam \
    Li_zIP1_2.uniq_sorted.bam \
    Li_zIP1_2.uniq_sorted.bam \
    Li_zIP1_2.uniq_sorted.bam

samtools merge -@ 10 -o Li_zIP1_3.cp.bam \
    Li_zIP1_3.uniq_sorted.bam \
    Li_zIP1_3.uniq_sorted.bam

samtools merge -@ 10 -o Li_zIP1_4.cp.bam \
    Li_zIP1_4.uniq_sorted.bam \
    Li_zIP1_4.uniq_sorted.bam


# 14.2. Downsample
# Li_zIP1_1  1077146/(370306*3)=0.9708
samtools view -b --subsample-seed 1 --subsample 0.9708 Li_zIP1_1.cp.bam > Li_zIP1_1.cp_ds.bam
# Li_zIP1_2 1077146/(528227*3)=0.6801
samtools view -b --subsample-seed 1 --subsample 0.6801 Li_zIP1_2.cp.bam > Li_zIP1_2.cp_ds.bam
# Li_zIP1_3 1077146/(637921*2)=0.8449
samtools view -b --subsample-seed 1 --subsample 0.8449 Li_zIP1_3.cp.bam > Li_zIP1_3.cp_ds.bam
# Li_zIP1_4 1077146/(1016270*2)=0.5308
samtools view -b --subsample-seed 1 --subsample 0.5308 Li_zIP1_4.cp.bam > Li_zIP1_4.cp_ds.bam
# Li_zIP1_5 -> No downsampling
cp Li_zIP1_5.uniq_sorted.bam Li_zIP1_5.cp_ds.bam


samtools stats Li_zIP1_1.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 1077867
samtools stats Li_zIP1_2.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 1079277
samtools stats Li_zIP1_3.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 1077300
samtools stats Li_zIP1_4.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 1078030
samtools stats Li_zIP1_5.cp_ds.bam | grep '^SN' | grep 'reads mapped:' # 1077146
'''


# 15. Call peaks
rule call_peaks_exomePeak:
    input:
        gtf = config["ZF_genome_anno"],
        INPUT = "other_study/Li/results/mapping/Li_zInput10.uniq_sorted.bam",
        IP = "other_study/Li/results/mapping/{sample}.cp_ds.bam"
    output:
        "results/peaks/{sample}/peaks.bed"
    log:
        "logs/exomePeak/{sample}.log"
    shell:
        """
        Rscript module/exomePeak2.R \
        {input.gtf} \
        {input.INPUT} \
        {input.IP} \
        {output} 2> {log}
        """