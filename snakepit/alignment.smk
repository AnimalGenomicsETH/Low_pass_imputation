rule all:
    input:
        multiqc = multiQC + "multiqc_report.html",
        stats_1 = expand(sorted_alignment + "{sample}/{sample}.stats",sample=sample),
        stats_2 = expand(dedup_alignment + "{sample}/{sample}.stats",sample=sample),
        metrics = expand(metrics + "{sample}/{sample}.alignment_summary_metrics",sample=sample),
        idxstats = expand(metrics + "{sample}/{sample}_idxstats",sample=sample),
        sam_flags = expand(metrics + "sam_flag_logs/{sample}.log",sample=sample),
        coverage = coverage + "summary_coverage.txt"

rule QC:
    input:
        R1 = raw_data + "{sample}_R1.fastq.gz",
        R2 = raw_data + "{sample}_R2.fastq.gz"
    output:
        R1 = QC + "{sample}/{sample}_R1.fastq.gz",
        R2 = QC + "{sample}/{sample}_R2.fastq.gz",
        R_HTML = QC + "{sample}/{sample}_fastp.html",
        R_JSON = QC + "{sample}/{sample}_fastp.json"
    params:
        "-q 15 -u 40 -g"
    shell:
        FASTP + " -i {input.R1} -o {output.R1} -I {input.R2} -O {output.R2} -h {output.R_HTML} -j {output.R_JSON} {params} >/dev/null"

rule multiqc:
    input:
        fastp_json = expand(rules.QC.output.R_JSON,sample=sample)
    output:
        multiQC + "multiqc_report.html"
    params:
        json = "-k json"
    shell:
        MULTIQC + " {params.json} " + QC + "* -o" + multiQC

checkpoint split_fastq:
    input:
        R1 = rules.QC.output.R1,
        R2 = rules.QC.output.R2
    output:
        split_fastqs = temp(directory(split_fastq + "{sample}"))
    params:
        prefix = split_fastq + "{sample}/{sample}_",
        folder = split_fastq + "{sample}"
    shell:
        "mkdir {params.folder} \n" +
        FASTQ_SPLITTER + " -o {params.prefix} {input.R1} {input.R2}"

rule alignment:
    input:
        R1 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R1.fq.gz",
        R2 = split_fastq + "{sample}/{sample}_{flowcell}_{lane}_R2.fq.gz"
    output:
        BAM = temp(alignment + "{sample}/{sample}_{flowcell}_{lane}.bam")
    params:
        rg = "@RG\\tID:{flowcell}.{lane}\\tCN:ETH\\tLB:{sample}\\tPL:illumina\\tPU:{flowcell}:{lane}\\tSM:{sample}",
        bwa_mem = "-M -t 12 -R",
        ref = ref
    shell:
        BWA_MEM + " mem {params.bwa_mem} '{params.rg}' {params.ref} {input.R1} {input.R2} | " + SAMBLASTER + " -M | " + SAMTOOLS + " view -Sb - > {output.BAM}"

rule sort:
    input:
        BAM = rules.alignment.output.BAM
    output:
        sorted_BAM = temp(sorted_alignment + "{sample}/{sample}_{flowcell}_{lane}.bam")
    shell:
        SAMBAMBA + " sort -t 6 --out {output.sorted_BAM} {input.BAM}"

def list_files(wildcards):
    checkpoint_output = checkpoints.split_fastq.get(sample = wildcards.sample).output.split_fastqs
    all_wildcards = glob_wildcards(os.path.join(checkpoint_output, "{sample}_{flowcell}_{lane}_R1.fq.gz"))
    all_files = []
    for sample, flowcell, lane in zip(all_wildcards.sample, all_wildcards.flowcell, all_wildcards.lane):
        all_files.append(f"{sorted_alignment}" + f"{sample}/{sample}_{flowcell}_{lane}.bam")
    return(all_files)

rule merge:
    input: list_files
    output: temp(sorted_alignment + "{sample}/{sample}.bam")
    run:
        if len(input) == 1:
            shell("mv {input} {output} \n" + "mv {input}.bai {output}.bai")
        else:
            shell(SAMBAMBA + " merge -t 6 {output} {input}")

rule flagstat_1:
    input:
        BAM = rules.merge.output
    output:
        stats = sorted_alignment + "{sample}/{sample}.stats"
    params:
        threads = "--nthreads 10"
    shell:
        SAMBAMBA + " flagstat -t 3 {params.threads} {input.BAM} > {output.stats}"

rule mark_duplicates:
    input:
        BAM = rules.merge.output
    output:
        BAM = dedup_alignment + "{sample}/{sample}.bam",
        metrics= dedup_alignment +"{sample}/{sample}.metrics.txt"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        PICARD_TOOLS + " MarkDuplicates I={input} O={output.BAM} M={output.metrics}"

rule build_index:
    input:
        BAM = rules.mark_duplicates.output.BAM
    output:
        BAI = dedup_alignment + "{sample}/{sample}.bam.bai"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        PICARD_TOOLS + " BuildBamIndex I={input.BAM} O={output.BAI}"

rule flagstat_2:
    input:
        BAM = rules.mark_duplicates.output.BAM
    output:
        stats = dedup_alignment + "{sample}/{sample}.stats"
    params:
        threads = "--nthreads 10"
    shell:
        SAMBAMBA + " flagstat -t 3 {params.threads} {input.BAM} > {output.stats}"

rule picard_metrics:
    input:
        BAM = rules.mark_duplicates.output.BAM,
        BAI = rules.build_index.output.BAI,
        ref = ref
    output:
        metrics = metrics + "{sample}/{sample}.alignment_summary_metrics"
    params:
        prefix = metrics + "{sample}/{sample}"
    shell:
        LOAD_PICARD_TOOLS +
        LOAD_JAVA +
        LOAD_R +
        PICARD_TOOLS + " CollectMultipleMetrics I={input.BAM} O={params.prefix} R={input.ref}"

rule samtools_idxstats:
    input:
        BAM = rules.mark_duplicates.output.BAM,
        BAI = rules.build_index.output.BAI
    output:
        stats = metrics + "{sample}/{sample}_idxstats"
    shell:
        SAMTOOLS + " idxstats {input.BAM} > {output.stats}"

rule sam_flags:
    input:
        BAM = rules.mark_duplicates.output.BAM,
        BAI = rules.build_index.output.BAI
    params:
        metrics + "sam_flag_reads.txt"
    output:
        metrics + "sam_flag_logs/{sample}.log"
    shell:
        SAM_FLAGS + "{input.BAM} {output} {params} {wildcards.sample}"

rule coverage:
    input:
        BAM = rules.mark_duplicates.output.BAM,
        BAI = rules.build_index.output.BAI
    output:
        coverage + "{sample}.coverage"
    params:
        output = coverage + "{sample}.output",
        bed = coverage + "{sample}.output.per-base.bed.gz"
    shell:
        GET_COVERAGE + "{input.BAM} {output} {params.output} {params.bed}"

rule summary_cov:
    input:
        file_list = expand(rules.coverage.output,sample=sample),
        info = info_file
    output:
        coverage = coverage + "summary_coverage.txt"
    script:
        AVERAGE_COVERAGE
