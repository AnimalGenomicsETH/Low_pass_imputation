rule all:
    input:
        R1 = expand(downsampling + "{coverage}/{sample}/{sample}_R1.fastq.gz",sample=sample,coverage=coverage),
        R2 = expand(downsampling + "{coverage}/{sample}/{sample}_R2.fastq.gz",sample=sample,coverage=coverage)"

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

rule downsampling:
    input:
        R1 = rules.QC.output.R1,
        R2 = rules.QC.output.R2
    output:
        R1 = downsampling + "{coverage}/{sample}/{sample}_R1.fastq.gz",
        R2 = downsampling + "{coverage}/{sample}/{sample}_R2.fastq.gz"
    params:
        cov = lambda wildcards: fold[wildcards.coverage]
    shell:
        SEQTK + " sample -s189 {input.R1} {params.cov} | gzip > {output.R1} \n" +
        SEQTK + " sample -s189 {input.R2} {params.cov} | gzip > {output.R2} \n"