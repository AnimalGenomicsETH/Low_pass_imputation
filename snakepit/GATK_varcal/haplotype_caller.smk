rule all:
    input:
        expand (hap_call + "gvcf/{chromosome}/{sample}_{chromosome}.g.vcf.gz", chromosome=chr_list, sample=sample)

rule recalibrator_creator:
    input:
        bam = dedup_bam + "{sample}/{sample}.bam",
        ref = assembly,
        db = BQSR
    output:
        temp(hap_call + "bam_recal/{sample}_recalibrator.table")
    params:
        chr = int_list
    shell:
        LOAD_JAVA +
        GATK4 +
        " BaseRecalibrator " +
        " -I {input.bam} " +
        " {params.chr} " +
        " -R {input.ref} " +
        " --known-sites {input.db} " +
        " -O {output}"

rule base_recalibrator:
    input:
        recal = rules.recalibrator_creator.output,
        bam = dedup_bam + "{sample}/{sample}.bam",
        ref = assembly
    output:
        temp(hap_call + "bam_recal/{sample}_recalibrated.bam")
    params:
        chr = int_list
    shell:
        LOAD_JAVA +
        GATK4 +
        " ApplyBQSR " +
        " -R {input.ref} " +
        " {params.chr} " +
        " -I {input.bam} " +
        " --bqsr-recal-file {input.recal} " +
        " -O {output} "

rule haplotype_caller:
    input:
        recal_bam = rules.base_recalibrator.output,
        ref = assembly
    output:
        hap_call + "gvcf/{chromosome}/{sample}_{chromosome}.g.vcf.gz"
    params:
        chr = "{chromosome}"
    shell:
        LOAD_JAVA +
        GATK4 +
        " HaplotypeCaller " +
        " -I {input.recal_bam} " +
        " -R {input.ref} " +
        " -L {params.chr} " +
        " -O {output} " +
        " --ERC GVCF "
