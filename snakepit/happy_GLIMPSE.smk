wildcard_constraints:
    sample = r"[A-Z0-9]+"

rule all:
    input:
        truth_concat_stats = truth_set + "autosomes_beagle.stats",
        truth_biSNPs_stats = truth_set + "autosomes_biSNPs_beagle.stats",
        truth_split_stats = expand(truth_set + "{sample}_autosomes_biSNPs_beagle.stats", sample = sample),
        query_concat_stats = expand(query_set + "{red}/{shuf}/{cov}/imputed/all_autosomes.stats", red = reds, shuf = shufs, cov = covs),
        query_split_stats = expand(query_set + "{red}/{shuf}/{cov}/imputed/{sample}_autosomes.stats", red = reds, shuf = shufs, cov = covs, sample = sample),
        happy_summary = expand(happy_out + "{red}/{shuf}/{cov}/{sample}.summary.csv", red = reds, shuf = shufs, cov = covs, sample = sample)

rule truth_set_concat:
    output:
        VCF = truth_set + "autosomes_beagle.vcf.gz",
        stats = truth_set + "autosomes_beagle.stats"
    params:
        list = "-f" + truth_set + "autosomal_list.txt",
        out = "-Oz -o"
    threads: 1
    resources:
        mem_mb = 600,
        walltime = '00:30'
    shell:
        '''
        bcftools concat {params.list} {params.out} {output.VCF} \n \
        stats {output.VCF} > {output.stats} \
        '''

rule truth_bi_SNPs:
    input:
        VCF = rules.truth_set_concat.output.VCF
    output:
        VCF = truth_set + "autosomes_biSNPs_beagle.vcf.gz",
        stats = truth_set + "autosomes_biSNPs_beagle.stats",
    params:
        view = "view -m2 -M2 -v snps",
        out = "-Oz -o",
        index = "-fp vcf"
    threads: 2
    resources:
        mem_mb = 500,
        walltime = '00:30'
    shell:
        '''
        bcftools {params.view} {input.VCF} --threads {threads} {params.out} {output.VCF} \n \
        tabix -fp vcf {output.VCF} \n \
        bcftools stats {output.VCF} > {output.stats} \n
        '''

rule truth_set_split:
    input:
        VCF = rules.truth_bi_SNPs.output.VCF
    output:
        VCF = truth_set + "{sample}_autosomes_biSNPs_beagle.vcf.gz",
        stats = truth_set + "{sample}_autosomes_biSNPs_beagle.stats"
    params:
        split_sample = "view -s {sample}",
        out = "| bgzip -c >"
    threads: 2
    resources:
        mem_mb = 600,
        walltime = '00:30'
    shell:
        '''
        bcftools {params.split_sample} {input.VCF} --threads {threads} {params.out} {output.VCF} \n \
        tabix -fp vcf {output.VCF} \n \
        bcftools stats {output.VCF} > {output.stats} \n
        '''

rule query_set_concat:
    output:
        VCF = query_set + "{red}/{shuf}/{cov}/imputed/all_autosomes.vcf.gz",
        stats = query_set + "{red}/{shuf}/{cov}/imputed/all_autosomes.stats"
    params:
        list = "-f" + query_set + "{red}/{shuf}/{cov}/imputed/autosomal_list.txt",
        out = "-Oz -o"
    threads: 1
    resources:
        mem_mb = 600,
        walltime = '00:30'
    shell:
        '''
        bcftools concat {params.list} {params.out} {output.VCF} \n \
        bcftools stats {output.VCF} > {output.stats} \
        '''

rule query_set_split:
    input:
        VCF = rules.query_set_concat.output.VCF
    output:
        VCF = query_set + "{red}/{shuf}/{cov}/imputed/{sample}_autosomes.vcf.gz",
        stats = query_set + "{red}/{shuf}/{cov}/imputed/{sample}_autosomes.stats"
    params:
        split_sample = "view -s {sample}",
        out = "| bgzip -c >"
    threads: 2
    resources:
        mem_mb = 600,
        walltime = '00:30'
    shell:
        '''
        bcftools {params.split_sample} {input.VCF} threads {threads} {params.out} {output.VCF} \n \
        bcftools stats {output.VCF} > {output.stats} \n
        '''

rule happy:
    input:
        truth_VCF = rules.truth_set_split.output.VCF,
        query_VCF = rules.query_set_split.output.VCF,
        ref = ref + "ARS_UCD.fasta"
    output:
        CSV = happy_out + "{red}/{shuf}/{cov}/{sample}.summary.csv",
        VCF = happy_out + "{red}/{shuf}/{cov}/{sample}.vcf.gz"
    params:
        singularity = "singularity exec",
        binding = "-B " + happy_in + ":" + happy_in + " -B " + truth_set + ":" + truth_set + " -B " + query_set + ":" + query_set + " -B" + happy_out + ":" + happy_out + " -B " + ref + ":" + ref + " -B " + happy_in + "temp/{red}/{shuf}/{cov}/{sample}:/tmp ",
        out = happy_out + "{red}/{shuf}/{cov}/{sample}",
        tmp = happy_in + "temp/{red}/{shuf}/{cov}/{sample}",
        happy = HAPPY,
    threads: 12
    resources:
        mem_mb = 500,
        walltime = '04:00'
    shell:
        '''
        mkdir -p {params.tmp} \n \
        export HGREF={input.ref} \n \
        {params.singularity} {params.binding} {params.happy} /opt/hap.py/bin/hap.py {input.truth_VCF} {input.query_VCF} -r {input.ref} -o {params.out} --threads {threads}
        '''
