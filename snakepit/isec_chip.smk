configfile: "config.yaml"

chroms = list(range(1,30))

wildcard_constraints:
    sample = "[^/]+"

rule all:
    input:
        seq_stats = seq_data + "autosomes_33.stats",
        SNP_stats = seq_data + "autosomes_33_biSNPs.stats",
        isec_stats = expand(seq_data + "query_isec/{isec}.stats", isec=[f'{i:04}' for i in range(4)]),
        sample_stats = expand(seq_data + "sample_vcfs/{sample}/{isec}.stats", sample = samples, isec=[f'{i:04}' for i in range(4)]),
        chr_stats = expand(seq_data + "chromosome_vcfs/{chr}.stats", chr = chroms),
        seq_bed = seq_data + "intersections.bed",
        seq_stats = seq_data + "intersections.txt"

rule seq_samples:
    input:
        vcf = seq + "autosomes_gatk4.vcf.gz"
    output:
        vcf = seq_data + "autosomes_33.vcf.gz",
        index = seq_data + "autosomes_33.vcf.gz.tbi",
        stats = seq_data + "autosomes_33.stats"
    params:
        view_sample = "view -Oz --threads 4 -S " + samps,
        index = "-fp vcf"
    threads: 4
    resources:
        mem_mb = 3000,
        walltime = '45'
    shell:
        '''
        bcftools {params.view_sample} {input.vcf} -o {output.vcf} \n \
        tabix {params.index} {output.vcf} \n \
        bcftools stats {output.vcf} > {output.stats} \n
        '''

rule select_SNP:
    input:
        vcf = rules.seq_samples.output.vcf,
        ref = assembly
    output:
        vcf = seq_data + "autosomes_33_biSNPs.vcf.gz",
        stats = seq_data + "autosomes_33_biSNPs.stats"
    params:
        biSNPs = "view -m2 -M2 -v snps",
        index = "-fp vcf"
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '45'
    shell:
        '''
        bcftools {params.biSNPs} --threads {threads} -Oz -o {output.vcf} {input.vcf} \n \
        tabix {params.index} {output.vcf} \n \
        bcftools stats {output.vcf} > {output.stats} \n
        '''

rule bcftools_isec:
    input:
        seq_vcf = rules.select_SNP.output.vcf,
        chip_vcf = chip + "chip_all_filtered.vcf.gz"
    output:
        isecs = (seq_data + f"query_isec/{ISEC:04}.vcf.gz" for ISEC in range(4))
    params:
        isec = "isec -c some -p query_isec"
    threads: 4
    resources:
        mem_mb = 3000,
        walltime = '30'
    shell:
        '''
        bcftools {params.isec} --threads {threads} -O z {input.chip_vcf} {input.seq_vcf} \n
        '''

rule manipulate_VCF:
    input:
        vcf = seq_data + "query_isec/{isec}.vcf.gz"
    output:
         index = seq_data + "query_isec/{isec}.vcf.gz.tbi",
         stats = seq_data + "query_isec/{isec}.stats"
    resources:
        mem_mb = 500,
        walltime = '30'
    shell:
        '''
        tabix -fp vcf {input.vcf}
        bcftools stats {input.vcf} > {output.stats}
        '''

rule seq_split_sample:
    input:
        seq_vcf = seq_data + "query_isec/0003.vcf.gz",
        chip_vcf = chip + "{sample}.vcf.gz"
    output:
        vcf = seq_data + "sample_vcfs/{sample}/{sample}.vcf.gz",
        index = seq_data + "sample_vcfs/{sample}/{sample}.vcf.gz.tbi",
        stats = seq_data + "sample_vcfs/{sample}/{sample}.stats",
        isecs = (seq_data + f"sample_vcfs/{{sample}}/{ISEC:04}.vcf.gz" for ISEC in range(4))
    params:
        split_sample = "view --threads 2 -Oz -s {sample}",
        isec = "isec -c all -p sample_vcfs/{sample}",
        index = "-fp vcf"
    threads: 4
    resources:
        mem_mb = 1500,
        walltime = '30'
    shell:
        '''
        bcftools {params.split_sample} {input.seq_vcf} -o {output.vcf} \n \
        tabix {params.index} {output.vcf} \n \
        bcftools stats {output.vcf} > {output.stats} \n \
        bcftools {params.isec} --threads {threads} -O z {input.chip_vcf} {output.vcf} \n
        '''

rule manipulate_split_VCF:
    input:
        vcf = seq_data + "sample_vcfs/{sample}/{isec}.vcf.gz"
    output:
        index = seq_data + "sample_vcfs/{sample}/{isec}.vcf.gz.tbi",
        stats = seq_data + "sample_vcfs/{sample}/{isec}.stats"
    resources:
        mem_mb = 500,
        walltime = '30'
    shell:
        '''
        tabix -fp vcf {input.vcf}
        bcftools stats {input.vcf} > {output.stats}
        '''

rule seq_split_chr:
    input:
        vcf = seq_data + "query_isec/0003.vcf.gz"
    output:
        vcf = seq_data + "chromosome_vcfs/{chr}.vcf.gz",
        index = seq_data + "chromosome_vcfs/{chr}.vcf.gz.tbi",
        stats = seq_data + "chromosome_vcfs/{chr}.stats"
    params:
        view_chr = "view -r",
        view_out = "-Oz -o",
        index = "-fp vcf"
    resources:
        mem_mb = 3000,
        walltime = '45'
    shell:
        '''
        bcftools {params.view_chr} {wildcards.chr} {input.vcf} {params.view_out} {output.vcf} \n \
        tabix {params.index} {output.vcf} \n \
        bcftools stats {output.vcf} > {output.stats} \n
        '''

rule bedtools_intersect:
    input:
        DV_vcf = DV_data + "autosomes_33_biSNPs.vcf.gz",
        GATK_vcf = GATK_data + "autosomes_33_biSNPs.vcf.gz",
        chip_vcf = chip + "chip_all_filtered.vcf.gz"
    output:
        stats = seq_data + "intersections.txt",
        bed = seq_data + "intersections.bed"
    threads: 2
    resources:
        mem_mb = 500,
        walltime = '02:00'
    shell:
        '''
        bedtools multiinter -cluster -i {input.DV_vcf} {input.GATK_vcf} {input.chip_vcf} | cut -f 5 | sort | uniq -c > {output.stats} \n \
        bedtools multiinter -cluster -i {input.DV_vcf} {input.GATK_vcf} {input.chip_vcf} > {output.bed} \n
        '''