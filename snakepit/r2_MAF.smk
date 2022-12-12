wildcard_constraints:
    pan = r'BSW|non_BSW',
    red = r"\d+_samples"

rule all:
    input:
        stats = expand(r2 + "{pan}_panel/{red}/{rep}/{cov}/{sample}_r2.stats", pan = panels, red = reductions, rep = replicas, cov = coverages, sample = sample),
        multi_stats = expand(r2 + "multibreed_panel/{per}/{rep}/{cov}/{sample}_r2.stats", per = percentages, rep = replicas, cov = coverages, sample = sample)

rule r2_MAF:
    input:
        truth = truth + "autosomes_biSNPs_beagle.vcf.gz",
        query = query + "{pan}_panel/{red}/{rep}/{cov}/imputed/all_autosomes.vcf.gz"
    output:
        stats = r2 + "{pan}_panel/{red}/{rep}/{cov}/{sample}_r2.stats"
    params:
        stats = "stats -s {sample}"
    threads: 2
    resources:
        mem_mb = 100,
        walltime = '00:15'
    shell:
        '''
        bcftools {params.stats} --threads {threads} {input.truth} {input.query} > {output.stats}
        '''

rule r2_MAF_multi:
    input:
        truth = truth + "autosomes_biSNPs_beagle.vcf.gz",
        query = query + "multibreed_panel/{per}/{rep}/{cov}/imputed/all_autosomes.vcf.gz"
    output:
        stats = r2 + "multibreed_panel/{per}/{rep}/{cov}/{sample}_r2.stats"
    params:
        stats = "stats -s {sample}"
    threads: 2
    resources:
        mem_mb = 100,
        walltime = '00:15'
    shell:
        '''
        bcftools {params.stats} --threads {threads} {input.truth} {input.query} > {output.stats}
        '''
