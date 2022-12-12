rule all:
    input:
        happy_summary = expand(happy_out + "{stage}/{varcal}/{sample}.summary.csv", stage = stages, varcal = varcals, sample = samples),
        happy_stats = expand(happy_out + "{stage}/{varcal}/{sample}.stats", stage = stages, varcal = varcals, sample = samples)

rule happy:
    input:
        true_VCF = chip_dir + "{sample}.vcf.gz",
        ref = ref + "ARS_UCD.fasta"
    output:
        CSV = happy_out + "{stage}/{varcal}/{sample}.summary.csv",
        VCF = happy_out + "{stage}/{varcal}/{sample}.vcf.gz"
    params:
        query_VCF = lambda wildcards: [config['files'][wildcards.stage][varcal] for varcal in varcals],
        singularity = "singularity exec",
        binding = "-B " + happy_in + ":" + happy_in + " -B " + chip_dir + ":" + chip_dir + " -B " + seq_dir + ":" + seq_dir + " -B " + happy_out + ":" + happy_out + " -B " + ref + ":" + ref + " -B " + happy_in + "temp/{stage}/{varcal}/{sample}:/tmp ",
        out = happy_out + "{stage}/{varcal}/{sample}",
        tmp = happy_in + "temp/{stage}/{varcal}/{sample}",
        happy = HAPPY,
        threads = "--threads 12"
    threads: 12
    resources:
        mem_mb = 750,
        walltime = '04:00'
    shell:
        '''
        mkdir -p {params.tmp} \n \
        export HGREF={input.ref} \n \
        {params.singularity} {params.binding} {params.happy} /opt/hap.py/bin/hap.py {input.true_VCF} {params.query_VCF}{wildcards.sample}/0003.vcf.gz -r {input.ref} -o {params.out} {params.threads}
        '''

rule stats:
    input:
        VCF = rules.happy.output.VCF
    output:
        stats = happy_out + "{stage}/{varcal}/{sample}.stats"
    threads: 1
    resources:
        mem_mb = 250,
        walltime = '00:45'
    shell:
        '''
        bcftools stats {input.VCF} > {output.stats}
        '''
