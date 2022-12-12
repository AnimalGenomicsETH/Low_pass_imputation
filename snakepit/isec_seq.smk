configfile: "config.yaml"

CALLERS=['DV','GATK']

rule all:
    input:
        stats = expand(intersect + "{stage}/{mode}/{isec}.stats", stage = config['resources'], isec=[f'{i:04}' for i in range(4)],mode=['all','none'])

rule bcftools_isec:
    output:
        isecs = (intersect + f"{{stage}}/{{mode}}/{ISEC:04}.vcf.gz" for ISEC in range(4))
    params:
        files = lambda wildcards: [config['resources'][wildcards.stage][caller] for caller in CALLERS]
    threads: 4
    resources:
        mem_mb = 1000,
        walltime = '00:30'
    shell:
        '''
        bcftools isec -c {wildcards.mode} -p {wildcards.stage}/{wildcards.mode} --threads {threads} -O z {params.files}
        '''

rule manipulate_VCF:
    input:
        intersect + "{stage}/{mode}/{isec}.vcf.gz"
    output:
        index = intersect + "{stage}/{mode}/{isec}.vcf.gz.tbi",
        stats = intersect + "{stage}/{mode}/{isec}.stats"
    resources:
        mem_mb = 500,
        walltime = '00:15'
    shell:
        '''
        tabix -fp vcf {input} \n \
        bcftools stats {input} > {output.stats}
        '''
