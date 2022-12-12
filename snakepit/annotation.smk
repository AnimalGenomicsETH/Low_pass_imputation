wildcard_constraints:
    varcal = r'GATK|DV'

rule all:
    input:
        vep_vcf_index = expand(vep_out + "{stage}/{varcal}/vep.vcf.gz.tbi", stage = stages, varcal = varcals)

rule VEP:
    input:
        gtf = gtf3 + "Bos_taurus.ARS-UCD1.2.104.chr.gtf3.gz",
        ref = assembly + "ARS_UCD.fasta"
    output:
        stats = vep_out + "{stage}/{varcal}/vep.stats.txt",
        vcf = vep_out + "{stage}/{varcal}/vep.vcf.gz",
        index = vep_out + "{stage}/{varcal}/vep.vcf.gz.tbi"
    params:
        VEP = config['tools']['VEP'],
        vcf = lambda wildcards: config["resources"][wildcards.stage][wildcards.varcal],
        flags = "--hgvs --symbol -o stdout --stats_text",
        compres = "| bgzip -c >",
        binding = "-B " + assembly + ":" + assembly + " -B " + gtf3 + ":" + gtf3 + " -B " + sing_dir + ":" + sing_dir + " -B " + vep_out + ":" + vep_out,
        index = "-fp vcf"
    threads: 4
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        '''
        singularity exec {params.binding} {params.VEP} /opt/vep/src/ensembl-vep/vep \
        -i {params.vcf} -gtf {input.gtf} -fasta {input.ref} {params.flags} --stats_file {output.stats} --vcf {params.compres} {output.vcf} \n \
        tabix {params.index} {output.vcf}
        '''