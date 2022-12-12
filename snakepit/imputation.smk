rule all:
    input:
        expand(varcal_path + "cohort.{chromosome}.stats", chromosome=chroms),
        expand(beagle_path + "{chromosome}_beagle.vcf.gz.tbi", chromosome=chroms),
        expand(beagle_path + "{chromosome}_beagle.stats", chromosome=chroms)

rule stats_vcf:
    input:
        varcal_path + "cohort.{chromosome}.vcf.gz"
    output:
        stats = varcal_path + "cohort.{chromosome}.stats"
    shell:
        "bcftools stats {input} > {output.stats} \n"

rule beagle_imputation:
    input:
        varcal_path + "cohort.{chromosome}.vcf.gz"
    output:
        beagle_path + "{chromosome}_beagle.vcf.gz"
    params:
        ch = "{chromosome}"
    shell:
        LOAD_JAVA +
        " java -Xss25m -Xmx40G -jar " + BEAGLE +
        " gl={input} " +
        " nthreads=18 " +
        " out=" + beagle_path + "{params.ch}_beagle "

rule finish_vcf:
    input:
        rules.beagle_imputation.output
    output:
        tbi = beagle_path + "{chromosome}_beagle.vcf.gz.tbi",
        stats = beagle_path + "{chromosome}_beagle.stats"
    shell:
        "tabix -fp vcf {input} \n" +
        "bcftools stats {input} > {output.stats} \n"
