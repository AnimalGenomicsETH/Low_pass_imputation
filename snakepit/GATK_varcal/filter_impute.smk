rule all:
    input:
        expand(vcf_beagle_path + "{chromosome}/beagle.vcf.gz.tbi", chromosome=chroms),
        expand(vcf_beagle_path + "{chromosome}/beagle.stats", chromosome=chroms)

rule select_SNP:
    input:
        vcf = vcf_genot_path + "{chromosome}_gatk4.vcf.gz",
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/snp.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants  " +
        " -R {input.ref} " +
        " -V {input.vcf}  " +
        " --select-type-to-include SNP  " +
        " --output {output} "

rule filter_SNP:
    input:
        vcf = rules.select_SNP.output,
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/snp_filtered.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        "  VariantFiltration " +
        "  -R {input.ref} " +
        "  -V {input.vcf} " +
        " --filter \"QD < 2.0\" --filter-name \"snp_QD2\" "
        " --filter \"QUAL < 30.0\" --filter-name \"snp_QUAL30\" "
        " --filter \"SOR > 3.0\" --filter-name \"snp_SOR3\" "
        " --filter \"FS > 60.0\" --filter-name \"snp_FS60\" "
        " --filter \"MQ < 40.0\" --filter-name \"snp_MQ40\" "
        " --filter \"MQRankSum < -12.5\" --filter-name \"snp_MQRankSum-12.5\" "
        " --filter \"ReadPosRankSum < -8.0\" --filter-name \"snp_ReadPosRankSum-8\" "
        " --missing-values-evaluate-as-failing true "
        " --output {output} "

rule select_indel:
    input:
        vcf = vcf_genot_path + "{chromosome}_gatk4.vcf.gz",
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/indel.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants  " +
        " -R {input.ref} " +
        " -V {input.vcf}  " +
        " --select-type-to-include INDEL " +
        " --output {output} "


rule filter_indel:
    input:
        vcf = rules.select_indel.output,
        ref = assembly
    output:
        vcf_beagle_path + "{chromosome}/indel_filtered.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        "  VariantFiltration " +
        "  -R {input.ref} " +
        "  -V {input.vcf} " +
        " --filter \"QD < 2.0\" --filter-name \"indel_QD2\" "
        " --filter \"QUAL < 30.0\" --filter-name \"indel_QUAL30\" "
        " --filter \"SOR > 10.0\" --filter-name \"indel_SOR10\" "
        " --filter \"FS > 200.0\" --filter-name \"indel_FS200\" "
        " --filter \"ReadPosRankSum < -20.0\" --filter-name \"indel_ReadPosRankSum-20\" "
        " --missing-values-evaluate-as-failing true "
        " --output {output} "

rule merge_variants:
    input:
        sic=rules.filter_SNP.output,
        idc=rules.filter_indel.output
    output:
        vcf_beagle_path + "{chromosome}/merge_var.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        " MergeVcfs " +
        " --INPUT {input.sic} " +
        " --INPUT {input.idc} " +
        " --OUTPUT {output} "

rule remove_filtered:
    input:
        var = rules.merge_variants.output
    output:
        vcf_beagle_path + "{chromosome}/prune_var.vcf.gz"
    shell:
        LOAD_JAVA +
        GATK4 +
        " SelectVariants "
        " -V {input.var} "
        " --exclude-filtered "
        "--output {output}"

rule beagle_imputation:
    input:
        rules.remove_filtered.output
    output:
        vcf_beagle_path + "{chromosome}/beagle.vcf.gz"
    params:
        ch = "{chromosome}"
    shell:
        LOAD_JAVA +
        " java -Xss25m -Xmx40G -jar " + BEAGLE +
        " gl={input} " +
        " nthreads=18 " +
        " out=" + vcf_beagle_path + "{params.ch}/beagle "

rule index_creation:
    input:
        select_SNP = rules.select_SNP.output,
        filter_SNP = rules.filter_SNP.output,
        select_indel = rules.select_indel.output,
        filter_indel = rules.filter_indel.output,
        merge_variants = rules.merge_variants.output,
        remove_filtered = rules.remove_filtered.output,
        beagle_imputation = rules.beagle_imputation.output
    output:
        select_SNP = vcf_beagle_path + "{chromosome}/snp.vcf.gz.tbi",
        filter_SNP = vcf_beagle_path + "{chromosome}/snp_filtered.vcf.gz.tbi",
        select_indel = vcf_beagle_path + "{chromosome}/indel.vcf.gz.tbi",
        filter_indel = vcf_beagle_path + "{chromosome}/indel_filtered.vcf.gz.tbi",
        merge_variants = vcf_beagle_path + "{chromosome}/merge_var.vcf.gz.tbi",
        remove_filtered = vcf_beagle_path + "{chromosome}/prune_var.vcf.gz.tbi",
        beagle_imputation = vcf_beagle_path + "{chromosome}/beagle.vcf.gz.tbi"
    shell:
        TABIX + " -fp vcf {input.select_SNP} \n" +
        TABIX + " -fp vcf {input.filter_SNP} \n" +
        TABIX + " -fp vcf {input.select_indel} \n" +
        TABIX + " -fp vcf {input.filter_indel} \n" +
        TABIX + " -fp vcf {input.merge_variants} \n" +
        TABIX + " -fp vcf {input.remove_filtered} \n" +
        TABIX + " -fp vcf {input.beagle_imputation} \n"

rule stats_creation:
    input:
        select_SNP = rules.select_SNP.output,
        filter_SNP = rules.filter_SNP.output,
        select_indel = rules.select_indel.output,
        filter_indel = rules.filter_indel.output,
        merge_variants = rules.merge_variants.output,
        remove_filtered = rules.remove_filtered.output,
        beagle_imputation = rules.beagle_imputation.output,
        beagle_index = rules.index_creation.output.beagle_imputation
    output:
        select_SNP = vcf_beagle_path + "{chromosome}/snp.stats",
        filter_SNP = vcf_beagle_path + "{chromosome}/snp_filtered.stats",
        select_indel = vcf_beagle_path + "{chromosome}/indel.stats",
        filter_indel = vcf_beagle_path + "{chromosome}/indel_filtered.stats",
        merge_variants = vcf_beagle_path + "{chromosome}/merge_var.stats",
        remove_filtered = vcf_beagle_path + "{chromosome}/prune_var.stats",
        beagle_imputation = vcf_beagle_path + "{chromosome}/beagle.stats"
    shell:
        BCFTOOLS + " stats {input.select_SNP} > {output.select_SNP} \n" +
        BCFTOOLS + " stats {input.filter_SNP} > {output.filter_SNP} \n" +
        BCFTOOLS + " stats {input.select_indel} > {output.select_indel} \n" +
        BCFTOOLS + " stats {input.filter_indel} > {output.filter_indel} \n" +
        BCFTOOLS + " stats {input.merge_variants} > {output.merge_variants} \n" +
        BCFTOOLS + " stats {input.remove_filtered} > {output.remove_filtered} \n" +
        BCFTOOLS + " stats {input.beagle_imputation} > {output.beagle_imputation} \n"
