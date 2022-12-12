chroms = list(range(1,30))

rule all:
    input:
        expand(GLIMPSE_path + "{coverage}/imputed/merged_recoded.bed", coverage=coverage)

rule prep_ref_panel:
    input:
        hap_panel + "{chromosome}_beagle.vcf.gz"
    output:
        bcf = GLIMPSE_path + "hap_panel/chr_{chromosome}.bcf",
        csi = GLIMPSE_path + "hap_panel/chr_{chromosome}.bcf.csi"
    params:
        bial_snps = "-m 2 -M 2 -v snps",
        threads = "--threads 4",
        out = "-Ob -o",
        index = "-f"
    shell:
        "bcftools view {params.bial_snps} {params.threads} {input} {params.out} {output.bcf} \n" +
        "bcftools index {params.index} {output.bcf} \n"

rule extract_var_sites:
    input:
        rules.prep_ref_panel.output.bcf
    output:
        sites = GLIMPSE_path + "hap_panel/chr_{chromosome}.sites.vcf.gz",
        sites_csi = GLIMPSE_path + "hap_panel/chr_{chromosome}.sites.vcf.gz.csi",
        sites_tsv = GLIMPSE_path + "hap_panel/chr_{chromosome}.sites.tsv.gz",
        sites_tbi = GLIMPSE_path + "hap_panel/chr_{chromosome}.sites.tsv.gz.tbi"
    params:
        bial_snps = "-G -m 2 -M 2 -v snps",
        out = "-Oz -o",
        index = "-f",
        tab = "-f'%CHROM\\t%POS\\t%REF,%ALT\\n'",
        bgzip = "-c",
        tabix = "-s1 -b2 -e2"
    shell:
        "bcftools view {params.bial_snps} {input} {params.out} {output.sites} \n" +
        "bcftools index {params.index} {output.sites} \n" +
        "bcftools query {params.tab} {output.sites} | bgzip {params.bgzip} > {output.sites_tsv} \n" +
        "tabix {params.tabix} {output.sites_tsv} \n"

rule compute_GLs:
    input:
        VCF = rules.extract_var_sites.output.sites,
        TSV = rules.extract_var_sites.output.sites_tsv,
        BAM = down_samples + "{coverage}/{sample}/{sample}.bam"
    output:
        VCF = GLIMPSE_path + "{coverage}/samples/{sample}/{sample}_{chromosome}.vcf.gz",
        CSI = GLIMPSE_path + "{coverage}/samples/{sample}/{sample}_{chromosome}.vcf.gz.csi"
    params:
        ref = "-f" + ref,
        format = "-I -E -a 'FORMAT/DP' -T",
        chr = "-r {chromosome}",
        pipe = "-Ou |",
        call = "-Aim -C alleles -T",
        out = "-Oz -o",
        index = "-f"
    shell:
        "bcftools mpileup {params.ref} {params.format} {input.VCF} {params.chr} {input.BAM} {params.pipe} bcftools call {params.call} {input.TSV} {params.out} {output.VCF} \n" +
        "bcftools index {params.index} {output.VCF} \n"

rule merge_samples:
    input:
        expand(GLIMPSE_path + "{{coverage}}/samples/{sample}/{sample}_{{chromosome}}.vcf.gz", sample=sample)
    output:
        list = GLIMPSE_path + "{coverage}/merged/list_{chromosome}.txt",
        VCF = GLIMPSE_path + "{coverage}/merged/{chromosome}.vcf.gz",
        csi = GLIMPSE_path + "{coverage}/merged/{chromosome}.vcf.gz.csi",
    params:
        cov = "{coverage}",
        chr = "{chromosome}",
        merge = "-m none -r {chromosome}",
        out = "-Oz -o",
        list = "-l",
        index = "-f"
    shell:
        "ls " + GLIMPSE_path + "{params.cov}/samples/*/*_{params.chr}.vcf.gz > {output.list} \n" +
        "bcftools merge {params.merge} {params.out} {output.VCF} {params.list} {output.list} \n" +
        "bcftools index {params.index} {output.VCF} \n"

rule chunk:
    input:
        VCF = rules.extract_var_sites.output.sites
    output:
        GLIMPSE_path + "hap_panel/chunks/chunk_chr_{chromosome}.txt"
    params:
        chr = "--region {chromosome}",
        sizes = "--window-size 2000000 --buffer-size 200000"
    shell:
        GLIMPSE + "GLIMPSE_chunk_static --input {input.VCF} {params.chr} {params.sizes} --output {output} \n"

rule genetic_map:
    input:
        rules.extract_var_sites.output.sites_tsv
    output:
        txt = GLIMPSE_path + "hap_panel/genetic_map/chr_{chromosome}.txt",
        gz = GLIMPSE_path + "hap_panel/genetic_map/chr_{chromosome}.txt.gz"
    params:
        awk = "'{print $2, 1, $2/1000000}'"
    shell:
        "echo pos chr cM > {output.txt} \n" +
        "zcat {input} | awk {params.awk} >> {output.txt} \n" +
        "cat {output.txt} | bgzip -c > {output.gz} \n"

rule impute_phase:
    input:
        chunk = rules.chunk.output,
        VCF = rules.merge_samples.output.VCF,
        ref = rules.prep_ref_panel.output.bcf,
        map = rules.genetic_map.output.gz
    output:
        BCF = GLIMPSE_path + "{coverage}/imputed/{chromosome}.bcf",
        bcsi = GLIMPSE_path + "{coverage}/imputed/{chromosome}.bcf.csi",
        VCF = GLIMPSE_path + "{coverage}/imputed/{chromosome}.vcf.gz",
        vcsi = GLIMPSE_path + "{coverage}/imputed/{chromosome}.vcf.gz.csi"
    params:
        #Use a * to indicate the {ID} value is an unknown wildcard
        chunk_file =  GLIMPSE_path + "{coverage}/merged/{chromosome}_*.bcf",
        chunk_files = GLIMPSE_path + "{coverage}/merged/{chromosome}_list",
        index = "index -f",
        view = "view",
        bgzip = "bgzip -c >"
    run:
        with open(params.chunk_files,'w') as fout, open(str(input.chunk),'r') as chunks_in:
            #enumerate, starting from, so ID counts as 1,2,3... as we go through the file
            for ID,line in enumerate(chunks_in,1):
                #split the line and take the IRG and ORG values
                (IRG,ORG) = line.split('\t')[2:4]
                #Now we create a new string, replacing the * wildcard with the ID of this iteration
                chunk = params.chunk_file.replace('*',str(ID))
                shell(GLIMPSE + "GLIMPSE_phase_static --input {input.VCF} --reference {input.ref} --map {input.map} --input-region " + IRG + " --output-region " + ORG + " --output " + chunk + " --thread 8")
                shell('bcftools index -f ' + chunk)
                #write out the file name for this ID to the list of all names
                fout.write(chunk+'\n')
        #now merge
        shell(GLIMPSE + "GLIMPSE_ligate_static --input {params.chunk_files} --output {output.BCF}")
        shell('bcftools {params.index} {output.BCF}')
        shell('bcftools {params.view} {output.BCF} | {params.bgzip} {output.VCF}')
        shell('bcftools {params.index} {output.VCF}')

rule merge_bcf:
    input:
        BCF = expand(GLIMPSE_path + "{{coverage}}/imputed/{chromosome}.bcf",chromosome=chroms)
    output:
        list = GLIMPSE_path + "{coverage}/imputed/imputed_list",
        VCF = GLIMPSE_path + "{coverage}/imputed/merged.vcf.gz",
        vcsi = GLIMPSE_path + "{coverage}/imputed/merged.vcf.gz.csi"
    params:
        concat = "concat -f",
        out = "-Oz -o",
        index = "index -f"
    shell:
        "ls {input.BCF} >> {output.list} \n" +
        "bcftools {params.concat} {output.list} {params.out} {output.VCF} \n" +
        "bcftools {params.index} {output.VCF} \n"

rule plink:
    input:
        rules.merge_bcf.output.VCF
    output:
        bim = GLIMPSE_path + "{coverage}/imputed/merged.bim",
        ids = GLIMPSE_path + "{coverage}/imputed/recode_SNP_ids",
        recode = GLIMPSE_path + "{coverage}/imputed/merged_recoded.bed"
    params:
        mem = "--memory 20000",
        par = "--set-missing-var-ids @_#\$1_\$2 --const-fid --allow-no-sex --make-bed --cow --allow-extra-chr --biallelic-only --snps-only --list-duplicate-vars suppress-first --chr 1-29 --keep-allele-order",
        PLINK = GLIMPSE_path + "{coverage}/imputed/merged",
        awk = """awk '{print $2, $1"_"$4}' >""",
        keep = "--keep-allele-order --cow --update-name",
        bed = "--make-bed --out",
        recode = GLIMPSE_path + "{coverage}/imputed/merged_recoded"
    shell:
        PLINK + "plink {params.mem} --vcf {input} {params.par} --out {params.PLINK} \n" +
        "cat {output.bim} | {params.awk} {output.ids} \n" +
        "plink {params.mem} --bfile {params.PLINK} {params.keep} {output.ids} {params.bed} {params.recode} \n"