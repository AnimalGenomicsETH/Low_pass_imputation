wildcard_constraints:
    red = r"\d+_samples"

rule all:
    input:
        pure = expand(accuracy + "{pure}/{red}/{rep}/{cov}/{sam}_60.summary", pure = pures, red = reductions, rep = replicas, cov = coverages, sam = sample),
        stats_multi = expand(accuracy + "multibreed_panel/{per}/{rep}/{cov}/{sam}_60.summary", pure = pures, per = percentages, rep = replicas, cov = coverages, sam = sample)

rule accuracy_pure:
    input:
        vcf = glimpse + "{pure}/{red}/{rep}/{cov}/imputed/{sam}_autosomes.vcf.gz"
    output:
        summary = accuracy + "{pure}/{red}/{rep}/{cov}/{sam}.summary",
        average = accuracy + "{pure}/{red}/{rep}/{cov}/{sam}.average"
    params:
        query = "query -f '%INFO/INFO\n'",
        awk_all = "| awk '{print n+=1, $1}'",
        awk_av = "| awk '{c+=$1;n+=1} END {print c/n}'"
    resources:
        mem_mb = 150,
        walltime = '01:00'
    shell:
        '''
        bcftools {params.query} {input.vcf} {params.awk_all} > {output.summary} \n \
        bcftools {params.query} {input.vcf} {params.awk_av} > {output.average} \n
        '''

rule threshold_pure:
    input:
        summary = rules.accuracy_pure.output.summary
    output:
        summary = accuracy + "{pure}/{red}/{rep}/{cov}/{sam}_60.summary"
    params:
        awk = "awk '$2 > 0.6 {print $1, $2}'"
    resources:
        mem_mb = 200,
        walltime = '00:15'
    shell:
        '''
        {params.awk} {input.summary} > {output.summary}
        '''

rule accuracy_multi:
    input:
        vcf = glimpse + "multibreed_panel/{per}/{rep}/{cov}/imputed/{sam}_autosomes.vcf.gz"
    output:
        summary = accuracy + "multibreed_panel/{per}/{rep}/{cov}/{sam}.summary",
        average = accuracy + "multibreed_panel/{per}/{rep}/{cov}/{sam}.average"
    params:
        query = "query -f '%INFO/INFO\n'",
        awk_all = "| awk '{print n+=1, $1}'",
        awk_av = "| awk '{c+=$1;n+=1} END {print c/n}'"
    resources:
        mem_mb = 150,
        walltime = '01:00'
    shell:
        '''
        bcftools {params.query} {input.vcf} {params.awk_all} > {output.summary} \n \
        bcftools {params.query} {input.vcf} {params.awk_av} > {output.average} \n
        '''

rule threshold_multi:
    input:
        summary = rules.accuracy_multi.output.summary
    output:
        summary = accuracy + "multibreed_panel/{per}/{rep}/{cov}/{sam}_60.summary"
    params:
        awk = "awk '$2 > 0.6 {print $1, $2}'"
    resources:
        mem_mb = 200,
        walltime = '00:15'
    shell:
        '''
        {params.awk} {input.summary} > {output.summary}
        '''
