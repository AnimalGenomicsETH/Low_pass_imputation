module load jdk

GATK4=gatk-4.2.0.0/gatk
REF=file.fasta
gvcf_dir=haplotype_caller/gvcf
Chr=MYCHR

cd $TMPDIR
while read sample
do
    cp ${gvcf_dir}/${Chr}/${sample}_${Chr}.g.vcf{.gz,.gz.tbi} .
    echo $sample
done < selected_samples.txt

export TILEDB_DISABLE_FILE_LOCKING=1

${GATK4} \
GenomicsDBImport \
    --sample-name-map /joint_genotyping/maps/${Chr}.map \
    --genomicsdb-workspace-path db_${Chr} \
    -L ${Chr} \
    --reader-threads 25 \
    --batch-size 50 \
    -R ${REF}

${GATK4} \
GenotypeGVCFs  \
-R ${REF}  \
-L ${Chr} \
-O ${Chr}_gatk4.vcf.gz \
-V gendb://db_${Chr}

if [[ -f ${Chr}_gatk4.vcf.gz.tbi ]]
then
    mv ${Chr}_gatk4.vcf.gz ${Chr}_gatk4.vcf.gz.tbi ${LS_SUBCWD}
fi
