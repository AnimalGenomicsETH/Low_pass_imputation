all_reads=`samtools view -c $1`
MQ_sam_flag=`samtools view -c -q 10 -F 1796 $1`
multi_sam_flag=`samtools view -c -F 256 $1`
echo -e $4'\t'${all_reads}'\t'${MQ_sam_flag}'\t'${multi_sam_flag} >> $3
echo 'Sample' $4 'has been successfully explored'  > $2