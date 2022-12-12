cd $TMPDIR

export LD_LIBRARY_PATH=/htslib/

for chr in {1..32}
do
if [[ $chr = "30" ]]; then
    chr=X
fi
if [[ $chr = "31" ]]; then
    chr=Y
fi
if [[ $chr = "32" ]]; then
    chr=MT
fi

mosdepth -c $chr $3 $1

if [[ $chr != "X" ]]; then
coverage_sum=`zcat $4 | awk '{ total += ($3-$2)*$4; count++ } END {print total}'`;length=`zcat $4 | tail -n 1 | awk '{print $3}'`; coverage=`printf %.2f $(echo $coverage_sum/$length | bc -l)`; echo $1 $chr $length $coverage >> $2
fi

if [[ $chr = "X" ]]; then
coverage_sum=`zcat $4 | awk '$2<133300518' | awk '{ total += ($3-$2)*$4; count++ } END {print total}'`;length=133300518; coverage=`printf %.2f $(echo $coverage_sum/$length | bc -l)`; echo $1 $chr 133300518 $coverage >> $2
coverage_sum=`zcat $4 | awk '$2>133300518' | awk '{ total += ($3-$2)*$4; count++ } END {print total}'`;length=5708626; coverage=`printf %.2f $(echo $coverage_sum/$length | bc -l)`; echo $1 PAR 5708626 $coverage >> $2
fi
done

rm $3*