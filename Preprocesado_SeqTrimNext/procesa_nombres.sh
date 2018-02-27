for i in `ls -d *` ;  do
pushd .
cd $i/output_files
echo `date` procesando $i
if ! test -e paired_1_.fastq_old ; then
time ~/scripts/seqtk seq -l0 paired_1_.fastq > paired_1_4lines.fastq
mv paired_1_.fastq paired_1_.fastq_old
mv paired_1_4lines.fastq paired_1_.fastq
sleep 1
time ~/scripts/seqtk seq -l0 paired_2_.fastq > paired_2_4lines.fastq
mv paired_2_.fastq paired_2_.fastq_old
mv paired_2_4lines.fastq paired_2_.fastq
ls -l paired*
sleep 1
fi
popd
done