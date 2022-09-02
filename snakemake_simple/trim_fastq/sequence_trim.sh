# Simple command to trim sequence in fastq
for sample in *.fastq
do  
   id=$(echo ${sample} | sed 's/.fastq//')  
   # 1. Trim fastq file
   seqtk trimfq -b 5 -e 10 $sample > ${id}.trimmed.fastq
   # 2. Compress fastq file
   gzip -c ${id}.trimmed.fastq > ${id}.trimmed.fastq.gz
   # 3. Remove intermediate files  
   rm ${id}.trimmed.fastq
done
