#!/bin/bash
ls /home/lsf/data/machilus34/raw_data/*gz | cut -d"/" -f 7 | cut -d"." -f 1 | uniq | while read id
do
time bwa mem -t 5 -R '@RG\tID:foo\tPL:illumina\tSM:GWHACDM' /home/lsf/data/machilus34/reference/GWHACDM00000000.genome.fasta \
/home/lsf/data/machilus34/raw_data/${id}.clean.1.fastq.gz /home/lsf/data/machilus34/raw_data/${id}.clean.2.fastq.gz | /home/lsf/bin/samtools-1.12/samtools view -Sb - > /home/lsf/data/machilus34/bam/${id}.bam  && echo "** bwa mapping to bam done **"
done
#mv ./*.bam /home/lsf/data/machilus34/bam
cd /home/lsf/data/machilus34/bam
ls /home/lsf/data/machilus34/bam/*.bam | cut -d"/" -f 8 | cut -d"." -f 1 | while read id
do
time /home/lsf/bin/samtools-1.12/samtools sort -@ 4 -m 4G -O bam -o /home/lsf/data/machilus34/sortbam/${id}.sorted.bam /home/lsf/data/machilus34/bam/${id}.bam && echo "** samtools sort bam done **"
done
#mv *.sort.bam /home/lsf/data/machilus34/raw_data/sortbam/
cd /home/lsf/data/machilus34/sortbam/
ls /home/lsf/data/machilus34/sortbam/*.sort.bam | cut -d"/" -f 7 | cut -d"_" -f 1 | while read id
do
time bcftools mpileup -Ou -I -f /home/lsf/data/machilus34/reference/GWHACDM00000000.genome.fasta ${id}.sorted.bam | bcftools call -c -Ov  vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${id}.fq.gz  && echo "** bcf calling done **"
done
mv *.fq.gz /home/lsf/data/machilus34/fqgz/
cd /home/lsf/data/machilus34/raw_data/fqgz/
ls /home/lsf/data/machilus34/raw_data/fqgz/*.fq.gz | cut -d "/" -f 7 | cut -d "_" -f 1 | while read id
do
time fq2psmcfa -q 20 ${id}.fq.gz > ${id}.psmcfa  && echo "** fq to psmcfa done **"
done

mv *.psmcfa /home/lsf/data/machilus34/psmcfa

