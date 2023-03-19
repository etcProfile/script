```bash
#!/bin/bash

# Define variables
raw_data_path=/home/lsf/data/machilus34/raw_data
reference_path=/home/lsf/data/machilus34/reference/GWHACDM00000000.genome.fasta
bam_path=/home/lsf/data/machilus34/bam
sortbam_path=/home/lsf/data/machilus34/sortbam
fqgz_path=/home/lsf/data/machilus34/fqgz
psmcfa_path=/home/lsf/data/machilus34/psmcfa

# Process fastq files
ls ${raw_data_path}/*gz | cut -d"/" -f 7 | cut -d"." -f 1 | uniq | while read id
do
    time bwa mem -t 5 -R '@RG\tID:foo\tPL:illumina\tSM:GWHACDM' ${reference_path} \
    ${raw_data_path}/${id}.clean.1.fastq.gz ${raw_data_path}/${id}.clean.2.fastq.gz | /home/lsf/bin/samtools-1.12/samtools view -Sb - > ${bam_path}/${id}.bam  && echo "** bwa mapping to bam done **"
done

# Process bam files
cd ${bam_path}
ls ${bam_path}/*.bam | cut -d"/" -f 8 | cut -d"." -f 1 | while read id
do
    time /home/lsf/bin/samtools-1.12/samtools sort -@ 4 -m 4G -O bam -o ${sortbam_path}/${id}.sorted.bam ${bam_path}/${id}.bam && echo "** samtools sort bam done **"
done

# Process sortbam files
cd ${sortbam_path}
ls ${sortbam_path}/*.sort.bam | cut -d"/" -f 7 | cut -d"_" -f 1 | while read id
do
    time bcftools mpileup -Ou -I -f ${reference_path} ${id}.sorted.bam | bcftools call -c -Ov  vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${fqgz_path}/${id}.fq.gz  && echo "** bcf calling done **"
done

# Process fqgz files
cd ${raw_data_path}/fqgz
ls ${raw_data_path}/fqgz/*.fq.gz | cut -d "/" -f 7 | cut -d "_" -f 1 | while read id
do
    time fq2psmcfa -q 20 ${id}.fq.gz > ${psmcfa_path}/${id}.psmcfa  && echo "** fq to psmcfa done **"
done

mv *.psmcfa ${psmcfa_path}
```
