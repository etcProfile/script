#!/bin/bash
#找同源基因的脚本，很简单，最后还建了个树，不过要跑很久 Author: Luo Shaofan

# 声明变量
database_file="database.fasta"
query_file="query.fasta"
result_file="result.out"

# 建库
makeblastdb -in $database_file -dbtype prot -out database
mmseqs createdb $database_file database
diamond makedb --in $database_file -d database
hmmbuild database.hmm $database_file

# 搜索
blastp -query $query_file -db database -outfmt 6 -evalue 1e-5 | awk '{print $2}' > blastp.out
mmseqs search $query_file database result tmp --threads 4 --format-output 'query,target' | awk '{print $2}' > mmseqs.out
diamond blastp -q $query_file -d database -o diamond.out -e 1e-5 --threads 4 | awk '{print $2}' > diamond.out
hmmsearch --tblout hmm.out --noali --notextw database.hmm $query_file | awk '{print $1}' > hmm.out

# 并集
cat blastp.out mmseqs.out diamond.out hmm.out | sort -u > $result_file

# 找出基因序列
less $query_file | seqkit grep -f $result_file > result.fasta

# 建树
mafft --maxiterate 1000 --localpair result.fasta > result.aln
trimal -in result.aln -out result.trimmed.aln -automated1
iqtree -s result.trimmed.aln -bb 1000 -m MFP -nt AUTO -bnni
