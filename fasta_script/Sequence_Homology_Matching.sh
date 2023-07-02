#!/bin/bash
#找同源基因的脚本，很简单，最后还建了个树，不过要跑很久
# 声明变量
database_file="database.fasta"
#如果有从pfam下载的种子文件，则声明文件即可，如没有，则为空
hmmsed="your.hmm file"
query_file="query.fasta"
result_file="result.out"
evalue=1e-5
union=T
intersection=T
#判断要用的软件是不是全部都安装了，以及输出软件的版本
software_list=(blast hmmsearch mmseq diamond mafft trimal iqtree)
for software in "${software_list[@]}"; do
    if ! command -v "$software" &> /dev/null; then
        echo "$software is not installed"
    else
        echo "$software version: $($software --version)"
    fi
done
# 建库
date | xargs echo "Building database... | " 
makeblastdb -in $database_file -dbtype prot -out ./blastpdb/database
mmseqs createdb $database_file ./mmseqsdb/database
diamond makedb --in $database_file -d ./dianmonddb/database
# hmm判断，有则match,无则build一个，然后match
if [ -n "$hmmsed" ]; then
    hmmsearch --tblout hmm.out --noali --notextw $hmmsed $query_file | awk '{print $1}' > hmm.out
else
    hmmbuild database.hmm $database_file
    hmmsearch --tblout hmm.out --noali --notextw database.hmm $query_file | awk '{print $1}' > hmm.out
fi
date | xargs echo "Complete database creation～！ | " 
date | xargs echo "Homologous sequences are being matched... ｜"
# 搜索
date | xargs echo "blastp is running... |"
blastp -query $query_file -db database -outfmt 6 -evalue $evalue | awk '{print $2}' > blastp.out
date | xargs echo "Blastp completed and the evalue is ${evalue} |"
date | xargs echo "mmseqs is running... |"
mmseqs search $query_file database result tmp --evalue $evalue --threads 4 --format-output 'query,target' | awk '{print $2}' > mmseqs.out
date | xargs echo "mmseq completed and the evalue is ${evalue} |"
date | xargs echo "diamond is running..."
diamond blastp -q $query_file -d database -o diamond.out -e $evalue --threads 4 | awk '{print $2}' > diamond.out
date | xargs echo "diamond completed and evalue is ${evalue} |"
# 判断用户输入：交集还是并集，如果是交集，则输出交集的list
if [ "$union" = "T" ]; then
    cat blastp.out mmseqs.out diamond.out hmm.out ｜ sort |  uniq -d > $result_file.intersection
    less $query_file | seqkit grep -f $result_file.intersection > result_seq_intersection.fasta
fi
# 如果是并集，则输出并集的list
if [ "$union" = "T" ]; then
    cat blastp.out mmseqs.out diamond.out hmm.out | sort -u > $result_file.union
    less $query_file | seqkit grep -f $result_file.union > result_seq_union.fasta
fi
