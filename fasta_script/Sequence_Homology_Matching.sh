#!/bin/bash
#找同源基因的脚本，很简单，最后还建了个树，不过要跑很久
# 声明变量
database_file="database.fasta"
#如果有从pfam下载的种子文件，则声明文件即可，如没有，则为空
hmmseed="your.hmm file"
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
    hmmsearch --tblout hmm.out --noali --notextw $hmmsed $query_file 
else
    hmmbuild database.hmm $database_file
    hmmsearch --tblout hmm.out --noali --notextw database.hmm $query_file 
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
#标准化序列
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${id}_cds.fa > ${id}_cds.tmp.fa
#删除终止密码子
sed -e "s/TGA$//" -e "s/TAA$//" -e "s/TAG$//" ${id}_cds.tmp.fa > ${id}_no_stop_codon.fa
#核酸序列翻译为蛋白序列
transeq -sequence ${id}_no_stop_codon.fa -outseq ${id}_protein.fa
#去掉header里面冗余的部分，使得核酸和蛋白的序列header一样
sed 's/_1//g' -i ${id}_protein.fa
#用mafft来序列对齐
# Perform a multiple sequence alignment of the protein sequences in the file ${dir}.pep.fasta using the MAFFT program
mafft --maxiterate 1000 --thread 40 --localpair ${dir}.pep.fasta > all_${dir}_pep.aln.fasta
# Convert the protein alignment to a codon alignment using the pal2nal.pl program
pal2nal.pl  all_${dir}_pep.aln.fasta ${dir}.cds.fasta -nogap -codontable 11 -output clustal > ${dir}_cds_pal2nal.aln
# Convert the codon alignment to the input format required by the KaKs_Calculator program
AXTConvertor ${dir}_cds_pal2nal.aln ${dir}_input.axt
# Calculate the Ka/Ks ratio for the gene using the KaKs_Calculator program
KaKs_Calculator -i ${dir}_input.axt -o ${dir}.result -c 11 -m YN 
# Extract the relevant information from the output file and save it to ${dir}.result.list
cat ${dir}.result | cut -f 1,2,3,4,5,6 | grep -v 'Sequence' | less -S  > ${dir}.result.list

