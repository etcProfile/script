#remove genefamily which is 0 in the target species 
cat Orthogroups.GeneCount.tsv | awk -F '\t' '$5>0{print$0}' > Orthogroups.GeneCount_smar0.tsv
#remove genefamily which is <0 in the other species 
awk '{count=0; for(i=1;i<=NF;i++) if($i=="0" && length($i)==1) count++; if(count=3) print $0}' Orthogroups.GeneCount_smar0.tsv
#get single copy genefamily in the Orthogroups 
awk '{for(i=2;i<=NF;i++)if($i!=1)next}1' Orthogroups.GeneCount.tsv  > single_copy.tsv
#get single copy genefamily list
cut -f 1 single_copy.tsv > single_copy.list
 #clean up the genome annotation data (make same header of one gene between protein.fasta and cds.fasta)
 #remove the context after the first space of the header 
sed '/^>/s/^>\([^ ]*\) .*/>\1 /' Tpra.cds.fa 
#remove the context after the second point (".") of the header 
sed '/^>/s/\.[^\.]*\.[^\.]*$//' Pvul_protein.fa
#get all nucl sequence
cat Aedg.cds.fa Aeve.cds.fa Ljap.cds.fa Pvul.cds.fa Smar_cds.fa Tpra.cds.fa >> all_cds.fa
 ls | while read id; do cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}; sed 's/.\([0-9]\
)\.p//g' -i ${id}.list; seqkit grep -f ${id}.list > ${id}_cds.fa; done
#change seqeuence header, same as the species name
sed '/^>/s/Aedg.*/Adg/; /^>/s/Ae.*/Aeve/; /^>/s/Lj.*/Ljap/; /^>/s/Tp.*/Tpra/; /^>/s/EVM.*/Smar/; /^>/s/Pv5.*/Pvul/; /^>/s/Adg.*/Aedg/' ${id}_cds.fa > ${id}_cds_clean.fa
#Translate nucleic acid sequences to aa sequneces
trans 

