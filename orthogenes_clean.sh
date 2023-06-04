#remove genefamily which is 0 in the target species 
cat Orthogroups.GeneCount.tsv | awk -F '\t' '$5>0{print$0}' > Orthogroups.GeneCount_smar0.tsv
#remove genefamily which is <0 in the other species 
awk '{count=0; for(i=1;i<=NF;i++) if($i=="0" && length($i)==1) count++; if(count=3) print $0}' Orthogroups.GeneCount_smar0.tsv
#get single copy genefamily in the Orthogroups 
awk '{for(i=2;i<=NF;i++)if($i!=1)next}1' Orthogroups.GeneCount.tsv  > single_copy.tsv
