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
 ls | while read id; do cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}; sed 's/.\([0-9]\)\.p//g' -i ${id}.list; seqkit grep -f ${id}.list > ${id}_cds.fa; done
#change seqeuence header, same as the species name
sed '/^>/s/Aedg.*/Adg/; /^>/s/Ae.*/Aeve/; /^>/s/Lj.*/Ljap/; /^>/s/Tp.*/Tpra/; /^>/s/EVM.*/Smar/; /^>/s/Pv5.*/Pvul/; /^>/s/Adg.*/Aedg/' ${id}_cds.fa > ${id}_cds_clean.fa
#
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${id}_cds_clean.fa > ${id}_cds.tmp.fa
#
sed -e "s/TGA$//" -e "s/TAA$//" -e "s/TAG$//" ${id}_cds.tmp.fa > ${id}_no_stop_codon.fa
#Translate nucleic acid sequences to aa sequneces
transeq - -sequence ${id}_no_stop_codon.fa -outseq ${id}_protein.fa
mafft --quiet --localpair --maxiterate 1000 --thread 40 ${id}_protein.fa >  ${id}_protein.aln.fa
pal2nal.pl ${id}_protein.aln.fa ${id}_no_stop_codon.fa -nogap -codontable 11 -output paml > cds.pal2nal
#find the wrong og 
for dir in  ~/data/selection_analysis/PSGanalysis/single_copy/*; do
    if [ -d "$dir" ]; then
        if [ -e "$dir/cds.pal2nal" ]; then
            if [ ! -s "$dir/cds.pal2nal" ]; then
                echo "$dir"
            fi
        fi
    fi
done | cut -d "/" -f 9 | sort | less
#use branch site model in the paml 
mkdir alrt_model null_model 
#creat codeml.ctl files (codeml.ctl and codeml_alrt.ctl)
#creat tree files 
ls | grep -v "script.sh" | grep -v "nohup.out" | while read id
do
        cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/alrt_model
        codeml codeml_alrt.ctl   
         grep -E "lnL|foreground w" mlc  | awk '{print$5}' | tr '\n' ' '  |  xargs echo -n ${id} >> ../${id}_result.out
        cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/null_model
        grep -E  "lnL|foreground w" mlc | awk '{print$5}' | tr '\n' ' ' | xargs echo " " >> ../${id}_result.out
        codeml codeml.ctl
done

cat *${id}_result.out >> all_result.out
awk '{print$0,$2-$4}' all_result.out | awk '{print$0,sqrt($6*$6)}' | awk '{$6=1;print$0}' | awk '{print$0,($7*2)}' >> all_result.tmp.out
less all_result.tmp.out | awk '{print$6,$8}' | while read id; do chi2 ${id}; done | grep -v ^$ | awk '{print$6}' > chi2_result.out
paste -d " " all_result.tmp.out chi2_result.out > final_paml.out

ls | grep -v "sh" | grep -v "out" | while read id
do
        cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/branch_model/M0_model
        grep "lnL" mlc | awk  -F : '{print$4}' | awk '{print$1}' | xargs echo -n ${id} >> ../${id}_branch_model_result.out
        grep "omega" mlc | awk '{print$4}' | xargs echo -n  " " >> ../${id}_branch_model_result.out
        cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/branch_model/M2_model
        grep "lnL" mlc | awk  -F : '{print$4}' | awk '{print$1}'  | xargs echo -n " " >> ../${id}_branch_model_result.out
        grep "(dN/dS)" mlc | awk '{print"w1="$5";""w2="$6}' | tr '\n' ' '  | xargs echo " " >> ../${id}_branch_model_result.out
done
cat *${id}_result.out >> all_result.out
awk '{print$0,$2-$4}' all_branch_model_result.out | awk '{print$0,sqrt($6*$6)}'  | awk '{print$0,($7*2)}' | awk '{$9=1;print$0}' >> all_result.tmp.out
less all_result.tmp.out | awk '{print$9,$8}' | while read id; do chi2 ${id}; done | grep -v ^$ | awk '{print$6}' > chi2_result.out
paste -d " " all_result.tmp.out chi2_result.out > final_paml.out
        
      
 



