#序列标准化
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${id}_cds.fa > ${id}_cds.tmp.fa
#删除终止密码子
sed -e "s/TGA$//" -e "s/TAA$//" -e "s/TAG$//" ${id}_cds.tmp.fa > ${id}_no_stop_codon.fa
#核酸序列翻译为蛋白序列
transeq -sequence ${id}_no_stop_codon.fa -outseq ${id}_protein.fa
#去掉header里面冗余的部分，使得核酸和蛋白的序列header一样
sed 's/_1//g' -i ${id}_protein.fa
#用mafft来序列对齐
mafft --quiet --localpair --maxiterate 1000 --thread 40 ${id}_protein.fa >  ${id}_protein.aln.fa
#用pal2nal 蛋白对齐文件根据密码子表来转换为核酸对齐的文件
pal2nal.pl ${id}_protein.aln.fa ${id}_no_stop_codon.fa -nogap -codontable 11 -output paml > cds.pal2nal
#找出哪些输出结果是有问题的
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
cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/alrt_model
odeml codeml_alrt.ctl   
grep "lnL" mlc | awk '{print$5}' | tr '\n' ' ' | xargs echo -n ${id} >> ../${id}_result.out
cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/null_model
codeml codeml.ctl
grep "lnL" mlc | awk '{print$5}' | tr '\n' ' ' | xargs echo " " >> ../${id}_result.out
cat *${id}_result.out >> all_result.out
awk '{print$0,$2-$4}' all_result.out | awk '{print$0,sqrt($6*$6)}' | awk '{$6=1;print$0}' | awk '{print$0,($7*2)}' >> all_result.tmp.out
less all_result.tmp.out | awk '{print$6,$8}' | while read id; do chi2 ${id}; done | grep -v ^$ | awk '{print$6}' > chi2_result.out
paste -d " " all_result.tmp.out chi2_result.out > final_paml.out
#branch model 
cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/branch_model/M0_model
grep "lnL" mlc | awk  -F : '{print$4}' | awk '{print$1}' | xargs echo -n ${id} >> ../${id}_branch_model_result.out
grep "omega" mlc | awk '{print$4}' | xargs echo -n  " " >> ../${id}_branch_model_result.out
cd ~/data/selection_analysis/PSGanalysis/single_copy/${id}/branch_model/M2_model
grep "lnL" mlc | awk  -F : '{print$4}' | awk '{print$1}'  | xargs echo -n " " >> ../${id}_branch_model_result.out
grep "(dN/dS)" mlc | awk '{print"w1="$5";""w2="$6}' | tr '\n' ' '  | xargs echo " " >> ../${id}_branch_model_result.out

cat ${id}_branch_model_result.out >> all_branch_model_result.out
awk '{print$0,$2-$4}' all_branch_model_result.out | awk '{print$0,sqrt($6*$6)}'  | awk '{print$0,($7*2)}' | awk '{$9=1;print$0}' >> all_result.tmp.out
less all_result.tmp.out | awk '{print$9,$8}' | while read id; do chi2 ${id}; done | grep -v ^$ | awk '{print$6}' > chi2_result.out
paste -d " " all_result.tmp.out chi2_result.out > final_paml.out
