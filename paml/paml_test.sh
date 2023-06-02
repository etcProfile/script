```
#!/bin/bash
# This script was written by Luo Shaofan on May 19, 2022
# Define variables
seagrass_path="/home/lsf/seagrass/chloroplast/New_cp/kaks"
paml_path="/home/lsf/seagrass/chloroplast/New_cp/kaks/branch_model"

# List all directories
ls -l | grep ^d | awk '{print$9}' | while read id
do
# Merge fasta files
#perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' ${gene}.fasta > ${gene}.tmp.fasta
#sed -e "s/TGA$//" -e "s/TAA$//" -e "s/TAG$//" ${gene}.tmp.fasta > ${id}_${gene}_no_stop_codon.fasta
#transeq -sequence ${id}_${gene}_no_stop_codon.fasta -outseq ${id}_${gene}_protein.fasta
#sed -i "s/\_1//g" ${id}_${gene}_protein.fasta

concatenate_fasta.py *no_stop_codon.fasta --separate -o all_${id}_cds.fasta 
concatenate_fasta.py *protein.fasta --separate -o all_${id}_pep.fasta 
# Perform multiple sequence alignment
mafft --retree 2 --maxiterate 1000 --thread 40 all_${id}_pep.fasta > all_${id}_pep.aln.fasta 
#pal2nal.pl all_${id}_pep.aln.fasta all_${id}_cds.fasta -nogap -codontable 11 -output paml > cds.pal2nal
# Create folders
mkdir -p ${paml_path}/alrt_model
mkdir -p ${paml_path}/null_model
# Copy files
cp ./*.pal2nal ${paml_path}/alrt_model
cp ./*.pal2nal ${paml_path}/null_model
cp ${seagrass_path}/all_cds.aln.fasta.treefile ${paml_path}/null_model
cp ${seagrass_path}/all_cds.aln.fasta.treefile ${paml_path}/alrt_model
cp ${seagrass_path}/codeml.ctl ${paml_path}/null_model
cp ${seagrass_path}/codeml_alrt.ctl ${paml_path}/alrt_model
# Enter folders
cd ${paml_path}/${id}/paml/null_model
# Run codeml
codeml codeml.ctl
# Enter folders
cd ${paml_path}/${id}/paml/alrt_model
# Run codeml
codeml codeml_alrt.ctl
done
```
