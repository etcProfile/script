```bash
#!/bin/bash

# This script was written by Luo Shaofan on May 12, 2022.
# List all directories in the current directory and loop over them
ls | while read dir 
do
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
done
```
