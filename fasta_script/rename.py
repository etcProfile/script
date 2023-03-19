from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Rename headers in a FASTA file.')
parser.add_argument('-in', dest='input_file', required=True, help='Input FASTA file')
parser.add_argument('-out', dest='output_file', required=True, help='Output FASTA file')
parser.add_argument('-rename', dest='rename_file', required=True, help='Rename list file. Each line should contain an old header and a new header separated by a space.')
parser.add_argument('-help', dest='help', action='store_true', help='Print usage')
args = parser.parse_args()

if args.help:
    print('Author: Luo Shaofan; if you have any questions please send email to lsf763557959@gmail.com')
    print('Usage: python script.py -in input.fasta -out output.fasta -rename rename.list')
    print('rename.list format: each line should contain an old header and a new header separated by a space.')
else:
    # Read the input FASTA file
    fasta_sequences = SeqIO.parse(open(args.input_file),'fasta')

    # Read the rename list file
    rename_dict = {}
    with open(args.rename_file, 'r') as f:
        for line in f:
            old_header, new_header = line.strip().split()
            if old_header in rename_dict:
                print(f"Error: {old_header} is duplicated in rename.list")
            rename_dict[old_header] = new_header

    # Check if the headers in the FASTA file are unique
    fasta_headers = set()
    for fasta in fasta_sequences:
        if fasta.id in fasta_headers:
            print(f"Error: {fasta.id} is duplicated in {fasta.id}")
        fasta_headers.add(fasta.id)

    # Iterate over each sequence in the FASTA file and replace the header of each sequence
    with open(args.output_file, 'w') as f:
        for fasta in SeqIO.parse(open(args.input_file),'fasta'):
            old_header = fasta.id
            if old_header in rename_dict:
                fasta.id = rename_dict[old_header]
                fasta.description = ''
            SeqIO.write(fasta, f, "fasta")
