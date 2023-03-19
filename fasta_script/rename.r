input_file <- 'input.fasta'
output_file <- 'output.fasta'

fasta <- readLines(input_file)
new_fasta <- character(length(fasta))

for (i in 1:length(fasta)) {
  if (substr(fasta[i], 1, 1) == '>') {
    header <- unlist(strsplit(fasta[i], split='_'))
    new_header <- paste(rev(header), collapse='_')
    new_fasta[i] <- paste0('>', new_header)
  } else {
    new_fasta[i] <- fasta[i]
  }
}

writeLines(new_fasta, output_file)
