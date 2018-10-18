#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

library(aphid)
library(seqinr)

HmmFile <- args[1]
FaFile <- args[2]
OutputPrefix <- args[3]

print('1) Reading HMM')
hmm_HD <- readPHMM(HmmFile)
print('2) Reading Seqs')
Seqs <- read.fasta(FaFile, set.attributes = FALSE, seqtype = 'AA')

#Function to Return Viterbi Matches
returnViterbiResult <- function(seq, hmm, alntype){
  v <- Viterbi(x = hmm, y = seq, type = alntype)
  #Parse Result to '|' delimited PFam Start | Seq Start | Path
  paste(v$start[1], v$start[2], paste(v$path, collapse = ','), sep = '|')
}
#returnViterbiResult(l[[1]], hmm_HD, 'local')

#Output Match Data
writeLikeFasta <- function(listlike, outputname){
  o <- c()
  for(name in names(listlike)){
    val <- listlike[[name]]
    header <- paste('>', name, sep = '')
    o <- c(o, header, val)
  }
  fileConn<-file(outputname)
  writeLines(c(o), fileConn)
  close(fileConn)
}

#Find Viterbi Matches
print('3) Finding Viterbi Matches')
print('Starting LOCAL search...')
hmm_matches <- lapply(Seqs, returnViterbiResult, hmm_HD, 'local')
writeLikeFasta(hmm_matches, paste(OutputPrefix, '.Viterbi_local', sep = ''))
print('!Done LOCAL search...')

print('Starting GLOBAL search...')
hmm_matches <- lapply(Seqs, returnViterbiResult, hmm_HD, 'global')
writeLikeFasta(hmm_matches, paste(OutputPrefix, '.Viterbi_global', sep = ''))
print('!Done GLOBAL search...')

print('Starting SEMIGLOBAL search...')
hmm_matches <- lapply(Seqs, returnViterbiResult, hmm_HD, 'semiglobal')
writeLikeFasta(hmm_matches, paste(OutputPrefix, '.Viterbi_semiglobal', sep = ''))
print('!Done SEMIGLOBAL search...')
