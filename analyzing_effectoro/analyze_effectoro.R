#!/usr/bin/env Rscript

#
# INFO:
#   This R script is used to analyze the output of the EffectorO script. This will
#   serve as a benchmark to test the efficiency of later EffectorO versions.
#
# Input: FASTA file from EffectorO output
# Outputs:
#   - PNGs of histograms for probabilities and sequence lengths
#   - frequency table files for probabilities and sequence lengths
#

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    stop("A fasta file input is required: ./analyze_effectoro.R <FASTA_FILE> [<MIN_PROB_CUTOFF>]")
}

# get FASTA file
path <- getwd()
fasta_file <- readLines(paste(path, "/", args[1], sep=""))

probabilities <- c()
length_of_seqs <- c()

for (line in 1:length(fasta_file)) {
  # extract probability values from headers
  if (line %% 2 == 1) {
    probability_str <- unlist(strsplit(fasta_file[line], " "))[3]
    probability_val <- as.numeric(unlist(strsplit(probability_str, "="))[2])
    probabilities <- c(probabilities, probability_val)
  }
  # extract length of sequences
  else {
    length_of_seqs <- c(length_of_seqs, nchar(fasta_file[line]) - 1) # doesn't count '*' at end of sequences
  }
}

# create histogram for probabilities
png("effector_probabilities_hist.png", width=800, height=400)
hist(probabilities, main="Histogram of Effector Probabilities", xlab="Probability", ylab="Frequency")
dev.off()

# create histogram for length of sequences
png("length_of_seqs_hist.png", width=800, height=400)
hist(length_of_seqs, main="Histogram of Sequence Lengths", xlab="Sequence Length", ylab="Frequency", breaks=20)
dev.off()

# create frequency text files of both histograms 
probabilities_freq <- table(probabilities)
sorted_probabilities_freq <- probabilities_freq[order(-probabilities_freq)]
write.table(sorted_probabilities_freq, file="effector_probabilities_frequencies.txt", sep='\t', row.names=FALSE)
length_of_seqs_freq <- table(length_of_seqs)
sorted_length_of_seqs_freq <- length_of_seqs_freq[order(-length_of_seqs_freq)]
write.table(length_of_seqs_freq, file="length_of_seqs_frequencies.txt", sep='\t', row.names=FALSE)

# index and create file consisting of only sequences w/ prob > 0.8
if (length(args) >= 2) {
  min_cutoff <- as.numeric(args[2])
  if (is.na(min_cutoff))
    stop("Minimum probability cutoff must be a number.")

  file_filtered_effectors <- "predicted_effectors_filtered.fasta"
  qualifying_prob_indices <- which(probabilities >= min_cutoff)
  lines_to_write <- c()

  for (line in qualifying_prob_indices * 2)
    lines_to_write <- c(lines_to_write, fasta_file[line - 1], fasta_file[line])
  writeLines(lines_to_write, file_filtered_effectors)
}
