#!/usr/bin/env Rscript

##################################################################################
#
# INFO:
#   This R script is used to analyze the output of the EffectorO script. This
#   will serve as a benchmark to test the efficiency of later EffectorO versions.
#
#   Warning: This script requires the "optparse" packages to be installed.
#
# Inputs:
#   - FASTA file from EffectorO output (required)
#   - minimum probability cutoff (optional)
#   - maximum probability cutoff (optional)
#
#   Note: if min & max cutoff are not provided then a filtered FASTA file will
#         not be created.
#
# Outputs:
#   - PNGs of histograms for probabilities and sequence lengths
#   - frequency table files for probabilities and sequence lengths
#   - if applicable, filtered FASTA file
#
##################################################################################

require("optparse")

option_list <- list(
  make_option(
    c("--input", "-i"),
    type = "character",
    help = "Input FASTA file path (required)."
  ),
  make_option(
    c("--min_cutoff", "-L"),
    type = "numeric",
    help = "Minimum effector probability cutoff.",
  ),
  make_option(
    c("--max_cutoff", "-H"),
    type = "numeric",
    help = "Maximum effector probability cutoff."
  )
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Input fasta file required.")
}

# get FASTA file
path <- getwd()
fasta_file <- readLines(paste(path, "/", opt$input, sep=""))

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

# create file consisting of only sequences within cutoff parameters
min_cutoff <- opt$min_cutoff
max_cutoff <- opt$max_cutoff

## don't create file if cutoffs aren't specified
if (is.null(min_cutoff) && is.null(max_cutoff))
  q(status=0)

filtered_fasta_file <- "predicted_effectors_filtered"
## apply min cutoff if applicable
if (!is.null(min_cutoff)) {
  qualifying_probability_indices <- which(probabilities >= min_cutoff)
  filtered_fasta_file <- paste(filtered_fasta_file, "_L", min_cutoff, sep='')
}
## apply max cutoff if applicable
if (!is.null(max_cutoff)) {
  qualifying_probability_indices <- which(probabilities <= max_cutoff)
  filtered_fasta_file <- paste(filtered_fasta_file, "_H", max_cutoff, sep='')
}
filtered_fasta_file <- paste(filtered_fasta_file, ".fasta", sep='')

## create file
lines_to_write <- c()
for (line in qualifying_probability_indices * 2)
  lines_to_write <- c(lines_to_write, fasta_file[line - 1], fasta_file[line])
writeLines(lines_to_write, filtered_fasta_file)
