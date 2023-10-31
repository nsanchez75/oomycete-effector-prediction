#!/usr/bin/env Rscript

##################################################################################
#
# INFO:
#   This R script is used to analyze the output of the EffectorO script. This
#   will serve as a benchmark to test the efficiency of later EffectorO versions.
#
#   Warning: This script requires the "optparse" packages to be installed.
#
#
# Inputs:
#   - FASTA file from EffectorO output (required)
#   - input file that lists sequences that will be focused on (optional)
#   - minimum probability cutoff (optional)
#   - maximum probability cutoff (optional)
#
#   Note: if either min or max cutoff or input file are not provided then a filtered FASTA file will
#         not be created.
#
#
# Outputs:
#   - PNGs of histograms for probabilities and sequence lengths
#   - frequency table files for probabilities and sequence lengths
#   - if applicable, filtered FASTA file
#
#   Note: if a filtered FASTA file is given, a directory based on the file's name is created.
#
##################################################################################

require("optparse")


# some defined functions:

create_freq_files <- function(data, out_filename) {
  freq_table <- table(data)
  sorted_freq_table <- freq_table[order(-freq_table)]
  write.table(sorted_freq_table, file=out_filename, sep='\t', row.names=FALSE)
}

is_seq_in_filter_list <- function(name, f_list) {
  for (valid_name in f_list)
    if (name == valid_name) return(TRUE)
  return(FALSE)
}


option_list <- list(
  make_option(
    c("--input", "-i"),
    type = "character",
    help = "Input FASTA file path (required)."
  ),
  make_option(
    c("--filter_list", "-f"),
    type = "character",
    help = "Input file that lists names of sequences to be focused on."
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
if (!file.exists(opt$input))
  stop(paste("Input FASTA file", opt$input, "does not exist."))
fasta_file <- readLines(opt$input)

# apply filtered list if applicable
filter_list_check <- !is.null(opt$filter_list)
if (filter_list_check) {
  if (!file.exists(opt$filter_list))
    stop(paste("Filter list file", opt$filter_list, "does not exist."))
  filter_list <- readLines(opt$filter_list)

  # create directory for filter list file
  filter_list_filename <- gsub("[.][a-z,A-Z,0-9]+", '', opt$filter_list)
  filter_list_filename <- paste(filter_list_filename, "_analysis_results", sep='')
  dir.create(filter_list_filename)
  # switch to new directory for all outputs
  setwd(filter_list_filename)
}

# init vectors
probabilities <- c()
length_of_seqs <- c()
filter_list_lines <- c()

for (line in seq(1, length(fasta_file), by=2)) {
  seq_header <- unlist(strsplit(fasta_file[line], " "))
  seq_name <- sub('>', '', seq_header[1])
  if (!filter_list_check || is_seq_in_filter_list(seq_name, filter_list)) {
    # remove found name from list
    filter_list <- filter_list[filter_list != seq_name]

    # extract probability values from headers
    probability_val <- as.numeric(unlist(strsplit(seq_header[3], "="))[2])
    probabilities <- c(probabilities, probability_val)
    # determine sequence length
    length_of_seqs <- c(length_of_seqs, nchar(fasta_file[line + 1]) - 1) # doesn't count '*' at end of sequences

    if (filter_list_check)
      filter_list_lines <- c(filter_list_lines, fasta_file[line], fasta_file[line + 1]) # used to create FASTA file for filter list
  }
}

# warn user if some filter names were not found in the file
if (filter_list_check) {
  warning_file <- file(paste(filter_list_filename, "_warnings.log", sep=''), 'w')
  warnings_list <- c()
  if (length(filter_list) > 0) {
    warnings_list <- c(warnings_list, "Warning: these sequence names are not found in the input FASTA file:")
    for (name in filter_list)
      warnings_list <- c(warnings_list, paste('\t', name, sep=''))
  }

  for (line in warnings_list)
    print(sub("\\t", '', line))

  writeLines(warnings_list, warning_file)
  close(warning_file)
}


# create histograms
png("effector_probabilities_hist.png", width=800, height=400)
hist(probabilities, main="Histogram of Effector Probabilities", xlab="Probability", xlim=c(0.0,1.0), ylab="Frequency")
dev.off()
png("length_of_seqs_hist.png", width=800, height=400)
hist(length_of_seqs, main="length_of_seqs_hist.png", xlab="Sequence Lengths", ylab="Frequency", breaks=20)
dev.off()


# create frequency text files of both histograms
create_freq_files(probabilities, "effector_probabilities_frequencies.txt")
create_freq_files(length_of_seqs, "length_of_seqs_frequencies.txt")


# create subfile of FASTA file consisting of filter list names
if (filter_list_check)
  writeLines(filter_list_lines, paste(filter_list_filename, "_predicted_effectors.fasta", sep=''))


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
