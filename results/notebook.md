# EffectorO Updates

## Running EffectorO

Created a conda environment `effectoro.yml` that can be used to run EffectorO. Here is the code to create conda environment: `conda env create -f [path/to/]effectoro.yml`

Ran EffectorO on `B_lac-SF5.protein.fasta` to test it and it works.

## Running EffectorO and SignalP

Created pipelines that would run EffectorO first and then SignalP and vice versa.

## Testing EffectorO Effectiveness

### Analyzing ORFs through EffectorO directly (w/o going through signalP first)

Developed an R script that provides analysis on the effector sequences (probability and sequence length distributions). These scripts can be found in `scripts/EffectorO_analysis`.

(10/24/2023)

- added a feature to produce fasta file of seqs filtered by probability
- used optparse for argument parsing (makes it that `analyze_effectoro.r` to require the optparse package)

(10/25/2023)

- simplified the R script by adding functions that remove redundant code
  - note: histogram creation function includes a bin_count parameter so made numbers for prob. and seq_length 22 and 20 respectively (22 was chosen just b/c Bremia's histograms by default made this & 20 was chosen b/c the sequence lengths were very much right-skewed)
  - create_table now uses "data" for left-hand column title instead of "probabilities" and "length_of_seqs"

(10/26/2023)

- removed function create_histogram because histograms of probabilities and sequence lengths have different parameters
- implemented option to input a file that lists sequence names a user would want to focus on in their inputted FASTA file
- made it so the probabilities histogram uses domain from [0.0 - 1.0]
- R script creates new directory if dealing with a list of sequence names to filter (all produced files will go into the newly made directory)

(10/27/2023)

- updated analysis_effector.r to log any sequence names not found in the input FASTA file

### Analyzing Effector Probability on Confirmed Effectors & Non-Effectors

(10/25/2023)

Downloadined datasets from Munir's directory. So far me and Kelsey have only found WY and CRN files displaying what sequences are classified as each.

TODO: figure out where Munir created these secreted proteins (maybe reread Kelsey's paper)

Datasets:

- Validated Effectors
  - WY Domain
  - CRN
  - RXLR (FIND)
  - Signal Peptides (FIND / maybe just use signalP on EffectorO output which I had already done?)
  - RNA-seq expressed (TODO: figure this out after above ones are done)
  - Lineage Specificity (TODO: figure this out after above ones are done)
- Validated Non-Effectors
  - TODO: figure this part out later

Update 10/26/2023: I analyzed WY and CRN with my revamped R script. I found that all of the sequences indicated by the filter list files can be found in the EffectorO-produced FASTA file. However, I don't notice any correlation between the fact that they are known effectors and their probabilities.

Update 10/30/2023: I fixed the R script to properly log what sequences from the filter list were not found in the input FASTA file. I found that there were some sequences from the known effectors lists that were not found. Thus, I need to talk to Kelsey about the effectiveness of EffectorO.
