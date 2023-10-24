# EffectorO Updates

## Running EffectorO

Created a conda environment `effectoro.yml` that can be used to run EffectorO

- Code to create conda environment: `conda env create -f [path/to/]effectoro.yml`

Ran EffectorO on `B_lac-SF5.protein.fasta` to test it and it works

## Running EffectorO Before SignalP

Created pipelines that would run EffectorO first and then SignalP and vice versa (probably not very useful but oh well)

## Testing EffectorO Effectiveness

Developed an R script that provides analysis on the effector sequences (probability and sequence length distributions). These scripts can be found in `scripts/EffectorO_analysis`.

Update 10/24/2023:

- added a feature to produce fasta file of seqs filtered by probability
- used optparse for argument parsing (makes it that `analyze_effectoro.r` to require the optparse package)

TODO: **Check in with Kelsey to see if I'm doing things properly (pipelines are in my frappe stuff)**
