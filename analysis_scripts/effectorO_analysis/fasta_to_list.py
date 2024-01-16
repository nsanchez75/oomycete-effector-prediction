import argparse
import os

def main():
  parser = argparse.ArgumentParser(prog='fasta_to_list.py', description="A python script that produces a list of all the sequence header names in a FASTA file.")
  parser.add_argument("--input", "-i", type=str, required=True, help="Input sequence for FASTA file.")
  args = parser.parse_args()
  
  # check if input exists
  FASTA_FILE = str(args.input)
  if not os.path.isfile(FASTA_FILE):
    exit(f"{FASTA_FILE} either is a file or does not exist.")

  # create output file
  ## determine name of FASTA_file
  after_dirnames = FASTA_FILE.rfind('/')
  ## add filterlist indicator as an extension so it gets removed in analyze_effectoro.r
  OUTFILENAME = FASTA_FILE[after_dirnames + 1:].split('.')[0] + ".filterlist.txt"

  # write inputs to output
  with open(FASTA_FILE, 'r') as fin, open(OUTFILENAME, 'w') as fout:
    for line in fin:
      if line[0] == '>':
        fout.write(line[1:].strip().split()[0] + '\n')


if __name__ == "__main__":
  main()
