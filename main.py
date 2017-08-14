# run_max.py
import subprocess

# Define command and arguments
command = 'Rscript'
path2script = '/home/garner1/Dropbox/pipelines/fastq2cloud/module_1/MC_model_from_fasta.R'

# Build subprocess command
cmd = [command, path2script]

# check_output will run the command and store to result
x = subprocess.check_output(cmd, universal_newlines=True)

print('The Markov-Chain model has been created')
