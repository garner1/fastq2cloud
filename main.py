from module_2.mod_create_sentences import *
import subprocess
import os
# # Define command and arguments
# command = 'Rscript'
# path2script = '/home/garner1/Dropbox/pipelines/fastq2cloud/module_1/MC_model_from_fasta.R'

# # Build subprocess command
# cmd = [command, path2script]

# # check_output will run the command and store to result
# x = subprocess.check_output(cmd, universal_newlines=True)

# print('The Markov-Chain model has been created')

#################################

# bashCommand = "bash module_2/parse_fastq.sh /home/garner1/Work/dataset/bliss/BB23_R1.fastq.gz"
# import subprocess
# process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
# output, error = process.communicate()

bashCommand = "parallel ./module_2/create_sentences.py {} /home/garner1/Work/dataset/fastq2cloud/transitionMatrix_3to3.csv ::: `ls /home/garner1/Work/dataset/fastq2cloud/reads_*`"
import subprocess
process = subprocess.Popen(bashCommand.split(), shell=True,stdout=subprocess.PIPE)
output, error = process.communicate()

# create_corpus('/home/garner1/Work/dataset/fastq2cloud/reads_aa','/home/garner1/Work/dataset/fastq2cloud/transitionMatrix_3to3.csv')
