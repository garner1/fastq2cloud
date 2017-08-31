library(Biostrings)
library(data.table)  

# Set the dimensions of the MC model 
in_dim = 3
out_dim = 3
  
# Specify the input fastq file
# filename = "/home/garner1/Work/dataset/fastq2cloud/test1K.fastq.gz"
args = commandArgs(trailingOnly=TRUE)
seq = readDNAStringSet(args[1], format="fastq")

# Build the MC homogeneous transition matrix
Tmat = oligonucleotideTransitions(seq, left = in_dim, right = out_dim, as.prob = FALSE)

# Define the information content
Tmat_norm <- -log(Tmat/rowSums(Tmat))

# Build the dataframe associating to a word, obtained joining input and output state, the information
key <- c()
value <- c()
i <- 0
for(row in 1:4**in_dim){
  for(col in 1:4**out_dim){
    i <- i+1
    key[[i]] <- paste(rownames(Tmat)[row],colnames(Tmat)[col],sep = '')
    value[[i]] <- Tmat_norm[row,col]
  }
}

df <- data.frame(key,value) #THERE MIGHT BE INFINITE VALUES DUE TO DIVISION BY 0
filename <- paste("/media/DS2415_seq/silvanog/Genome_segmentation/transitionMatrix_fromFastq_",as.character(in_dim),"to",as.character(out_dim),".csv", sep = '')
fwrite(x=df, file=filename, append=FALSE, row.names=TRUE)



