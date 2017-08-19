
#THIS SCRIPT CALCULATE THE INFORMATION OF THE PROBABILITY OF THE TRANSITION MATRIX FROM A GIVEN LIST OF FASTA SEQUENCES

## try http:// if https:// URLs are not supported
#install.packages("XML",repos = "http://cran.us.r-project.org")
#source("https://bioconductor.org/biocLite.R")
#biocLite('BSgenome.Hsapiens.UCSC.hg19')
require('BSgenome.Hsapiens.UCSC.hg19')
library(Biostrings)
library(data.table)  

revstring <- function(s) paste(rev(strsplit(s,"")[[1]]),collapse="") # reverse the string
# revstring(compl_rows,2)

# chromosome names
cnames = seqnames(Hsapiens)

# loop of the dimensions of the markov chain, with in_dim from 3 to 8 and out_dim from 6 to 1, respectivelly 
# In the end you have the transition matrix summed over forward and reverse genome and all chromosomes (including the non-standard ones)
# The transition matrix is evaluated only from 5' to 3', since this is what you get in the fastq files (read1\2 are always from 5' to 3')
# The matrix for the reverse genome is obtained analytically from the forward matrix considering the transposition and rev-compl transformation needed to define it
in_dim = 3
out_dim = 3
Tmat = matrix(0, nrow = 4**in_dim, ncol = 4**out_dim)
for(name in cnames){
  print(name)
  seqChr = unmasked(Hsapiens[[name]]) # get the chromosome sequence
  Tmat_plus = oligonucleotideTransitions(seqChr, left = in_dim, right = out_dim, as.prob = FALSE)
  Tmat_minus = Tmat_plus
  for(row in rownames(Tmat)){
    for(col in colnames(Tmat)){
      newcol = revstring(chartr("ATGC","TACG",row)) # complement
      newrow = revstring(chartr("ATGC","TACG",col)) # complement
      Tmat_minus[newrow,newcol] = Tmat_plus[row,col]
      }
  }
  Tmat = Tmat + Tmat_plus + Tmat_minus
}

Tmat_norm <- -log(Tmat/rowSums(Tmat)) # evaluate the information of the probability

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
df <- data.frame(key,value)
filename <- paste("/media/DS2415_seq/silvanog/Genome_segmentation/transitionMatrix_",as.character(in_dim),"to",as.character(out_dim),".csv", sep = '')
fwrite(x=df, file=filename, append=FALSE, row.names=TRUE)
