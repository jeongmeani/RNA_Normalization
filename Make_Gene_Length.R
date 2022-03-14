#!/usr/bin/Rscript
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

# it refers to https://github.com/dpryan79/Answers/blob/master/SEQanswers_42420/GTF2LengthGC.R
## elementLengths() was deprecated some time ago (2016 I think) and replaced with elementNROWS(). The generic is in the S4Vectors package.
## So Line 17 -> modify elementLengths(grl) to elementNROWS(grl)

GTFfile = "/data/ref/gencode/v24/gencode.v15.annotation.gtf"
FASTAfile = "./hg19.fa"

#Load the annotation and reduce it
GTF <- import.gff(GTFfile, format="gtf", genome="hg19", feature.type="exon")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

#Open the fasta file
FASTA <- FaFile(FASTAfile)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
output <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(output) <- c("Length", "GC")

write.table(output, file="GC_lengths.tsv", sep="\t", quote = FALSE)
