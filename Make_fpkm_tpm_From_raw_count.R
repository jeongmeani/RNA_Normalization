## THIS script is to calculate FPKM and TPM
## AUTHOR = KJM
## Calculation formular of this refers to TPM -> https://statquest.org/rpkm-fpkm-and-tpm-clearly-explained/

#### FPKM ####
# 1. Count up the total reads in a sample and divide that number by 1,000,000 – this is our “per million” scaling factor.
# 2. Divide the read counts by the “per million” scaling factor. This normalizes for sequencing depth, giving you reads per million (RPM)
# 3. Divide the RPM values by the length of the gene, in kilobases(length/1,000). This gives you RPKM. RPKM = FPKM

#### TPM ####

# 1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
# 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# 3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.

library(dplyr)


args <- commandArgs(trailingOnly=TRUE)

count_matrix <- read.table(args[1], header=T)
gene_length <- read.table(args[2], header=T)
annotation<-read.table('/data/ref/hg19/annotation.hg19_for_Normal.csv',header=T, sep='\t',quote="",row.names = NULL)
head(annotation)

join_matrix=full_join(count_matrix, gene_length, by = c("Gene_name" = "Gene_name"))



Make_FPKM <-function(matrix_join){
	col_len=length(colnames(matrix_join))-2
	for (x in 2:col_len){
		S_name=paste(colnames(matrix_join)[x],'.FPKM', sep="")
		read_sum_per_sample=sum(matrix_join[,x])/1000000000 # THIS IS Step 1 of FPKM
#		matrix_join$FPKM=matrix_join[,x]/read_sum_per_sample # this version 1 to step2 of FPKM
		matrix_join[,S_name]=matrix_join[,x]/read_sum_per_sample # this version 1 to step2 of FPKM
		matrix_join[,S_name]=matrix_join[,S_name] / matrix_join$Length
}
	return(list("matrix" = matrix_join, "col_len"=col_len))
}

Make_TPM <- function(matrix_join, col_len){
	for (x in 2:col_len){
		S_name=paste(colnames(matrix_join)[x],'.TPM', sep="")
		matrix_join[,S_name]=matrix_join[,x] / matrix_join$Length # Step 1 of TPM => Divide the read counts by the length of each gene in kilobases.
		RPK_value_sum=sum(matrix_join[,S_name])/1000000 # Step 2 => Count up all the RPK values in a sample and divide this number by 1,000,000
		matrix_join[,S_name]=matrix_join[,S_name]/RPK_value_sum #Step 3 => Divide the RPK values by the “per million” scaling factor.
}
	matrix_join2=left_join(matrix_join,annotation,by='Gene_name')
	write.table(matrix_join2, file='Normalization.tsv', sep="\t", quote=FALSE, row.names=FALSE)
}

matrix_FPKM_and_len_list<-Make_FPKM(join_matrix)
Make_TPM(matrix_FPKM_and_len_list$matrix, matrix_FPKM_and_len_list$col_len)
