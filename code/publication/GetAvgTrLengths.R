#Here, we download and calculate average transcript lengths for the mouse genes, which is used in
#gene length normalization

#
#BiocManager::install("DESeq2")
library(tidyverse)

dataFolderModels = "../data/Models/"


#get average transcript length for each gene
#####################################

#BiocManager::install("GenomicFeatures")
library(GenomicFeatures)
library(biomaRt)

ensembl_us_west = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="uswest.ensembl.org")
listDatasets(ensembl_us_west)
ensembl = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")

geneConvTable <- getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart = ensembl)#fails sometimes, just try again!

listMarts() #gave us version 100

#Annoyingly, this fails somtimes, just try again, and again, and again... Sometimes you need to try 20 times....
#txdbM = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl", transcript_ids=NULL, circ_seqs=DEFAULT_CIRC_SEQS, filter=NULL, id_prefix="ensembl_", host="www.ensembl.org", port=80, taxonomyId=NA, miRBaseBuild=NA)
txdbM = makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

tlM = transcriptLengths(txdbM, with.cds_len=FALSE)

#

#convert the gene in transcript lengths to symbol id

tlmGenes = tlM$gene_id
ind = match(tlmGenes,geneConvTable$ensembl_gene_id)
newTlmGenes = geneConvTable$mgi_symbol[ind]

trLengths = as_tibble(tlM)
trLengths = trLengths %>% add_column(gene_symbol=newTlmGenes)
#skip all lines with no gene symbol
trLengths = trLengths[trLengths$gene_symbol != "",]
#now get the mean of all genes:
avgTrLengths = trLengths %>% group_by(gene_symbol, gene_id) %>% summarize(txlen = mean(tx_len))
colnames(avgTrLengths) = c("geneSymbol","ensemblGeneId","trLength")

#write to file
save(avgTrLengths, file=paste0(dataFolderModels, "AvgTrLengths.RData"))

# save as TSV file including the column of Ensembl ids
write.table(avgTrLengths, "avgTrLenMouse.tsv", sep = "\t", dec = ".",
            quote = FALSE, row.names = F)



