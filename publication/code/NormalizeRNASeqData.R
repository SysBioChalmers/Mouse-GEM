#Here we do the following:
#1. We gather all mouse model RNA-Seq data into common matrices (metadata and counts)
#2. We normalize the counts using the following procedure:
#	a. we calculate size factors for samples using DeSeq2
#	b. we divide the counts by average transcript length
#	c. we scale the counts of each sample with the size factor 
#	d. we scale all counts so a sample on average has 10^6 reads
#	The goal of this procedure is to get something similar to TPM (which is suitable for the modeling), 
#	but where the samples are normalized properly against each other.
#3. We save the matrices. The normalization is done on inner joined data (i.e. only using genes common to all datasets),
#   while the final data retains all genes.
#4. We group all samples belonging to the same condition and take the average expression for each condition. These new
#	matrices are saved as "grouped". The files are also saved as .mat files to simplify import into matlab.
#BiocManager::install("DESeq2")
library(DESeq2)
library(tidyverse)

dataFolderModels = "../data/Models"

#Help function that handles a list with 3 matrices: list(Metadata, fullJoinData, innerJoinData)
#The function adds an additional dataset to the matrices and returns a new list.
AddDataset <- function(dataList, pathCounts, pathMeta, dsid) {
  met = read_delim(paste0(dataFolderModels, pathMeta), delim=',')
  met = met %>% add_column(Dataset=dsid)
  #rename the sampleID column to "Sample ID" for some samples, so all are the same
  if (assertthat::has_name(met, "Sample ID")) {
    colnames(met)[which(colnames(met) == "Sample ID")] = "sampleID"
  }
  cnt = read_delim(paste0(dataFolderModels, pathCounts), delim=',') 
  
  if (length(dataList) == 0) { #first call
    return(list(met,cnt,cnt))
  } else {
    newMet = bind_rows(dataList[[1]], met)
    newCnt = full_join(dataList[[2]], cnt, by="genes");
    newCntInner = inner_join(dataList[[3]], cnt, by="genes");
    
    return(list(newMet,newCnt,newCntInner))
  }
   
}


###############################################
#load metadata and counts for all datasets, and build matrices with all data
###############################################

dataList = list();
#APP23
dataList = AddDataset(dataList, "APP23/GSE80465_counts.csv", "APP23/GSE80465_metadata.csv", "GSE80465")
dataList = AddDataset(dataList, "APP23/GSE104630_counts.csv", "APP23/GSE104630_metadata.csv", "GSE104630")
#APPPS1
dataList = AddDataset(dataList, "APPPS1/GSE63943_counts.csv", "APPPS1/GSE63943_metadata.csv", "GSE63943")
dataList = AddDataset(dataList, "APPPS1/GSE100070_counts.csv", "APPPS1/GSE100070_metadata.csv", "GSE100070")
dataList = AddDataset(dataList, "APPPS1/GSE110741_counts.csv", "APPPS1/GSE110741_metadata.csv", "GSE110741")
#APPswe-PSEN1dE9-C57BL6
dataList = AddDataset(dataList, "APPswe-PSEN1dE9-C57BL6/GSE131659_counts.csv", "APPswe-PSEN1dE9-C57BL6/GSE131659_metadata.csv", "GSE131659")
#APPswe-PSEN1dE9-C57BL6
dataList = AddDataset(dataList, "APPswe-PSEN1dE9-line85/GSE93678_counts.csv", "APPswe-PSEN1dE9-line85/GSE93678_metadata.csv", "GSE93678")
dataList = AddDataset(dataList, "APPswe-PSEN1dE9-line85/GSE104424_counts.csv", "APPswe-PSEN1dE9-line85/GSE104424_metadata.csv", "GSE104424")
dataList = AddDataset(dataList, "APPswe-PSEN1dE9-line85/GSE136861_counts.csv", "APPswe-PSEN1dE9-line85/GSE136861_metadata.csv", "GSE136861")
#Trem2KO-JAX
dataList = AddDataset(dataList, "Trem2KO-JAX/GSE124266_counts.csv", "Trem2KO-JAX/GSE124266_metadata.csv", "GSE124266")
dataList = AddDataset(dataList, "Trem2KO-JAX/GSE134031_counts.csv", "Trem2KO-JAX/GSE134031_metadata.csv", "GSE134031")
#Trem2KO-KOMP-APPPS1
dataList = AddDataset(dataList, "Trem2KO-KOMP-APPPS1/GSE104381_counts.csv", "Trem2KO-KOMP-APPPS1/GSE104381_metadata.csv", "GSE104381")

dim(dataList[[2]])#see how many genes we have
dim(dataList[[3]])#see how many genes we have
length(unique(dataList[[1]]$Dataset))#should be the number of datasets

#use the inner_join data for the normalization (to be fair across samples), and remove the gene column
countsForNorm = dataList[[3]][,-1]

#use DESeq2 to normalize. The design matrix is not used for normalization, so we don't care about it.
dds <- DESeqDataSetFromMatrix(countData = countsForNorm, colData = dataList[[1]], design = ~ 1)
dds <- estimateSizeFactors(dds)

sizeFactors = dds$sizeFactor

#now add transcript length to the data
load(file=paste0(dataFolderModels, "AvgTrLengths.RData"))

withTxLenFull = inner_join(avgTrLengths, dataList[[2]], by="genes")


#now we do three things: (we now work on the full_join matrix)
#1. we divide by average transcript length
#2. We scale each sample with the size factor 
#3. we scale all counts so a sample on average has 10^6 reads
tpmish = withTxLenFull

#1
tpmish[,c(-1,-2)] = tpmish[,c(-1,-2)]/tpmish$txLen #So, division divides each col with the vector, i.e. each row is divided by the same number

tpmish = tpmish[,-2] #we don't need the transcript length anymore

#test that the counts actually got divided by transcript length
x1 = tpmish[tpmish$genes == "Dlat",6, drop=T] #1.15
x2 = avgTrLengths$txLen[avgTrLengths$genes == "Dlat"]
xcmp = withTxLenFull[withTxLenFull$genes == "Dlat", 7]/x2
x1 == xcmp#ok

#2
forSum = tpmish[,-1]
forSum[is.na(forSum)] = 0
pseudoCountsPerSamp = colSums(forSum)

normPseudoCountsPerSamp = pseudoCountsPerSamp/sizeFactors 

#3  
globalScaleFactor = 10^6/mean(normPseudoCountsPerSamp)

finalCountsPerSamp = normPseudoCountsPerSamp * globalScaleFactor
#now scale the original full_join matrix containing NA

finalTPMish = tpmish
for (i in 1:length(finalCountsPerSamp)) {
  s = sum(finalTPMish[,i+1], na.rm=T)
  finalTPMish[,i+1] = finalTPMish[,i+1]*finalCountsPerSamp[i]/s
}

#test, should be 1M
mean(colSums(finalTPMish[,-1], na.rm=T)) #looks good!

#write files with all samples and all metadata
write_tsv(finalTPMish, paste0(dataFolderModels, "AllTPM.txt"))
write_tsv(dataList[[1]], paste0(dataFolderModels, "AllMeta.txt"))

metadata = dataList[[1]]


#pathMeta = "APP23/GSE104630_metadata.csv"
#pathCounts = "APP23/GSE104630_counts.csv"
#dsid = "GSE104630"

#pathMeta = "APP23/GSE80465_metadata.csv"
#pathCounts = "APP23/GSE80465_counts.csv"
#dsid = "GSE80465"

#############################################
#now group the data and calculate the mean expression for each gene within each group
#always group on dataset
#############################################

#datasets = unique(dataList[[1]]$Dataset)

#datasets = c("GSE80465",  
#             "GSE104630", 
#             "GSE63943",  
#             "GSE100070", 
#             "GSE110741", 
#             "GSE131659", 
#             "GSE93678",  
#             "GSE104424", 
#             "GSE136861", 
#             "GSE124266", 
#             "GSE134031", 
#             "GSE104381")

#merge_cols = c(c("Genotype", "Tissue", "Celltype", "Gender", "Age"), #GSE80465
#               c("Genotype", "Tissue", "Celltype", "Gender", "Age"))

groups = unique(metadata[,c(-1,-2)])
print(groups, n=10000)
write_tsv(groups, paste0(dataFolderModels, "GroupMetaData.txt"))


#filter away all but one dataset for now
#groups = groups[groups$Dataset == "GSE80465",]
#groups = groups[1:8,]

numGroups = dim(groups)[1]

#the easiest thing is to just loop
expr = matrix(0,nrow = dim(finalTPMish)[1], ncol = numGroups)

metdatcmp = paste(metadata$Dataset, metadata$Model, metadata$Genotype, metadata$Tissue, metadata$Celltype, metadata$Gender, metadata$Age)
groupcmp = paste(groups$Dataset, groups$Model, groups$Genotype, groups$Tissue, groups$Celltype, groups$Gender, groups$Age)

#calculate the mean expression for each 
for (i in 1:numGroups) {
  expr[,i] = rowMeans(finalTPMish[,c(F, metdatcmp == groupcmp[[i]]) ])#skip the gene column
}
#test
(25.676527337533116 + 23.353251711062786 + 27.020447520480584)/3#from the allTPMs file, the first 3 samples, 3rd gene
expr[3,1]#these should be the same, which they are
(4.20841269365805 + 4.746841768855952 + 5.32535397060764)/3#from the allTPMs file, sample 7-9, 5th gene
expr[5,3]#should be the same, which they are, good


groupTPMs = bind_cols(finalTPMish[,1], as_tibble(expr)) #gives a warning about column names, can be ignored, we're setting them below
colnames(groupTPMs)[-1] = groupcmp

write_tsv(groupTPMs, paste0(dataFolderModels, "TPMsPerCond.txt"))

#also write to matlab file and a text file with the genes only
write_tsv(groupTPMs[,1], path=paste0(dataFolderModels, "TPMsPerCondGenes.txt"), col_names = F)

library("R.matlab")

toExp = as.matrix(groupTPMs[,-1])

writeMat(paste0(dataFolderModels, "TPMsPerCond.mat"), TPMs = toExp, sampleNames = colnames(groupTPMs)[-1])

#also create 3 chunks of the models, to make it finish faster
chunk1 = groupTPMs[, c(1,2:41)]
chunk2 = groupTPMs[, c(1,42:81)]
chunk3 = groupTPMs[, c(1,82:(dim(groupTPMs)[2]))]

expData <- function(d, suffix) {
  toExp = as.matrix(d[,-1])
  
  writeMat(paste0(dataFolderModels, "TPMsPerCond", suffix, ".mat"), TPMs = toExp, sampleNames = colnames(d)[-1])
  
}

expData(chunk1, "_chunk1")
expData(chunk2, "_chunk2")
expData(chunk3, "_chunk3")
