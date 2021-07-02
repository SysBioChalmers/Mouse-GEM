
## DE analysis of datasets from transgenic AD mice


# Load packages
library(DESeq2)


## AD mouse models:
# 1. APPPS1
# 2. APPswe/PSEN1dE9 (line 85)
# 3. APPswe/PSEN1dE9 (C57BL6)
# 4. Trem2 KO (JAX)
# 5. Trem2 KO (KOMP) x APPPS1
# 6. APP23


## APPPS1
########################################################################
#GSE100070	48	APPPS1	
#GSE110741	96	APPPS1	consider about THY-Tau22
#GSE63943	  69	APPPS1	treatment study
sampleDir = "../data/APPPS1"
DEdir     = "../data/APPPS1"


## GSE100070 - APPPS1
setwd(sampleDir)
countData <- read.csv('GSE100070_counts.csv')
metaData  <- read.csv('GSE100070_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Age))


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","APPPS14","WT4"))
GSE100070_APPPS1_DentateGyrusMixedMale4 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE100070_APPPS1_DentateGyrusMixedMale4, "GSE100070_APPPS1_DentateGyrusMixedMale4.txt",
            append = FALSE, sep = "\t", dec = ".")


## GSE110741	96	APPPS1	consider about THY-Tau22
setwd(sampleDir)
countData <- read.csv('GSE110741_counts.csv')
metaData  <- read.csv('GSE110741_metaData.csv')
metaData$group <- factor(paste0(metaData$Model, metaData$Genotype, metaData$Age))

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","APPPS1Transgenic10","APPPS1Wildtype10"))
GSE110741_APPPS1_HippocampusMixedMale10 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1Transgenic4","APPPS1Wildtype4"))
GSE110741_APPPS1_HippocampusMixedMale4 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE110741_APPPS1_HippocampusMixedMale10, "GSE110741_APPPS1_HippocampusMixedMale10.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE110741_APPPS1_HippocampusMixedMale4, "GSE110741_APPPS1_HippocampusMixedMale4.txt",
            append = FALSE, sep = "\t", dec = ".")


## GSE63943	  69	APPPS1	treatment study
setwd(sampleDir)
countData <- read.csv('GSE63943_counts.csv')
metaData  <- read.csv('GSE63943_metaData.csv')
#metaData$group <- factor(paste0(metaData$Genotype, metaData$Tissue, metaData$Celltype, metaData$Gender, metaData$Age))
metaData$group <- factor(paste0(metaData$Genotype, metaData$Tissue, metaData$Age))
metaData

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","TransgenicvehicleHippocampalCA110","WTvehicleHippocampalCA110"))
GSE63943_APPPS1_HippocampalCA1MixedNA10 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE63943_APPPS1_HippocampalCA1MixedNA10, "GSE63943_APPPS1_HippocampalCA1MixedNA10.txt",
            append = FALSE, sep = "\t", dec = ".")



## APPswe/PSEN1dE9 (line 85)
########################################################################
#GSE104424	38	APPswe/PSEN1dE9 (line 85)	
#GSE136861	60	APPswe/PSEN1dE9 (line 85)	
#GSE93678 	13	APPswe/PSEN1dE9 (line 85)	
sampleDir = "../data/APPswe:PSEN1dE9-line85"
DEdir     = "../data/APPswe:PSEN1dE9-line85"


# GSE104424	38	APPswe/PSEN1dE9 (line 85)	
setwd(sampleDir)
countData <- read.csv('GSE104424_counts.csv')
metaData  <- read.csv('GSE104424_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Age))


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","APPdelta9Placebo10","WildtypePlacebo10"))
GSE104424_APPPS1line85_PlaceboHippocampusMixedMale10 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPdelta9Placebo6","WildtypePlacebo6"))
GSE104424_APPPS1line85_PlaceboHippocampusMixedMale6 <- res[ ! (is.na(res$pvalue)), ]

# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE104424_APPPS1line85_PlaceboHippocampusMixedMale10, "GSE104424_APPPS1line85_PlaceboHippocampusMixedMale10.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE104424_APPPS1line85_PlaceboHippocampusMixedMale6, "GSE104424_APPPS1line85_PlaceboHippocampusMixedMale6.txt",
            append = FALSE, sep = "\t", dec = ".")


# GSE136861	60	APPswe/PSEN1dE9 (line 85)
setwd(sampleDir)
countData <- read.csv('GSE136861_counts.csv', header = TRUE, sep = ",")
metaData  <- read.csv('GSE136861_metaData.csv', header = TRUE, sep = ",")
metaData$group <- factor(paste0(metaData$Genotype, metaData$Age))


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","Transgenic2","WT2"))
GSE136861_APPPS1line85_BrainMixedFemale2 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","Transgenic4","WT4"))
GSE136861_APPPS1line85_BrainMixedFemale4 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","Transgenic5","WT5"))
GSE136861_APPPS1line85_BrainMixedFemale5 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","Transgenic6","WT6"))
GSE136861_APPPS1line85_BrainMixedFemale6 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE136861_APPPS1line85_BrainMixedFemale2, "GSE136861_APPPS1line85_BrainMixedFemale2.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE136861_APPPS1line85_BrainMixedFemale4, "GSE136861_APPPS1line85_BrainMixedFemale4.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE136861_APPPS1line85_BrainMixedFemale5, "GSE136861_APPPS1line85_BrainMixedFemale5.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE136861_APPPS1line85_BrainMixedFemale6, "GSE136861_APPPS1line85_BrainMixedFemale6.txt", 
            append = FALSE, sep = "\t", dec = ".")


# GSE93678 	13	APPswe/PSEN1dE9 (line 85)	
setwd(sampleDir)
countData <- read.csv('GSE93678_counts.csv')
metaData  <- read.csv('GSE93678_metaData.csv')


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~Genotype, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("Genotype","Transgenicvehicle","WTvehicle"))
GSE93678_APPPS1line85_HippocampusMixedFemale13 <- res[ ! (is.na(res$pvalue)), ]

# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE93678_APPPS1line85_HippocampusMixedFemale13, "GSE93678_APPPS1line85_HippocampusMixedFemale13.txt",
            append = FALSE, sep = "\t", dec = ".")



## APPswe/PSEN1dE9 (C57BL6)
########################################################################
# GSE131659	93	APPswe/PSEN1dE9 (C57BL6)	
sampleDir = "../data/APPswe:PSEN1dE9-C57BL6"
DEdir     = "../data/APPswe:PSEN1dE9-C57BL6"

setwd(sampleDir)
countData <- read.csv('GSE131659_counts.csv')
#removeDuplicate = countData %>% group_by(genes) %>% summarise_all(funs(sum))
#countData <- as.data.frame(removeDuplicate)
metaData  <- read.csv('GSE131659_metaData.csv')
#metaData$group <- factor(paste0(metaData$Genotype, metaData$Tissue, metaData$Celltype, metaData$Gender, metaData$Age))
metaData$group <- factor(paste0(metaData$Genotype, metaData$Gender))
metaData

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)


res <- results(DEdds, contrast=c("group","TransgenicFemale","WTFemale"))
GSE131659_APPPS1C57BL6_BrainMixedFemale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TransgenicMale","WTMale"))
GSE131659_APPPS1C57BL6_BrainMixedMale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1WSBFemale","WTWSBFemale"))
GSE131659_APPPS1WSB_BrainMixedFemale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1WSBMale","WTWSBMale"))
GSE131659_APPPS1WSB_BrainMixedMale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1PWKFemale","WTPWKFemale"))
GSE131659_APPPS1PWK_BrainMixedFemale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1PWKMale","WTPWKMale"))
GSE131659_APPPS1PWK_BrainMixedMale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1CASTFemale","WTCASTFemale"))
GSE131659_APPPS1CAST_BrainMixedFemale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","APPPS1CASTMale","WTCASTMale"))
GSE131659_APPPS1CAST_BrainMixedMale8 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE131659_APPPS1C57BL6_BrainMixedFemale8, "GSE131659_APPPS1C57BL6_BrainMixedFemale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1C57BL6_BrainMixedMale8, "GSE131659_APPPS1C57BL6_BrainMixedMale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1WSB_BrainMixedFemale8, "GSE131659_APPPS1WSB_BrainMixedFemale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1WSB_BrainMixedMale8, "GSE131659_APPPS1WSB_BrainMixedMale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1PWK_BrainMixedFemale8, "GSE131659_APPPS1PWK_BrainMixedFemale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1PWK_BrainMixedMale8, "GSE131659_APPPS1PWK_BrainMixedMale8.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1CAST_BrainMixedFemale8, "GSE131659_APPPS1CAST_BrainMixedFemale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE131659_APPPS1CAST_BrainMixedMale8, "GSE131659_APPPS1CAST_BrainMixedMale8.txt", 
            append = FALSE, sep = "\t", dec = ".")



## Trem2 KO (JAX)
########################################################################
#GSE134031	56	Trem2 KO (JAX)	
#GSE124266	74	Trem2 KO (JAX)	
sampleDir = "../data/Trem2KO-JAX"
DEdir     = "../data/Trem2KO-JAX"


# GSE134031	56	Trem2 KO (JAX)
setwd(sampleDir)
countData <- read.csv('GSE134031_counts.csv')
removeDuplicate = countData %>% group_by(genes) %>% summarise_all(funs(sum))
countData <- as.data.frame(removeDuplicate)

metaData  <- read.csv('GSE134031_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Celltype, metaData$Age))
metaData

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)


res <- results(DEdds, contrast=c("group","TREM2KOAstrocyte16","WTAstrocyte16"))
GSE134031_TREM2KOJAX_BrainAstrocyteFemale16 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOAstrocyte2","WTAstrocyte2"))
GSE134031_TREM2KOJAX_BrainAstrocyteFemale2 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOMicroglia16","WTMicroglia16"))
GSE134031_TREM2KOJAX_BrainMicrogliaFemale16 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOMicroglia2","WTMicroglia2"))
GSE134031_TREM2KOJAX_BrainMicrogliaFemale2 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE134031_TREM2KOJAX_BrainAstrocyteFemale16, "GSE134031_TREM2KOJAX_BrainAstrocyteFemale16.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE134031_TREM2KOJAX_BrainAstrocyteFemale2, "GSE134031_TREM2KOJAX_BrainAstrocyteFemale2.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE134031_TREM2KOJAX_BrainMicrogliaFemale16, "GSE134031_TREM2KOJAX_BrainMicrogliaFemale16.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE134031_TREM2KOJAX_BrainMicrogliaFemale2, "GSE134031_TREM2KOJAX_BrainMicrogliaFemale2.txt", 
            append = FALSE, sep = "\t", dec = ".")


# GSE124266	74	Trem2 KO (JAX) "ignore Trem2KOHet for now"
setwd(sampleDir)
countData <- read.csv('GSE124266_counts.csv')
metaData  <- read.csv('GSE124266_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Celltype, metaData$Age))


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","controlDietTREM2KOMicroglia13","controlDietWTMicroglia13"))
GSE124266_TREM2KOJAX_BrainMicrogliaMixed13 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","controlDietTREM2KOMicroglia11","controlDietWTMicroglia11"))
GSE124266_TREM2KOJAX_BrainMicrogliaMixed11 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","controlDietTREM2KONonMicroglia13","controlDietWTNonMicroglia13"))
GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed13 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","controlDietTREM2KONonMicroglia11","controlDietWTNonMicroglia11"))
GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed11 <- res[ ! (is.na(res$pvalue)), ]


# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE124266_TREM2KOJAX_BrainMicrogliaMixed13, "GSE124266_TREM2KOJAX_BrainMicrogliaMixed13.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE124266_TREM2KOJAX_BrainMicrogliaMixed11, "GSE124266_TREM2KOJAX_BrainMicrogliaMixed11.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed13, "GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed13.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed11, "GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed11.txt", 
            append = FALSE, sep = "\t", dec = ".")



## Trem2 KO (KOMP) x APPPS1
########################################################################
# GSE104381	128	Trem2 KO (KOMP) x APPPS1	
sampleDir = "../data/Trem2KO-KOMP:APPPS1"
DEdir     = "../data/Trem2KO-KOMP:APPPS1"


setwd(sampleDir)
countData <- read.csv('GSE104381_counts.csv')
metaData  <- read.csv('GSE104381_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Tissue, metaData$Age))
metaData

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~ group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)


res <- results(DEdds, contrast=c("group","TREM2KOHippocampus8","WTHippocampus8"))
GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOHippocampus4","WTHippocampus4"))
GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale4 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOCortex8","WTCortex8"))
GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale8 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","TREM2KOCortex4","WTCortex4"))
GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale4 <- res[ ! (is.na(res$pvalue)), ]

# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale8, "GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale4, "GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale4.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale8, "GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale8.txt", 
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale4, "GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale4.txt", 
            append = FALSE, sep = "\t", dec = ".")


## APP23
########################################################################
#GSE80465	  24	APP23	
#GSE104630	39	APP23	lipopolysaccharides treatment
sampleDir = "../data/APP23"
DEdir     = "../data/APP23"


# GSE80465	  24	APP23	
setwd(sampleDir)
countData <- read.csv('GSE80465_counts.csv')
metaData  <- read.csv('GSE80465_metaData.csv')
metaData$group <- factor(paste0(metaData$Genotype, metaData$Age))


# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~group, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("group","HETAPP232","WT2"))
GSE80465_APP23_ForebrainMixedMale2 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","HETAPP236","WT6"))
GSE80465_APP23_ForebrainMixedMale6 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","HETAPP2318","WT18"))
GSE80465_APP23_ForebrainMixedMale18 <- res[ ! (is.na(res$pvalue)), ]

res <- results(DEdds, contrast=c("group","HETAPP2324","WT24"))
GSE80465_APP23_ForebrainMixedMale24 <- res[ ! (is.na(res$pvalue)), ]

# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE80465_APP23_ForebrainMixedMale2, "GSE80465_APP23_ForebrainMixedMale2.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE80465_APP23_ForebrainMixedMale6, "GSE80465_APP23_ForebrainMixedMale6.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE80465_APP23_ForebrainMixedMale18, "GSE80465_APP23_ForebrainMixedMale18.txt",
            append = FALSE, sep = "\t", dec = ".")

write.table(GSE80465_APP23_ForebrainMixedMale24, "GSE80465_APP23_ForebrainMixedMale24.txt",
            append = FALSE, sep = "\t", dec = ".")


# GSE104630	39	APP23	lipopolysaccharides treatment
# gene count data is unavailable
setwd(sampleDir)
countData <- read.csv('GSE104630_counts.csv')
metaData  <- read.csv('GSE104630_metaData.csv')
#metaData$group <- factor(paste0(metaData$Genotype, metaData$Tissue, metaData$Celltype, metaData$Gender, metaData$Age))

# Construct DESEQDataSet Object
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design= ~Genotype, tidy = TRUE)

# Differential expression analysis
DEdds <- DESeq(dds)

res <- results(DEdds, contrast=c("Genotype","APP23PBS","WTPBS"))
GSE104630_APP23_BrainMicrogliaFemale9 <- res[ ! (is.na(res$pvalue)), ]

# output DE results for gene-set analysis
setwd(DEdir)
write.table(GSE104630_APP23_BrainMicrogliaFemale9, "GSE104630_APP23_BrainMicrogliaFemale9.txt",
            append = FALSE, sep = "\t", dec = ".")


