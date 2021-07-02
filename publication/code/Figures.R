

library(pheatmap)
library(viridisLite)
library(RColorBrewer)
library(grid)
library(dplyr)
library(tidyr)
library(fmsb)
library(readxl)
library(ggpubr)



# Specify Box folder code directory
setwd('./')

##############################
### Fig 2B-2F. Radar plots ###
##############################


####
mouseGEMs <- read_excel("Dataset_S01.xlsx", 
                     sheet = "Table1",range = "C3:I5",
                     col_names = FALSE)
legend = c("Mouse1", "iMM1865", "MMR")

####
ratGEMs <- read_excel("Dataset_S01.xlsx", 
                      sheet = "Table1",range = "C7:I8",
                      col_names = FALSE)
legend = c("Rat1", "iRno")

####
zebrafishGEMs <- read_excel("Dataset_S01.xlsx", 
                      sheet = "Table1",range = "C10:I12",
                      col_names = FALSE)
legend = c("Zebrafish1", "ZebraGEM2", "ZebraGEM")

####
flyGEMs <- read_excel("Dataset_S01.xlsx",
                      sheet = "Table1",range = "C15:I17",
                      col_names = FALSE)
legend = c("Fruitfly1", "BMID000000141998", "FlySilico")

####
wormGEMs <- read_excel("Dataset_S01.xlsx", 
                      sheet = "Table1",range = "C19:I21",
                      col_names = FALSE)
legend = c("Worm1", "iCEL1314", "iCEL1273")


# Create data from GEMs
data <- mouseGEMs
data <- rbind(c(4000, 15000, 10000, 10000, 5000, 800, 80), c(1000, 5000, 5000, 4000, 2000, 0, 0), data)

# Create data from GEMs
data <- ratGEMs
data <- rbind(c(4000, 15000, 10000, 10000, 5000, 800, 80), c(1000, 5000, 4000, 4000, 2000, 0, 0), data)

# Create data from GEMs
data <- zebrafishGEMs
data <- rbind(c(4000, 15000, 10000, 10000, 5000, 800, 80), c(1000, 2000, 2000, 2000, 1000, 0, 0), data)

# Create data from GEMs
data <- flyGEMs
data <- rbind(c(5000, 15000, 8000, 10000, 8000, 600, 80), c(0, 200, 0, 200, 0, 0, 0), data)

# Create data from GEMs
data <- wormGEMs
data <- rbind(c(2000, 15000, 10000, 10000, 5000, 800, 80), c(1000, 1000, 1000, 1000, 0, 0, 0), data)    


## Common part ------------------
# Add colnames
colnames(data) <- c("Genes", "Reactions", "Enzymatic reactions",  "Metabolites", "Unique metabolites", "Enzyme complex", "Memote score")

# Prepare color
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9), rgb(0.6,0.2,0.5,0.9), rgb(0.4,0.2,0.5,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4), rgb(0.4,0.5,0.5,0.4), rgb(0.4,0.5,0.5,0.4)  )

# Custom the radarChart
radarchart( data, axistype=2,
            #custom polygon
            #pcol=colors_border , pfcol=colors_in , plwd=4, plty=1 , 
            pcol=colors_border , plwd=4, plty=1 , 
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", cglwd=1.1,
            caxislabels = seq(1000,4000,750),
            #custom labels
            vlcex=1.1
)

# Legend
legend(x=0.6, y=1.2, legend, bty = "n", pch=20 , col=colors_border , text.col = "black", cex=1.4, pt.cex=1.6)



#########################
### GENE SET HEATMAPs ###
#########################


if (!dir.exists('FigParts')) {dir.create('FigParts')}

GStype = 'ReMet'  # 'Combine', 'ReMet', or 'Subsys'

# Fig S6
if (GStype == 'Subsys') {
  filter_type = 'top each'  # 'top each' or 'pval'
  filter_threshold = 4
  heatmap_height = 6
  heatmap_width = 8.6

# Fig 3D
} else if (GStype == 'ReMet') {
  filter_type = 'top each'
  filter_threshold = 6
  heatmap_height = 6
  heatmap_width = 7.9

# Fig S5
} else if (GStype == 'Combine') {
  filter_type = 'top each'
  filter_threshold = 3
  heatmap_height = 8
  heatmap_width = 10
}


# load data
p_up <- read.delim(paste0('FigData/GSAres', GStype, '_pUp.txt'), row.names=1)
p_down <- read.delim(paste0('FigData/GSAres', GStype, '_pDown.txt'), row.names=1)
if (!all(colnames(p_up) == colnames(p_down)) || !all(rownames(p_up) == rownames(p_down))) {
  stop('*** DATA MATRICES ARE NOT ALIGNED! ***')
}

# extract information from column names
gse <- unlist(lapply(colnames(p_up), function(x) unlist(strsplit(x, '_'))[1]))
coldata <- data.frame(Model = unlist(lapply(colnames(p_up), function(x) unlist(strsplit(x, '_'))[2])), row.names=colnames(p_up))
other <- unlist(lapply(colnames(p_up), function(x) unlist(strsplit(x, '_'))[3]))
coldata$Month <- as.factor(as.numeric(regmatches(other, regexpr('\\d+$', other))))

coldata$Tissue <- NA
coldata$Tissue[grepl('Hippocamp', other)] <- 'Hippocampus'
coldata$Tissue[grepl('BrainMixed', other) | grepl('Forebrain', other)] <- 'Brain'
coldata$Tissue[grepl('Cortex', other)] <- 'Cortex'
coldata$Tissue[grepl('DentateGyrus', other)] <- 'Dentate Gyrus'
coldata$Tissue[grepl('Microglia', other)] <- 'Microglia'
coldata$Tissue[grepl('Astrocyte', other)] <- 'Astrocyte'
coldata$Tissue <- as.factor(coldata$Tissue)

modelnames <- data.frame(short=c('Trem2KOKOMPAPPPS1','APPPS1','APP23','APPPS1line85','TREM2KOJAX','APPPS1C57BL6'),
                         nice=c('Trem2KO (KOMP) x APPPS1','APPPS1','APP23','APPswe/PSEN1dE9','Trem2KO (JAX)','APPswe/PSEN1dE9'))
coldata$Model <- as.factor(modelnames$nice[match(coldata$Model, modelnames$short)])

# coldata$Stage <- rep('Unspecified', length(coldata$Month))
# coldata$Stage[coldata$Model == 'APP23' & coldata$Month %in% c('2','6')] <- 'Early'
# coldata$Stage[coldata$Model == 'APP23' & coldata$Month %in% c('18','24')] <- 'Late'
# coldata$Stage[coldata$Model == 'APPPS1' & gse %in% c('GSE63943', 'GSE110741') & coldata$Month %in% c('3','4')] <- 'Early'
# coldata$Stage[coldata$Model == 'APPPS1' & gse %in% c('GSE63943', 'GSE110741') & coldata$Month %in% c('10','20')] <- 'Late'
# coldata$Stage[coldata$Model == 'APPswe/PSEN1dE9' & gse %in% c('GSE104424', 'GSE136861') & coldata$Month %in% c('2','4','5','6')] <- 'Early'
# coldata$Stage[coldata$Model == 'APPswe/PSEN1dE9' & gse %in% c('GSE104424', 'GSE136861') & coldata$Month %in% c('10')] <- 'Late'
# coldata$Stage[coldata$Model == 'Trem2KO (JAX)' & coldata$Month %in% c('2','11')] <- 'Early'
# coldata$Stage[coldata$Model == 'Trem2KO (JAX)' & coldata$Month %in% c('13','16')] <- 'Late'
# coldata$Stage[coldata$Model == 'Trem2KO (KOMP) x APPPS1' & coldata$Month %in% c('4')] <- 'Early'
# coldata$Stage[coldata$Model == 'Trem2KO (KOMP) x APPPS1' & coldata$Month %in% c('8')] <- 'Late'
# coldata$Stage <- as.factor(coldata$Stage)

coldata$deposit <- rep('No', length(coldata$Month))
coldata$deposit[coldata$Model == 'APP23' & coldata$Month %in% c('2','6')] <- 'Before'
coldata$deposit[coldata$Model == 'APP23' & coldata$Month %in% c('18','24')] <- 'After'
coldata$deposit[coldata$Model == 'APPPS1'] <- 'After'
coldata$deposit[coldata$Model == 'APPswe/PSEN1dE9' & coldata$Month %in% c('2','4','5','6')] <- 'Before'
coldata$deposit[coldata$Model == 'APPswe/PSEN1dE9' & coldata$Month %in% c('8','10','13')] <- 'After'
coldata$deposit[coldata$Model == 'Trem2KO (KOMP) x APPPS1'] <- 'After'
coldata$deposit <- factor(coldata$deposit, levels=c('No', 'Before', 'After'))

coldata <- select(coldata, !Month)  # remove Month from coldata


# extract information from rows
if (GStype == 'Combine') {
  gs_lib <- unlist(lapply(rownames(p_up), function(x) unlist(strsplit(x, '-'))[1]))
  gs <- tolower(unlist(lapply(rownames(p_up), function(x) paste(unlist(strsplit(x, '-'))[-1], collapse=' '))))
  gs <- gsub('(^.)','\\U\\1', gs, perl=T)  # capitalize first letter of every name
  long_indx <- unlist(lapply(gs, nchar)) > 40
  gs <- substr(gs, 1, 40)
  gs[long_indx] <- paste0(gs[long_indx], '...')
  
  rowdata <- data.frame(Database=as.factor(gs_lib), row.names=rownames(p_up))
  
  # manually reformat some gene set names
  gs <- gsub('Dna ', 'DNA ', gs)
  gs[gs == 'Glycosaminoglycan biosynthesis chondroit...'] <- 'GAG biosynthesis: chondroitin sulfate'
  gs[gs == 'Rho gtpases activate wasps and waves'] <- 'RHO GTPases Activate WASPs and WAVEs'
  gs[gs == 'Slc mediated transmembrane transport'] <- 'SLC-mediated transmembrane transport'
  gs[gs == 'Interaction between l1 and ankyrins'] <- 'Interaction between L1 and Ankyrins'
  gs[gs == 'G2m checkpoint'] <- 'G2M checkpoint'
  gs[gs == 'Activation of the pre replicative comple...'] <- 'Activation of pre-replicative complex'
  gs[gs == 'Citric acid cycle tca cycle'] <- 'TCA cycle'
  gs[gs == 'Mtorc1 signaling'] <- 'mTORC1 signaling'
  gs[gs == 'Il2 stat5 signaling'] <- 'IL2 STAT5 signaling'
  gs[gs == 'Signaling by notch'] <- 'Signaling by Notch'
  gs[gs == 'Hiv infection'] <- 'HIV infection'
  gs[gs == 'Cytokine cytokine receptor interaction'] <- 'Cytokine-cytokine receptor interaction'
  gs[gs == 'Metabolism of rna'] <- 'Metabolism of RNA'
  gs[gs == 'Uv response dn'] <- 'UV response: Down'
  gs[gs == 'Asparagine n linked glycosylation'] <- 'Asparagine N-linked glycosylation'
  gs[gs == 'Notch4 intracellular domain regulates tr...'] <- 'NOTCH4 Intracell. Domain Reg. Transcr.'
  gs[gs == 'Tnfa signaling via nfkb'] <- 'TNFa signaling via NF-kB'
  gs[gs == 'E2f targets'] <- 'E2F targets'
  gs[gs == 'Mapk signaling pathway'] <- 'MAPK signaling pathway'

} else {
  gs <- rownames(p_up)
  
  # manually reformat some gene set names
  gs[gs == 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'] <- 'TCA cycle and dicarboxylate metab.'
  gs[gs == 'Glycosphingolipid biosynthesis-lacto and neolacto series'] <- 'Glycosphingolipid synthesis-lacto and neolacto'
  gs[gs == 'Phenylalanine, tyrosine and tryptophan biosynthesis'] <- 'Phe, Tyr and Trp biosynthesis'
  gs[gs == 'Valine, leucine, and isoleucine metabolism'] <- 'Val, Leu, and Ile metabolism'
  gs[gs == 'Amino sugar and nucleotide sugar metabolism'] <- 'Amino & nucleotide sugar metabolism'
  gs[gs == 'Chondroitin / heparan sulfate biosynthesis'] <- 'Chondroitin / heparan sulfate biosynth.'
  
  gs[gs == '(6Z,9Z,12Z,15Z)-octadecatetraenoyl-CoA[c]'] <- 'octadecatetraenoyl-CoA[c]'
  gs[gs == '(5Z,8Z,11Z,14Z,17Z)-eicosapentaenoyl-CoA[c]'] <- 'eicosapentaenoyl-CoA[c]'
  gs[gs == '(7Z)-octadecenoyl-CoA[c]'] <- 'octadecenoyl-CoA[c]'
  gs[gs == '(6Z,9Z)-octadecadienoyl-CoA[c]'] <- 'octadecadienoyl-CoA[c]'
  
  rowdata <- NA
}

# merge directional pvalues and convert to signed log10 pvalues
p_dir <- p_up
p_dir[p_down < p_up] <- p_down[p_down < p_up]
log_p_dir <- -log10(p_up)
log_p_dir[p_down < p_up] <- log10(p_down[p_down < p_up])

# filter and prepare data for plotting
if (casefold(filter_type) == 'top each') {
  keep <- unique(as.vector(apply(p_dir, 2, function(x) order(x))[1:filter_threshold, ]))
} else if (casefold(filter_type) == 'pval') {
  keep <- rowSums(p_dir < filter_threshold) > 0
}
dat <- log_p_dir[keep, ]
gs <- gs[keep]


# set up color palette functions
pal_RdYlBu <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))
pal_YlOrBr <- colorRampPalette(rev(brewer.pal(9, 'YlOrBr')))
pal_YlGnBu <- colorRampPalette(rev(brewer.pal(9, 'YlGnBu')))
pal_purp <- colorRampPalette(rev(brewer.pal(9, 'Purples')))

# specify annotation colors
annColors <- list(Model = rev(magma(nlevels(coldata$Model))) %>% setNames(levels(coldata$Model)),
                  #Model = rev(viridis(nlevels(coldata$Model))) %>% setNames(levels(coldata$Model)),
                  #Model = pal_YlOrBr(nlevels(coldata$Model)) %>% setNames(levels(coldata$Model)),
                  Tissue = brewer.pal(nlevels(coldata$Tissue), 'Set1') %>% setNames(levels(coldata$Tissue)),
                  #Stage = c('white','gray70','black') %>% setNames(c('Unspecified','Early','Late')),
                  `AB Deposition` = c('white','gray70','black') %>% setNames(c('No','Before','After')))
if (GStype == 'Combine') {
  annColors$Database <- brewer.pal(nlevels(rowdata$Database), 'Dark2') %>% setNames(levels(rowdata$Database))
}

colnames(coldata)[colnames(coldata) == 'deposit'] <- 'AB Deposition'

# generate heatmap
pdf(file=paste0('FigParts/GSAheatmap_', GStype, '.pdf'), width=heatmap_width, height=heatmap_height, onefile=F)
pheatmap(dat,
         scale='none',
         color=pal_RdYlBu(100), #brewer.pal(11,'RdYlBu'),
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         clustering_method='ward.D2',  # ward.D, ward.D2, single, complete, average, mcquitty, median, centroid
         breaks=seq(-4,4,len=100),
         # cutree_rows=4,
         # cutree_cols=4,
         # labels_col=rep('',ncol(dat)),  # remove column labels
         labels_col=gse,
         labels_row=gs,
         # angle_col=0,
         annotation_col=coldata,
         annotation_row=rowdata,
         annotation_colors=annColors,
         annotation_names_col=T)

grid_linewidth <- 0.5
grid_color <- 'white'
grid.ls(grid.force())
if (GStype == 'Combine') {
  grid.gedit('row_annotation[.]', gp=gpar(col=grid_color, lwd=grid_linewidth))
}
grid.gedit('matrix::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('annotation_legend::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('col_annotation[.]', gp=gpar(col=grid_color, lwd=grid_linewidth))
dev.off()



##############
### Fig 5A ###
##############

# define annotations
annotation <- c(
  'GSE110741/Hippocampus/Male/4',      # APPPS1
  'GSE110741/Hippocampus/Male/10',     # APPPS1
  'GSE100070/DentateGyrus/Male/4',     # APPPS1
  'GSE63943/HippocampalCA1/NA/10',     # APPPS1
  'GSE80465/Forebrain/Male/2',         # APP23
  'GSE80465/Forebrain/Male/6',         # APP23
  'GSE80465/Forebrain/Male/18',        # APP23
  'GSE80465/Forebrain/Male/24',        # APP23
  'GSE131659/Brain/Female/8',          # APPPS1C57BL6
  'GSE131659/Brain/Male/8',            # APPPS1C57BL6
  'GSE93678/Hippocampus/Female/13',    # APPPS1line85
  'GSE104424/Hippocampus/Male/6',      # APPPS1line85
  'GSE104424/Hippocampus/Male/10',     # APPPS1line85
  'GSE136861/Brain/Female/2',          # APPPS1line85
  'GSE136861/Brain/Female/4',          # APPPS1line85
  'GSE136861/Brain/Female/5',          # APPPS1line85
  'GSE136861/Brain/Female/6',          # APPPS1line85
  'GSE104381/Cortex/Male/4',           # Trem2KOKOMPAPPPS1
  'GSE104381/Cortex/Male/8',           # Trem2KOKOMPAPPPS1
  'GSE104381/Hippocampus/Male/4',      # Trem2KOKOMPAPPPS1
  'GSE104381/Hippocampus/Male/8',      # Trem2KOKOMPAPPPS1
  'GSE124266/Microglia/Mixed/11',      # TREM2KOJAX
  'GSE124266/Microglia/Mixed/13',      # TREM2KOJAX
  'GSE134031/Astrocyte/Female/2',      # TREM2KOJAX
  'GSE134031/Astrocyte/Female/16',     # TREM2KOJAX
  'GSE134031/Microglia/Female/2',      # TREM2KOJAX
  'GSE134031/Microglia/Female/16'      # TREM2KOJAX
  #'GSE124266/NonMicroglia/Mixed/11',  # TREM2KOJAX
  #'GSE124266/NonMicroglia/Mixed/13'   # TREM2KOJAX
)

# define datatypes
data_types <- c(
  'GSE110741_APPPS1_HippocampusMixedMale4',
  'GSE110741_APPPS1_HippocampusMixedMale10',
  'GSE100070_APPPS1_DentateGyrusMixedMale4',
  'GSE63943_APPPS1_HippocampalCA1MixedNA10',
  'GSE80465_APP23_ForebrainMixedMale2',
  'GSE80465_APP23_ForebrainMixedMale6',
  'GSE80465_APP23_ForebrainMixedMale18',
  'GSE80465_APP23_ForebrainMixedMale24',
  'GSE131659_APPPS1C57BL6_BrainMixedFemale8',
  'GSE131659_APPPS1C57BL6_BrainMixedMale8',
  'GSE93678_APPPS1line85_HippocampusMixedFemale13',
  'GSE104424_APPPS1line85_HippocampusMixedMale6',
  'GSE104424_APPPS1line85_HippocampusMixedMale10',
  'GSE136861_APPPS1line85_BrainMixedFemale2',
  'GSE136861_APPPS1line85_BrainMixedFemale4',
  'GSE136861_APPPS1line85_BrainMixedFemale5',
  'GSE136861_APPPS1line85_BrainMixedFemale6',
  'GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale4',
  'GSE104381_Trem2KOKOMPAPPPS1_CortexMixedMale8',
  'GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale4',
  'GSE104381_Trem2KOKOMPAPPPS1_HippocampusMixedMale8',
  'GSE124266_TREM2KOJAX_BrainMicrogliaMixed11',
  'GSE124266_TREM2KOJAX_BrainMicrogliaMixed13',
  'GSE134031_TREM2KOJAX_BrainAstrocyteFemale2',
  'GSE134031_TREM2KOJAX_BrainAstrocyteFemale16',
  'GSE134031_TREM2KOJAX_BrainMicrogliaFemale2',
  'GSE134031_TREM2KOJAX_BrainMicrogliaFemale16'
  #'GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed11',
  #'GSE124266_TREM2KOJAX_BrainNonMicrogliaMixed13'
)

targetGenes <- c('Gm2a','Hexa','Hexb','Ctsb','Ctsc','Ctsd','Ctsf','Ctss','Ctsz','Nme4')

fnames <- paste0(data_types, '.txt')
foldChange <- as.data.frame(matrix(NA, nrow=length(targetGenes), ncol=length(fnames)), row.names=targetGenes) %>% setNames(annotation)
pValues <- as.data.frame(matrix(NA, nrow=length(targetGenes), ncol=length(fnames)), row.names=targetGenes) %>% setNames(annotation)

# load each dataset individually
for (i in seq(length(fnames))) {
  DEdata <- read.delim(file.path('ADmodelDE', fnames[i]), sep='\t')
  geneIDs <- rownames(DEdata)     # gene abbreviations
  log2FC <- DEdata$log2FoldChange     # log2(fold-changes) (disease vs. control)
  pvals <- DEdata$padj;     # fold-change significance (p-values)
  foldChange[, i] <- log2FC[match(targetGenes, geneIDs)]
}

# convert annotations into dataframe
gse <- unlist(lapply(annotation, function(x) unlist(strsplit(x, '/'))[1]))
tissue <- unlist(lapply(annotation, function(x) unlist(strsplit(x, '/'))[2]))
tissue_short <- tissue
# tissue_short[tissue_short %in% 'Forebrain'] <- 'FBrain'
tissue_short[tissue_short %in% 'DentateGyrus'] <- 'Dentate'
# tissue_short[tissue_short %in% 'Astrocyte'] <- 'Astro.'
# tissue_short[tissue_short %in% 'Microglia'] <- 'Microg.'
tissue_short[tissue_short %in% c('Hippocampus','HippocampalCA1')] <- 'Hippo.'
sex <- unlist(lapply(annotation, function(x) unlist(strsplit(x, '/'))[3]))
sex_short <- sex
sex_short[sex_short %in% 'Male'] <- 'M'
sex_short[sex_short %in% 'Female'] <- 'F'
sex_short[sex_short %in% 'Mixed'] <- 'Mix'
sex_short[sex_short %in% 'NA'] <- ''
month <- unlist(lapply(annotation, function(x) unlist(strsplit(x, '/'))[4]))
month <- as.numeric(month)
model <- unlist(lapply(data_types, function(x) unlist(strsplit(x, '_'))[2]))
model[model %in% c('APPPS1C57BL6', 'APPPS1line85')] <- 'APPswe/PSEN1dE9'
model[model %in% c('TREM2KOJAX')] <- 'Trem2KO(JAX)'
model[model %in% c('Trem2KOKOMPAPPPS1')] <- 'Trem2KO(KOMP) x APPPS1'
model <- factor(model, levels=c('APP23', 'APPPS1', 'APPswe/PSEN1dE9', 'Trem2KO(KOMP) x APPPS1', 'Trem2KO(JAX)'))
shortname <- paste(gse, tissue_short, sex_short, month)
annot_data <- data.frame(shortname, gse, tissue, sex, month, model, row.names=annotation)

annot_data$x_ordering <- (month+10)*10 + as.numeric(factor(tissue_short))

# prepare data for boxplot
dat <- pivot_longer(foldChange %>% rownames_to_column('Gene'), cols=c(colnames(foldChange)), names_to='annot', values_to='Log2FC')
dat <- cbind(dat, annot_data)

# generate boxplot
pdf(file='FigParts/lysosomal_gene_exp_boxplot.pdf', width=8.5, height=4)
ggplot(dat, aes(x=reorder(shortname, x_ordering), y=Log2FC, fill=model)) +
  # geom_boxplot() +   # keep outliers
  geom_boxplot(outlier.shape=NA) +  # remove outliers
  facet_grid(. ~ model, scales='free_x', space='free_x', labeller=label_wrap_gen(width=10)) +
  # theme_bw() +
  theme(axis.text=element_text(size=10, color='black'),
        axis.text.x=element_text(size=10, angle=45, hjust=1),
        plot.margin=margin(t=10, r=10, l=50, b=10),
        legend.position='none') +
  coord_cartesian(ylim=c(-1.5,2.5)) +   # rescale if outliers removed
  xlab('') +
  ylab('Log2 Fold Change')
invisible(dev.off())

# title('Expression pattern of genes involved in GM2 and peptide degradation pathways');


# Perform Kruskal-Wallis test (non-parametric ANOVA)
# ==================================================

# perform KW test on each mouse model type
models <- levels(dat$model)
kw_pvals <- NULL
for (m in models) {
  
  keep <- dat$model %in% m
  x <- data.frame(sample=dat$shortname[keep], log2FC=dat$Log2FC[keep])
  
  res <- kruskal.test(log2FC ~ sample, data=x)
  kw_pvals <- append(kw_pvals, res$p.value)
  
  # optional: test pair-wise differences
  # pair_res <- pairwise.wilcox.test(x$log2FC, x$sample, p.adjust.method = "BH")
}
names(kw_pvals) <- models



########################################
### TINIT MODEL STRUCTURE COMPARISON ###
########################################

# load data
ham <- read.delim('FigData/init_models_hamming_similarity.txt', row.names=1)
subsys <- read.delim('FigData/init_models_subsystem_coverage.txt', row.names=1)

# extract information from rownames
rnames <- rownames(ham)

gse <- unlist(lapply(rnames, function(x) unlist(strsplit(x, ' '))[1]))

model <- rep('APP23', length(gse))
model[grepl('APPPS1', rnames) & !grepl('Trem2', rnames)] <- 'APPPS1'
model[grepl('APPswe/PSEN1dE9', rnames)] <- 'APPswe/PSEN1dE9'
model[grepl('TREM2KO', rnames) & !grepl('KOMP', rnames)] <- 'Trem2KO(JAX)'
model[grepl('Trem2', rnames)] <- 'Trem2KO(KOMP) x APPPS1'
model <- factor(model, levels=c('APP23', 'APPPS1', 'APPswe/PSEN1dE9', 'Trem2KO(KOMP) x APPPS1', 'Trem2KO(JAX)'))

sex <- rep('', length(gse))
sex[grepl('Male', rnames)] <- 'M'
sex[grepl('Female', rnames)] <- 'F'

age <- unlist(lapply(rnames, function(x) tail(unlist(strsplit(x, ' ')),1)))

tissue <- rep('Brain', length(gse))
tissue[grepl('Hippocamp', rnames)] <- 'Hippo.'
tissue[grepl('Forebrain', rnames)] <- 'Forebrain'
tissue[grepl('Dentate gyrus', rnames)] <- 'Dentate'
tissue[grepl('Microglia', rnames)] <- 'Microglia'
tissue[grepl('Astrocyte', rnames)] <- 'Astrocyte'
tissue[grepl('Cortex', rnames)] <- 'Cortex'
tissue <- as.factor(tissue)

deposit <- rep('No', length(gse))
deposit[model == 'APP23' & age %in% c('2','6')] <- 'Before'
deposit[model == 'APP23' & age %in% c('18','24')] <- 'After'
deposit[model == 'APPPS1'] <- 'After'
deposit[model == 'APPswe/PSEN1dE9' & age %in% c('2','4','5','6')] <- 'Before'
deposit[model == 'APPswe/PSEN1dE9' & age %in% c('8','10','13')] <- 'After'
deposit[model == 'Trem2KO(KOMP) x APPPS1'] <- 'After'
deposit <- factor(deposit, levels=c('No', 'Before', 'After'))

# define labels and update row/column names
name_labels <- paste(gse, tissue, sex, age)
name_labels <- gsub('  ', ' ', name_labels)

rownames(ham) <- name_labels
colnames(ham) <- name_labels
colnames(subsys) <- name_labels


# calculate differences in subsystem coverage (% deviation from mean)
sub_diff <- (subsys - rowMeans(subsys)) / rowMeans(subsys) * 100


# set up color palette functions
pal_RdYlBu <- colorRampPalette(rev(brewer.pal(11, 'RdYlBu')))
pal_YlOrBr <- colorRampPalette(rev(brewer.pal(9, 'YlOrBr')))
pal_YlGnBu <- colorRampPalette(rev(brewer.pal(9, 'YlGnBu')))
pal_purp <- colorRampPalette(rev(brewer.pal(9, 'Purples')))

# specify annotation colors
annColors <- list(Model = rev(magma(nlevels(model))) %>% setNames(levels(model)),
                  Tissue = brewer.pal(nlevels(tissue), 'Set1') %>% setNames(levels(tissue)),
                  `AB Deposition` = c('white','gray70','black') %>% setNames(c('No','Before','After')))
coldata <- data.frame(Model=model, Tissue=tissue, Deposition=deposit)
colnames(coldata)[3] <- 'AB Deposition'
rownames(coldata) <- name_labels


# Fig S2. generate Hamming similarity heatmap
# ===================================
pdf(file='FigParts/tINIT_Ham_Heatmap.pdf', width=10, height=7.5, onefile=F)
pheatmap(ham,
         scale='none',
         color=viridis(100),
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         clustering_method='ward.D2',
         breaks=seq(0.8,1,len=100),
         labels_col=name_labels,
         labels_row=name_labels,
         annotation_col=coldata,
         #annotation_row=coldata,
         annotation_colors=annColors,
         angle_col=90)

grid_linewidth <- 0.5
grid_color <- 'white'
grid.ls(grid.force())
grid.gedit('row_annotation[.]', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('matrix::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('annotation_legend::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('col_annotation[.]', gp=gpar(col=grid_color, lwd=grid_linewidth))
dev.off()


# Fig S3. generate subsystem coverage heatmap
# ===================================

# filter data (remove uninteresting subsystems)
maxdiff <- apply(abs(sub_diff), 1, max)
# dat <- sub_diff[maxdiff > 10, ]  # difference threshold
dat <- sub_diff[order(maxdiff, decreasing=T)[1:30], ]  # top N most different subsystems

# abbreviate some rownames
rownames(dat) <- sub('Beta oxidation', 'B-oxid.', rownames(dat))
rownames(dat) <- sub('\\(peroxisomal\\)', '[p]', rownames(dat))
rownames(dat) <- sub('Glycosphingolipid', 'GSL', rownames(dat))
rownames(dat) <- sub('lacto and neolacto series', 'lacto and neolacto', rownames(dat))
rownames(dat) <- sub('tyrosine and tryptophan biosynthesis', 'tyrosine, tryptophan synthesis', rownames(dat))


# generate heatmap
pdf(file='FigParts/tINIT_subsys_heatmap.pdf', width=11, height=8.3, onefile=F)  # height=6.8 for 20 rows, 8.3 for 30, 9.8 for 40
pheatmap(dat,
         scale='none',
         color=pal_RdYlBu(100),
         clustering_distance_rows='euclidean',
         clustering_distance_cols='euclidean',
         clustering_method='ward.D2',
         breaks=seq(-25,25,len=100),
         labels_col=name_labels,
         annotation_col=coldata,
         annotation_colors=annColors,
         annotation_names_col=T,
         angle_col=90)

grid_linewidth <- 0.5
grid_color <- 'white'
grid.ls(grid.force())
grid.gedit('matrix::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('annotation_legend::GRID.rect', gp=gpar(col=grid_color, lwd=grid_linewidth))
grid.gedit('col_annotation[.]', gp=gpar(col=grid_color, lwd=grid_linewidth))
dev.off()


##############
### Fig S9 ###
##############

# Reanalyze CSF peptides detected from Alzheimers Res Ther. 2019 11:82
# https://alzres.biomedcentral.com/articles/10.1186/s13195-019-0533-9

# Load data from original paper
pilotStudy <- read_excel("Dataset_S04.xlsx", 
                         sheet = "
                         ",
                         range = cell_rows(1:87))
pilotStudy$type <- factor(paste0(pilotStudy$Group, pilotStudy$Normalizeda))


study2 <- read_excel("Dataset_S04.xlsx", 
                     sheet = "CSFClinicalStudyII",
                     range = cell_rows(1:92))



# boxplots for pilot study and clinical study II, using selected peptides

# pilot study: CTSB_80-87, CTSB_210-220, GM2A_89-96, GM2A_170-179 
ADCTSB8087 <- subset(pilotStudy$`CTSB_80-87 (L/H)c`, pilotStudy$type == "ADYes")
ControlCTSB8087 <- subset(pilotStudy$`CTSB_80-87 (L/H)c`, pilotStudy$type == "ControlYes")
ADCTSB210220 <- subset(pilotStudy$`CTSB_210-220 (L/H)c`, pilotStudy$type == "ADYes")
ControlCTSB210220 <- subset(pilotStudy$`CTSB_210-220 (L/H)c`, pilotStudy$type == "ControlYes")
ADGM2A8996 <- subset(pilotStudy$`GM2A_89-96 (L/H)c`, pilotStudy$type == "ADYes")
ControlGM2A8996 <- subset(pilotStudy$`GM2A_89-96 (L/H)c`, pilotStudy$type == "ControlYes")
ADGM2A170179 <- subset(pilotStudy$`GM2A_170-179 (L/H)c`, pilotStudy$type == "ADYes")
ControlGM2A170179 <- subset(pilotStudy$`GM2A_170-179 (L/H)c`, pilotStudy$type == "ControlYes")

# clinical study II: CTSD_80-87, CTSd_210-220
ADCTSD5572 <- subset(study2$`CTSD_55-72 (L/H)a`, study2$Group == "AD")
ControlCTSD5572 <- subset(study2$`CTSD_55-72 (L/H)a`, study2$Group == "Control")
ADCTSD112122 <- subset(study2$`CTSD_112-122 (L/H)a`, study2$Group == "AD")
ControlCTSD112122 <- subset(study2$`CTSD_112-122 (L/H)a`, study2$Group == "Control")


selectedCases <- data.frame(
  Condition = c(rep(c("Control"),length(ControlCTSB8087)), rep(c("AD"),length(ADCTSB8087)), 
                rep(c("Control"),length(ControlCTSB210220)), rep(c("AD"),length(ADCTSB210220)), 
                rep(c("Control"),length(ControlCTSD5572)), rep(c("AD"),length(ADCTSD5572)), 
                rep(c("Control"),length(ControlCTSD112122)), rep(c("AD"),length(ADCTSD112122)),
                rep(c("Control"),length(ControlGM2A8996)), rep(c("AD"),length(ADGM2A8996)),
                rep(c("Control"),length(ControlGM2A170179)), rep(c("AD"),length(ADGM2A170179))),
  LH_ratio = c(ControlCTSB8087, ADCTSB8087, ControlCTSB210220, ADCTSB210220,
               ControlCTSD5572, ADCTSD5572, ControlCTSD112122, ADCTSD112122, 
               ControlGM2A8996, ADGM2A8996, ControlGM2A170179, ADGM2A170179),
  peptide = c(rep(c("CTSB_80-87"), length(ControlCTSB8087) + length(ADCTSB8087)),
              rep(c("CTSB_210-220"), length(ControlCTSB210220) + length(ADCTSB210220)),
              rep(c("CTSD_55-72"), length(ControlCTSD5572) + length(ADCTSD5572)),
              rep(c("CTSD_112-122"), length(ControlCTSD112122) + length(ADCTSD112122)),
              rep(c("GM2A_89-96"), length(ControlGM2A8996) + length(ADGM2A8996)),
              rep(c("GM2A_170-179"), length(ControlGM2A170179) + length(ADGM2A170179)))
)

# Box plot facetted by peptide group
p <- ggboxplot(selectedCases, x = "Condition", y = "LH_ratio",
               color = "Condition", palette = "jco",
               add = "jitter", short.panel.labs = TRUE)

p + facet_wrap(~peptide, scales="free") +
  stat_compare_means(label = "p.format", label.x.npc = 'right', label.y.npc = 'top') + 
  theme(text = element_text(size = 12, face="bold")) +
  theme(panel.border = element_rect(fill=NA,colour="black", size=1))






