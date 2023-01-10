# ROC for each simulations 
# Use simulation 1 as example and try the others
# Results generated for the simulations of expecting detection of HVGs
# Is these HVGs ranked highly in a method

library(Seurat)
library(scran)
library(dplyr)
library(patchwork)
# high + low BCV
# Low BCV + High BCV
dat1 <- read.csv("~/simulation_highBCV_dropout.csv", stringsAsFactors = 1, row.names = 1)# 2000 genes with 400 cells

dat2 <- read.csv("~/simulation_lowBCV_dropout.csv",stringsAsFactors = 1, row.names = 1) # 10000 genes with 400 cells

dat3 <- read.csv("~/simulation_moderateBCV_dropout.csv",stringsAsFactors = 1, row.names = 1)# 2000 genes with 400 cells



# For smaller dataset, random 200 high bcv and 800 low bcv
set.seed(123)
dat3.index <- sample(2000,200)
dat3 <- dat3[dat3.index,]
dat2.index <- sample(10000, 800)
dat2 <- dat2[dat2.index,]

# Change gene names for first dataset
rownames(dat3) <- paste("Gene",c(801:1000),sep = "")
rownames(dat2) <- paste("Gene",c(1:800),sep = "")

dat <- rbind(dat2,dat3)

#Sim1 <- read.csv("D:/PhD/Trip to China/Simulation/Simulation_NB.csv", row.names = 1)

# Increase the dropout rate by round up 
# sim1.1 <- round(Sim1*0.15) == 0
# Sim1[round(Sim1*0.15) == 0] <- 0
Sim1.obj <- CreateSeuratObject(dat)
Sim1.obj <- NormalizeData(Sim1.obj)
Sim1.norm <- GetAssayData(object = Sim1.obj, slot = "data")
Sim1.norm <- as.data.frame(as.matrix(Sim1.norm))
Sim1.norm <- Sim1.norm[rowSums(Sim1.norm)!=0,]
sum(Sim1.norm==0)/(988*400) # 68.6% 
 
# Uses what TMS did for basic statistics 
# mean 
Sim1_mean <- rowMeans(Sim1.norm)

# SD
Sim1_sd <- apply(Sim1.norm,1,function(x) sd(x))

#IQR -- Interquartile range 
Sim1_iqr <- apply(Sim1.norm,1,function(x) IQR(x))

#MAD -- Median absolute deviation 
Sim1_mad <- apply(Sim1.norm,1,function(x) mad(x))

# CV -- coefficient of variation 
Sim1_CV <- Sim1_sd/Sim1_mean

# var
Sim1_var <- apply(Sim1.norm,1,function(x) var(x)) 

# Fano Factor
Sim1_FF <- Sim1_var/Sim1_mean

# edgeR
Sim1.raw <- as.data.frame(as.matrix(dat))
Sim1.raw <- Sim1.raw[rowSums(Sim1.raw)!=0,] 
library(edgeR)
y <- DGEList(Sim1.raw)
y <- estimateDisp(y,robust=TRUE)
plotBCV(y)
Sim1_bcv <- sqrt(y$tagwise.dispersion)
plot(y$AveLogCPM,Sim1_bcv)
names(Sim1_bcv) <- rownames(Sim1.raw)

library("DESeq2")
# Normal way of doing it
coldata <- data.frame(name = rep("sim1",400))
dds_Sim1 <- DESeqDataSetFromMatrix(dat, colData = coldata, design = ~1)
dds_Sim1 <- estimateSizeFactors(dds_Sim1, type = 'iterate')
dds_Sim1 <- estimateDispersions(dds_Sim1)
plotDispEsts(dds_Sim1)
hist(log1p(dispersions(dds_Sim1)[!is.na(dispersions(dds_Sim1))]))
DESeq2_glm <- log1p(dispersions(dds_Sim1))
names(DESeq2_glm) <- rownames(dat)
# One specifically design for sc RNA-seq
library("glmGamPoi")
dds_Sim1_glmgam <- estimateDispersions(dds_Sim1, fitType = "glmGamPoi")
plotDispEsts(dds_Sim1_glmgam)
plot(log1p(mcols(dds_Sim1_glmgam)$baseMean),log1p(mcols(dds_Sim1_glmgam)$dispGeneEst))
DESeq2_glmgam <- log1p(mcols(dds_Sim1_glmgam)$dispGeneEst)
names(DESeq2_glmgam) <- rownames(dat)

Sim1_lcv.df <- cbind(Sim1_mean,Sim1_CV)
colnames(Sim1_lcv.df) <- c("Mean","CV")
# ordered cv by mean expression
Sim1_lcv.df <- Sim1_lcv.df[order(Sim1_lcv.df[,1]),]
# an example based on paper suggested parameter - follow python package # with 
LCV_Sim1 <- c()
window_size = 50
for (k in 1:dim(Sim1_lcv.df)){
  if (k <= window_size/2) 
    pos <-  match(Sim1_lcv.df[k,2],Sim1_lcv.df[(1:window_size),2])
  else if (k >window_size/2 & k< dim(Sim1_lcv.df)[1]- window_size/2){
    window_order <- Sim1_lcv.df[c((k-window_size/2):(k+window_size/2-1)),2]
    window_order <- window_order[order(window_order)]
    
    pos <- match(Sim1_lcv.df[k,2], window_order)
  }
  else pos <- match(Sim1_lcv.df[k,2],Sim1_lcv.df[c((dim(Sim1_lcv.df)[1]-window_size+1):(dim(Sim1_lcv.df)[1])),2])
  lcv_value <- pos/(window_size)*100
  names(lcv_value) <- rownames(Sim1_lcv.df)[k]  
  LCV_Sim1 <- c(LCV_Sim1, lcv_value)
  
}
hist(LCV_Sim1)
LCV_Sim1 <- LCV_Sim1[match(names(Sim1_sd),names(LCV_Sim1))]

Sim1_seurat <- FindVariableFeatures(Sim1.obj,selection.method =  "mvp")
Sim1_seurat_mvp <- Sim1_seurat@assays[["RNA"]]@meta.features[["mvp.dispersion"]]
Sim1_seurat_mvp_mean <- Sim1_seurat@assays[["RNA"]]@meta.features[["mvp.mean"]]

Sim1_seurat.1 <- FindVariableFeatures(Sim1.obj,selection.method =  "vst",nfeatures =20138 )
Sim1_seurat_vst <- Sim1_seurat.1@assays[["RNA"]]@meta.features[["vst.variance.standardized"]]
Sim1_seurat_vst_mean <- Sim1_seurat.1@assays[["RNA"]]@meta.features[["vst.mean"]]
plot(Sim1_seurat_mvp_mean,Sim1_seurat_mvp)
plot(log1p(Sim1_seurat_vst_mean),Sim1_seurat_vst)
# Extract variable feature information
Seurat_mvp <- HVFInfo(Sim1_seurat,selection.method =  'mvp')$dispersion.scaled
names(Seurat_mvp) <- rownames(HVFInfo(Sim1_seurat,selection.method =  'mvp'))

Seurat_vst <- HVFInfo(Sim1_seurat.1,selection.method =  'vst')$variance.standardized
names(Seurat_vst) <- rownames(HVFInfo(Sim1_seurat.1,selection.method =  'vst'))
# remove completely no expression genes - same as others
Seurat_mvp <- Seurat_mvp[HVFInfo(Sim1_seurat,selection.method =  'mvp')$mean > 0 ]
# remove completely no expression genes - same as others
Seurat_vst <- Seurat_vst[HVFInfo(Sim1_seurat.1,selection.method =  'vst')$mean > 0 ]

library(scran)
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays = list(counts = dat))
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
dec <- modelGeneVar(sce)
plot(dec$mean, dec$bio, xlab="Mean log-expression", ylab="Variance")
curve(metadata(dec)$trend(x), col="blue", add=TRUE)
Sim1_scran <- dec$bio
names(Sim1_scran) <- dec@rownames

Sim1_DM <- NA
names(Sim1_scran) <- str_replace(names(Sim1_scran),pattern="-",replacement="-")


Sim1_measurements <- list(Sim1_sd, Sim1_iqr, Sim1_mad,Sim1_CV, Sim1_FF, Sim1_bcv,DESeq2_glm,DESeq2_glmgam,Sim1_DM,LCV_Sim1,Seurat_mvp,Seurat_vst,Sim1_scran)
names(Sim1_measurements) <- c("sd","iqr","mad","CV","FF","edgeR","DESeq2_glm","DESeq2_glmgam",'DM',"LCV","Seurat_mvp","Seurat_vst","scran")


# ROC curve 
library("stringr")
library(caret)
HVG_labels <- rownames(Sim1[801:1000,])
HVG_labels <- str_replace(HVG_labels,pattern="_",replacement="")
Accuracy1
Accuracy2 
Accuracy3
Accuracy4 
acc <- rbind(Accuracy1,Accuracy2,Accuracy3,Accuracy4)
acc <- as.data.frame(acc)
acc <- acc[,-9]
boxplot(acc)
sapply(Sim1_measurements,function(metric)length(intersect(names(metric[order(metric)[801:1000]]),HVG_labels)))

library(reshape2)
library(ggplot2)
acc <- melt(acc)
acc$variable <- as.factor(acc$variable)
levels(acc$variable) <- c('SD',"IQR",'MAD',"CV","FF","edgeR","DESeq2","glmGamPoi","LCV","Seurat_mvp","Seurat_vst",'scran')
acc$variable <- factor(acc$variable, levels = c('SD',"IQR",'MAD',"CV","FF","LCV","DESeq2","edgeR","glmGamPoi",'scran',"Seurat_mvp","Seurat_vst"))
saveRDS(acc, "Performance_on_simulation.rds")
ggplot(acc, aes(x = variable, y = value/200, fill= category)) + 
  geom_boxplot(outlier.shape = NA) +scale_fill_manual(values=c("#9DC3E6", "#FF8282","#F5CDFF"))+
  geom_jitter()+ theme_classic()+ylab("Rediscovery rate")+xlab("")
