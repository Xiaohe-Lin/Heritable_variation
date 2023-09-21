### transcriptom analysis in R                                                                                ###                                                       
###---------------------------------------------------------------------------------------------------------- ### 
#       Input files(sample information including experimental design and expression matrix : .csv);           ###      
### R packages required: DESeq2, GenomicRanges, tidyverse, clusterProfiler, grid, patchwork                   ###                                        
###-----------------------------------------------------------------------------------------------------------###

# 1. Analysis of differential expressed genes (DEGs) using DEseq2----
library(DESeq2)
ecotype <- c('Abd','Ang','Tol','TRE')
generation <- c('F1','F2','F3','F4')
for (i in ecotype) {
  for (j in generation) {
    cts <- as.matrix(read.csv(paste0(i,'/',i,'.counts.',j,'.csv'),row.names=1,header=T))
    coldata <- read.csv(paste0(i,'/',i,'.info.',j,'.csv'), row.names=1,header = T)
    coldata$treat <- factor(coldata$treat,levels = c('control','Cd1','Cd2','Drought'))
    coldata$batch <- factor(coldata$batch)
    coldata <- coldata[order(coldata$treat),]
    cts <- cts[,rownames(coldata)]
    cpm <- apply(cts,2, function(x) { x/sum(x)*1000000 })
    cts.retain <- cts[rowSums(cpm[,c(1:6)]) > 1 | rowSums(cpm[,c(1:3,7:9)]) > 1 | rowSums(cpm[,c(1:3,10:12)]) > 1, ]
    dim(cts.retain)
    dds <- DESeqDataSetFromMatrix(countData = cts.retain,
                                  colData = coldata,
                                  design = ~ treat+batch)
    
    dds$treat <- relevel(dds$treat, ref ="control")
    dim(dds)
    suppressMessages(dds2 <- DESeq(dds)) 
    dds.name <- c(resultsNames(dds2)[2:4])
    dds.name
    for (m in dds.name) {
      dds.result <- results(dds2, name=m)
      resOrdered <- dds.result[order(dds.result$pvalue),]
      resSig <- subset(resOrdered,pvalue<0.05 & abs(log2FoldChange) > 1)
      k=gsub('treat_','',m)
      k=gsub('_vs_control','',k)
      write.csv(as.data.frame(resOrdered),file=paste0(i,'/DEG/',i,'.',j,'.',k,'.all.csv'))
      write.csv(as.data.frame(resSig),file=paste0(i,'/DEG/',i,'.',j,'.',k,'.remove.csv'))
    }
  }
}

# 2. GO enrichment ----
library(clusterProfiler)
library(tidyverse)
ecotype <- c('Abd','Ang','Tol','TRE')
treatment <- c('Cd1','Cd2','Drought')
for (i in ecotype) {
  for (j in treatment) {
    F1 <- read.csv(paste0(i,'/DEG/',i,'.F1.',j,'.remove.csv'),header=T)
    F2 <- read.csv(paste0(i,'/DEG/',i,'.F2.',j,'.remove.csv'),header=T)
    F3 <- read.csv(paste0(i,'/DEG/',i,'.F3.',j,'.remove.csv'),header=T)
    F4 <- read.csv(paste0(i,'/DEG/',i,'.F4.',j,'.remove.csv'),header=T)
    all <- rbind(F1,F2,F3,F4)
    all.sub <- as.data.frame(table(all$X))
    all.sub <- all.sub[all.sub$Freq>1,] # extracting genes at least expressed in two generations
    assign(paste0(i,'.',j),unique(all.sub$Var1))  
  }
}

for (i in ecotype) {
  Cd1 <- get(paste0(i,'.Cd1')) 
  Cd2 <- get(paste0(i,'.Cd2'))
  Drought <- get(paste0(i,'.Drought'))
  gene.list = list(Cd1 = Cd1, Cd2 = Cd2, Drought = Drought)  
  gene.GO <- compareCluster(gene.list, fun="enrichGO",ont = "BP", OrgDb = "org.At.tair.db",
                            keyType = "TAIR", pvalueCutoff=1,pAdjustMethod = "BH",qvalueCutoff = 1) 
  gene.GO.cluster <- clusterProfiler::simplify(gene.GO,cutoff=0.7,by='pvalue',select_fun=min)
  gene.GO.filter <- gene.GO.cluster[gene.GO.cluster@compareClusterResult$pvalue<0.05,]
  write.csv(gene.GO.filter,file = paste0('GO/',i,'GO.union.csv'))  
}

# 3. Transposable element sampling ----
library(GenomicRanges)
library(grid)
library(patchwork)

# Obtainning heritable DEGs
ecotype <- c('Abd','Ang','Tol','TRE')
generation <- c('F1','F2','F3','F4')
treatment <- c('Cd1','Cd2','Drought')
DEG.all <- c()
for (i in ecotype) {
  for (j in generation) {
    for (k in treatment) {
      a <- read.csv(paste0(i,'/DEG/',i,'.',j,'.',k,'.remove.csv'),header=T)
      a$ecotype <- i
      a$Generation <- j
      a$treatment <- k
      DEG.all <- rbind(DEG.all,a)
    }
  }
}

heritable_DEG <- DEG.all %>%
  group_by(ecotype, treatment, X) %>%
  summarise(DEG_number = n()) %>%
  filter(DEG_number>1) 

heritable_overlap_DEGs <- as.data.frame(table(heritable_DEG$X))
heritable_overlap_DEGs <- heritable_overlap_DEGs[heritable_overlap_DEGs$Freq>1,]

# Read and preprocess TE classification data
TE.anno <- read.table('TAIR10_Transposable_Elements.txt',header = T)
TE.anno <- TE.anno[,c(1,5,6)]
TE.anno[grep('LINE',TE.anno$Transposon_Super_Family),]$Transposon_Super_Family <- 'LINE'
TE.anno[grep('Rath',TE.anno$Transposon_Super_Family),]$Transposon_Super_Family <- 'SINE'
table(TE.anno$Transposon_Super_Family)
TE.1 <- data.frame(Var1=unique(TE.anno$Transposon_Super_Family))

# Read and preprocess TE and gene annotation gff
arab.gff<-read.table("TAIR10_GFF3_genes_transposons.gff", sep="\t",header=F)
names(arab.gff)<-c("seqnames", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")
arab.gff <- arab.gff[arab.gff$seqnames%in%c('Chr1','Chr2','Chr3','Chr4','Chr5'),]
gr<-makeGRangesFromDataFrame(arab.gff,keep.extra.columns=T)
gr<-renameSeqlevels(gr,c("1","2","3","4","5"))
gr<-gr[strand(gr)!="*"]
genome(gr)<-"TAIR10"
gr.TE<-gr[gr$feature=='transposable_element']
gr.TE$attributes <- str_split_fixed(gr.TE$attributes,';',3)[,1]
gr.TE$attributes <- gsub('ID=','',gr.TE$attributes)

arab.gene.gff <- arab.gff[arab.gff$feature=='gene',]
arab.gene.gff$start <- arab.gene.gff$start-2000
gr.gene <-makeGRangesFromDataFrame(arab.gene.gff,keep.extra.columns=T)
gr.gene<-renameSeqlevels(gr.gene,c("1","2","3","4","5"))
gr.gene<-gr.gene[strand(gr.gene)!="*"]
genome(gr.gene)<-"TAIR10"
gr.gene <- gr.gene[gr.gene$feature=='gene']
gr.gene$attributes <- str_split_fixed(gr.gene$attributes,';',3)[,1]
gr.gene$attributes <- gsub('ID=','',gr.gene$attributes)

rm(gr)
rm(arab.gff)
rm(arab.gene.gff)

#  Data subset selection
heritable_overlap_DEGs.uniq <- unique(heritable_overlap_DEGs$Var1)
length(heritable_overlap_DEGs.uniq)
gr.data <- gr.gene[gr.gene$attributes %in% heritable_overlap_DEGs.uniq,]
gr.data.anno.TE <- findOverlaps(gr.data,gr.TE)
gr.data.TE.sub <- data.frame(gr.TE[gr.data.anno.TE@to,])
gr.data.TE.sub <- gr.data.TE.sub[,c(1:5,7,10)]
colnames(gr.data.TE.sub)[7] <- colnames(TE.anno)[1]
TE.sub <- left_join(gr.data.TE.sub,TE.anno,by=colnames(TE.anno)[1])
TE.sub.1 <- data.frame(table(TE.sub$Transposon_Super_Family))
TE.sub.1 <- left_join(TE.1,TE.sub.1,by='Var1')
TE.sub.1[is.na(TE.sub.1)] <- 0


# Extracting genes with expression levels in treatments of low cadmium 
df_all <- c()
ecotype <- c('Abd','Ang','Tol','TRE')
generation <- c('F1','F2','F3','F4')
for (i in ecotype) {
  for (j in generation) {
    df <- read.csv(paste0(i,'/DEG/',i,'.',j,'.Cd1.all.csv'),header=T,row.names = 1)
    df$gene <- rownames(df)
    df <- df[,c('gene','pvalue')]
    df_all <- rbind(df_all,df)
  }
}  

df.stat <- as.data.frame(table(df_all$gene))
df.stat <- df.stat[df.stat$Freq==16,]
gene.all <- unique(df.stat$Var1)
rm(df_all)
rm(df)
rm(df.stat)

# Gene sampling and TE analysis
TE.result <- c()
for (i in 1:10000) {
  gene.sub <- sample(gene.all,size = 2669)
  gr.data <- gr.gene[gr.gene$attributes %in% gene.sub,]
  gr.data.anno.TE <- findOverlaps(gr.data,gr.TE)
  gr.data.TE.sub <- data.frame(gr.TE[gr.data.anno.TE@to,])
  gr.data.TE.sub <- gr.data.TE.sub[,c(1:5,7,10)]
  colnames(gr.data.TE.sub)[7] <- colnames(TE.anno)[1]
  TE.aa <- left_join(gr.data.TE.sub,TE.anno,by=colnames(TE.anno)[1])
  TE.aa.1 <- as.data.frame(table(TE.aa$Transposon_Super_Family))
  TE.aa.1 <- left_join(TE.1,TE.aa.1,by='Var1')
  if (i==1) {
    TE.result <- TE.aa.1
  }else{
    TE.result <- cbind(TE.result,TE.aa.1[,-1])
  }
}

colnames(TE.result)[2:10001] <- 1:10000
TE.result[is.na(TE.result)] <- 0
TE.var <- TE.result$Var1
TE.var==TE.sub.1$Var1

# Calculation of probability
for (i in TE.var) {
  df <- TE.result[TE.result$Var1 == i,]
  s1 <- as.data.frame(as.numeric(df[1,2:10001]))
  head(s1)
  colnames(s1) <- 'a'
  observed.value <- TE.sub.1[TE.sub.1$Var1==i,2]
  s2 <- s1[s1$a > observed.value,]
  pvalue <- length(s2)/10000
  grob <- grobTree(textGrob(paste0('P = ',pvalue),
                            x=0.6,y=0.9, hjust=0,
                            gp=gpar(col="black", fontsize=8)))
  assign(paste0('p_',i),ggplot(data = s1,aes(x=a))+
           geom_density(aes(y=..density..),bw = 1,linetype=2,color='black',fill='#fcaf7c',alpha=0.8,size=1.2)+
           coord_cartesian(xlim=c(0,NA))+
           annotation_custom(grob)+
           geom_vline(xintercept = as.numeric(observed.value),color='#04686b')+
           xlab(i)+
           theme_test())
}
