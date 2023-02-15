if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("arrayQualityMetrics") 

install.packages(c("ggplot", "illuminaHumanv3.db", "genefilter", "limma","RColorBrewer","enrichR","pathview","data.table"))



library(arrayQualityMetrics)
library(illuminaHumanv3.db)
library(genefilter)
library(limma)
library(gplots)
library(RColorBrewer)
library(enrichR)
library(pathview)
library(data.table)


# read data from illumina idat files

idatfiles <- list.files(,pattern="idat")
bgxfile <- list.files(,pattern="bgx")
raw.data<-read.idat(idatfiles, bgxfile)
raw.data$other$Detection<-detectionPValues(raw.data)


# normalization of Illumina data

norm.data<-neqc(raw.data)
i<-which(duplicated(norm.data$genes$Probe_Id))
norm.data<-norm.data[-i,]
rownames(norm.data$E)<-norm.data$genes$Probe_Id
tt<-norm.data$E

# QC report

eset<-ExpressionSet(assayData=tt,annotation="illuminaHumanv3")
QC_report_nor<-arrayQualityMetrics(eset)

# removing of non-expressed illumina probes

expressed <- rowSums(norm.data$other$Detection < 0.05) >= 1
norm.data <- norm.data[expressed,]
dim(norm.data)
tt<-norm.data$E

# additional filtering of Illumina probes

eset<-ExpressionSet(assayData=tt,annotation="illuminaHumanv3")
eset<-nsFilter(eset,require.entrez=T,
               remove.dupEntrez=F,
               var.filter=F)
eset<-eset$eset
tt<-exprs(eset)
dim(tt)

fwrite(as.data.frame(tt),file="matrix.csv",row.names=T)




# identification of DEGs (case of two groups)

eset<-fread("matrix.csv",data.table=F)
rownames(eset)<-eset$V1
eset<-eset[,-1]

samples<-colnames(eset)
samples<-sapply(1:length(samples), function(x) {
  strsplit(samples[x],split="_",fixed=T)[[1]][1]
})
colnames(eset)<-samples

groups <- read.delim("Phenotype_2.txt")
i<-match(colnames(eset),groups$ID)
groups<-groups[i,]
groups<-groups$three
table(groups)
groups <- factor(groups,levels=c("less_50","more_50"))

?model.matrix
design <- model.matrix(~groups,eset)
colnames(design)<-c("less_50","more_50vsless_50")
fit <- lmFit(eset,design)
efit <- eBayes(fit)
tt <- topTable(efit, coef="more_50vsless_50", n=500000, sort.by="p")

tt <- cbind(Symbol = unlist(mget(rownames(tt), illuminaHumanv3SYMBOL)), EntrezID = unlist(mget(rownames(tt), illuminaHumanv3ENTREZID)), Description = unlist(mget(rownames(tt), illuminaHumanv3GENENAME)),tt)
tt<-na.omit(tt)


write.table(tt,file="DEG limma.txt", row.names=T,sep='\t',quote=FALSE)


################################################

# pathway enrichment analysis
?listEnrichrDbs
dbs <- listEnrichrDbs()

#deg3<-unique(tt$Symbol[tt$logFC>0.7 & tt$adj.P.Val<0.1])
#deg4<-unique(tt$Symbol[tt$logFC<(-0.7) & tt$adj.P.Val<0.1])

deg<-unique(tt$Symbol[tt$logFC>0.5 & tt$P.Value<0.05]) 
deg2<-unique(tt$Symbol[tt$logFC<(-0.5) & tt$P.Value<0.05])

write.table(deg,file="giper.txt", row.names=T,sep='\t',quote=FALSE)
write.table(deg2,file="gipo.txt", row.names=T,sep='\t',quote=FALSE)

enriched <- enrichr(deg, "GO_Biological_Process_2021")
enriched2 <- enrichr(deg2, "GO_Biological_Process_2021")

dep_bio<-enriched[[1]]
dep_bio<-dep_bio[dep_bio$Adjusted.P.value<0.1,] 

dep_bio2<-enriched2[[1]]
dep_bio2<-dep_bio2[dep_bio2$Adjusted.P.value<0.1,]

write.table(dep_bio,file="dep_bio.txt", row.names=T,sep='\t',quote=FALSE)
write.table(dep_bio2,file="dep_bio2.txt", row.names=T,sep='\t',quote=FALSE)


enriched_Kegg <- enrichr(deg, "KEGG_2021_Human")
enriched_Kegg2 <- enrichr(deg2, "KEGG_2021_Human")

dep<-enriched_Kegg[[1]]
dep<-dep[dep$Adjusted.P.value<0.1,] 

#dep2<-enriched_Kegg2[[1]]
#dep2<-dep2[dep2$Adjusted.P.value<0.1,]


#creation of KEGG pathway view

fold<-tt[tt$P.Value<0.05 & abs(tt$logFC)>0.7,4]
names(fold)<-tt[tt$P.Value<0.05 & abs(tt$logFC)>0.7,2]

#a

pv.out <- pathview(gene.data = fold, pathway.id = "04612",
                   species = "hsa", out.suffix = "Antigen processing and presentation ", kegg.native = T)

fold<-tt$logFC
names(fold)<-tt$EntrezID
pv.out <- pathview(gene.data = fold, pathway.id = "04612",
                   species = "hsa", out.suffix = "Antigen processing and presentation 1 ", kegg.native = T)

#b
pv.out <- pathview(gene.data = fold, pathway.id = "04350",
                   species = "hsa", out.suffix = "TGF-beta signaling pathway", kegg.native = T)

fold<-tt$logFC
names(fold)<-tt$EntrezID
pv.out <- pathview(gene.data = fold, pathway.id = "04350",
                   species = "hsa", out.suffix = "TGF-beta signaling pathway 1 ", kegg.native = T)

#c

pv.out <- pathview(gene.data = fold, pathway.id = "05332",
                   species = "hsa", out.suffix = "Graft-versus-host disease", kegg.native = T)

fold<-tt$logFC
names(fold)<-tt$EntrezID
pv.out <- pathview(gene.data = fold, pathway.id = "05332",
                   species = "hsa", out.suffix = "Graft-versus-host disease 1 ", kegg.native = T)


#d

pv.out <- pathview(gene.data = fold, pathway.id = "04650",
                   species = "hsa", out.suffix = "Natural killer cell mediated cytotoxicity", kegg.native = T)

fold<-tt$logFC
names(fold)<-tt$EntrezID
pv.out <- pathview(gene.data = fold, pathway.id = "04650",
                   species = "hsa", out.suffix = "Natural killer cell mediated cytotoxicity 1 ", kegg.native = T)


#2
deg2<-unique(tt$Symbol[tt$logFC<(-0.5) & tt$P.Value<0.05])
enriched_Kegg <- enrichr(deg2, "KEGG_2021_Human")

dep<-enriched_Kegg[[1]]
dep<-dep[dep$Adjusted.P.value<0.1,] 

pv.out <- pathview(gene.data = fold, pathway.id = "00830",
                   species = "hsa", out.suffix = "Retinol metabolism ", kegg.native = T)

fold<-tt$logFC
names(fold)<-tt$EntrezID
pv.out <- pathview(gene.data = fold, pathway.id = "00830",
                   species = "hsa", out.suffix = "Retinol metabolism 1 ", kegg.native = T)

#creation of heatmap


#dat <- eset[rownames(tt)[tt$adj.P.Val<0.1 & abs(tt$logFC)>0.7],]
dat <- eset[rownames(tt)[tt$P.Value<0.05 & abs(tt$logFC)>0.5],]

rbpal <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
rbpal <- rev(rbpal)

pat.col<-as.character(groups)
pat.col[pat.col=="more_50"] <- "red"
pat.col[pat.col=="less_50"] <- "cyan"

hclustfun<-function(x) hclust(x, method="ward.D2")
distfun<-function(x) as.dist((1-cor(t(x)))/2)

#nm<-tt$Symbol[tt$adj.P.Val<0.1 & abs(tt$logFC)>0.7]
nm<-tt$Symbol[tt$P.Value<0.05 & abs(tt$logFC)>0.5]

heatmap.2(data.matrix(dat), col=rbpal, 
          trace="none",density.info="none",
          ColSideColors=pat.col,
          labRow=nm,
          scale="row",
          hclustfun=hclustfun,
          distfun=distfun,
          margins = c(5, 5),
          cexCol=0.5,
          cexRow = 0.4)

