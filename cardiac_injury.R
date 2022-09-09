library(DESeq2)
setwd("C:\\Users\\user\\Desktop\\covid card\\1")
rt = read.table("covid_cardgse184715.txt",header = T,row.names = 1)

colData <- data.frame(row.names = colnames(rt),
                      condition =
                        factor(c("Mock","Mock","Mock","Infection","Infection","Infection"),
                               levels = c("Mock","Infection")))
colData
dds <- DESeqDataSetFromMatrix(countData = round(rt), colData = colData,
                              design = ~ condition)


nrow(dds)
dds <- dds[rowSums(counts(dds))>1,]
vsd <- vst(dds, blind = FALSE)
sampleDists <- dist(t(assay(vsd)))
hc <- hclust(sampleDists, method = "ward.D2")
plot(hc, hang = -1)

dds <- DESeq(dds)

normalized_counts <- as.data.frame(counts(dds, normalized=TRUE))


contrast <- c("condition", "Infection", "Mock")
dd1 <- results(dds, contrast=contrast, alpha = 0.05)
plotMA(dd1, ylim=c(-2,2))

dd2 <- lfcShrink(dds, contrast=contrast, res=dd1)
plotMA(dd2, ylim=c(-2,2))

resApe <-lfcShrink(dds, coef=2,type="apeglm")
plotMA(resApe, ylim = c(-10,10))

summary(resApe, alpha = 0.05)
library(dplyr)
library(tibble)
res <- resApe %>% 
  data.frame() %>% 
  rownames_to_column("gene_id")
library(AnnotationDbi)
library(org.Hs.eg.db)

res$entrez <- mapIds(org.Hs.eg.db,
                     keys=res$gene_id,
                     column="ENTREZID",
                     keytype="SYMBOL",
                     multiVals="first")

library(dplyr)
gene_df <- res %>% 
  dplyr::select(gene_id,log2FoldChange,entrez) %>% 
  
  filter(entrez!="NA") %>% 
  
  distinct(entrez,.keep_all = T)

geneList <- gene_df$log2FoldChange
names(geneList) = gene_df$entrez
geneList = sort(geneList, decreasing = TRUE)
head(geneList)
library(clusterProfiler)
library(stats)
gseaKEGG <- gseKEGG(geneList     = geneList,
                    organism     = 'hsa',
                    nPerm        = 1000,
                    minGSSize    = 20,
                    pvalueCutoff = 0.05,
                    verbose      = FALSE)
library(ggplot2)
pdf("gsea.pdf",height = 10,width = 14)
dotplot(gseaKEGG,showCategory=10,split=".sign")+facet_grid(~.sign)
dev.off()



sizeFactors(dds)

res <- results(dds)
head(res)
class(res)
res <- as.data.frame(res)
res <- cbind(rownames(res), res)

colnames(res) <- c("gene_id", "baseMean", "log2FoldChange", "lfcSE", "stat",
                   "pval", "padj")

write.table(res, "COVID-19gse184715.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE,
            na = "")


res = read.table("deglist.txt",sep = "\t",header = T,row.names = 1)

resSig <- res[which(res$padj < 0.001 & abs(res$log2FoldChange) > 2),]
resSig[which(resSig$log2FoldChange > 0), "up_down"] <- "Up"
resSig[which(resSig$log2FoldChange < 0), "up_down"] <- "Down"
gene_dd <- resSig %>% 
  dplyr::select(gene_id,baseMean,log2FoldChange,lfcSE,pvalue,padj,up_down,entrez) %>% 
  
  filter(entrez!="NA") %>% 
 
  distinct(entrez,.keep_all = T)

write.table(resSig, "Deg.genefc2wuNA.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")


dr =res

res = na.omit(res)

 res$type = ifelse(res$log2FoldChange >= 2 & res$padj <=0.001,"up-regulated", ifelse(res$log2FoldChange<= -2 & res$padj <=0.001,"down-regulated","none siginificant"))

 table(res$type)

 
 title = paste0('The number of up gene is ', nrow(res[res$type == 'up-regulated',]),'\nThe number of down gene is ', nrow(res[res$type == 'down-regulated',]))

 res = as.data.frame(res)
library(ggplot2)
 pdf("vol.pdf",width = 10 ,height = 8)
 ggplot(res, aes(x = log2FoldChange, y = -log10(padj)))+geom_point(aes(color = type)) + scale_color_manual(values = c("red","grey","green"),limits = c('up-regulated','none siginificant',"down-regulated")) + theme_bw(base_size = 20) + theme(plot.title = element_text(size=15,hjust = 0.5) + theme_classic())
dev.off()
 p = ggplot(res, aes(log2FoldChange, -log10(padj),color= type))
 
 p +geom_point()
 
  x_lim = max(res$log2FoldChange, - res$log2FoldChange)
 
  gg = p + geom_point(size =1) + xlim(-x_lim, x_lim) +labs(x= "log2FoldChange", y = "-log10(padj)") +scale_color_manual(values = c("#A52A2A","grey","#f8766d")) + geom_hline(aes(yintercept=-1*log10(0.05)),colour="black", linetype="dashed") +     geom_vline(xintercept=c(-2,2),colour="black", linetype="dashed")
 
  print(gg)
 
  
  library(stringr)
  library(enrichplot)
  library(clusterProfiler)
  library(GOplot)
  library(DOSE)
  library(ggnewscale)
  library(topGO)
  library(circlize)
  library(ComplexHeatmap)
  library(enrichplot)
  library(GOSemSim)
  library(DOSE)
  

  GO_database <- 'org.Hs.eg.db' 
  
  
  gene <- bitr(resSig$gene_id,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)  
  
  rownames(gene) = gene$SYMBOL
  ht = intersect(rownames(gene),rownames(resSig))
  ty = resSig[ht,]
  gene = gene[ht,]
  yy = cbind(ty,gene)
  write.table(yy,"3097genes.txt",sep = '\t',row.names = T,col.names = T)
  yy = read.table("3097genes.txt",sep = "\t",header = T,row.names = 1)
  
  GSEA_input <- yy$log2FoldChange
  names(GSEA_input) = yy$ENTREZID
  GSEA_input = sort(GSEA_input, decreasing = TRUE)
  GSEA_KEGG <- gseKEGG(GSEA_input, organism = KEGG_database, pvalueCutoff = 0.05)#GSEA富集分析
  
  gene_dd <- resSig %>% 
    dplyr::select(gene_id,baseMean,log2FoldChange,lfcSE,pvalue,padj,entrez) %>% 
    
    filter(entrez!="NA") %>% 
    
    distinct(entrez,.keep_all = T)
  
  
  
GO<-enrichGO( gene_dd$entrez,
                OrgDb = GO_database,
                keyType = "ENTREZID",
                ont = "ALL",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = T)

KEGG<-enrichKEGG(gene_dd$entrez,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)

pdf("Gov2.pdf",width = 12,height = 18)
dotplot(GO,showCategory = 15,split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
dev.off()

pdf("Keg1.pdf",width = 10,height = 8)
dotplot(KEGG)
dev.off()
enrichplot::cnetplot(GO,circular=FALSE,colorEdge = TRUE)#基因-通路关联网络图
enrichplot::cnetplot(KEGG,circular=FALSE,colorEdge = TRUE)#circlua


GO2 <- pairwise_termsim(GO)
KEGG2 <- pairwise_termsim(KEGG)
pdf("Goplot1.pdf",width = 20,height = 14)
enrichplot::emapplot(GO2,showCategory = 50, color = "p.adjust", layout = "kk")#通路间关联网络图
dev.off()
pdf("Kegg1.pdf",width = 20,height = 14)
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")
dev.off()
gseaplot2(GSEA_KEGG,1:26)


cp = as.data.frame(GO@result)
write.table(cp, "GOresult.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")


ct = as.data.frame(KEGG@result)
write.table(ct, "KEGGresult.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")

cl = as.data.frame(gseaKEGG@result)
write.table(cl, "GSEAresult.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")


write.table(resSig, "DEGSresult2705.xls",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE, na = "")
