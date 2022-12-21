#Pathway analysis SHAM + Zymozan vs VNS + Zymozan
intGroups <- c("GROUP_B2", "Reference_2")
intGroups <- paste(intGroups[1], "vs", intGroups[2], sep = " ")

#GO
library(clusterProfiler)
library(org.Mm.eg.db)

#Select genes with lFC > 1
vnsGeneIds <- res[[intGroups]] %>% as.data.frame() %>% 
  .[.$padj < 0.05 & .$log2FoldChange >= 1,] %>% .[!is.na(.$padj),] %>% rownames(.) %>%
  AnnotationDbi::select(org.Mm.eg.db, keys = ., columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENSEMBL")
vnsGeneIds$LFC <- res[[intGroups]] %>% as.data.frame() %>% .[match(vnsGeneIds$ENSEMBL, rownames(.)), "log2FoldChange"]

vnsGO <- enrichGO(gene          = unique(vnsGeneIds$ENTREZID),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "ALL",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)

pdf(file.path(figDir, "VNS_GOenrichmentBar.pdf"), width = 8, height = 6)
barplot(vnsGO)
dev.off()

#Gene Set Enrichment Analysis of Gene Ontology
vnsGeneData <- as.data.frame(res[[intGroups]])
vnsGeneData$ENSEMBL <- rownames(vnsGeneData)

vnsGeneData <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(vnsGeneData), columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENSEMBL") %>% 
  dplyr::left_join(vnsGeneData, ., by = "ENSEMBL") 

vnsIntGeneLFC <- vnsGeneData %>% .[, c("SYMBOL", "ENTREZID", "log2FoldChange")]

vnsGeneLFC <- vnsIntGeneLFC$log2FoldChange
names(vnsGeneLFC) <- vnsIntGeneLFC$ENTREZID
vnsGeneLFC %<>% .[!duplicated(names(.))]
vnsGeneLFC %<>%  sort(., TRUE)
vnsGeneLFC %<>% .[!is.na(.)]

vnsEnr <- gseGO(geneList     = vnsGeneLFC,
                OrgDb        = org.Mm.eg.db,
                ont          = "ALL",
                minGSSize    = 100,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE)

pdf(file.path(figDir, "VNS_GSE-GO_map.pdf"), width = 15, height = 15)
goplot(vnsEnr)
dev.off()

pdf(file.path(figDir, "VNS_GSE-GO_dot.pdf"), width = 8, height = 15)
dotplot(vnsEnr, showCategory = 30)
dev.off()

#Get gene names
library(pheatmap)
vnsGSEAgenes <- vnsEnr@result$core_enrichment %>% unlist() %>% str_split("/") %>% unlist() %>%
  AnnotationDbi::select(org.Mm.eg.db, keys = ., columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENTREZID")

#Filter changed genes
changedGenes <- res[[intGroups]] %>% as.data.frame() %>% .[.$log2FoldChange >= 1 & .$pvalue <= 0.05,] %>% rownames(.)

vnsGSEAgenes %<>% .[.$ENSEMBL %in% changedGenes,]

vnsTPM <- tpms %>% .[names(.) %in% vnsGSEAgenes$SYMBOL] %>% 
  lapply(., function(gene) {
    gene %<>% .[.$Group %in% c("Reference_2", "GROUP_A2"),]
    
    TPMs <- gene$TPMs
    TPMs <- TPMs/mean(TPMs, na.rm = TRUE)
    names(TPMs) <- gene$Group
    
    return(TPMs)
  }) %>% do.call("rbind", .)

codeDf <- data.frame(Group = rep(c("Reference_2", "GROUP_B2"), each = 4),
                     Code = paste(rep(c("SHAM + Z", "VNS + Z"), each = 4), c(1:4), sep = " ")
)

colnames(vnsTPM) <- codeDf$Code

pheatmap(vnsTPM, metadata = codeDf, cluster_cols = FALSE)


#Pathway analysis SHAM + Zymozan vs VX + Zymozan
intGroups <- c("GROUP_A2", "Reference_2")
intGroups <- paste(intGroups[1], "vs", intGroups[2], sep = " ")

#GO
#Select genes with lFC > 1
vxGeneIds <- res[[intGroups]] %>% as.data.frame() %>% 
  .[.$padj < 0.05 & .$log2FoldChange >= 1,] %>% .[!is.na(.$padj),] %>% rownames(.) %>%
  AnnotationDbi::select(org.Mm.eg.db, keys = ., columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENSEMBL")
vxGeneIds$LFC <- res[[intGroups]] %>% as.data.frame() %>% .[match(vxGeneIds$ENSEMBL, rownames(.)), "log2FoldChange"]

vxGO <- enrichGO(gene          = unique(vxGeneIds$ENTREZID),
                 OrgDb         = org.Mm.eg.db,
                 ont           = "ALL",
                 pAdjustMethod = "fdr",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)

pdf(file.path(figDir, "VX_GOenrichmentBar.pdf"), width = 8, height = 6)
barplot(vxGO)
dev.off()

#Gene Set Enrichment Analysis of Gene Ontology
vxGeneData <- as.data.frame(res[[intGroups]])
vxGeneData$ENSEMBL <- rownames(vxGeneData)

vxGeneData <- AnnotationDbi::select(org.Mm.eg.db, keys = rownames(vxGeneData), columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENSEMBL") %>% 
  dplyr::left_join(vxGeneData, ., by = "ENSEMBL") 

vxIntGeneLFC <- vxGeneData %>% .[, c("SYMBOL", "ENTREZID", "log2FoldChange")]

vxGeneLFC <- vxIntGeneLFC$log2FoldChange
names(vxGeneLFC ) <- vxIntGeneLFC$ENTREZID
vxGeneLFC  %<>% .[!duplicated(names(.))]
vxGeneLFC  %<>%  sort(., TRUE)
vxGeneLFC  %<>% .[!is.na(.)]

vxEnr <- gseGO(geneList     = vxGeneLFC ,
               OrgDb        = org.Mm.eg.db,
               ont          = "ALL",
               minGSSize    = 100,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

pdf(file.path(figDir, "VX_GSE-GO_map.pdf"), width = 15, height = 15)
goplot(vxEnr)
dev.off()

#Get gene names
library(pheatmap)
vxGSEAgenes <- vxEnr@result$core_enrichment %>% unlist() %>% str_split("/") %>% unlist() %>%
  AnnotationDbi::select(org.Mm.eg.db, keys = ., columns = c("SYMBOL", "ENTREZID", "ENSEMBL"), keytype = "ENTREZID")

#Filter changed genes
changedGenes <- res[[intGroups]] %>% as.data.frame() %>% .[.$log2FoldChange >= 1 & .$pvalue <= 0.05,] %>% rownames(.)

vxGSEAgenes %<>% .[.$ENSEMBL %in% changedGenes,]

vxTPM <- tpms %>% .[names(.) %in% vxGSEAgenes$SYMBOL] %>% 
  lapply(., function(gene) {
    gene %<>% .[.$Group %in% c("Reference_2", "GROUP_A2"),]
    
    TPMs <- gene$TPMs
    TPMs <- TPMs/mean(TPMs, na.rm = TRUE)
    names(TPMs) <- gene$Group
    
    return(TPMs)
    }) %>% do.call("rbind", .)

codeDf <- data.frame(Group = rep(c("Reference_2", "GROUP_A2"), each = 4),
                     Code = paste(rep(c("SHAM + Z", "VX + Z"), each = 4), c(1:4), sep = " ")
                     )

colnames(vxTPM) <- codeDf$Code

pheatmap(vxTPM, metadata = codeDf, cluster_cols = FALSE)



