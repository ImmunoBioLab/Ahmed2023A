library(ggplot2)
library(rstatix)
library(ggrepel)
library(openxlsx)

#Figures for SHAM + Zymozan vs VNS + Zymozan
intGenes <- c("Pnpla3", "Csf2", "Cxcl1", "Il6", "Ccl2", "Mif", "Serpine1", "Ccl5", "Des", "Acta2", "Syp", "Igfbp3", "Gfap", "Mcam")
intGroups <- c("GROUP_B2", "Reference_2")

#Volcano
bestGenes <- res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>% dplyr::arrange(., desc(abs(.$log2FoldChange)))
bestGenes$MGI <- ensemblData %>% .[match(rownames(bestGenes), .$ensembl_gene_id), "mgi_symbol"]
bestGenes %<>% .[.$MGI %in% intGenes & abs(.$log2FoldChange) > 1 & .$padj < 0.05,]

pdf(file.path(figDir, "VNS_Volcano.pdf"), width = 3.2, height = 3.2)
res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>%
  ggplot(data = .) +
  geom_point(aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.2, fill = "gray80", color = "gray30") +
  geom_point(data = . %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),],
             aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.7, color = "tomato2", fill = "firebrick3", size = 2, shape = 21) +
  geom_text_repel(data = bestGenes, aes(x = log2FoldChange, y = log(padj, 2)*-1, label = MGI), 
                  max.overlaps = 30, size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  geom_vline(xintercept = c(-1, 1), size = 1, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = log(0.05, 2)*-1, size = 1, alpha = 0.5, linetype = "dashed") +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = expression(paste(log[2]," FC", sep = "")), y = expression(paste(-log[2], " adjusted p value", sep = ""))) +
  theme_classic() +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14)
  )
dev.off()

#Export  
wb <- createWorkbook()
addWorksheet(wb, "DE")
writeData(wb, 1,  as.data.frame(res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]]), rowNames = TRUE)
saveWorkbook(wb, file.path(dataDir, "VNS_DE.xlsx"), overwrite = TRUE)


#Boxplots
boxData <- tpms %>% .[names(.) %in% intGenes] %>% do.call("rbind", .)
codeDf <- data.frame(Group = c(intGroups),
                     Code = c("VNS + Z", "SHAM + Z"))
boxData %<>% dplyr::left_join(., codeDf, by = "Group") %>%
  .[.$Group %in% intGroups,]
boxData$Code %<>% factor(., levels = c("SHAM + Z", "VNS + Z"))

statTPMs <- boxData %>%
  dplyr::group_by_at(c("mgi_symbol")) %>%
  t_test(TPMs ~ Code) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
statTPMs %<>% add_xy_position(x = "Code", dodge = 0.8)
statTPMs <- boxData %>% dplyr::group_by(mgi_symbol) %>% 
  dplyr::summarize(newY = max(TPMs)*1.1) %>% 
  dplyr::left_join(statTPMs, ., "mgi_symbol") 

tiff(file.path(figDir, "VNS_boxplot.tiff"), width = 8, height = 5, units = "in", res = 600)
ggplot(data = boxData, aes(x = Code, y = TPMs)) +
  geom_boxplot(fill = rep(c("white", "black"), times = length(unique(boxData$mgi_symbol)))) +
  facet_wrap(~mgi_symbol, scale = "free_y", nrow = 2) +
  ggpubr::stat_pvalue_manual(statTPMs,  label = "p", tip.length = 0.03, size = 3, bracket.size = 0.6, y.position = "newY") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12), axis.title.x = element_blank())
dev.off()

#Export  
wb <- createWorkbook()
addWorksheet(wb, "TPMs")
writeData(wb, 1, boxData)
addWorksheet(wb, "Stats")
writeData(wb, 2, statTPMs)
saveWorkbook(wb, file.path(dataDir, "VNS_boxplot.xlsx"), overwrite = TRUE)



#Figures for SHAM + Zymozan vs VX + Zymozan
intGroups <- c("GROUP_A2", "Reference_2")

#Volcano
bestGenes <- res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>% dplyr::arrange(., desc(abs(.$log2FoldChange)))
bestGenes$MGI <- ensemblData %>% .[match(rownames(bestGenes), .$ensembl_gene_id), "mgi_symbol"]

bestGenes %<>% .[.$MGI %in% intGenes & abs(.$log2FoldChange) > 1 & .$padj < 0.05,]

pdf(file.path(figDir, "VX_Volcano.pdf"), width = 3.2, height = 3.2)
res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>%
  ggplot(data = .) +
  geom_point(aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.2, fill = "gray80", color = "gray30") +
  geom_point(data = . %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),],
             aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.7, color = "tomato2", fill = "firebrick3", size = 2, shape = 21) +
  geom_text_repel(data = bestGenes, aes(x = log2FoldChange, y = log(padj, 2)*-1, label = MGI), 
                  max.overlaps = 30, size = 3, box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  geom_vline(xintercept = c(-1, 1), size = 1, alpha = 0.5, linetype = "dashed") +
  geom_hline(yintercept = log(0.05, 2)*-1, size = 1, alpha = 0.5, linetype = "dashed") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(limits = c(-4, 4)) +
  labs(x = expression(paste(log[2]," FC", sep = "")), y = expression(paste(-log[2], " adjusted p value", sep = ""))) +
  theme_classic() +
  theme(axis.line = element_line(size = 1), axis.ticks = element_line(size = 1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 14)
  )
dev.off()

#Export  
wb <- createWorkbook()
addWorksheet(wb, "DE")
writeData(wb, 1,  as.data.frame(res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]]), rowNames = TRUE)
saveWorkbook(wb, file.path(dataDir, "VX_DE.xlsx"), overwrite = TRUE)


#Boxplots
boxData <- tpms %>% .[names(.) %in% intGenes] %>% do.call("rbind", .)
codeDf <- data.frame(Group = c(intGroups),
                     Code = c("VX + Z", "SHAM + Z"))
boxData %<>% dplyr::left_join(., codeDf, by = "Group") %>%
  .[.$Group %in% intGroups,]
boxData$Code %<>% factor(., levels = c("SHAM + Z", "VX + Z"))

statTPMs <- boxData %>%
  dplyr::group_by_at(c("mgi_symbol")) %>%
  t_test(TPMs ~ Code) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
statTPMs %<>% add_xy_position(x = "Code", dodge = 0.8)
statTPMs <- boxData %>% dplyr::group_by(mgi_symbol) %>% 
  dplyr::summarize(newY = max(TPMs)*1.1) %>% 
  dplyr::left_join(statTPMs, ., "mgi_symbol") 

tiff(file.path(figDir, "VX_boxplot.tiff"), width = 6, height = 5, units = "in", res = 600)
ggplot(data = boxData, aes(x = Code, y = TPMs)) +
  geom_boxplot(fill = rep(c("white", "black"), times = length(unique(boxData$mgi_symbol)))) +
  facet_wrap(~mgi_symbol, scale = "free_y", nrow = 2) +
  ggpubr::stat_pvalue_manual(statTPMs,  label = "p", tip.length = 0.03, size = 3, bracket.size = 0.6, y.position = "newY") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 12), axis.title.x = element_blank())
dev.off()

#Export  
wb <- createWorkbook()
addWorksheet(wb, "TPMs")
writeData(wb, 1, boxData)
addWorksheet(wb, "Stats")
writeData(wb, 2, statTPMs)
saveWorkbook(wb, file.path(dataDir, "VX_boxplot.xlsx"), overwrite = TRUE)


#Cross-FC plot
vnsGenes <- res[["GROUP_B2 vs Reference_2"]] %>% as.data.frame() %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),] %>% rownames()
vxGenes <- res[["GROUP_A2 vs Reference_2"]] %>% as.data.frame() %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),] %>% rownames()
chGenes <- ensemblData %>% .[match(intGenes, .$mgi_symbol), "ensembl_gene_id"] %>% c(., vnsGenes, vxGenes) %>% unique()

geneTable <- data.frame(MGI = ensemblData %>% .[match(chGenes, .$ensembl_gene_id), "mgi_symbol"],
                        ENSEMBL = chGenes
                        )

vnsFC <- res[["GROUP_B2 vs Reference_2"]] %>% as.data.frame() %>% .[rownames(.) %in% geneTable$ENSEMBL, c("log2FoldChange", "padj")]
vnsFC %<>% dplyr::mutate(., ENSEMBL = rownames(.))
colnames(vnsFC) <- c("VNS_FC", "VNS_padj", "ENSEMBL")

vxFC <- res[["GROUP_A2 vs Reference_2"]] %>% as.data.frame() %>% .[rownames(.) %in% geneTable$ENSEMBL, c("log2FoldChange", "padj")]
vxFC %<>% dplyr::mutate(., ENSEMBL = rownames(.))
colnames(vxFC) <- c("VX_FC", "VX_padj", "ENSEMBL")

geneTable %<>% dplyr::left_join(., vxFC, "ENSEMBL") 
geneTable %<>% dplyr::left_join(., vnsFC, "ENSEMBL")

geneTable %<>% .[!is.na(.$VX_FC),]

shapiro.test(geneTable$VX_FC)
shapiro.test(geneTable$VNS_FC)

corData <- cor.test(geneTable$VX_FC, geneTable$VNS_FC, method = "spearman")
corData <- data.frame(p = formatC(corData$p.value, format = "e", digits = 2),
                      rho = round(corData$estimate, 2))

xLimits <- c(min(geneTable$VX_FC)-0.5, max(geneTable$VX_FC)+0.5)

pdf(file.path(figDir, "VNS vs VX FC.pdf"), width = 4, height = 4)
ggpubr::ggarrange(plotlist  = list(
  ggplot() +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_point(geneTable[geneTable$VNS_FC > 10,], mapping = aes(x = VX_FC, y = VNS_FC), color = "grey80", fill = "gray90", alpha = 0.9, size = 3, shape = 21) +
        scale_y_continuous(limits = c(22.5, 25), breaks = c(22.5, 25)) +
    scale_x_continuous(limits = xLimits) +
    theme_classic() +
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.title = element_blank(), axis.ticks.x = element_blank()),
  ggplot() +
    geom_hline(yintercept = 0, color = "gray50", linetype = "dashed") +
    geom_vline(xintercept = 0, color = "gray50", linetype = "dashed") +
    geom_point(geneTable[geneTable$VNS_FC < 10,], mapping = aes(x = VX_FC, y = VNS_FC), color = "grey80", fill = "gray90", alpha = 0.9, size = 3, shape = 21) +
    geom_point(geneTable[geneTable$MGI %in% intGenes,], mapping = aes(x = VX_FC, y = VNS_FC), color = "tomato", alpha = 0.7, size = 3) +
    geom_text_repel(data = geneTable[geneTable$MGI %in% intGenes,], aes(x = VX_FC, y = VNS_FC, label = MGI), max.overlaps = 30, size = 3) +
    geom_text(data = corData, mapping = aes(x = Inf, y = Inf, hjust = 1, vjust = 1,
                                            label = paste0("rho = ", rho)), inherit.aes = FALSE, parse = FALSE) +
    geom_text(data = corData, mapping = aes(x = Inf, y = Inf, hjust = 1, vjust = 2.85,
                                            label = paste0("p = ", p)), inherit.aes = FALSE) +
    scale_y_continuous(limits = c(-5, 10), breaks = c(-5, -2,5, 0, 2,5, 5, 7.5, 10)) +
    scale_x_continuous(limits = xLimits) +
    labs(x = expression(paste("VX vs Sham ", log[2]," FC", sep = "")), y = expression(paste("VNS vs Sham ", log[2]," FC", sep = ""))) +
    theme_classic() +
    theme(axis.title.y = element_text(hjust = 0.8))
  ),
  ncol = 1, nrow = 2, align = "v", heights = c(0.15, 1)
)
dev.off()