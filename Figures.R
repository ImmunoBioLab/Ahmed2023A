library(ggplot2)
library(rstatix)
library(ggrepel)

#Figures for SHAM + Zymozan vs VNS + Zymozan
intGenes <- c("Pnpla3", "Csf2", "Cxcl1", "Il6", "Ccl2", "Mif", "Serpine11", "Ccl5")
intGroups <- c("GROUP_B2", "Reference_2")

#Volcano
bestGenes <- res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>% dplyr::arrange(., desc(abs(.$log2FoldChange)))
bestGenes$MGI <- ensemblData %>% .[which(.$ensembl_gene_id %in% rownames(bestGenes)), c("ensembl_gene_id", "mgi_symbol")] %>% unique() %>% 
  .[match(rownames(bestGenes), .$ensembl_gene_id), "mgi_symbol"]

bestGenes %<>% .[.$MGI %in% intGenes,]

tiff(file.path(figDir, "VNS_Volcano.tiff"), width = 5, height = 10, res = 600, units = "in")
res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>%
  ggplot(data = .) +
  geom_point(aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.2, fill = "gray80", color = "gray30") +
  geom_point(data = . %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),],
             aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.35, fill = "tomato", color = "red") +
  geom_text_repel(data = bestGenes, aes(x = log2FoldChange, y = log(padj, 2)*-1, label = MGI), max.overlaps = 30, size = 3) +
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
saveWorkbook(wb, file.path(dataDir, "VNS_DE.xlsx"))


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

tiff(file.path(figDir, "VNS_boxplot.tiff"), width = 6, height = 5, units = "in", res = 600)
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
saveWorkbook(wb, file.path(dataDir, "VNS_boxplot.xlsx"))



#Figures for SHAM + Zymozan vs VX + Zymozan
intGenes <- c("Pnpla3", "Csf2", "Cxcl1", "Il6", "Ccl2", "Mif", "Serpine11", "Ccl5")
intGroups <- c("GROUP_A2", "Reference_2")

#Volcano
bestGenes <- res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>% dplyr::arrange(., desc(abs(.$log2FoldChange)))
bestGenes$MGI <- ensemblData %>% .[which(.$ensembl_gene_id %in% rownames(bestGenes)), c("ensembl_gene_id", "mgi_symbol")] %>% unique() %>% 
  .[match(rownames(bestGenes), .$ensembl_gene_id), "mgi_symbol"]

bestGenes %<>% .[.$MGI %in% intGenes,]

tiff(file.path(figDir, "VX_Volcano.tiff"), width = 5, height = 10, res = 600, units = "in")
res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>%
  ggplot(data = .) +
  geom_point(aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.2, fill = "gray80", color = "gray30") +
  geom_point(data = . %>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05 & !is.na(.$padj),],
             aes(x = log2FoldChange, y = log(padj, 2)*-1), alpha = 0.35, fill = "tomato", color = "red") +
  geom_text_repel(data = bestGenes, aes(x = log2FoldChange, y = log(padj, 2)*-1, label = MGI), max.overlaps = 30, size = 3) +
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
saveWorkbook(wb, file.path(dataDir, "VX_DE.xlsx"))


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
saveWorkbook(wb, file.path(dataDir, "VX_boxplot.xlsx"))
