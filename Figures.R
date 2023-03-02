library(ggplot2)
library(rstatix)
library(ggrepel)
library(openxlsx)

#Figures for SHAM + Zymozan vs VX + Zymozan
intGenes <- c("Pnpla3", "Csf2", "Cxcl1", "Il6", "Ccl2", "Mif", "Serpine1", "Ccl5", "Des", "Acta2", "Syp", "Igfbp3", "Gfap", "Mcam")
intGroups <- c("GROUP_A2", "Reference_2")

#Fig. 4A: Volcano
#Select genes with log2 FC >1 and adjusted p value < 0.05
sigGenes <- res[[paste(intGroups[1], "vs", intGroups[2], sep = " ")]] %>% as.data.frame() %>% dplyr::arrange(., desc(abs(.$log2FoldChange)))
sigGenes %<>% .[abs(.$log2FoldChange) > 1 & .$padj < 0.05,] %>% .[!is.na(.$baseMean),]
sigGenes$MGI <- ensemblData %>% .[match(rownames(sigGenes), .$ensembl_gene_id), "mgi_symbol"]

bestGenes <- sigGenes %>% dplyr::arrange(., desc(abs(.$log2FoldChange))) %>% .[c(1:10),]

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


#Figure 4B: Boxplots
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

pdf(file.path(figDir, "VX_boxplot.pdf"), width = 6, height = 5)
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


#Figure S4: Boxplots
#Remove Pnpla3 (plotted in Fig. 4) and order
upGenes <- sigGenes %>% .[.$log2FoldChange > 0, "MGI"]
downGenes <- sigGenes %>% dplyr::arrange(., abs(.$log2FoldChange)) %>% .[.$log2FoldChange < 0, "MGI"]

sigNames <- c(upGenes, downGenes) %>% .[!. %in% "Pnpla3"] 

#Get data
sigBox <- tpms %>% .[names(.) %in% sigNames] %>% do.call("rbind", .)
sigBox %<>% dplyr::left_join(., codeDf, by = "Group") %>%
  .[.$Group %in% intGroups,]
sigBox$Code %<>% factor(., levels = c("SHAM + Z", "VX + Z"))
sigBox$mgi_symbol %<>% factor(., levels = sigNames)

#Calculate stats
statSig <- sigBox %>%
  dplyr::group_by_at(c("mgi_symbol")) %>%
  t_test(TPMs ~ Code) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj")
statSig %<>% add_xy_position(x = "Code", dodge = 0.8)
statSig <- sigBox %>% dplyr::group_by(mgi_symbol) %>% 
  dplyr::summarize(newY = max(TPMs)*1.1) %>% 
  dplyr::left_join(statSig, ., "mgi_symbol") 
statSig$mgi_symbol %<>% factor(., levels = sigNames)

#Plot
pdf(file.path(figDir, "VX_boxplot_sigGenes.pdf"), width = 10, height = 5)
ggplot(data = sigBox, aes(x = Code, y = TPMs)) +
  geom_boxplot(fill = rep(c("white", "black"), times = length(unique(sigBox$mgi_symbol)))) +
  facet_wrap(~mgi_symbol, scale = "free_y", nrow = 2) +
  ggpubr::stat_pvalue_manual(statSig,  label = "p", tip.length = 0.03, size = 3, bracket.size = 0.6, y.position = "newY") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12), axis.title.x = element_blank())
dev.off()

#Export  
wb <- createWorkbook()
addWorksheet(wb, "TPMs")
writeData(wb, 1, sigBox)
addWorksheet(wb, "Stats")
writeData(wb, 2, statSig)
saveWorkbook(wb, file.path(dataDir, "VX_boxplot_sigGenes.xlsx"), overwrite = TRUE)