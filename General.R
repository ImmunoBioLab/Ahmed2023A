library(stringr)
library(magrittr)

setwd("Z:/Groups/Peder Olofsson/VSS/Bioconductor/OA_Liver_2022-05-25/")

dir.create(file.path(getwd(), "Figures"))
figDir <- file.path(getwd(), "Figures")
dir.create(file.path(getwd(), "Data"))
dataDir <- file.path(getwd(), "Data")


#Get sample names
samples <- file.path(getwd(), "quant_data") %>% dir() %>% str_remove(., "_quants") %>% str_trim() %>% factor()

#Get sample info
sampleInfo <- file.path(getwd(), "sampleInfo.csv") %>% readr::read_csv() %>% as.data.frame()
colnames(sampleInfo)[colnames(sampleInfo) == "Sample Name"] <- "names"
sampleInfo$names %<>% factor(., levels = samples)
sampleInfo %<>% dplyr::arrange(., names)
sampleInfo$files <- file.path(getwd(), "quant_data", str_c(samples, "quants", sep = "_"), "quant.sf")


#Get groups
library(openxlsx)
groupDf <- read.xlsx("20220531 Blinding groups for Vladimir.xlsx")

sampleInfo %<>% dplyr::left_join(., groupDf[, c("names", "Groups", "Treatment", "Surgery")], by = "names") %>% as.data.frame()
rownames(sampleInfo) <- sampleInfo$names 


#Get MGI symbols
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

geneNames <- lapply(sampleInfo$files, function(quant) {
  dataDf <- readr::read_tsv(quant)
  
  dataDf$Name
}) %>% unlist() %>% unique()

ensemblData <- getBM(attributes = c('ensembl_transcript_id_version', 'ensembl_gene_id', "mgi_symbol"),
                     filters = 'ensembl_transcript_id_version',
                     values = geneNames, 
                     mart = ensembl)

all(geneNames %in% unique(ensemblData$ensembl_transcript_id_version))

unique(ensemblData$ensembl_gene_id) %>% length()
unique(ensemblData$mgi_symbol) %>% length()


#Get TPMs
tpms <- lapply(sampleInfo$files, function(quant) {
  dataDf <- readr::read_tsv(quant)
})
names(tpms) <- sampleInfo$names

tpms <-lapply(seq_along(tpms), function(index) {
  tpms[[index]] %<>% dplyr::mutate(Sample = names(tpms)[index])
}) %>% do.call("rbind", .)
colnames(tpms)[1] <- "ensembl_transcript_id_version"

tpms %<>% dplyr::left_join(., ensemblData, by = "ensembl_transcript_id_version")

tpms <- lapply(unique(tpms$mgi_symbol), function(gene) {
  gene <- tpms %>% .[.$mgi_symbol == gene,]
  
  if(sum(gene$TPM) > 0) {
    tpmDf <- dplyr::group_by(gene, Sample) %>%
      dplyr::summarise(TPMs = sum(TPM))
    
    tpmDf$mgi_symbol <- unique(gene$mgi_symbol)
    
    return(tpmDf)
  } else {
    cat("No reads detected! Removing", unique(gene$mgi_symbol), "\n")
  }
})
names(tpms) <- sapply(tpms, function(gene) unique(gene$mgi_symbol))

tpms %<>% .[sapply(., function(gene) !is.null(gene))]

#Add groups
tpms %<>% lapply(., function(gene) {
  gene$Group <- sampleInfo %>% .[match(.$names, gene$Sample), "Groups"]
  
  return(gene)
})

sumTpms <- lapply(tpms, function(gene) {
  sumDf <- dplyr::group_by(gene, Group) %>%
    dplyr::summarize(uMean = mean(TPMs),
                     SD = sd(TPMs),
                     Var = var(TPMs),
                     Count = sum(complete.cases(TPMs)))
  
  sumDf %<>% dplyr::mutate(uSEM = .$SD/sqrt(.$Count),
                           Mean = .$uMean/.$uMean[[nrow(sumDf)]])
  sumDf %<>% dplyr::mutate(SEM = .$uSEM/.$uMean[[nrow(sumDf)]])
  
  return(sumDf)
})
names(sumTpms) <- names(tpms)


#Make DESeq object
library(tximport)
filePath <- sampleInfo$files
names(filePath) <- rownames(sampleInfo)

txi <- tximport(filePath, type = "salmon", tx2gene = ensemblData)

library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = sampleInfo, design = ~Groups)

#Filtering genes which have less then 10 counts in all samples combined
sum(rowSums(counts(dds)) < 10)
dds <- dds[rowSums(counts(dds)) >= 10,]

dds %<>% DESeq()


#Establish comparisons
contrasts <- list(c("GROUP_A1", "Reference_1"),
                  c("GROUP_A2", "Reference_2"),
                  c("GROUP_B2", "Reference_2"),
                  c("GROUP_A2", "GROUP_A1"),
                  c("GROUP_B2", "GROUP_A1"),
                  c("GROUP_B2", "GROUP_A2")
)

#Calculate DE per comparison
res <- list()
for(i in seq_along(contrasts)) {
  res[[i]] <- results(dds, contrast = c("Groups", contrasts[[i]][1], contrasts[[i]][2]))
}
names(res) <- sapply(contrasts, function(group) paste(group[1], "vs", group[2]))