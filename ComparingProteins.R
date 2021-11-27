#load in data
library(readr)
library(stringr)
GWAS_MS <- read_csv("projectCode/data/GWAS/MS.csv")
GWAS_Hashimoto <- read_csv("projectCode/data/GWAS/Hashimotos.csv")
GWAS_Celiac <- read_csv("projectCode/data/GWAS/Celiac.csv")
GWAS_RArthritis <- read_csv("projectCode/data/GWAS/RArthritis.csv")
GWAS_SLupusE <- read_csv("projectCode/data/GWAS/SystemicLupusE.csv")
GWAS_T1Diabetes <- read_csv("projectCode/data/GWAS/Type1Diabetes.csv")
GWAS_Asthma <- read_csv("projectCode/data/GWAS/Asthma.csv")
GWAS_Graves <- read_csv("projectCode/data/GWAS/Graves.csv")
GWAS_IBD <- read_csv("projectCode/data/GWAS/IBD.csv")
GWAS_POS <- read_csv("projectCode/data/GWAS/POS.csv")

#separate columns Celiac
C_gene <- do.call(rbind, lapply(1: NROW(GWAS_Celiac), function(i)
  setNames(data.frame(GWAS_Celiac$`Variant and risk allele`[i],
                      GWAS_Celiac$`P-value`[i],
                      GWAS_Celiac$`P-value annotation`[i],
                      GWAS_Celiac$RAF[i],
                      GWAS_Celiac$OR[i],
                      GWAS_Celiac$Beta[i],
                      GWAS_Celiac$CI[i],
                      GWAS_Celiac$`Mapped gene`[i],
                      GWAS_Celiac$`Reported trait`[i],
                      GWAS_Celiac$`Trait(s)`[i],
                      GWAS_Celiac$`Background trait(s)`[i],
                      GWAS_Celiac$`Study accession`[i],
                      GWAS_Celiac$Location[i],
                      unlist(strsplit(as.character(GWAS_Celiac$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_Celiac$`Mapped gene`[i]), ", ")),
                                                  "\\d{4}")),
           c(names(GWAS_Celiac), "Gene"))))

#Separate Columns MS
MS_gene <- do.call(rbind, lapply(1: NROW(GWAS_MS), function(i)
  setNames(data.frame(GWAS_MS$`Variant and risk allele`[i],
                      GWAS_MS$`P-value`[i],
                      GWAS_MS$`P-value annotation`[i],
                      GWAS_MS$RAF[i],
                      GWAS_MS$OR[i],
                      GWAS_MS$Beta[i],
                      GWAS_MS$CI[i],
                      GWAS_MS$`Mapped gene`[i],
                      GWAS_MS$`Reported trait`[i],
                      GWAS_MS$`Trait(s)`[i],
                      GWAS_MS$`Background trait(s)`[i],
                      GWAS_MS$`Study accession`[i],
                      GWAS_MS$Location[i],
                      unlist(strsplit(as.character(GWAS_MS$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_MS$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_MS), "Gene"))))

#Separate Columns Hashimoto
H_gene <- do.call(rbind, lapply(1: NROW(GWAS_Hashimoto), function(i)
  setNames(data.frame(GWAS_Hashimoto$`Variant and risk allele`[i],
                      GWAS_Hashimoto$`P-value`[i],
                      GWAS_Hashimoto$`P-value annotation`[i],
                      GWAS_Hashimoto$RAF[i],
                      GWAS_Hashimoto$OR[i],
                      GWAS_Hashimoto$Beta[i],
                      GWAS_Hashimoto$CI[i],
                      GWAS_Hashimoto$`Mapped gene`[i],
                      GWAS_Hashimoto$`Reported trait`[i],
                      GWAS_Hashimoto$`Trait(s)`[i],
                      GWAS_Hashimoto$`Background trait(s)`[i],
                      GWAS_Hashimoto$`Study accession`[i],
                      GWAS_Hashimoto$Location[i],
                      unlist(strsplit(as.character(GWAS_Hashimoto$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_Hashimoto$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_Hashimoto), "Gene"))))

#Separate Columns Arthritis
A_gene <- do.call(rbind, lapply(1: NROW(GWAS_RArthritis), function(i)
  setNames(data.frame(GWAS_RArthritis$`Variant and risk allele`[i],
                      GWAS_RArthritis$`P-value`[i],
                      GWAS_RArthritis$`P-value annotation`[i],
                      GWAS_RArthritis$RAF[i],
                      GWAS_RArthritis$OR[i],
                      GWAS_RArthritis$Beta[i],
                      GWAS_RArthritis$CI[i],
                      GWAS_RArthritis$`Mapped gene`[i],
                      GWAS_RArthritis$`Reported trait`[i],
                      GWAS_RArthritis$`Trait(s)`[i],
                      GWAS_RArthritis$`Background trait(s)`[i],
                      GWAS_RArthritis$`Study accession`[i],
                      GWAS_RArthritis$Location[i],
                      unlist(strsplit(as.character(GWAS_RArthritis$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_RArthritis$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_RArthritis), "Gene"))))


#Separate Columns Lupus
L_gene <- do.call(rbind, lapply(1: NROW(GWAS_SLupusE), function(i)
  setNames(data.frame(GWAS_SLupusE$`Variant and risk allele`[i],
                      GWAS_SLupusE$`P-value`[i],
                      GWAS_SLupusE$`P-value annotation`[i],
                      GWAS_SLupusE$RAF[i],
                      GWAS_SLupusE$OR[i],
                      GWAS_SLupusE$Beta[i],
                      GWAS_SLupusE$CI[i],
                      GWAS_SLupusE$`Mapped gene`[i],
                      GWAS_SLupusE$`Reported trait`[i],
                      GWAS_SLupusE$`Trait(s)`[i],
                      GWAS_SLupusE$`Background trait(s)`[i],
                      GWAS_SLupusE$`Study accession`[i],
                      GWAS_SLupusE$Location[i],
                      unlist(strsplit(as.character(GWAS_SLupusE$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_SLupusE$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_SLupusE), "Gene"))))


#Separate Columns Type 1 Diabetes
D_gene <- do.call(rbind, lapply(1: NROW(GWAS_T1Diabetes), function(i)
  setNames(data.frame(GWAS_T1Diabetes$`Variant and risk allele`[i],
                      GWAS_T1Diabetes$`P-value`[i],
                      GWAS_T1Diabetes$`P-value annotation`[i],
                      GWAS_T1Diabetes$RAF[i],
                      GWAS_T1Diabetes$OR[i],
                      GWAS_T1Diabetes$Beta[i],
                      GWAS_T1Diabetes$CI[i],
                      GWAS_T1Diabetes$`Mapped gene`[i],
                      GWAS_T1Diabetes$`Reported trait`[i],
                      GWAS_T1Diabetes$`Trait(s)`[i],
                      GWAS_T1Diabetes$`Background trait(s)`[i],
                      GWAS_T1Diabetes$`Study accession`[i],
                      GWAS_T1Diabetes$Location[i],
                      unlist(strsplit(as.character(GWAS_T1Diabetes$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_T1Diabetes$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_T1Diabetes), "Gene"))))


#separate columns IBD
IBD_gene <- do.call(rbind, lapply(1: NROW(GWAS_IBD), function(i)
  setNames(data.frame(GWAS_IBD$`Variant and risk allele`[i],
                      GWAS_IBD$`P-value`[i],
                      GWAS_IBD$`P-value annotation`[i],
                      GWAS_IBD$RAF[i],
                      GWAS_IBD$OR[i],
                      GWAS_IBD$Beta[i],
                      GWAS_IBD$CI[i],
                      GWAS_IBD$`Mapped gene`[i],
                      GWAS_IBD$`Reported trait`[i],
                      GWAS_IBD$`Trait(s)`[i],
                      GWAS_IBD$`Background trait(s)`[i],
                      GWAS_IBD$`Study accession`[i],
                      GWAS_IBD$Location[i],
                      unlist(strsplit(as.character(GWAS_IBD$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_IBD$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_IBD), "Gene"))))

#separate columns Graves
G_gene <- do.call(rbind, lapply(1: NROW(GWAS_Graves), function(i)
  setNames(data.frame(GWAS_Graves$`Variant and risk allele`[i],
                      GWAS_Graves$`P-value`[i],
                      GWAS_Graves$`P-value annotation`[i],
                      GWAS_Graves$RAF[i],
                      GWAS_Graves$OR[i],
                      GWAS_Graves$Beta[i],
                      GWAS_Graves$CI[i],
                      GWAS_Graves$`Mapped gene`[i],
                      GWAS_Graves$`Reported trait`[i],
                      GWAS_Graves$`Trait(s)`[i],
                      GWAS_Graves$`Background trait(s)`[i],
                      GWAS_Graves$`Study accession`[i],
                      GWAS_Graves$Location[i],
                      unlist(strsplit(as.character(GWAS_Graves$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_Graves$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_Graves), "Gene"))))

#separate columns Asthma
AS_gene <- do.call(rbind, lapply(1: NROW(GWAS_Asthma), function(i)
  setNames(data.frame(GWAS_Asthma$`Variant and risk allele`[i],
                      GWAS_Hashimoto$`P-value`[i],
                      GWAS_Hashimoto$`P-value annotation`[i],
                      GWAS_Hashimoto$RAF[i],
                      GWAS_Hashimoto$OR[i],
                      GWAS_Hashimoto$Beta[i],
                      GWAS_Hashimoto$CI[i],
                      GWAS_Hashimoto$`Mapped gene`[i],
                      GWAS_Hashimoto$`Reported trait`[i],
                      GWAS_Hashimoto$`Trait(s)`[i],
                      GWAS_Hashimoto$`Background trait(s)`[i],
                      GWAS_Hashimoto$`Study accession`[i],
                      GWAS_Hashimoto$Location[i],
                      unlist(strsplit(as.character(GWAS_Hashimoto$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_Hashimoto$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_Hashimoto), "Gene"))))

#separate columns POS
POS_gene <- do.call(rbind, lapply(1: NROW(GWAS_POS), function(i)
  setNames(data.frame(GWAS_POS$`Variant and risk allele`[i],
                      GWAS_POS$`P-value`[i],
                      GWAS_POS$`P-value annotation`[i],
                      GWAS_POS$RAF[i],
                      GWAS_POS$OR[i],
                      GWAS_POS$Beta[i],
                      GWAS_POS$CI[i],
                      GWAS_POS$`Mapped gene`[i],
                      GWAS_POS$`Reported trait`[i],
                      GWAS_POS$`Trait(s)`[i],
                      GWAS_POS$`Background trait(s)`[i],
                      GWAS_POS$`Study accession`[i],
                      GWAS_POS$Location[i],
                      unlist(strsplit(as.character(GWAS_POS$`Mapped gene`[i]), ", ")),
                      str_extract(unlist(strsplit(as.character(GWAS_POS$`Mapped gene`[i]), ", ")),
                                  "\\d{4}")),
           c(names(GWAS_POS), "Gene"))))


CeliacGene <- na.omit(unique(C_gene$Gene))
MSGene <- na.omit(unique(MS_gene$Gene))
HashimotoGene <- na.omit(unique(H_gene$Gene))
RArthritisGene <- na.omit(unique(A_gene$Gene))
LupusGene <- na.omit(unique(L_gene$Gene))
T1DiabetesGene <- na.omit(unique(D_gene$Gene))
IBDGene <- na.omit(unique(IBD_gene$Gene))
GravesGene <- na.omit(unique(G_gene$Gene))
AsthmaGene <- na.omit(unique(AS_gene$Gene))
POSGene <- na.omit(unique(POS_gene$Gene))

library(VennDiagram)
library(RColorBrewer)

# Chart
myCol <- brewer.pal(4, "Spectral")

venn.diagram(
  x = list(MSGene, RArthritisGene, LupusGene, T1DiabetesGene),
  category.names = c("MS" , "RArthritis" , "Lupus", "T1Diabetes"),
  filename = 'Genes_venn_diagramm.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

myCol2 <- brewer.pal(5, "Spectral")

# Compare MS, Arthritis, Lupus, T1Diabetes, Asthma

venn.diagram(
  x = list(MSGene, RArthritisGene, LupusGene, T1DiabetesGene, AsthmaGene),
  category.names = c("MS" , "RArthritis" , "Lupus", "T1Diabetes", "Asthma"),
  filename = 'Genes_venn_diagramm4.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol2,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# Asthma has good amount in common with Arthritis and T1Diabetes, but not together.

# Compare IBD, Graves, Celiac, Asthma

venn.diagram(
  x = list(IBDGene, GravesGene, CeliacGene, AsthmaGene),
  category.names = c("IBD", "Graves", "Celiac", "Asthma"),
  filename = 'Genes_venn_diagramm2.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# IBD and asthma most in common, but some with IBD, Celiac and asthma, and just celiac and asthma

# Compare IBD, Graves, Hashimoto's, Celiac, Asthma

venn.diagram(
  x = list(IBDGene, GravesGene, HashimotoGene, CeliacGene, AsthmaGene),
  category.names = c("IBD", "Graves", "Hashimoto's", "Celiac", "Asthma"),
  filename = 'Genes_venn_diagramm3.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol2,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# Nothing really of interest on this one

# Comparing MS, Arthritis, Lupus, and Asthma

venn.diagram(
  x = list(MSGene, RArthritisGene, LupusGene, AsthmaGene),
  category.names = c("MS", "RArthritis", "Lupus", "Asthma"),
  filename = 'Genes_venn_diagramm5.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# Lupus and Arthritis most similar, Asthma/ arthritis, arthritis/ MS, MS/ Lupus

# Compare IBD, T1Diabetes, Celiac, Asthma

venn.diagram(
  x = list(IBDGene, T1DiabetesGene, CeliacGene, AsthmaGene),
  category.names = c("IBD", "T1Diabetes", "Celiac", "Asthma"),
  filename = 'Genes_venn_diagramm6.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

# 4 common Genes, many when excluding asthma, asthma/IBD, IBD/T1Diabetes, T1Diabetes/asthma

# Compare IBD, Graves, Hashimoto's, Celiac

venn.diagram(
  x = list(IBDGene, GravesGene, HashimotoGene, CeliacGene, AsthmaGene),
  category.names = c("IBD", "Graves", "Hashimoto's", "Celiac", "Asthma"),
  filename = 'Genes_venn_diagramm8.png',
  output=TRUE,
  lwd = 1.5,
  lty = 'blank',
  fill = myCol2,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  cat.cex = .4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.fontfamily = "sans"
)

filter <- AsthmaGene[AsthmaGene %in% IBDGene]
filter2 <- AsthmaGene[AsthmaGene %in% T1DiabetesGene]
filter3 <- AsthmaGene[AsthmaGene %in% RArthritisGene]
filter4 <- AsthmaGene[AsthmaGene %in% LupusGene]
filter5 <- AsthmaGene[AsthmaGene %in% POSGene]
filter6 <- AsthmaGene[AsthmaGene %in% CeliacGene]
filter7 <- AsthmaGene[AsthmaGene %in% MSGene]
filter8 <- AsthmaGene[AsthmaGene %in% GravesGene]
filter9 <- AsthmaGene[AsthmaGene %in% HashimotoGene]

tbl_a <- AsthmaGene[AsthmaGene %in% filter]
tbl_a2 <- AsthmaGene[AsthmaGene %in% filter2]
tbl_a3 <- AsthmaGene[AsthmaGene %in% filter3]
tbl_a4 <- AsthmaGene[AsthmaGene %in% filter4]
tbl_a5 <- AsthmaGene[AsthmaGene %in% filter5]
tbl_a6 <- AsthmaGene[AsthmaGene %in% filter6]
tbl_a7 <- AsthmaGene[AsthmaGene %in% filter7]
tbl_a8 <- AsthmaGene[AsthmaGene %in% filter8]
tbl_a9 <- AsthmaGene[AsthmaGene %in% filter9]

tbl_d <- T1DiabetesGene[T1DiabetesGene %in% filter2]
tbl_ar <- RArthritisGene[RArthritisGene %in% filter3]
tbl_l <- LupusGene[LupusGene %in% filter4]
tbl_POS <- POSGene[POSGene %in% filter5]
tbl_IBD <- IBDGene[IBDGene %in% filter]
tbl_c <- CeliacGene[CeliacGene %in% filter6]
tbl_MS <- MSGene[MSGene %in% filter7]
tbl_g <- GravesGene[GravesGene %in% filter8]
tbl_h <- HashimotoGene[HashimotoGene %in% filter9]

data.frame(tbl_a3, tbl_ar)
data.frame(tbl_a6, tbl_c)
data.frame(tbl_a2, tbl_d)
data.frame(tbl_a8, tbl_g)
data.frame(tbl_a9, tbl_h)
data.frame(tbl_a, tbl_IBD)
data.frame(tbl_a4, tbl_l)
data.frame(tbl_a7, tbl_MS)
data.frame(tbl_a5, tbl_POS)

data.frame(tbl_a, tbl_d, tbl_ar, tbl_l, tbl_POS, tbl_IBD, tbl_c, tbl_MS, tbl_g, tbl_h)

a_list <- AsthmaGene
d_list <- c()
ar_list <- c()
l_list <- c()
POS_list <- c()
IBD_list <- c()
c_list <- c()
MS_list <- c()
g_list <- c()
h_list <- c()

for (i in AsthmaGene) {
  c_genes <- paste(unique(CeliacGene[CeliacGene == i]), collapse = "")
  IBD_genes <- paste(unique(IBDGene[IBDGene == i], 1), collapse = "")
  d_genes <- paste(unique(T1DiabetesGene[T1DiabetesGene == i], 1), collapse = "")
  ar_genes <- paste(unique(RArthritisGene[RArthritisGene == i], 1), collapse = "")
  l_genes <- paste(unique(LupusGene[LupusGene == i], 1), collapse = "")
  MS_genes <- paste(unique(MSGene[MSGene == i], 1), collapse = "")
  g_genes <- paste(unique(GravesGene[GravesGene == i], 1), collapse = "")
  h_genes <- paste(unique(HashimotoGene[HashimotoGene == i], 1), collapse = "")
  
  d_list <- append(d_list, d_genes)
  ar_list <- append(ar_list, ar_genes)
  l_list <- append(l_list, l_genes)
  IBD_list <- append(IBD_list, IBD_genes)
  c_list <- append(c_list, c_genes)
  MS_list <- append(MS_list, MS_genes)
  g_list <- append(g_list, g_genes)
  h_list <- append(h_list, h_genes)
  
}

data.frame(a_list, c_list, IBD_list, d_list, ar_list, l_list, MS_list, g_list, h_list)
