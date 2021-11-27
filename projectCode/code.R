#load in data
library(readr)
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
GWAS_Celiac$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_Celiac$`Variant and risk allele`)
GWAS_Celiac$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_Celiac$`Variant and risk allele`)

#Separate Columns MS
GWAS_MS$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_MS$`Variant and risk allele`)
GWAS_MS$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_MS$`Variant and risk allele`)

#Separate Columns Hashimoto
GWAS_Hashimoto$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_Hashimoto$`Variant and risk allele`)
GWAS_Hashimoto$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_Hashimoto$`Variant and risk allele`)

#Separate Columns Arthritis
GWAS_RArthritis$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_RArthritis$`Variant and risk allele`)
GWAS_RArthritis$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_RArthritis$`Variant and risk allele`)

#Separate Columns Lupus
GWAS_SLupusE$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_SLupusE$`Variant and risk allele`)
GWAS_SLupusE$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_SLupusE$`Variant and risk allele`)

#Separate Columns Type 1 Diabetes
GWAS_T1Diabetes$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_T1Diabetes$`Variant and risk allele`)
GWAS_T1Diabetes$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_T1Diabetes$`Variant and risk allele`)

#separate columns IBD
GWAS_IBD$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_IBD$`Variant and risk allele`)
GWAS_IBD$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_IBD$`Variant and risk allele`)

#separate columns Graves
GWAS_Graves$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_Graves$`Variant and risk allele`)
GWAS_Graves$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_Graves$`Variant and risk allele`)

#separate columns Asthma
GWAS_Asthma$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_Asthma$`Variant and risk allele`)
GWAS_Asthma$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_Asthma$`Variant and risk allele`)

#separate columns POS
GWAS_POS$ID <- sub("*\\-<b>[A,T,G,C,?]</b>", "", GWAS_POS$`Variant and risk allele`)
GWAS_POS$allele <- gsub(".*-<b>(.+)</b>.*", "\\1", GWAS_POS$`Variant and risk allele`)



#See if any SNPs Match up

CeliacSNP <- unique(GWAS_Celiac$ID)
# Pathway Name Entities Found Entities ratio Entities p-value Entities FRD Reactions Found Reactions Total Reactions Ratio Species Name
# RUNX1 and FOXP3 control the development of regulatory T lymphocytes (Tregs) 6 17 0.001 2.5E-7 8E-5 9 20 0.001 Homo sapiens
# Costimulation by the CD28 family 7 97 0.007 6.35E-4 1.01E-1 18 35 0.003 Homo sapiens
# Interleukin-37 signaling 4 36 0.003 2.07E-3 2.2E-1 8 14 0.001 Homo sapiens
# RUNX3 regulates BCL2L11 (BIM) transcription 2 6 0 3.63E-3 2.86E-1 2 2 0 Homo sapiens
# RUNX3 regulates CDKN1A transcription 2 8 0.001 6.32E-3 3.98E-1 6 6 0 Homo sapiens
# Nef and signal transduction 2 9 0.001 7.93E-3 4.15E-1 2 5 0 Homo sapiens
# The role of Nef in HIV-1 replication and disease pathogenesis 3 29 0.002 9.22E-3 4.15E-1 5 20 0.001 Homo sapiens
# Interleukin-23 signaling 2 11 0.001 1.16E-2 4.53E-1 10 12 0.001 Homo sapiens
# Insulin-like Growth Factor-2 mRNA Binding Proteins (IGF2BPs/IMPs/VICKZs) bind RNA 2 13 0.001 1.59E-2 5.41E-1 3 3 0 Homo sapiens
# Interleukin-2 signaling 2 14 0.001 1.83E-2 5.41E-1 12 19 0.001 Homo sapiens
# Translocation of ZAP-70 to Immunological synapse 3 42 0.003 2.44E-2 5.41E-1 4 4 0 Homo sapiens
# Phosphorylation of CD3 and TCR zeta chains 3 45 0.003 2.91E-2 5.41E-1 5 7 0.001 Homo sapiens
# PD-1 signaling 3 45 0.003 2.91E-2 5.41E-1 1 5 0 Homo sapiens
# Interferon gamma signaling 8 250 0.018 3.22E-2 5.41E-1 3 16 0.001 Homo sapiens
# Interleukin-2 family signaling 3 47 0.003 3.25E-2 5.41E-1 20 59 0.004 Homo sapiens
# Keratan sulfate degradation 2 22 0.002 4.19E-2 5.41E-1 1 7 0.001 Homo sapiens
# Nef mediated downregulation of CD28 cell surface expression 1 3 0 4.3E-2 5.41E-1 3 3 0 Homo sapiens
# SMAD2/3 MH2 Domain Mutants in Cancer 1 3 0 4.3E-2 5.41E-1 1 1 0 Homo sapiens
# SMAD4 MH2 Domain Mutants in Cancer 1 3 0 4.3E-2 5.41E-1 1 1 0 Homo sapiens
# Loss of Function of SMAD4 in Cancer 1 3 0 4.3E-2 5.41E-1 1 1 0 Homo sapiens

MSSNP <- unique(GWAS_MS$ID)
HashimotoSNP <- unique(GWAS_Hashimoto$ID)
RArthritisSNP <- unique(GWAS_RArthritis$ID)
LupusSNP <- unique(GWAS_SLupusE$ID)
T1DiabetesSNP <- unique(GWAS_T1Diabetes$ID)
IBDSNP <- unique(GWAS_IBD$ID)
GravesSNP <- unique(GWAS_Graves$ID)
AsthmaSNP <- unique(GWAS_Asthma$ID)
POSSNP <- unique(GWAS_POS$ID)


#AD Celiac

library(VennDiagram)
library(RColorBrewer)

# Chart
myCol <- brewer.pal(4, "Spectral")

# Compare MS, Arthritis, Lupus, T1Diabetes
venn.diagram(
  x = list(MSSNP, RArthritisSNP, LupusSNP, T1DiabetesSNP),
  category.names = c("MS" , "RArthritis" , "Lupus", "T1Diabetes"),
  filename = 'SNPs_venn_diagramm.png',
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

# Two SNPs in common for all. Lupus and RArthritis many in common and those two + T1Diabetes

myCol2 <- brewer.pal(5, "Spectral")

# Compare MS, Arthritis, Lupus, T1Diabetes, Asthma

venn.diagram(
  x = list(MSSNP, RArthritisSNP, LupusSNP, T1DiabetesSNP, AsthmaSNP),
  category.names = c("MS" , "RArthritis" , "Lupus", "T1Diabetes", "Asthma"),
  filename = 'SNPs_venn_diagramm4.png',
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
  x = list(IBDSNP, GravesSNP, CeliacSNP, AsthmaSNP),
  category.names = c("IBD", "Graves", "Celiac", "Asthma"),
  filename = 'SNPs_venn_diagramm2.png',
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
  x = list(IBDSNP, GravesSNP, HashimotoSNP, CeliacSNP, AsthmaSNP),
  category.names = c("IBD", "Graves", "Hashimoto's", "Celiac", "Asthma"),
  filename = 'SNPs_venn_diagramm3.png',
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
  x = list(MSSNP, RArthritisSNP, LupusSNP, AsthmaSNP),
  category.names = c("MS", "RArthritis", "Lupus", "Asthma"),
  filename = 'SNPs_venn_diagramm5.png',
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
  x = list(IBDSNP, T1DiabetesSNP, CeliacSNP, AsthmaSNP),
  category.names = c("IBD", "T1Diabetes", "Celiac", "Asthma"),
  filename = 'SNPs_venn_diagramm6.png',
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

# 4 common SNPs, many when excluding asthma, asthma/IBD, IBD/T1Diabetes, T1Diabetes/asthma

# Compare IBD, Graves, Hashimoto's, Celiac

venn.diagram(
  x = list(IBDSNP, GravesSNP, HashimotoSNP, CeliacSNP),
  category.names = c("IBD", "Graves", "Hashimoto's", "Celiac"),
  filename = 'SNPs_venn_diagramm7.png',
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

# One with all, Hashimoto's and Graves similar

#Looking at the overlap between Arthritis, lupus, asthma

filter <- RArthritisSNP[RArthritisSNP %in% LupusSNP]


tbl_RA <- GWAS_RArthritis[GWAS_RArthritis$ID %in% filter, c("ID", "allele", "P-value")]
tbl_l <- GWAS_SLupusE[GWAS_SLupusE$ID %in% filter, c("ID", "allele", "P-value")]
tbl_a <- GWAS_Asthma[GWAS_Asthma$ID %in% filter, c("ID", "allele", "P-value")]

RA_list <- c()
l_list <- c()
a_list <- c()

for (i in unique(tbl_a$ID)) {
  RA_alleles <- paste(unique(tbl_RA[tbl_RA$ID == i,2][[1]]), collapse = "")
  l_alleles <- paste(unique(tbl_l[tbl_l$ID == i,2][[1]]), collapse = "")
  a_alleles <- paste(unique(tbl_a[tbl_a$ID == i,2][[1]]), collapse = "")
  
  RA_list <- append(RA_list, RA_alleles)
  l_list <- append(l_list, l_alleles)
  a_list <- append(a_list, a_alleles)
  
}

RA_list
l_list
a_list    
SNP<- unique(tbl_a$ID)


#Data Frame for mix of CE C and AD SNPS
data.frame(SNP, l_list, RA_list, a_list)

# Result table
#         SNP l_list RA_list a_list
# 1 rs62324212      A       A      ?
# 2  rs1689510      C       C     ?C
# 3   rs911263      ?      ?A      C
# 4 rs61839660      ?       ?     ?T
# 5  rs7725052      C       C      T

# Looking at overlap between IBD, T1Diabetes, Asthma
filter2 <- IBDSNP[IBDSNP %in% AsthmaSNP]


tbl_IBD <- GWAS_IBD[GWAS_IBD$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_d <- GWAS_T1Diabetes[GWAS_T1Diabetes$ID %in% filter2, c("ID", "allele", "P-value")]
#tbl_c <- GWAS_Celiac[GWAS_Celiac$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_a2 <- GWAS_Asthma[GWAS_Asthma$ID %in% filter2, c("ID", "allele", "P-value")]


IBD_list <- c()
d_list <- c()
#c_list <- c()
a2_list <- c()

# Comparing all the diseases


for (i in unique(tbl_a2$ID)) {
  IBD_alleles <- paste(unique(tbl_IBD[tbl_IBD$ID == i,2][[1]]), collapse = "")
  d_alleles <- paste(unique(tbl_d[tbl_d$ID == i,2][[1]]), collapse = "")
  #c_alleles <- paste(unique(tbl_c[tbl_c$ID == i,2][[1]]), collapse = "")
  a2_alleles <- paste(unique(tbl_a2[tbl_a2$ID == i,2][[1]]), collapse = "")
  
  IBD_list <- append(IBD_list, IBD_alleles)
  d_list <- append(d_list, d_alleles)
  #c_list <- append(c_list, d_alleles)
  a2_list <- append(a2_list, a2_alleles)
  
}

IBD_list
d_list
#c_list
a2_list    
SNP<- unique(tbl_a2$ID)


#Data Frame for mix of CE C and AD SNPS
data.frame(SNP, d_list, IBD_list, c_list, a2_list)


# result table
#           SNP d_list IBD_list c_list a2_list
# 1   rs4845604              ?G             ?G
# 2  rs17293632             ?AT             ?T
# 3   rs7936312               G             ?T
# 4  rs62324212      A        A      A       ?
# 5  rs17622378              ?G              A
# 6    rs174535               ?              T
# 7  rs11236797               ?             ?A
# 8   rs2066844               ?             ?T
# 9   rs1811711              ?C              ?  
# 10  rs1799964               C              ?
# 11  rs2305480               T              G
# 12   rs449454               ?            ?AG
# 13  rs4129267               ?             CT
# 14  rs7927894               T              T
# 15  rs2155219              TA              T
# 16  rs6894249               ?              G
# 17  rs4795397               ?              A
# 18 rs76181804               ?              A
# 19  rs1689510      C        C      C      ?C
# 20   rs734999               C              C
# 21 rs17264332               A              A
# 22     rs1535               ?              A
# 23   rs653178      C       G?      C       T
# 24 rs12946510             A?T              C
# 25 rs61839660     C?        ?     C?      ?T
# 26  rs7725052      C        C      C       T


# Only the all in common
#           SNP d_list IBD_list c_list a2_list
# 4  rs62324212      A        A      A       ?
# 19  rs1689510      C        C      C      ?C
# 23   rs653178      C       G?      C       T
# 25 rs61839660     C?        ?     C?      ?T
# 26  rs7725052      C        C      C       T

# Comparing all the diseases
filter2 <- IBDSNP[IBDSNP %in% AsthmaSNP]

tbl_IBD <- GWAS_IBD[GWAS_IBD$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_d <- GWAS_T1Diabetes[GWAS_T1Diabetes$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_c <- GWAS_Celiac[GWAS_Celiac$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_a3 <- GWAS_Asthma[GWAS_Asthma$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_RA <- GWAS_RArthritis[GWAS_RArthritis$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_l <- GWAS_SLupusE[GWAS_SLupusE$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_MS <- GWAS_MS[GWAS_MS$ID %in% filter2, c("ID", "allele", "P-value")]
tbl_POS <- GWAS_POS[GWAS_POS$ID %in% filter2, c("ID", "allele", "P-value")]

RA_list <- c()
l_list <- c()
IBD_list <- c()
d_list <- c()
c_list <- c()
a3_list <- c()
MS_list <- c()
POS_list <- c()


for (i in unique(tbl_a3$ID)) {
  IBD_alleles <- paste(unique(tbl_IBD[tbl_IBD$ID == i,2][[1]]), collapse = "")
  d_alleles <- paste(unique(tbl_d[tbl_d$ID == i,2][[1]]), collapse = "")
  c_alleles <- paste(unique(tbl_c[tbl_c$ID == i,2][[1]]), collapse = "")
  a3_alleles <- paste(unique(tbl_a3[tbl_a3$ID == i,2][[1]]), collapse = "")
  RA_alleles <- paste(unique(tbl_RA[tbl_RA$ID == i,2][[1]]), collapse = "")
  l_alleles <- paste(unique(tbl_l[tbl_l$ID == i,2][[1]]), collapse = "")
  MS_alleles <- paste(unique(tbl_MS[tbl_MS$ID == i,2][[1]]), collapse = "")
  POS_alleles <- paste(unique(tbl_POS[tbl_POS$ID == i,2][[1]]), collapse = "")
  
  MS_list <- append(MS_list, MS_alleles)
  RA_list <- append(RA_list, RA_alleles)
  l_list <- append(l_list, l_alleles)
  IBD_list <- append(IBD_list, IBD_alleles)
  d_list <- append(d_list, d_alleles)
  c_list <- append(c_list, d_alleles)
  a3_list <- append(a3_list, a3_alleles)
  POS_list <- append(POS_list, POS_alleles)
  
}

MS_list
a3_list
d_list
c_list
IBD_list
RA_list
l_list
POS_list
SNP<- unique(tbl_a3$ID)

#Data Frame for mix of all the diseases
data.frame(SNP, MS_list, d_list, IBD_list, c_list, RA_list, l_list, a3_list, POS_list)

