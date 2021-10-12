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



#See if any SNPs Match up

CeliacSNP <- unique(GWAS_Celiac$ID)
MSSNP <- unique(GWAS_MS$ID)
HashimotoSNP <- unique(GWAS_Hashimoto$ID)
RArthritisSNP <- unique(GWAS_RArthritis$ID)
LupusSNP <- unique(GWAS_SLupusE$ID)
T1DiabetesSNP <- unique(GWAS_T1Diabetes$ID)
IBDSNP <- unique(GWAS_IBD$ID)
GravesSNP <- unique(GWAS_Graves$ID)
AsthmaSNP <- unique(GWAS_Asthma$ID)


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

#Looking at the overlap between AD C and CE

#filter <- CeliacSNP[CeliacSNP %in% EczemaSNP]


#tbl_ce <- GWAS_Celiac[GWAS_Celiac$ID %in% filter, c("ID", "allele", "P-value")]
#tbl_c <- GWAS_Crohns[GWAS_Crohns$ID %in% filter, c("ID", "allele", "P-value")]
#tbl_ad <- GWAS_Eczema[GWAS_Eczema$ID %in% filter, c("ID", "allele", "P-value")]

#ce_list <- c()
#c_list <- c()
#ad_list <- c()

#for (i in unique(tbl_ad$ID)) {
 # ce_alleles <- paste(unique(tbl_ce[tbl_ce$ID == i,2][[1]]), collapse = "")
#  c_alleles <- paste(unique(tbl_c[tbl_ce$ID == i,2][[1]]), collapse = "")
 # ad_alleles <- paste(unique(tbl_ad[tbl_ad$ID == i,2][[1]]), collapse = "")
  
  #ce_list <- append(ce_list, ce_alleles)
  #c_list <- append(c_list, c_alleles)
  #ad_list <- append(ad_list, ad_alleles)
  
#}

#ce_list
#c_list
#ad_list    
#SNP<- unique(tbl_ad$ID)


#Data Frame for mix of CE C and AD SNPS
#data.frame(SNP, c_list, ce_list, ad_list)
