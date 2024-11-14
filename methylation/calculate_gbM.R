###################
# calculate_gbM.R #
# Jessie Pelosi   #
# 2024            #
#                 #
###################

library(dplyr)

### LYMI G1 ###

# CG 

LYMI_G1_genes.CG <- read.delim("v1.1.9/LYMI-G1_methyl.CG.genes.bedgraph", header = F)

LYMI_G1_filtGenes <- LYMI_G1_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G1_filtGenes_summary <- LYMI_G1_filtGenes %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.773){binom.test(x=a, n=b, p=0.733, alternative = "greater")$p.value}

LYMI_G1_filtGenes_summary$pVal <- mapply(bt, LYMI_G1_filtGenes_summary$methylated_cytosines, 
                                         LYMI_G1_filtGenes_summary$total_cytosines)

LYMI_G1_filtGenes_summary$BFpVal <- p.adjust(LYMI_G1_filtGenes_summary$pVal, method = "bonferroni")

LYMIG1_BM_genes <- filter(LYMI_G1_filtGenes_summary, BFpVal < 0.05) 

LYMIG1_IM_genes <- filter(LYMI_G1_filtGenes_summary, BFpVal < 0.95 & BFpVal > 0.05)

LYMIG1_UM_genes <- filter(LYMI_G1_filtGenes_summary, BFpVal > 0.95) 

write.csv(file="LYMI-G1.UM.csv", LYMIG1_UM_genes, quote = F)
  
# CHG 

LYMI_G1_genes.CHG <- read.delim("v1.1.9/LYMI-G1_methyl.CHG.genes.bedgraph", header = F)

LYMI_G1_filtGenes.CHG <- LYMI_G1_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G1_filtGenes.CHG_summary <- LYMI_G1_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.645){binom.test(x=a, n=b, p=0.645, alternative = "greater")$p.value}

LYMI_G1_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_G1_filtGenes.CHG_summary$methylated_cytosines, 
                                         LYMI_G1_filtGenes.CHG_summary$total_cytosines)

LYMI_G1_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_G1_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIG1.CHG_BM_genes <- filter(LYMI_G1_filtGenes.CHG_summary, BFpVal < 0.05) #80 genes 

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIG1.CG_BM.tmp <- anti_join(LYMIG1_BM_genes, LYMIG1.CHG_BM_genes, by = c("V15"))

# CHH 

LYMI_G1_genes.CHH <- read.delim("v1.1.9/LYMI-G1_methyl.CHH.genes.bedgraph", header = F)

LYMI_G1_filtGenes.CHH <- LYMI_G1_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G1_filtGenes.CHH_summary <- LYMI_G1_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.019){binom.test(x=a, n=b, p=0.019, alternative = "greater")$p.value}

LYMI_G1_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_G1_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_G1_filtGenes.CHH_summary$total_cytosines)

LYMI_G1_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_G1_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIG1.CHH_BM_genes <- filter(LYMI_G1_filtGenes.CHH_summary, BFpVal < 0.05) #2460 genes 

# Remove BM genes that have adjusted PCHH < 0.05 

LYMIG1.CG_BM <- anti_join(LYMIG1.CG_BM.tmp, LYMIG1.CHH_BM_genes, by = c("V15")) #3858 genes

write.csv(file="LYMI-G1.BM.csv", LYMIG1.CG_BM, quote = F)


### LYMI-G2 ###

# CG 

LYMI_G2_genes.CG <- read.delim("v1.1.9/LYMI-G2_methyl.CG.genes.bedgraph", header = F)

LYMI_G2_filtGenes <- LYMI_G2_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G2_filtGenes_summary <- LYMI_G2_filtGenes %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.583){binom.test(x=a, n=b, p=0.583, alternative = "greater")$p.value}

LYMI_G2_filtGenes_summary$pVal <- mapply(bt, LYMI_G2_filtGenes_summary$methylated_cytosines, 
                                         LYMI_G2_filtGenes_summary$total_cytosines)

LYMI_G2_filtGenes_summary$BFpVal <- p.adjust(LYMI_G2_filtGenes_summary$pVal, method = "bonferroni")

LYMIG2_BM_genes <- filter(LYMI_G2_filtGenes_summary, BFpVal < 0.05) 

LYMIG2_IM_genes <- filter(LYMI_G2_filtGenes_summary, BFpVal < 0.95 & BFpVal > 0.05)

LYMIG2_UM_genes <- filter(LYMI_G2_filtGenes_summary, BFpVal > 0.95) 

write.csv(file="LYMI-G2.UM.csv", LYMIG2_UM_genes, quote = F)

# CHG 

LYMI_G2_genes.CHG <- read.delim("v1.1.9/LYMI-G2_methyl.CHG.genes.bedgraph", header = F)

LYMI_G2_filtGenes.CHG <- LYMI_G2_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G2_filtGenes.CHG_summary <- LYMI_G2_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.568){binom.test(x=a, n=b, p=0.568, alternative = "greater")$p.value}

LYMI_G2_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_G2_filtGenes.CHG_summary$methylated_cytosines, 
                                             LYMI_G2_filtGenes.CHG_summary$total_cytosines)

LYMI_G2_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_G2_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIG2.CHG_BM_genes <- filter(LYMI_G2_filtGenes.CHG_summary, BFpVal < 0.05)  #157 genes

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIG2.CG_BM.tmp <- anti_join(LYMIG2_BM_genes, LYMIG2.CHG_BM_genes, by = c("V15"))

# CHH 

LYMI_G2_genes.CHH <- read.delim("v1.1.9/LYMI-G2_methyl.CHH.genes.bedgraph", header = F)

LYMI_G2_filtGenes.CHH <- LYMI_G2_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G2_filtGenes.CHH_summary <- LYMI_G2_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.017){binom.test(x=a, n=b, p=0.017, alternative = "greater")$p.value}

LYMI_G2_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_G2_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_G2_filtGenes.CHH_summary$total_cytosines)

LYMI_G2_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_G2_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIG2.CHH_BM_genes <- filter(LYMI_G2_filtGenes.CHH_summary, BFpVal < 0.05) #1863 genes

# Remove BM genes that have adjusted PCHH < 0.05 

LYMIG2.CG_BM <- anti_join(LYMIG2.CG_BM.tmp, LYMIG2.CHH_BM_genes, by = c("V15")) #5195 genes

write.csv(file="LYMI-G2.BM.csv", LYMIG2.CG_BM, quote = F)

### LYMI-G4 ###

# CG 

LYMI_G4_genes.CG <- read.delim("v1.1.9/LYMI-G4_methyl.CG.genes.bedgraph", header = F)

LYMI_G4_filtGenes <- LYMI_G4_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G4_filtGenes_summary <- LYMI_G4_filtGenes %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.823){binom.test(x=a, n=b, p=0.823, alternative = "greater")$p.value}

LYMI_G4_filtGenes_summary$pVal <- mapply(bt, LYMI_G4_filtGenes_summary$methylated_cytosines, 
                                         LYMI_G4_filtGenes_summary$total_cytosines)

LYMI_G4_filtGenes_summary$BFpVal <- p.adjust(LYMI_G4_filtGenes_summary$pVal, method = "bonferroni")

LYMIG4_BM_genes <- filter(LYMI_G4_filtGenes_summary, BFpVal < 0.05) #4600

LYMIG4_IM_genes <- filter(LYMI_G4_filtGenes_summary, BFpVal < 0.95 & BFpVal > 0.05) #264

LYMIG4_UM_genes <- filter(LYMI_G4_filtGenes_summary, BFpVal > 0.95) #13820

write.csv(file="LYMIG4_UM_genes.csv", LYMIG4_UM_genes, quote = F)
write.csv(file="LYMIG4_BM_genes-tmp.csv", LYMIG4_BM_genes, quote = F) 

# CHG 

LYMI_G4_genes.CHG <- read.delim("v1.1.9/LYMI-G4_methyl.CHG.genes.bedgraph", header = F)

LYMI_G4_filtGenes.CHG <- LYMI_G4_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G4_filtGenes.CHG_summary <- LYMI_G4_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.668){binom.test(x=a, n=b, p=0.668, alternative = "greater")$p.value}

LYMI_G4_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_G4_filtGenes.CHG_summary$methylated_cytosines, 
                                             LYMI_G4_filtGenes.CHG_summary$total_cytosines)

LYMI_G4_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_G4_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIG4.CHG_BM_genes <- filter(LYMI_G4_filtGenes.CHG_summary, BFpVal < 0.05) #87 genes  

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIG4.CG_BM.tmp <- anti_join(LYMIG4_BM_genes, LYMIG4.CHG_BM_genes, by = c("V15")) #4536

# CHH 

LYMI_G4_genes.CHH <- read.delim("v1.1.9/LYMI-G4_methyl.CHH.genes.bedgraph", header = F)

LYMI_G4_filtGenes.CHH <- LYMI_G4_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_G4_filtGenes.CHH_summary <- LYMI_G4_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.021){binom.test(x=a, n=b, p=0.021, alternative = "greater")$p.value}

LYMI_G4_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_G4_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_G4_filtGenes.CHH_summary$total_cytosines)

LYMI_G4_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_G4_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIG4.CHH_BM_genes <- filter(LYMI_G4_filtGenes.CHH_summary, BFpVal < 0.05) #3345

write.csv(file="v1.1.9/LYMIG4.CHH_BM.genes",LYMIG4.CHH_BM_genes, quote = F)
LYMIG4.CHH_BM_genes <- read.csv("v1.1.9/LYMIG4.CHH_BM.genes", row.names = 1)

# Remove BM genes that have adjusted PCHH < 0.05 

LYMIG4.CG_BM <- anti_join(LYMIG4.CG_BM.tmp, LYMIG4.CHH_BM_genes, by = c("V15")) #2067

write.csv(file="LYMI-G4.BM.csv", LYMIG4.CG_BM, quote = F)

## SPOROPHYTE LEAF TISSUES 

# CG

LYMI_S1_genes.CG <- read.delim("v1.1.9/LYMI-S1_methyl.CG.genes.bedgraph", header = F)

LYMI_S1_filtGenes.CG <- LYMI_S1_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S1_filtGenes.CG_summary <- LYMI_S1_filtGenes.CG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.937){binom.test(x=a, n=b, p=0.937, alternative = "greater")$p.value}

LYMI_S1_filtGenes.CG_summary$pVal <- mapply(bt, LYMI_S1_filtGenes.CG_summary$methylated_cytosines, 
                                             LYMI_S1_filtGenes.CG_summary$total_cytosines)

LYMI_S1_filtGenes.CG_summary$BFpVal <- p.adjust(LYMI_S1_filtGenes.CG_summary$pVal, method = "bonferroni")

LYMIS1.CG_BM_genes <- filter(LYMI_S1_filtGenes.CG_summary, BFpVal < 0.05) #2004 

LYMIS1.CG_IM_genes <- filter(LYMI_S1_filtGenes.CG_summary, BFpVal > 0.05 & BFpVal < 0.95)

LYMIS1.CG_UM_genes <- filter(LYMI_S1_filtGenes.CG_summary, BFpVal > 0.95) #17095 

write.csv(file = "LYMI-S1.BM.CG.csv", LYMIS1.CG_BM_genes, quote = F)
write.csv(file = "LYMI-S1.UM.csv", LYMIS1.CG_UM_genes, quote = F)

# CHG 

LYMI_S1_genes.CHG <- read.delim("v1.1.9/LYMI-S1_methyl.CHG.genes.bedgraph", header = F)

LYMI_S1_filtGenes.CHG <- LYMI_S1_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S1_filtGenes.CHG_summary <- LYMI_S1_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.694){binom.test(x=a, n=b, p=0.694, alternative = "greater")$p.value}

LYMI_S1_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_S1_filtGenes.CHG_summary$methylated_cytosines, 
                                             LYMI_S1_filtGenes.CHG_summary$total_cytosines)

LYMI_S1_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_S1_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIS1.CHG_BM_genes <- filter(LYMI_S1_filtGenes.CHG_summary, BFpVal < 0.05) #358 genes  

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIS1.CG_BM.tmp <- anti_join(LYMIS1.CG_BM_genes, LYMIS1.CHG_BM_genes, by = c("V15")) #1856

write.csv(file= "LYMI-S1.CHG.BM.csv", LYMIS1.CG_BM.tmp, quote = F)

# CHH 

LYMI_S1_genes.CHH <- read.delim("v1.1.9/LYMI-S1_methyl.CHH.genes.bedgraph", header = F)

LYMI_S1_filtGenes.CHH <- LYMI_S1_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S1_filtGenes.CHH_summary <- LYMI_S1_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.013){binom.test(x=a, n=b, p=0.013, alternative = "greater")$p.value}

LYMI_S1_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_S1_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_S1_filtGenes.CHH_summary$total_cytosines)

LYMI_S1_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_S1_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIS1.CHH_BM_genes <- filter(LYMI_S1_filtGenes.CHH_summary, BFpVal < 0.05) #5237

write.csv(file="v1.1.9/LYMIS1.CHH_BM.genes",LYMIS1.CHH_BM_genes, quote = F)

# Remove BM genes that have adjusted PCHH < 0.05 

LYMIS1.CG_BM <- anti_join(LYMIS1.CG_BM.tmp, LYMIS1.CHH_BM_genes, by = c("V15")) #278

write.csv(file="LYMI-S1.BM.csv", LYMIS1.CG_BM, quote = F)

### LYMI S2

# CG

LYMI_S2_genes.CG <- read.delim("v1.1.9/LYMI-S2_methyl.CG.genes.bedgraph", header = F)

LYMI_S2_filtGenes.CG <- LYMI_S2_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S2_filtGenes.CG_summary <- LYMI_S2_filtGenes.CG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.936){binom.test(x=a, n=b, p=0.936, alternative = "greater")$p.value}

LYMI_S2_filtGenes.CG_summary$pVal <- mapply(bt, LYMI_S2_filtGenes.CG_summary$methylated_cytosines, 
                                            LYMI_S2_filtGenes.CG_summary$total_cytosines)

LYMI_S2_filtGenes.CG_summary$BFpVal <- p.adjust(LYMI_S2_filtGenes.CG_summary$pVal, method = "bonferroni")

LYMIS2.CG_BM_genes <- filter(LYMI_S2_filtGenes.CG_summary, BFpVal < 0.05) #856 

LYMIS2.CG_IM_genes <- filter(LYMI_S2_filtGenes.CG_summary, BFpVal > 0.05 & BFpVal < 0.95)

LYMIS2.CG_UM_genes <- filter(LYMI_S2_filtGenes.CG_summary, BFpVal > 0.95) #17928 

write.csv(file = "LYMI-S2.BM.CG.csv", LYMIS2.CG_BM_genes, quote = F)
write.csv(file = "LYMI-S2.UM.csv", LYMIS2.CG_UM_genes, quote = F)

# CHG 

LYMI_S2_genes.CHG <- read.delim("v1.1.9/LYMI-S2_methyl.CHG.genes.bedgraph", header = F)

LYMI_S2_filtGenes.CHG <- LYMI_S2_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S2_filtGenes.CHG_summary <- LYMI_S2_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.690){binom.test(x=a, n=b, p=0.690, alternative = "greater")$p.value}

LYMI_S2_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_S2_filtGenes.CHG_summary$methylated_cytosines, 
                                             LYMI_S2_filtGenes.CHG_summary$total_cytosines)

LYMI_S2_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_S2_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIS2.CHG_BM_genes <- filter(LYMI_S2_filtGenes.CHG_summary, BFpVal < 0.05) #73  

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIS2.CG_BM.tmp <- anti_join(LYMIS2.CG_BM_genes, LYMIS2.CHG_BM_genes, by = c("V15")) #844

write.csv(file= "LYMI-S2.CHG.BM.csv", LYMIS1.CG_BM.tmp, quote = F)

# CHH 

LYMI_S2_genes.CHH <- read.delim("v1.1.9/LYMI-S2_methyl.CHH.genes.bedgraph", header = F)

LYMI_S2_filtGenes.CHH <- LYMI_S2_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S2_filtGenes.CHH_summary <- LYMI_S2_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.014){binom.test(x=a, n=b, p=0.014, alternative = "greater")$p.value}

LYMI_S2_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_S2_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_S2_filtGenes.CHH_summary$total_cytosines)

LYMI_S2_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_S2_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIS2.CHH_BM_genes <- filter(LYMI_S2_filtGenes.CHH_summary, BFpVal < 0.05) #4379

write.csv(file="v1.1.9/LYMIS2.CHH_BM.genes",LYMIS2.CHH_BM_genes, quote = F)

# Remove BM genes that have adjusted PCHH < 0.05 
# I had to read-in this file b/c CHH crashed R and lost environment 
LYMIS2.CG_BM.tmp <- read.csv("LYMI-S2.CHG.BM.csv")

LYMIS2.CG_BM <- anti_join(LYMIS2.CG_BM.tmp, LYMIS2.CHH_BM_genes, by = c("V15")) #509

write.csv(file="LYMI-S2.BM.csv", LYMIS2.CG_BM, quote = F)

### LYMI S2

# CG

LYMI_S3_genes.CG <- read.delim("v1.1.9/LYMI-S3_methyl.CG.genes.bedgraph", header = F)

LYMI_S3_filtGenes.CG <- LYMI_S3_genes.CG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S3_filtGenes.CG_summary <- LYMI_S3_filtGenes.CG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.937){binom.test(x=a, n=b, p=0.937, alternative = "greater")$p.value}

LYMI_S3_filtGenes.CG_summary$pVal <- mapply(bt, LYMI_S3_filtGenes.CG_summary$methylated_cytosines, 
                                            LYMI_S3_filtGenes.CG_summary$total_cytosines)

LYMI_S3_filtGenes.CG_summary$BFpVal <- p.adjust(LYMI_S3_filtGenes.CG_summary$pVal, method = "bonferroni")

LYMIS3.CG_BM_genes <- filter(LYMI_S3_filtGenes.CG_summary, BFpVal < 0.05) #1697

LYMIS3.CG_IM_genes <- filter(LYMI_S3_filtGenes.CG_summary, BFpVal > 0.05 & BFpVal < 0.95)

LYMIS3.CG_UM_genes <- filter(LYMI_S3_filtGenes.CG_summary, BFpVal > 0.95) #17232

write.csv(file = "LYMI-S3.BM.CG.csv", LYMIS3.CG_BM_genes, quote = F)
write.csv(file = "LYMI-S3.UM.csv", LYMIS3.CG_UM_genes, quote = F)

# CHG 

LYMI_S3_genes.CHG <- read.delim("v1.1.9/LYMI-S3_methyl.CHG.genes.bedgraph", header = F)

LYMI_S3_filtGenes.CHG <- LYMI_S3_genes.CHG %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S3_filtGenes.CHG_summary <- LYMI_S3_filtGenes.CHG %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.699){binom.test(x=a, n=b, p=0.699, alternative = "greater")$p.value}

LYMI_S3_filtGenes.CHG_summary$pVal <- mapply(bt, LYMI_S3_filtGenes.CHG_summary$methylated_cytosines, 
                                             LYMI_S3_filtGenes.CHG_summary$total_cytosines)

LYMI_S3_filtGenes.CHG_summary$BFpVal <- p.adjust(LYMI_S3_filtGenes.CHG_summary$pVal, method = "bonferroni")

LYMIS3.CHG_BM_genes <- filter(LYMI_S3_filtGenes.CHG_summary, BFpVal < 0.05) #346  

# Remove BM genes that have adjusted PCGH < 0.05 

LYMIS3.CG_BM.tmp <- anti_join(LYMIS3.CG_BM_genes, LYMIS3.CHG_BM_genes, by = c("V15")) #1564

write.csv(file= "LYMI-S2.CHG.BM.csv", LYMIS1.CG_BM.tmp, quote = F)

# CHH 

LYMI_S3_genes.CHH <- read.delim("v1.1.9/LYMI-S3_methyl.CHH.genes.bedgraph", header = F)

LYMI_S3_filtGenes.CHH <- LYMI_S3_genes.CHH %>% 
  group_by(V2, V3) %>%
  mutate(allReads = sum(V5, V6)) %>%
  filter(allReads >= 2)

LYMI_S3_filtGenes.CHH_summary <- LYMI_S3_filtGenes.CHH %>%
  group_by(V15) %>%
  summarize(total_cytosines = n(), methylated_cytosines = sum(V5 > 1)) %>%
  filter(total_cytosines >= 20)

# Calculate p-value using binomial approach from Takuno and Gaut (2012)

bt <- function(a,b,p=0.017){binom.test(x=a, n=b, p=0.017, alternative = "greater")$p.value}

LYMI_S3_filtGenes.CHH_summary$pVal <- mapply(bt, LYMI_S3_filtGenes.CHH_summary$methylated_cytosines, 
                                             LYMI_S3_filtGenes.CHH_summary$total_cytosines)

LYMI_S3_filtGenes.CHH_summary$BFpVal <- p.adjust(LYMI_S3_filtGenes.CHH_summary$pVal, method = "bonferroni")

LYMIS3.CHH_BM_genes <- filter(LYMI_S3_filtGenes.CHH_summary, BFpVal < 0.05) #5399

write.csv(file="v1.1.9/LYMIS3.CHH_BM.genes",LYMIS3.CHH_BM_genes, quote = F)

# Remove BM genes that have adjusted PCHH < 0.05 
# I had to read-in this file b/c CHH crashed R and lost environment 
LYMIS3.CHH_BM_genes <- read.csv("v1.1.9/LYMIS3.CHH_BM.genes")

LYMIS3.CG_BM <- anti_join(LYMIS3.CG_BM.tmp, LYMIS3.CHH_BM_genes, by = c("V15")) #186

write.csv(file="LYMI-S3.BM.csv", LYMIS3.CG_BM, quote = F)
