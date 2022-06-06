#### This is the code of interpreting the result of RNA seq analysis by limma and edgeR 
#### Need to run after "1.Normalize.R" and "2.FitModel.R" 

###### import packages ####
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)

###### Write our analysis result of each comparison into data frame ####
##### not using any criteria yet ####
##### n = Inf to keep all the genes in each comparison 
BATvsBAT_Cold. <- topTreat(efit, coef=1, n=Inf)
WATvsWAT_Cold. <- topTreat(efit, coef=2, n=Inf)
WATvsBAT. <- topTreat(efit, coef=3, n=Inf)
WAT_ColdvsBAT_Cold. <- topTreat(efit, coef=4, n=Inf)

## There are some genes didn't get annotated earlier
## get rid of genes without name, not useful for our interpretation
BATvsBAT_Cold <- BATvsBAT_Cold.[!is.na(BATvsBAT_Cold.$ENTREZID),]
WATvsWAT_Cold <- WATvsWAT_Cold.[!is.na(WATvsWAT_Cold.$ENTREZID),]
WATvsBAT <- WATvsBAT.[!is.na(WATvsBAT.$ENTREZID),]
WAT_ColdvsBAT_Cold <- WAT_ColdvsBAT_Cold.[!is.na(WAT_ColdvsBAT_Cold.$ENTREZID),]

## Label genes based on our predered criteria 
## Here I used logFC >1 and adj.P-value < 0.05 as differential expressed 
BATvsBAT_Cold <- BATvsBAT_Cold %>%
  mutate(gene_type = case_when(logFC > 1 & adj.P.Val <= 0.05 ~ "up",
                               logFC < -1 & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
WATvsWAT_Cold <- WATvsWAT_Cold %>%
  mutate(gene_type = case_when(logFC > 1 & adj.P.Val <= 0.05 ~ "up",
                               logFC < -1 & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
WATvsBAT <- WATvsBAT %>%
  mutate(gene_type = case_when(logFC > 1 & adj.P.Val <= 0.05 ~ "up",
                               logFC < -1 & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))   
WAT_ColdvsBAT_Cold <- WAT_ColdvsBAT_Cold %>%
  mutate(gene_type = case_when(logFC > 1 & adj.P.Val <= 0.05 ~ "up",
                               logFC < -1 & adj.P.Val <= 0.05 ~ "down",
                               TRUE ~ "ns"))

## cout how many dfGenes we have in each comparison 
BATvsBAT_Cold %>% count(gene_type)
WATvsWAT_Cold %>% count(gene_type)
WATvsBAT %>% count(gene_type)
WAT_ColdvsBAT_Cold %>% count(gene_type)

