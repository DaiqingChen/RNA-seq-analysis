###### This is the code for making plots of RNA seq analysis ####
###### Need to run after "1.Normalize.R", "2.FitModel.R" and "3.Interpret.R"

###### import packages ####
library(ggplot2)
library(ggrepel)
library(tidyr)
library(dplyr)

########### Volcano plot #######
## Add colour, size and alpha (transparency) to volcano plot 
cols <- c("up" = "red", "down" = "blue", "ns" = "black") ##ffad73 #26b3ff
## plot by ggplot
pBATvsBAT_Cold <- ggplot(data=BATvsBAT_Cold, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(aes(colour = gene_type),size = 0.5) + scale_colour_manual(values = cols) + geom_vline(xintercept=c(-0.58, 0.58), col="black",linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red",linetype = "dashed")

pWATvsWAT_Cold <- ggplot(data=WATvsWAT_Cold, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(aes(colour = gene_type),size = 0.5) + scale_colour_manual(values = cols) + geom_vline(xintercept=c(-0.58, 0.58), col="black",linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red",linetype = "dashed")

pWATvsBAT <- ggplot(data=WATvsBAT, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(aes(colour = gene_type),size = 0.5) + scale_colour_manual(values = cols) + geom_vline(xintercept=c(-0.58, 0.58), col="black",linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red",linetype = "dashed")

pWAT_ColdvsBAT_Cold <- ggplot(data=WAT_ColdvsBAT_Cold, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(aes(colour = gene_type),size = 0.5) + scale_colour_manual(values = cols) + geom_vline(xintercept=c(-0.58, 0.58), col="black",linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="red",linetype = "dashed")

cowplot::plot_grid(pBATvsBAT_Cold, pWATvsWAT_Cold, pWATvsBAT, pWAT_ColdvsBAT_Cold, ncol=2, labels=LETTERS[1:4])

