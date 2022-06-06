###### This is code for read in data, setting up the analysis and normalize data ###
###### we are going to use "limma" and "edgeR" for our analysis ####
###### you can change to other packages like deseq ###

######### Read in data #############
setwd("/Users/lasia/Documents/Cold study") # set working direction
raw_data<-read.table("data_raw/GSE112582_raw_count.txt") # read in raw_data

######## Set the data into matrix format ###########
df<-matrix(as.numeric(unlist(raw_data[-1,c(2:17)])),nrow=46078,ncol=16) 
rownames(df) <- raw_data[-1,1] # set the names of rows(gene name)
colnames(df) <- raw_data[1,-1] # set the names of cols(sample group)
df<-df[,c(5,6,1,2,7,8,3,4,13:16,9:12)] # re-order the columns

###### Import the packages ###########
library(limma) 
library(edgeR)
library(Glimma)

###### Set as an edgeR object #########
dgeObj <- DGEList(df) # Creates a DGEList object
group <-c(rep("BAT",4),rep("BAT_Cold",4),rep("WAT",4),rep("WAT_Cold",4)) # annotate the group 
dgeObj$samples$group <- group # set the group in our object

######### Remove low count rows by package default setting ########
keep.exprs <- filterByExpr(dgeObj, group=group) # Identify genes to remove by filterByExpr function of edgeR 
dgeObj <- dgeObj[keep.exprs,, keep.lib.sizes=FALSE] # keep the ones we filtered
dim(dgeObj) # look at how many rows left
## 19481 our of 46078 rows left

dgeObj$samples$lib.size # library size left in each column
# [1] 61657313 67117261 58469005 55738918 54258389 65961936 66102217 64870310
# [9] 66201305 61552351 61041280 57448584 66638143 62757911 64846758 71268555

head(dgeObj$counts) # take a look of what left

######### Annotate the gene names ################
### there are multiple packages that can do this 
### we are using org.Mm.eg.db today, but showing the code of using clusterProfiler and biomaRt

BiocManager::install("org.Mm.eg.db")
library("org.Mm.eg.db")
geneList <- select(org.Mm.eg.db, 
                   keys = c(rownames(dgeObj$counts)),
                   columns = c("ENTREZID", "SYMBOL","GENENAME"),
                   keytype = "ENSEMBL") # ENSEMBL

## example code of using clusterProfiler ##
# library(clusterProfiler)
# geneList2<-bitr(rownames(dgeObj$counts), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")

## example code of using biomaRt
# library( "biomaRt" ) 
# ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
# geneList3 <- getBM( attributes = c("ensembl_gene_id", "mgi_symbol",'entrezgene_id','external_gene_name'), filters = "ensembl_gene_id",values = rownames(dgeObj$counts), mart = ensembl )


dgeObj$genes <- geneList # assign the gene names to our object
dgeObj # take a look at each categories in our object

########## Normalize ##############
dgeObj<-calcNormFactors(dgeObj,method="TMM") # type ?calcNormFactors to see the mthod available
normalized_count<-data.frame(dgeObj$counts) # keep the normalized count as a data frame
write.csv(normalized_count,"data_clean/normalized_count.csv") # write to file





## End of file ##



