##### This is the code of fitting model in edgeR and limma #####
##### Run after "1.Normalize.R" file #####

######### Take a look of how our data grouped ####
plotMDS(dgeObj$counts,col=rep(1:4, each=4)) # 4 colors for 4 groups
# we can see our data is well seperated by group 

######### Set design matrix ######
design <- model.matrix(~ 0+group)
rownames(design) <- colnames(dgeObj$counts)
design

######## Set contract matrix ####
contr.matrix <- makeContrasts(
  BATvsBAT_Cold = groupBAT_Cold - groupBAT, 
  WATvsWAT_Cold = groupWAT_Cold - groupWAT,
  WATvsBAT = groupBAT-groupWAT,
  WAT_ColdvsBAT_Cold = groupBAT_Cold - groupWAT_Cold,
  levels = colnames(design))
contr.matrix

######## Estimating the dispersion ######
# Removing heteroscedascity from count data
# The common dispersion estimates the overall BCV of the dataset, averaged over all genes:
dgeObj <- estimateDisp(dgeObj,design = design)
plotBCV(dgeObj) # Plot the BCV
# The BCV (Biological Coefficient of Variation) plot is a way to measure the biological variation within a particular condition. 
# A common dispersion (i.e. red line on the BCV plot) between 0.2 and 0.4 is usually considered reasonable and hence could detect more DE genes.

v <- voom(dgeObj, design, plot=TRUE) # look at the mean-variance trend

####### Fit into model ##########
vfit <- lmFit(v, design) # fit into lm model
vfit <- contrasts.fit(vfit, contrasts=contr.matrix) # by our design and contract matrix
efit <- eBayes(vfit) # fit ebayes model for better outcome
plotSA(efit)

####### Look at how many deGenes we got ####
de <- decideTests(efit)
summary(de) # the criteria is 5% FDR

tfit <- treat(efit, lfc=2) # change the criteria of significance to logFC=2
dt <- decideTests(tfit)
summary(dt)

# we can chose the criteria we want base on how many genes you want to focus on 



### End of file ###

