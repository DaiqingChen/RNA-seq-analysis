# RNA-seq-analysis
When we are on our way of conducting bulk RNA-seq analysis, there are plenty of tutorials telling about the standard workflow of getting a differential expressed gene list. Mainly with 2 conditions, a treatmemt and a baseline group, the differential expressed gene list gives you some idea of what genes are expressed differentially during the treatment. 
But what if we have more than 2 conditions? What should we do after getting the differential expressed gene list? I mean, mostly people won't recognize all the genes showing up in the list, how to we interepret the result?
Today I want to share my workflow of conducting RNAseq analysis, using the dataset from Cheng et al.(https://doi.org/10.1016/j.celrep.2018.05.021
 ), a 4 conditions bulk RNA-seq dataset. You can download from GSE112582(https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE112582).

Dataset desciption:
The purpose of this study is to figure out the key difference between white and brown adipose tissue that cause the different capacity of browning during cold exposure. There are two influences in this study: white or brown adipose tissue, also room temperature and cold exposure. Thus, there are 4 groups of samples: WAT(white adipose tissue) at room temperature, WAT at cold, BAT(brown adipose tissue) at room temperature and BAT at cold.

Before we start, I assume you already have some idea of conducting RNA-seq analysis. If you still having problems setting up R or downloading packages, I highly recommed you to go through this training: https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/#part-i .

I prefer edgeR over deSeq, especially for more than 2 conditions, for you have more freedom to customize your design matrix. 

