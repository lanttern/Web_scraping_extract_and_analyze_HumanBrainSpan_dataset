#--------------------------------------------------------------------------------
#Please load library before analysis
#--------------------------------------------------------------------------------
library(rChoiceDialogs)
library(colorRamps)
library(devtools)
library(psych)
library(coin)
library(gplots)

#--------------------------------------------------------------------------------
#Define 3 functions: 1) calculate correlation, 2) creat.heatmap, 3) calculate p value
#Please run these 3 function before doing mouse or human analysis
#--------------------------------------------------------------------------------

cal.corr = function(data, index1, index2, method){
  #
  cor_matrix_all = cor(t(data[, index1:index2]), method = method);
  return (cor_matrix_all)
}

creat.heatmap = function (data, rownames, scale, pwd.image){
  #rename row names
  row.names(data) = rownames;
  #transform data to be matrix
  datamatrix = data.matrix(data);
  #setup tiff file
  tiff(filename = pwd.image, res = 300, height = 10, width = 10, units = "in");
  #creat expression heatmap
  heatmap.2(
    datamatrix, 
    Rowv=FALSE, 
    Colv=FALSE,
    revC = FALSE,
    scale = scale,
    col = colorRampPalette(c("blue", "white", "red"))(4000),
    dendrogram="none",
    key = TRUE,
    trace = "none",
    keysize = 1,
    density.info="none",
    key.title = "",
    key.xlab = "",
    main = "");
  dev.off();
  print ("Heatmap created!")
}

cal.pvalue = function(data, index){
  print ("p value of all MET-candidates correlated to MET as against high expressed genes correlated to MET");
  all_p = wilcox.test(data[2:index], data[index+1:length(data)], correct = TRUE);
  print (all_p)
}


#--------------------------------------------------------------------------------
#Analysis for  human data
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Read human data
#--------------------------------------------------------------------------------
#get file directory
getwd()
#select folder, run the following two lines before running other command line
setwd("/Volumes")
pwd.h = rchoose.dir(default = getwd(), caption = "Select Directory") 
#choose human data
pwd.data.h = paste(pwd.h, "data/expression data/clean data/gene expression of MET-interacting candidates.csv", sep = "/")
pwd.rna.seq.data.h = paste(pwd.h, "data/expression data/clean data/expression_matrix_add_Row_and_Column_4mon to 48 mon_STC-ITC-V1C.csv", sep = "/")
pwd.image.dir.h = paste(pwd.h, "result", sep = "/")
#read human data
data.h = read.csv(pwd.data.h, header = TRUE, sep = ",")
rna.seq.data.h = read.csv(pwd.rna.seq.data.h, header = TRUE, sep = ",")
#check data
dim(data.h)
dim(rna.seq.data.h)

#--------------------------------------------------------------------------------
#Calculate correlation for human data
#--------------------------------------------------------------------------------
#calculate correlation for human data
row.names(data.h) = data.h$Symbol
cor.exprs.h = cal.corr(data.h, index1 = 2, index2 = dim(data.h)[2], method = "spearman")
cor.ITC.h = cal.corr(data.h, index1 = 2, index2 = 6, method = "spearman")
cor.STC.h = cal.corr(data.h, index1 = 7, index2 = 11, method = "spearman")
cor.V1C.h = cal.corr(data.h, index1 = 12, index2 = 16, method = "spearman")
cor.MET.h = cbind(cor.ITC.h[, 1], cor.STC.h[, 1], cor.V1C.h[, 1])

#--------------------------------------------------------------------------------
#Creat heatmap for human data
#--------------------------------------------------------------------------------
data.exprs.h = data.h[,2:dim(data.h)[2]]
pwd.image.exprs.h = paste(pwd.image.dir.h, "Heatmap of gene expression of 72 candidates in human cortex.tiff", 
                          sep = "/")
creat.heatmap(data.exprs.h, rownames(data.exprs.h), scale = "row", pwd.image.exprs.h)
#correlation heatmap
pwd.image.cor.exprs.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex.tiff", 
                              sep = "/")
pwd.image.cor.ITC.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex ITC.tiff", 
                            sep = "/")
pwd.image.cor.STC.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex STC.tiff", 
                            sep = "/")
pwd.image.cor.V1C.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex V1C.tiff", 
                            sep = "/")
pwd.image.cor.MET.h = paste(pwd.image.dir.h, "Heatmap of expression correlation to MET of 72 candidates in human cortex.tiff", 
                            sep = "/")
creat.heatmap(cor.exprs.h, rownames(cor.exprs.h), scale = "none", pwd.image.cor.exprs.h)
creat.heatmap(cor.ITC.h, rownames(cor.ITC.h), scale = "none", pwd.image.cor.ITC.h)
creat.heatmap(cor.STC.h, rownames(cor.STC.h), scale = "none", pwd.image.cor.STC.h)
creat.heatmap(cor.V1C.h, rownames(cor.V1C.h), scale = "none", pwd.image.cor.V1C.h)
creat.heatmap(cor.MET.h, rownames(cor.MET.h), scale = "none", pwd.image.cor.MET.h)


#--------------------------------------------------------------------------------
#randomly select data with mean expression over a cutoff(200) to create a mixed dataset
#--------------------------------------------------------------------------------
#create mixed data
cutoff = 200
dup.h = subset(data.h, rowMeans(data.h[,2:dim(data.h)[2]])>cutoff)
rna.seq.data.highexpr.h = subset(rna.seq.data.h, rowMeans(rna.seq.data.h[,2:dim(rna.seq.data.h)[2]])>cutoff)
row.names(rna.seq.data.highexpr.h) = rna.seq.data.highexpr.h$Symbol
dim(rna.seq.data.highexpr.h)
rna.seq.data.highexpr.h = subset(rna.seq.data.highexpr.h, !rownames(rna.seq.data.highexpr.h) %in% rownames(dup.h))
print ("Number and percentage of selected high expressed genes")
dim(rna.seq.data.highexpr.h)
dim(rna.seq.data.highexpr.h)[1]*100/dim(rna.seq.data.h)[1]
mix.data.h = rbind(data.h, rna.seq.data.highexpr.h)
#check mixed data
dim(mix.data.h)
#calculate correlation using mixed data
cor.mix.data.h = cal.corr(mix.data.h, index1 = 2, index2 = dim(mix.data.h)[2], method = "spearman")
cor.mix.ITC.h = cal.corr(mix.data.h, index1 = 2, index2 = 6, method = "spearman")
cor.mix.STC.h = cal.corr(mix.data.h, index1 = 7, index2 = 11, method = "spearman")
cor.mix.V1C.h = cal.corr(mix.data.h, index1 = 12, index2 = 16, method = "spearman")
cor.mix.MET.h = cbind(cor.mix.ITC.h[, 1], cor.mix.STC.h[, 1], cor.mix.V1C.h[, 1])

#permutation test
count = 0
num = 10000000
test0 = wilcox.test(cor.mix.data.h[2:73,1], cor.mix.data.h[74:dim(cor.mix.data.h)[1],1])
x = c(1:num)
for (i in 1:num){
  cor.mix.data.h.df = sample(cor.mix.data.h[2:dim(cor.mix.data.h)[1]], replace = TRUE)
  test = wilcox.test(cor.mix.data.h.df[1:72], cor.mix.data.h.df[73:length(cor.mix.data.h.df)])
  #test = wilcox.test(cor.mix.data.h.df[1:8], cor.mix.data.h.df[73:length(cor.mix.data.h.df)])
  x[i] = test$statistic
  if (test$statistic > test0$statistic) {count = count +1}
}
print(count); (count/num);
hist(x, breaks = c(min(x):test0$statistic, by = 1000))
abline(v=test0$statistic, lwd = 2, col = "red")


#create heatmap for mix data
pwd.image.mix.data.h = paste(pwd.image.dir.h, "Heatmap of expression of 72 candidates and random selected genes in human cortex.tiff", 
                             sep = "/")
creat.heatmap(mix.data.h[, 2:dim(mix.data.h)[2]], rownames(mix.data.h), scale = "row", pwd.image.mix.data.h)
pwd.image.cor.mix.data.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates and random selected genes in human cortex.tiff", 
                                 sep = "/")
creat.heatmap(cor.mix.data.h, rownames(cor.mix.data.h), scale = "none", pwd.image.cor.mix.data.h)
pwd.image.cor.mix.MET.h = paste(pwd.image.dir.h, "Heatmap of expression correlation to MET of 72 candidates and random selected genes in human cortex ITC, STC and V1C.tiff", 
                                sep = "/")
creat.heatmap(cor.mix.MET.h, rownames(cor.mix.MET.h), scale = "none", pwd.image.cor.mix.MET.h)


#--------------------------------------------------------------------------------
#Statistics for human data
#--------------------------------------------------------------------------------

#wilcoxon-mann-whitney test - compare MET binding NDD candidates and non-NDD candidates
print ("p value - compare if ranked correlation differently between NDD candidates and non-NDD candidates")
wilcox.test(cor.exprs.h[2:9,1], cor.exprs.h[10:dim(cor.exprs.h)[1],1])

#wilcoxon-mann-whitney test - compare brain structures
print ("p value - compare if ranked correlation differently between ITC and STC")
wilcox.test(cor.MET.h[2:dim(cor.MET.h)[1],1], cor.MET.h[2:dim(cor.MET.h)[1],2])
print ("p value - compare if ranked correlation differently between ITC and V1C")
wilcox.test(cor.MET.h[2:dim(cor.MET.h)[1],1], cor.MET.h[2:dim(cor.MET.h)[1],3])
print ("p value - compare if ranked correlation differently between STC and V1C")
wilcox.test(cor.MET.h[2:dim(cor.MET.h)[1],2], cor.MET.h[2:dim(cor.MET.h)[1],3])


#wilcoxon-mann-whitney test - compare MET binding candidates and non-binding candidates
print ("p value - compare if ranked correlation differently between MET binding candidates and non-binding candidates")
cal.pvalue (cor.mix.data.h[,1], dim(data.h)[1])

#wilcoxon-mann-whitney test - compare MET binding NDD candidates and non-binding candidates
print ("p value - compare if ranked correlation differently between MET binding NDD candidates and non-binding candidates")
wilcox.test(cor.mix.data.h[2:9,1], cor.mix.data.h[74:dim(cor.mix.data.h)[1],1])

#--------------------------------------------------------------------------------
#Versions of tools
#--------------------------------------------------------------------------------
#session info
devtools::session_info()
