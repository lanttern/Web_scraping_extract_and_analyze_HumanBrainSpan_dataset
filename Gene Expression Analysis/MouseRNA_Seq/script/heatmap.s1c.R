#--------------------------------------------------------------------------------
#Load library
#--------------------------------------------------------------------------------
library(rChoiceDialogs)
library(colorRamps)
library(devtools)
library(psych)
library(coin)
library(gplots)

#--------------------------------------------------------------------------------
#Define 3 functions: 1) calculate correlation, 2) creat.heatmap, 3) calculate p value
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
#Analysis for  mouse data
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#Read mouse data
#--------------------------------------------------------------------------------
#get file directory
getwd()
#select folder, run the following two lines before running other command line
setwd("/Volumes")
pwd = rchoose.dir(default = getwd(), caption = "Select Directory")
#choose mouse data
pwd.data.m = paste(pwd, "data/clean data/Expression of 72 candidates in mouse cortex without adult.csv", sep = "/")
pwd.rna.seq.data.m = paste(pwd, "data/clean data/mouse_rna_seq_rename_mean_without adult.csv", sep = "/")
pwd.image.dir.m = paste(pwd, "result", sep = "/")

#read mouse data
data.m = read.csv(pwd.data.m, header = TRUE, sep = ",")
rna.seq.data.m = read.csv(pwd.rna.seq.data.m, header = TRUE, sep = ",")

#check data
dim(data.m)
dim(rna.seq.data.m)

#--------------------------------------------------------------------------------
#Calculate correlation for mouse data
#--------------------------------------------------------------------------------
row.names(data.m) = data.m$Geneid
cor.exprs.m = cal.corr(data.m, index1 = 2, index2 = dim(data.m)[2], method = "spearman")
cor.SgL.m = cal.corr(data.m, index1 = 2, index2 = 6, method = "spearman")
cor.IgL.m = cal.corr(data.m, index1 = 7, index2 = dim(data.m)[2], method = "spearman")
cor.MET.m = cbind(cor.SgL.m[, 1], cor.IgL.m[, 1])

#--------------------------------------------------------------------------------
#Creat heatmap for mouse data
#--------------------------------------------------------------------------------
#expression heatmap
data.exprs.m = data.m[,2:dim(data.m)[2]]
pwd.image.exprs.m = paste(pwd.image.dir.m, "Heatmap of gene expression of 72 candidates in mouse cortex.tiff", 
                        sep = "/")
creat.heatmap(data.exprs.m, rownames(data.exprs.m), scale = "row", pwd.image.exprs.m)
#correlation heatmap
pwd.image.cor.exprs.m = paste(pwd.image.dir.m, "Heatmap of expression correlation of 72 candidates in mouse cortex.tiff", 
                              sep = "/")
pwd.image.cor.SgL.m = paste(pwd.image.dir.m, "Heatmap of expression correlation of 72 candidates in mouse cortex SgL.tiff", 
                            sep = "/")
pwd.image.cor.IgL.m = paste(pwd.image.dir.m, "Heatmap of expression correlation of 72 candidates in mouse cortex IgL.tiff", 
                            sep = "/")
pwd.image.cor.MET.m = paste(pwd.image.dir.m, "Heatmap of expression correlation to MET of 72 candidates in mouse cortex SgL and IgL.tiff", 
                            sep = "/")
creat.heatmap(cor.exprs.m, rownames(cor.exprs.m), scale = "none", pwd.image.cor.exprs.m)
creat.heatmap(cor.SgL.m, rownames(cor.SgL.m), scale = "none", pwd.image.cor.SgL.m)
creat.heatmap(cor.IgL.m, rownames(cor.IgL.m), scale = "none", pwd.image.cor.IgL.m)
creat.heatmap(cor.MET.m, rownames(cor.MET.m), scale = "none", pwd.image.cor.MET.m)


#--------------------------------------------------------------------------------
#Statistics for mouse data
#--------------------------------------------------------------------------------
#wilcoxon-mann-whitney test - compare MET binding candidates and non-binding candidates
cutoff = 100
dup.m = subset(data.m, rowMeans(data.m[,2:dim(data.m)[2]])>cutoff)
rna.seq.data.highexpr.m = subset(rna.seq.data.m, rowMeans(rna.seq.data.m[,2:dim(rna.seq.data.m)[2]])>cutoff)
row.names(rna.seq.data.highexpr.m) = rna.seq.data.highexpr.m$Geneid
dim(rna.seq.data.highexpr.m)
rna.seq.data.highexpr.m = subset(rna.seq.data.highexpr.m, !rownames(rna.seq.data.highexpr.m) %in% rownames(dup.m))
dim(rna.seq.data.highexpr.m)
dim(rna.seq.data.highexpr.m)[1]*100/dim(rna.seq.data.m)[1]
mix.data.m = rbind(data.m, rna.seq.data.highexpr.m)
dim(mix.data.m)
cor.mix.data.m = cal.corr(mix.data.m, index1 = 2, index2 = dim(mix.data.m)[2], method = "spearman")
cal.pvalue (cor.mix.data.m[,1], dim(data.m)[1])

#create heatmap for mix data
pwd.image.mix.data.m = paste(pwd.image.dir.m, "Heatmap of expression of 72 candidates and random selected genes in mouse cortex.tiff", 
                             sep = "/")
creat.heatmap(mix.data.m[, 2:dim(mix.data.m)[2]], rownames(mix.data.m), scale = "row", pwd.image.mix.data.m)
pwd.image.cor.mix.data.m = paste(pwd.image.dir.m, "Heatmap of expression correlation of 72 candidates and random selected genes in human cortex.tiff", 
                                 sep = "/")
creat.heatmap(cor.mix.data.m, rownames(cor.mix.data.m), scale = "none", pwd.image.cor.mix.data.m)

#wilcoxon-mann-whitney test - compare structures
print ("p value - compare if ranked correlation differently between brain regions")
wilcox.test(cor.MET.m[2:dim(data.m)[1],1], cor.MET.m[2:dim(data.m)[1],2])


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
pwd.data.h = paste(pwd.h, "data/expression data/clean data/gene expression of MET-interacting candidates-S1C.csv", sep = "/")
pwd.rna.seq.data.h = paste(pwd.h, "data/expression data/clean data/expression_matrix_add_Row_and_Column_4mon to 48 mon_STC-ITC-V1C-S1C.csv", sep = "/")
pwd.image.dir.h = paste(pwd.h, "result", sep = "/")
#read human data
data.h = read.csv(pwd.data.h, header = TRUE, sep = ",")
rna.seq.data.h = read.csv(pwd.rna.seq.data.h, header = TRUE, sep = ",")
#check data
dim(data.h)
dim(rna.seq.data.h)

d = t(data.h[,2:dim(data.h)[2]])
colnames(d) = data.h$Symbol
plot(d)

#--------------------------------------------------------------------------------
#Calculate correlation for human data
#--------------------------------------------------------------------------------
#calculate correlation for human data
row.names(data.h) = data.h$Symbol
cor.exprs.h = cal.corr(data.h, index1 = 2, index2 = dim(data.h)[2], method = "spearman")
cor.ITC.h = cal.corr(data.h, index1 = 2, index2 = 6, method = "spearman")
cor.STC.h = cal.corr(data.h, index1 = 7, index2 = 11, method = "spearman")
cor.V1C.h = cal.corr(data.h, index1 = 12, index2 = 16, method = "spearman")
cor.S1C.h = cal.corr(data.h, index1 = 17, index2 = 21, method = "spearman")
cor.MET.h = cbind(cor.ITC.h[, 1], cor.STC.h[, 1], cor.V1C.h[, 1], cor.S1C.h[, 1])

#--------------------------------------------------------------------------------
#Creat heatmap for human data
#--------------------------------------------------------------------------------
data.exprs.h = data.h[,2:dim(data.h)[2]]
pwd.image.exprs.h = paste(pwd.image.dir.h, "Heatmap of gene expression of 72 candidates in human cortex-S1C.tiff", 
                          sep = "/")
creat.heatmap(data.exprs.h, rownames(data.exprs.h), scale = "row", pwd.image.exprs.h)
#correlation heatmap
pwd.image.cor.exprs.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex.tiff", 
                              sep = "/")
pwd.image.cor.ITC.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex ITC-S1C.tiff", 
                            sep = "/")
pwd.image.cor.STC.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex STC-S1C.tiff", 
                            sep = "/")
pwd.image.cor.V1C.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex V1C-S1C.tiff", 
                            sep = "/")
pwd.image.cor.V1C.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates in human cortex S1C-S1C.tiff", 
                            sep = "/")
pwd.image.cor.MET.h = paste(pwd.image.dir.h, "Heatmap of expression correlation to MET of 72 candidates in human cortex-S1C.tiff", 
                            sep = "/")
creat.heatmap(cor.exprs.h, rownames(cor.exprs.h), scale = "none", pwd.image.cor.exprs.h)
creat.heatmap(cor.ITC.h, rownames(cor.ITC.h), scale = "none", pwd.image.cor.ITC.h)
creat.heatmap(cor.STC.h, rownames(cor.STC.h), scale = "none", pwd.image.cor.STC.h)
creat.heatmap(cor.V1C.h, rownames(cor.V1C.h), scale = "none", pwd.image.cor.V1C.h)
creat.heatmap(cor.MET.h, rownames(cor.MET.h), scale = "none", pwd.image.cor.MET.h)

#--------------------------------------------------------------------------------
#Statistics for human data
#--------------------------------------------------------------------------------
#wilcoxon-mann-whitney test - compare MET binding candidates and non-binding candidates
cutoff = 100
dup.h = subset(data.h, rowMeans(data.h[,2:dim(data.h)[2]])>cutoff)
rna.seq.data.highexpr.h = subset(rna.seq.data.h, rowMeans(rna.seq.data.h[,2:dim(rna.seq.data.h)[2]])>cutoff)
row.names(rna.seq.data.highexpr.h) = rna.seq.data.highexpr.h$Symbol
dim(rna.seq.data.highexpr.h)
rna.seq.data.highexpr.h = subset(rna.seq.data.highexpr.h, !rownames(rna.seq.data.highexpr.h) %in% rownames(dup.h))
dim(rna.seq.data.highexpr.h)
dim(rna.seq.data.highexpr.h)[1]*100/dim(rna.seq.data.h)[1]
mix.data.h = rbind(data.h, rna.seq.data.highexpr.h)
dim(mix.data.h)
cor.mix.data.h = cal.corr(mix.data.h, index1 = 2, index2 = dim(mix.data.h)[2], method = "spearman")
#permutation test
count = 0
num = 100
p0 = wilcox.test(cor.mix.data.h[2:73,], cor.mix.data.h[74:dim(cor.mix.data.h),])[3]
print (p0)
for (i in 1:num){
  cor.mix.data.h.df = cor.mix.data.h[sample(nrow(cor.mix.data.h[2:dim(cor.mix.data.h),])),]
  test = wilcox.test(cor.mix.data.h.df[1:72,], cor.mix.data.h.df[73:dim(cor.mix.data.h.df),])
  p = test[3]
  if (p <= 2.391514e-58) {count = count +1}
}
print(count); (count/num);




cal.pvalue (cor.mix.data.h[,1], dim(data.h)[1])

#create heatmap for mix data
pwd.image.mix.data.h = paste(pwd.image.dir.h, "Heatmap of expression of 72 candidates and random selected genes in human cortex-S1C.tiff", 
                            sep = "/")
creat.heatmap(mix.data.h[, 2:dim(mix.data.h)[2]], rownames(mix.data.h), scale = "row", pwd.image.mix.data.h)
pwd.image.cor.mix.data.h = paste(pwd.image.dir.h, "Heatmap of expression correlation of 72 candidates and random selected genes in human cortex-S1C.tiff", 
                             sep = "/")
creat.heatmap(cor.mix.data.h, rownames(cor.mix.data.h), scale = "none", pwd.image.cor.mix.data.h)

#wilcoxon-mann-whitney test - compare structures
print ("p value - compare if ranked correlation differently between ITC and STC")
wilcox.test(cor.MET.h[2:dim(data.h)[1],1], cor.MET.h[2:dim(data.h)[1],2])
print ("p value - compare if ranked correlation differently between ITC and V1C")
wilcox.test(cor.MET.h[2:dim(data.h)[1],1], cor.MET.h[2:dim(data.h)[1],3])
print ("p value - compare if ranked correlation differently between STC and V1C")
wilcox.test(cor.MET.h[2:dim(data.h)[1],2], cor.MET.h[2:dim(data.h)[1],3])

wilcox.test(cor.MET.h[2:dim(data.h)[1],1], cor.MET.h[2:dim(data.h)[1],4])
wilcox.test(cor.MET.h[2:dim(data.h)[1],2], cor.MET.h[2:dim(data.h)[1],4])
wilcox.test(cor.MET.h[2:dim(data.h)[1],3], cor.MET.h[2:dim(data.h)[1],4])

#--------------------------------------------------------------------------------
#Versions of tools
#--------------------------------------------------------------------------------
#session info
devtools::session_info()