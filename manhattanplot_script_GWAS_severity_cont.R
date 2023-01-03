library(dplyr)
library(qqman)
library(ggrepel)
library(QCEWAS)

data_files <- list.files("results_file/severity_cont/")

l=lapply(1:length(data_files), function(i){
  assign(paste0("res_", i),                                  
         read.delim(paste0("results_file/severity_cont/", data_files[i]), 
                    sep=" ", header = T, na.strings = "NA"))
})

p=lapply(1:length(l), function(i){
  df=l[[i]][which(l[[i]]$frequentist_add_pvalue < 0.05),]
  res=df[,c(1:4)]
  
  return(list(df=as.data.frame(df), 
              res=as.data.frame(res)))
})

res=do.call(rbind,lapply(p, "[[",2))
head(res)

CHR= as.numeric(res$chromosome)
BP=res$position
SNP=res$rsid
P=res$frequentist_add_pvalue

gwasResults=data.frame(CHR, BP, SNP, P)
gwasResults=na.omit(gwasResults)

manhattan(gwasResults, annotatePval = 0.00001, col=c("green","black","blue","orange"))

gwasResults_order <- gwasResults[order(P),]
top20=head(gwasResults_order, 20)
write.csv(top20, "severity_cont_top20.csv", row.names=FALSE, quote=FALSE)
region=rep("region", 20)
chr=top20[,1]
startpos=top20[,2]
endpos=top20[,2]
snpnexus=data.frame(region,chr,startpos,endpos)
write.csv(snpnexus, "severity_cont_NexuSNP.csv", row.names=FALSE, quote=FALSE) 

####Q-Qplot + lambda
getwd()

data_files <- list.files("results_file/onlypvalue/")

l=lapply(1:length(data_files), function(i){
  assign(paste0("res_", i),                                  
         read.delim(paste0("results_file/onlypvalue/", data_files[i]), 
                    sep=" ", header = T, na.strings = "NA"))
})

p=lapply(1:length(l), function(i){
  df=l[[i]][which(l[[i]]$frequentist_add_pvalue < 2),]
  res=df[,c(6,7)]
  
  return(list(df=as.data.frame(df), 
              res=as.data.frame(res)))
})

res=do.call(rbind,lapply(p, "[[",2))

P_lambda(res$frequentist_add_pvalue)

p_values <- as.data.frame(res$frequentist_add_pvalue)
colnames(p_values) <- "pvalue"
chisq <- qchisq(1-p_values$pvalue,1)
lambda_value <- median(chisq)/qchisq(0.5,1)
lambda_value <- round(lambda_value, digits = 4)
print(lambda_value)

file_name=c("severity_cont")
qq_filename <- paste("qqplot-", file_name, ".jpeg", sep="")

#qqplot:
jpeg(filename = qq_filename, width = 500, height = 500)
qq(p_values$pvalue, cex.axis = 1.5)
dev.off()







