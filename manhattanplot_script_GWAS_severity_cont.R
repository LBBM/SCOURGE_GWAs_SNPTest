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

qq(gwasResults$P)
P_lambda(gwasResults$P)

gwasResults_order <- gwasResults[order(P),]
top20=head(gwasResults_order, 20)
write.csv(top20, "severity_cont_top20.csv", row.names=FALSE, quote=FALSE)
region=rep("region", 20)
chr=top20[,1]
startpos=top20[,2]
endpos=top20[,2]
snpnexus=data.frame(region,chr,startpos,endpos)
write.csv(snpnexus, "severity_cont_NexuSNP.csv", row.names=FALSE, quote=FALSE) 
