library(dplyr)
library(qqman)
library(ggrepel)

setwd("D:/Git/Rstudio/Admixture")

data_files <- list.files("/OneDrive/Documentos/Latinos/snptest_files/results_file/severity_cont/")

l=lapply(1:length(data_files), function(i){
  assign(paste0("res_", i),                                  
         read.delim(paste0("/OneDrive/Documentos/Latinos/snptest_files/results_file/severity_cont/", data_files[i]), 
                    sep="\t", header = T, na.strings = "NA", skip = 11))
})

p=lapply(1:length(l), function(i){
  df=l[[i]][which(l[[i]]$frequentist_add_pvalue < 0.01),]
  res=df[,c(2:4,42,44:45)]
  
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

snps=c("1-9067157","1-64948270","1-77501822","1-155066988","1-155175305","1-155197995","1-155278322","2-60480453","3-45796521","3-45818159","3-45873093","3-101737122","3-146517122","4-25446871","4-167824478","5-132441275","6-31153455","6-31513129","6-32702531","6-41522644","7-75623396","7-100032719","8-60532539","9-21206606","9-33425186","9-133271182","10-79946568","11-1219991","11-34482745","12-112919637","12-132481571","13-112881427","16-89196249","17-40003082","17-46085231","17-49863303","19-4717660","19-10352442","19-10414696","19-48702888","19-50374423","21-33229937","21-33237639","21-33287378","21-33949755","21-41479527","23-15523993")

manhattan(gwasResults, annotatePval = 0.00001)




