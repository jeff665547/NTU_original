#Cancer Data Management(for Breast, Colorectal, Leukemia, Lung)
#RMA
library("affy")
#Quantile Normalization
#source("http://bioconductor.org/biocLite.R")
library("preprocessCore")

#Normalization function: Nor
Normalization = function(x){
  require(affy)
  require(preprocessCore)
  
  cat("Read Raw Data(ReadAffy(affy))\n")
  Data <<- ReadAffy()
  cat("Robust Multi-array Average(rma(affy))\n")
  DataRMA <<- rma(Data, normalize = TRUE)
  write.exprs(DataRMA, file = (paste(as.character(x), "txt", sep = ".")))
  
  Data = as.matrix(read.table((paste(as.character(x), "txt", sep = ".")), sep = "\t", header = T, row.names = 1))
  cat("Quantile Normalization(normalize.quantiles(preprocessCore))\n")
  matrixData = normalize.quantiles(Data)
  rownames(matrixData) = rownames(Data)
  colnames(matrixData) = gsub(".CEL.gz","",colnames(Data))
  DataQuantNorm <<- matrixData
  write.csv(DataQuantNorm, file = (paste(as.character(x), "_processed.csv", sep = "")), quote = F, row.names = T)
  
  cat("\n")
  cat("Probes: ", nrow(DataQuantNorm))
  cat("\nSamples: ", ncol(DataQuantNorm))
  cat("\nData: \n")
  print(head(DataQuantNorm[,1:3]))
  return(DataQuantNorm)
}
#############################################################################################################################################################################################################################
#Genes and Samples Selection(KEGG Pathway with /get/hsa.....)
select = function(DataQuantNorm, gene, samples, name, uniqueprobe = FALSE, switchnormal = FALSE){
  genelist <<- gene$V2
  #Gene tiltle name for GPL570
  GPL570 = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133_Plus_2(GPL570)/Candidate/small-cell lung cancer/GDS4794(part of GSE43346)/GDS4794array.csv")
  GPL570_gene_title = as.character(GPL570$IDENTIFIER)
  #Gene tiltle name for GPL96
  GPL96 = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133A(GPL96)/Candidate/Smoking/GDS3257(GSE10072)Cigarette smoking effect on lung adenocarcinoma/GSE10072array.csv")
  GPL96_gene_title = as.character(GPL96$IDENTIFIER)
  
  if(dim(DataQuantNorm)[1] == length(GPL570_gene_title)){
    IDENTIFIER <<- GPL570_gene_title
    cat("GPL570\n")
  }else if(dim(DataQuantNorm)[1] == length(GPL96_gene_title)){
    IDENTIFIER <<- GPL96_gene_title
    cat("GPL96\n")
  }
  
  geneflag <<- IDENTIFIER %in% genelist
  cat(paste0("\n Genes: ", length(genelist)))
  cat(paste0("\n Probes: ", sum(geneflag)))
  Normal <<- as.character(samples[agrepl("normal", as.character(samples$V1)),2])
  cat(paste0("\n Normal: ", length(Normal)))
  Cancer <<- as.character(samples[agrepl("tumor", as.character(samples$V1)),2])
  cat(paste0("\n Cancer: ", length(Cancer)))
  
  if(uniqueprobe){
    name = paste0(name, "_uniqueprobe")
    temp <<- DataQuantNorm[geneflag, ]
    probelist <<- IDENTIFIER[geneflag] ##實際拿到的基因probe
    cat(paste0("\n Unique Gene Probes: ", length(unique(probelist))))
    n = 0
    uniqueflag <<- c()
    alldata <<- c()
    for(i in 1: length(unique(probelist))){
      gene = unique(probelist)[i]
      if(sum(probelist %in% gene) == 1){
        n = n + 1
        genedata = temp[probelist %in% gene, ]
        alldata <<- rbind(alldata, genedata)
        uniqueflag <<- c(uniqueflag, TRUE)
      }else{
        n = n + sum(probelist %in% gene)
        genedata = temp[probelist %in% gene, ]
        CV = apply(genedata, 1, function(x){ sd(x)/mean(x) })
        alldata <<- rbind(alldata , genedata[which.max(CV), ])
        tempflag = rep(FALSE, sum(probelist %in% gene))
        tempflag[which.max(CV)] = TRUE
        uniqueflag <<- c(uniqueflag, tempflag)
      }
    }
    sum(uniqueflag)
    rownames(alldata) = (rownames(temp)[uniqueflag])
  }else{
    alldata <<- DataQuantNorm[geneflag, ]
  }
  
  DataQuantNorm_Normal <<- alldata[ , colnames(DataQuantNorm) %in% Normal]
  DataQuantNorm_Cancer <<- alldata[ , colnames(DataQuantNorm) %in% Cancer]

  setwd("C:/Users/jeff/Desktop/Study/Data/Final")
  if(switchnormal){
    name = paste0(name, "_NC")
    final = cbind(DataQuantNorm_Normal, DataQuantNorm_Cancer)   #先Normal後Cancer
  }else{
    name = paste0(name, "_CN")
    final = cbind(DataQuantNorm_Cancer, DataQuantNorm_Normal)   #先Cancer後Normal
  }
  write.csv(final, file = (paste(as.character(name), ".csv", sep = "")), quote = F)
  cat("\nFile saved as:")
  cat(paste(as.character(name), ".csv", sep = ""))
  return(final)
}
##############################################################################################################################################################################################################################

#Breast Cancer Data(GPL96: GSE15852)
setwd("C:/Users/jeff/Desktop/Study/Data/Breast cancer/Platform_HG-U133A(GPL96)/GSE15852/raw")
BreastDataA = Normalization("GSE15852")

gene_PI3K_Akt = read.csv("C:/Users/jeff/Desktop/Study/Data/Breast cancer/gene set/PI3K-Akt Pathway.csv", header = FALSE)
gene_Wnt = read.csv("C:/Users/jeff/Desktop/Study/Data/Breast cancer/gene set/Wnt signaling pathway.csv", header = FALSE)

samples = read.csv("C:/Users/jeff/Desktop/Study/Data/Breast cancer/Platform_HG-U133A(GPL96)/Probes and samples/GSE15852_selected.csv", header = FALSE)

#Cancer Normal
BreastSelectedA_PI3K_Akt_CN = select(BreastDataA, gene_PI3K_Akt, samples, "Breast_GPL96_PI3K-Akt")
BreastSelectedA_Wnt_CN = select(BreastDataA, gene_Wnt, samples, "Breast_GPL96_Wnt")
#unique probe
BreastSelectedA_PI3K_Akt_unique_CN = select(BreastDataA, gene_PI3K_Akt, samples, "Breast_GPL96_PI3K-Akt", uniqueprobe = TRUE)
BreastSelectedA_Wnt_unique_CN = select(BreastDataA, gene_Wnt, samples, "Breast_GPL96_Wnt", uniqueprobe = TRUE)

#Normal Cancer
BreastSelectedA_PI3K_Akt_NC = select(BreastDataA, gene_PI3K_Akt, samples, "Breast_GPL96_PI3K-Akt", switchnormal = TRUE)
BreastSelectedA_Wnt_NC = select(BreastDataA, gene_Wnt, samples, "Breast_GPL96_Wnt", switchnormal = TRUE)
#unique probe
BreastSelectedA_PI3K_Akt_unique_NC = select(BreastDataA, gene_PI3K_Akt, samples, "Breast_GPL96_PI3K-Akt", uniqueprobe = TRUE, switchnormal = TRUE)
BreastSelectedA_Wnt_unique_NC = select(BreastDataA, gene_Wnt, samples, "Breast_GPL96_Wnt", uniqueprobe = TRUE, switchnormal = TRUE)



##Breast Cancer Data(GPL570: GSE10810)
#setwd("C:/Users/jeff/Desktop/Study/Data/Breast cancer/Platform_HG-U133_Plus_2(GPL570)/GSE10810/raw")
#BreastDataB = Normalization("GSE10810")


#Colorectal Cancer Data(GPL570: GSE8671)
setwd("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/Platform_HG-U133_Plus_2(GPL570)/GSE8671/raw")
ColorectalDataB = Normalization("GSE8671")

gene_p53 = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/gene set/p53 signaling pathway.csv", header = FALSE)
gene_PI3K_Akt = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/gene set/PI3K-Akt Pathway.csv", header = FALSE)
gene_Ras = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/gene set/Ras signaling pathway.csv", header = FALSE)
gene_TGF_beta = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/gene set/TGF-beta signaling pathway.csv", header = FALSE)
gene_Wnt = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/gene set/Wnt signaling pathway.csv", header = FALSE)

samples = read.csv("C:/Users/jeff/Desktop/Study/Data/Colorectal cancer/Platform_HG-U133_Plus_2(GPL570)/Probes and samples/GSE8671_selected.csv", header = FALSE)

#Cancer Normal
ColorectalSelectedB_p53_CN = select(ColorectalDataB, gene_p53, samples, "Colorectal_GPL570_p53")
ColorectalSelectedB_PI3K_Akt_CN = select(ColorectalDataB, gene_PI3K_Akt, samples, "Colorectal_GPL570_PI3K-Akt")
ColorectalSelectedB_Ras_CN = select(ColorectalDataB, gene_Ras, samples, "Colorectal_GPL570_Ras")
ColorectalSelectedB_TGF_beta_CN = select(ColorectalDataB, gene_TGF_beta, samples, "Colorectal_GPL570_TGF-beta")
ColorectalSelectedB_Wnt_CN = select(ColorectalDataB, gene_Wnt, samples, "Colorectal_GPL570_Wnt")
#unique probe
ColorectalSelectedB_p53_unique_CN = select(ColorectalDataB, gene_p53, samples, "Colorectal_GPL570_p53", uniqueprobe = TRUE)
ColorectalSelectedB_PI3K_Akt_unique_CN = select(ColorectalDataB, gene_PI3K_Akt, samples, "Colorectal_GPL570_PI3K-Akt", uniqueprobe = TRUE)
ColorectalSelectedB_Ras_unique_CN = select(ColorectalDataB, gene_Ras, samples, "Colorectal_GPL570_Ras", uniqueprobe = TRUE)
ColorectalSelectedB_TGF_beta_unique_CN = select(ColorectalDataB, gene_TGF_beta, samples, "Colorectal_GPL570_TGF-beta", uniqueprobe = TRUE)
ColorectalSelectedB_Wnt_unique_CN = select(ColorectalDataB, gene_Wnt, samples, "Colorectal_GPL570_Wnt", uniqueprobe = TRUE)


#Normal Cancer
ColorectalSelectedB_p53_NC = select(ColorectalDataB, gene_p53, samples, "Colorectal_GPL570_p53", switchnormal = TRUE)
ColorectalSelectedB_PI3K_Akt_NC = select(ColorectalDataB, gene_PI3K_Akt, samples, "Colorectal_GPL570_PI3K-Akt", switchnormal = TRUE)
ColorectalSelectedB_Ras_NC = select(ColorectalDataB, gene_Ras, samples, "Colorectal_GPL570_Ras", switchnormal = TRUE)
ColorectalSelectedB_TGF_beta_NC = select(ColorectalDataB, gene_TGF_beta, samples, "Colorectal_GPL570_TGF-beta", switchnormal = TRUE)
ColorectalSelectedB_Wnt_NC = select(ColorectalDataB, gene_Wnt, samples, "Colorectal_GPL570_Wnt", switchnormal = TRUE)
#unique probe
ColorectalSelectedB_p53_unique_NC = select(ColorectalDataB, gene_p53, samples, "Colorectal_GPL570_p53", uniqueprobe = TRUE, switchnormal = TRUE)
ColorectalSelectedB_PI3K_Akt_unique_NC = select(ColorectalDataB, gene_PI3K_Akt, samples, "Colorectal_GPL570_PI3K-Akt", uniqueprobe = TRUE, switchnormal = TRUE)
ColorectalSelectedB_Ras_unique_NC = select(ColorectalDataB, gene_Ras, samples, "Colorectal_GPL570_Ras", uniqueprobe = TRUE, switchnormal = TRUE)
ColorectalSelectedB_TGF_beta_unique_NC = select(ColorectalDataB, gene_TGF_beta, samples, "Colorectal_GPL570_TGF-beta", uniqueprobe = TRUE, switchnormal = TRUE)
ColorectalSelectedB_Wnt_unique_NC = select(ColorectalDataB, gene_Wnt, samples, "Colorectal_GPL570_Wnt", uniqueprobe = TRUE, switchnormal = TRUE)


#Leukemia(GPL96: GSE9476)
setwd("C:/Users/jeff/Desktop/Study/Data/Leukemia/Platform_HG-U133A(GPL96)/GSE9476/raw")
LeukemiaDataA = Normalization("GSE9476")

gene_Jak_STAT = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/gene set/Jak-STAT signaling pathway.csv", header = FALSE)
gene_mTOR = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/gene set/mTOR signaling pathway.csv", header = FALSE)
gene_PI3K_Akt = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/gene set/PI3K-Akt Pathway.csv", header = FALSE)
gene_Wnt = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/gene set/Wnt signaling pathway.csv", header = FALSE)

samples = read.csv("C:/Users/jeff/Desktop/Study/Data/Leukemia/Platform_HG-U133A(GPL96)/Probes and samples/GSE9476_selected.csv", header = FALSE)

#Cancer Normal
LeukemiaSelectedA_Jak_STAT_CN = select(LeukemiaDataA, gene_Jak_STAT, samples, "Leukemia_GPL96_Jak-STAT")
LeukemiaSelectedA_mTOR_CN = select(LeukemiaDataA, gene_mTOR, samples, "Leukemia_GPL96_mTOR")
LeukemiaSelectedA_PI3K_Akt_CN = select(LeukemiaDataA, gene_PI3K_Akt, samples, "Leukemia_GPL96_PI3K-Akt")
LeukemiaSelectedA_Wnt_CN = select(LeukemiaDataA, gene_Wnt, samples, "Leukemia_GPL96_Wnt")
#unique probe
LeukemiaSelectedA_Jak_STAT_unique_CN = select(LeukemiaDataA, gene_Jak_STAT, samples, "Leukemia_GPL96_Jak-STAT", uniqueprobe = TRUE)
LeukemiaSelectedA_mTOR_unique_CN = select(LeukemiaDataA, gene_mTOR, samples, "Leukemia_GPL96_mTOR", uniqueprobe = TRUE)
LeukemiaSelectedA_PI3K_Akt_unique_CN = select(LeukemiaDataA, gene_PI3K_Akt, samples, "Leukemia_GPL96_PI3K-Akt", uniqueprobe = TRUE)
LeukemiaSelectedA_Wnt_unique_CN = select(LeukemiaDataA, gene_Wnt, samples, "Leukemia_GPL96_Wnt", uniqueprobe = TRUE)

#Normal Cancer
LeukemiaSelectedA_Jak_STAT_NC = select(LeukemiaDataA, gene_Jak_STAT, samples, "Leukemia_GPL96_Jak-STAT", switchnormal = TRUE)
LeukemiaSelectedA_mTOR_NC = select(LeukemiaDataA, gene_mTOR, samples, "Leukemia_GPL96_mTOR", switchnormal = TRUE)
LeukemiaSelectedA_PI3K_Akt_NC = select(LeukemiaDataA, gene_PI3K_Akt, samples, "Leukemia_GPL96_PI3K-Akt", switchnormal = TRUE)
LeukemiaSelectedA_Wnt_NC = select(LeukemiaDataA, gene_Wnt, samples, "Leukemia_GPL96_Wnt", switchnormal = TRUE)
#unique probe
LeukemiaSelectedA_Jak_STAT_unique_NC = select(LeukemiaDataA, gene_Jak_STAT, samples, "Leukemia_GPL96_Jak-STAT", uniqueprobe = TRUE, switchnormal = TRUE)
LeukemiaSelectedA_mTOR_unique_NC = select(LeukemiaDataA, gene_mTOR, samples, "Leukemia_GPL96_mTOR", uniqueprobe = TRUE, switchnormal = TRUE)
LeukemiaSelectedA_PI3K_Akt_unique_NC = select(LeukemiaDataA, gene_PI3K_Akt, samples, "Leukemia_GPL96_PI3K-Akt", uniqueprobe = TRUE, switchnormal = TRUE)
LeukemiaSelectedA_Wnt_unique_NC = select(LeukemiaDataA, gene_Wnt, samples, "Leukemia_GPL96_Wnt", uniqueprobe = TRUE, switchnormal = TRUE)


#Lung Cancer Data(GPL96: GSE7670)
setwd("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133A(GPL96)/GSE7670/raw")
LungDataA = Normalization("GSE7670")

gene_Cell_cycle = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Cell cycle pathway.csv", header = FALSE)
gene_EGFR = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/EGFR tyrosine kinase inhibitor resistance.csv", header = FALSE)
gene_Jak_STAT = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Jak-STAT signaling pathway.csv", header = FALSE)
gene_mapk = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/mapk signaling pathway.csv", header = FALSE)
gene_mTOR = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/mTOR signaling pathway.csv", header = FALSE)
gene_p53 = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/p53 signaling pathway.csv", header = FALSE)
gene_Wnt = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Wnt signaling pathway.csv", header = FALSE)

samples = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133A(GPL96)/Probes and samples/GSE7670_selected.csv", header = FALSE)

#Cancer Normal
LungSelectedA_Cell_cycle_CN = select(LungDataA, gene_Cell_cycle, samples, "Lung_GPL96_Cell cycle")
LungSelectedA_EGFR_CN = select(LungDataA, gene_EGFR, samples, "Lung_GPL96_EGFR")
LungSelectedA_Jak_STAT_CN = select(LungDataA, gene_Jak_STAT, samples, "Lung_GPL96_Jak-STAT")
LungSelectedA_mapk_CN = select(LungDataA, gene_mapk, samples, "Lung_GPL96_mapk")
LungSelectedA_mTOR_CN = select(LungDataA, gene_mTOR, samples, "Lung_GPL96_mTOR")
LungSelectedA_p53_CN = select(LungDataA, gene_p53, samples, "Lung_GPL96_p53")
LungSelectedA_Wnt_CN = select(LungDataA, gene_Wnt, samples, "Lung_GPL96_Wnt")
#unique probe
LungSelectedA_Cell_cycle_unique_CN = select(LungDataA, gene_Cell_cycle, samples, "Lung_GPL96_Cell cycle", uniqueprobe = TRUE)
LungSelectedA_EGFR_unique_CN = select(LungDataA, gene_EGFR, samples, "Lung_GPL96_EGFR", uniqueprobe = TRUE)
LungSelectedA_Jak_STAT_unique_CN = select(LungDataA, gene_Jak_STAT, samples, "Lung_GPL96_Jak-STAT", uniqueprobe = TRUE)
LungSelectedA_mapk_unique_CN = select(LungDataA, gene_mapk, samples, "Lung_GPL96_mapk", uniqueprobe = TRUE)
LungSelectedA_mTOR_unique_CN = select(LungDataA, gene_mTOR, samples, "Lung_GPL96_mTOR", uniqueprobe = TRUE)
LungSelectedA_p53_unique_CN = select(LungDataA, gene_p53, samples, "Lung_GPL96_p53", uniqueprobe = TRUE)
LungSelectedA_Wnt_unique_CN = select(LungDataA, gene_Wnt, samples, "Lung_GPL96_Wnt", uniqueprobe = TRUE)


#Normal Cancer
LungSelectedA_Cell_cycle_NC = select(LungDataA, gene_Cell_cycle, samples, "Lung_GPL96_Cell cycle", switchnormal = TRUE)
LungSelectedA_EGFR_NC = select(LungDataA, gene_EGFR, samples, "Lung_GPL96_EGFR", switchnormal = TRUE)
LungSelectedA_Jak_STAT_NC = select(LungDataA, gene_Jak_STAT, samples, "Lung_GPL96_Jak-STAT", switchnormal = TRUE)
LungSelectedA_mapk_NC = select(LungDataA, gene_mapk, samples, "Lung_GPL96_mapk", switchnormal = TRUE)
LungSelectedA_mTOR_NC = select(LungDataA, gene_mTOR, samples, "Lung_GPL96_mTOR", switchnormal = TRUE)
LungSelectedA_p53_NC = select(LungDataA, gene_p53, samples, "Lung_GPL96_p53", switchnormal = TRUE)
LungSelectedA_Wnt_NC = select(LungDataA, gene_Wnt, samples, "Lung_GPL96_Wnt", switchnormal = TRUE)
#unique probe
LungSelectedA_Cell_cycle_unique_NC = select(LungDataA, gene_Cell_cycle, samples, "Lung_GPL96_Cell cycle", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_EGFR_unique_NC = select(LungDataA, gene_EGFR, samples, "Lung_GPL96_EGFR", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_Jak_STAT_unique_NC = select(LungDataA, gene_Jak_STAT, samples, "Lung_GPL96_Jak-STAT", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_mapk_unique_NC = select(LungDataA, gene_mapk, samples, "Lung_GPL96_mapk", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_mTOR_unique_NC = select(LungDataA, gene_mTOR, samples, "Lung_GPL96_mTOR", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_p53_unique_NC = select(LungDataA, gene_p53, samples, "Lung_GPL96_p53", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedA_Wnt_unique_NC = select(LungDataA, gene_Wnt, samples, "Lung_GPL96_Wnt", uniqueprobe = TRUE, switchnormal = TRUE)

#Lung Cancer Data(GPL570: GSE19804)
setwd("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133_Plus_2(GPL570)/GDS3837 (GSE19804)/raw")
LungDataB = Normalization("GSE19804")

gene_Cell_cycle = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Cell cycle pathway.csv", header = FALSE)
gene_EGFR = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/EGFR tyrosine kinase inhibitor resistance.csv", header = FALSE)
gene_Jak_STAT = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Jak-STAT signaling pathway.csv", header = FALSE)
gene_mapk = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/mapk signaling pathway.csv", header = FALSE)
gene_mTOR = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/mTOR signaling pathway.csv", header = FALSE)
gene_p53 = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/p53 signaling pathway.csv", header = FALSE)
gene_Wnt = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/gene set/Wnt signaling pathway.csv", header = FALSE)

samples = read.csv("C:/Users/jeff/Desktop/Study/Data/Lung cancer/Platform_HG-U133_Plus_2(GPL570)/Probes and samples/GSE19804_selected.csv", header = FALSE)

#Cancer Normal
LungSelectedB_Cell_cycle_CN = select(LungDataB, gene_Cell_cycle, samples, "Lung_GPL570_Cell cycle")
LungSelectedB_EGFR_CN = select(LungDataB, gene_EGFR, samples, "Lung_GPL570_EGFR")
LungSelectedB_Jak_STAT_CN = select(LungDataB, gene_Jak_STAT, samples, "Lung_GPL570_Jak-STAT")
LungSelectedB_mapk_CN = select(LungDataB, gene_mapk, samples, "Lung_GPL570_mapk")
LungSelectedB_mTOR_CN = select(LungDataB, gene_mTOR, samples, "Lung_GPL570_mTOR")
LungSelectedB_p53_CN = select(LungDataB, gene_p53, samples, "Lung_GPL570_p53")
LungSelectedB_Wnt_CN = select(LungDataB, gene_Wnt, samples, "Lung_GPL570_Wnt")
#unique probe
LungSelectedB_Cell_cycle_unique_CN = select(LungDataB, gene_Cell_cycle, samples, "Lung_GPL570_Cell cycle", uniqueprobe = TRUE)
LungSelectedB_EGFR_unique_CN = select(LungDataB, gene_EGFR, samples, "Lung_GPL570_EGFR", uniqueprobe = TRUE)
LungSelectedB_Jak_STAT_unique_CN = select(LungDataB, gene_Jak_STAT, samples, "Lung_GPL570_Jak-STAT", uniqueprobe = TRUE)
LungSelectedB_mapk_unique_CN = select(LungDataB, gene_mapk, samples, "Lung_GPL570_mapk", uniqueprobe = TRUE)
LungSelectedB_mTOR_unique_CN = select(LungDataB, gene_mTOR, samples, "Lung_GPL570_mTOR", uniqueprobe = TRUE)
LungSelectedB_p53_unique_CN = select(LungDataB, gene_p53, samples, "Lung_GPL570_p53", uniqueprobe = TRUE)
LungSelectedB_Wnt_unique_CN = select(LungDataB, gene_Wnt, samples, "Lung_GPL570_Wnt", uniqueprobe = TRUE)


#Normal Cancer
LungSelectedB_Cell_cycle_NC = select(LungDataB, gene_Cell_cycle, samples, "Lung_GPL570_Cell cycle", switchnormal = TRUE)
LungSelectedB_EGFR_NC = select(LungDataB, gene_EGFR, samples, "Lung_GPL570_EGFR", switchnormal = TRUE)
LungSelectedB_Jak_STAT_NC = select(LungDataB, gene_Jak_STAT, samples, "Lung_GPL570_Jak-STAT", switchnormal = TRUE)
LungSelectedB_mapk_NC = select(LungDataB, gene_mapk, samples, "Lung_GPL570_mapk", switchnormal = TRUE)
LungSelectedB_mTOR_NC = select(LungDataB, gene_mTOR, samples, "Lung_GPL570_mTOR", switchnormal = TRUE)
LungSelectedB_p53_NC = select(LungDataB, gene_p53, samples, "Lung_GPL570_p53", switchnormal = TRUE)
LungSelectedB_Wnt_NC = select(LungDataB, gene_Wnt, samples, "Lung_GPL570_Wnt", switchnormal = TRUE)
#unique probe
LungSelectedB_Cell_cycle_unique_NC = select(LungDataB, gene_Cell_cycle, samples, "Lung_GPL570_Cell cycle", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_EGFR_unique_NC = select(LungDataB, gene_EGFR, samples, "Lung_GPL570_EGFR", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_Jak_STAT_unique_NC = select(LungDataB, gene_Jak_STAT, samples, "Lung_GPL570_Jak-STAT", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_mapk_unique_NC = select(LungDataB, gene_mapk, samples, "Lung_GPL570_mapk", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_mTOR_unique_NC = select(LungDataB, gene_mTOR, samples, "Lung_GPL570_mTOR", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_p53_unique_NC = select(LungDataB, gene_p53, samples, "Lung_GPL570_p53", uniqueprobe = TRUE, switchnormal = TRUE)
LungSelectedB_Wnt_unique_NC = select(LungDataB, gene_Wnt, samples, "Lung_GPL570_Wnt", uniqueprobe = TRUE, switchnormal = TRUE)




#All cancer data finish(for Breast, Colorectal, Leukemia, Lung), but not include candidate dataset.
#Remark: Ovarian cancer data information is not available at this file,
#Instaed, it has specific file for it(Data Management for tcga-ov.R).
#=============================================================================================================================================================================================================================
