#Cancer Data Management(for Ovarian.)
setwd("C:/Users/jeff/Desktop/Study/Data/Ovarian cancer/Platform_HG-U133A(GPL96)/")
Ovariandata = read.csv("tcga-ov/E377_13799.csv")

GPL96 = read.csv("gene annotatioin file(for GPL96).csv", header = TRUE)
head(Ovariandata[,1:6])

IDEN = GPL96$IDENTIFIER
AA = as.numeric(as.character(GPL96$Entrez.Gene))
BB = as.numeric(as.character(GPL96$X))
BB1 = as.numeric(as.character(GPL96$X.1))
BB2 = as.numeric(as.character(GPL96$X.2))
BB3 = as.numeric(as.character(GPL96$X.3))
BB4 = as.numeric(as.character(GPL96$X.4))
BB5 = as.numeric(as.character(GPL96$X.5))
BB6 = as.numeric(as.character(GPL96$X.6))
BB7 = as.numeric(as.character(GPL96$X.7))
BB8 = as.numeric(as.character(GPL96$X.8))
BB9 = as.numeric(as.character(GPL96$X.9))
BB10 = as.numeric(as.character(GPL96$X.10))
BB11 = as.numeric(as.character(GPL96$X.11))
BB12 = as.numeric(as.character(GPL96$X.12))
BB13 = as.numeric(as.character(GPL96$X.13))
BB14 = as.numeric(as.character(GPL96$X.14))
BB15 = as.numeric(as.character(GPL96$X.15))
BB16 = as.numeric(as.character(GPL96$X.16))
BB17 = as.numeric(as.character(GPL96$X.17))
BB18 = as.numeric(as.character(GPL96$X.18))
BB19 = as.numeric(as.character(GPL96$X.19))
BB20 = as.numeric(as.character(GPL96$X.20))


AllIDEN = c()
for (i in 1:length(Ovariandata$geneid)) {
  print(i)
  add = as.character(GPL96$IDENTIFIER[which(AA == Ovariandata$geneid[i])])
  if(length(add) == 0){
    add = as.character(GPL96$IDENTIFIER[which(BB == Ovariandata$geneid[i])])
    if(length(add) == 0){
      add = as.character(GPL96$IDENTIFIER[which(BB1 == Ovariandata$geneid[i])])
      if(length(add) == 0){
        add = as.character(GPL96$IDENTIFIER[which(BB2 == Ovariandata$geneid[i])])
        if(length(add) == 0){
          add = as.character(GPL96$IDENTIFIER[which(BB3 == Ovariandata$geneid[i])])
          if(length(add) == 0){
            add = as.character(GPL96$IDENTIFIER[which(BB4 == Ovariandata$geneid[i])])
            if(length(add) == 0){
              add = as.character(GPL96$IDENTIFIER[which(BB5 == Ovariandata$geneid[i])])
              if(length(add) == 0){
                add = as.character(GPL96$IDENTIFIER[which(BB6 == Ovariandata$geneid[i])])
                if(length(add) == 0){
                  add = as.character(GPL96$IDENTIFIER[which(BB7 == Ovariandata$geneid[i])])
                  if(length(add) == 0){
                    add = as.character(GPL96$IDENTIFIER[which(BB8 == Ovariandata$geneid[i])])
                    if(length(add) == 0){
                      add = as.character(GPL96$IDENTIFIER[which(BB9 == Ovariandata$geneid[i])])
                      if(length(add) == 0){
                        add = as.character(GPL96$IDENTIFIER[which(BB10 == Ovariandata$geneid[i])])
                        if(length(add) == 0){
                          add = as.character(GPL96$IDENTIFIER[which(BB11 == Ovariandata$geneid[i])])
                          if(length(add) == 0){
                            add = as.character(GPL96$IDENTIFIER[which(BB12 == Ovariandata$geneid[i])])
                            if(length(add) == 0){
                              add = as.character(GPL96$IDENTIFIER[which(BB13 == Ovariandata$geneid[i])])
                              if(length(add) == 0){
                                add = as.character(GPL96$IDENTIFIER[which(BB14 == Ovariandata$geneid[i])])
                                if(length(add) == 0){
                                  add = as.character(GPL96$IDENTIFIER[which(BB15 == Ovariandata$geneid[i])])
                                  if(length(add) == 0){
                                    add = as.character(GPL96$IDENTIFIER[which(BB16 == Ovariandata$geneid[i])])
                                    if(length(add) == 0){
                                      add = as.character(GPL96$IDENTIFIER[which(BB17 == Ovariandata$geneid[i])])
                                      if(length(add) == 0){
                                        add = as.character(GPL96$IDENTIFIER[which(BB18 == Ovariandata$geneid[i])])
                                        if(length(add) == 0){
                                          add = as.character(GPL96$IDENTIFIER[which(BB19 == Ovariandata$geneid[i])])
                                          if(length(add) == 0){
                                            add = as.character(GPL96$IDENTIFIER[which(BB20 == Ovariandata$geneid[i])])
                                            if(length(add) == 0){
                                              add = NA
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  if(length(add) > 1){
    add = add[1]
  }
  AllIDEN = c(AllIDEN, add)
}

Ovariandata$geneid = AllIDEN

OvarianDataA = Ovariandata

#Genes and Samples Selection(KEGG Pathway with /get/hsa.....)
select = function(DataQuantNorm, gene, samples, name, uniqueprobe = FALSE, switchnormal = FALSE){
  genelist <<- gene$V2
  
  #Gene tiltle name for GPL570
  IDENTIFIER <<- Ovariandata$geneid
  
  geneflag <<- IDENTIFIER %in% genelist
  cat(paste0("\n Genes: ", length(genelist)))
  cat(paste0("\n Probes: ", sum(geneflag)))
  Normal <<- as.character(samples[agrepl("0", as.character(samples$V1), max.distance = 0),2])
  cat(paste0("\n Normal: ", length(Normal)))
  Cancer <<- as.character(samples[agrepl("1", as.character(samples$V1), max.distance = 0),2])
  cat(paste0("\n Cancer: ", length(Cancer)))
  
  name = paste0(name, "_mean_uniqueprobe")
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
      genedata <<- temp[probelist %in% gene, ]
      CV = apply(genedata[,-1], 1, function(x){ sd(x)/mean(x) })
      alldata <<- rbind(alldata , genedata[which.max(CV), ])
      tempflag = rep(FALSE, sum(probelist %in% gene))
      tempflag[which.max(CV)] = TRUE
      uniqueflag <<- c(uniqueflag, tempflag)
    }
  }
  sum(uniqueflag)
  rownames(alldata) = (rownames(temp)[uniqueflag])
  
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

##########################################################################################
setwd("C:/Users/jeff/Desktop/Study/Data/Ovarian cancer/gene set/")
gene_mTOR = read.csv("mTOR signaling pathway.csv", header = FALSE)
gene_Jak_STAT = read.csv("Jak-STAT signaling pathway.csv", header = FALSE)
gene_p53 = read.csv("p53 signaling pathway.csv", header = FALSE)
gene_PI3K_Akt = read.csv("PI3K-Akt Pathway.csv", header = FALSE)

setwd("C:/Users/jeff/Desktop/Study/Data/Ovarian cancer/Platform_HG-U133A(GPL96)/Probes and samples")
samples = read.csv("tcga_ov_preselected.csv", header = FALSE)

#Cancer Normal, Unique mean gene
OvarianSelectedA_mTOR_mean_unique_CN = select(OvarianDataA, gene_mTOR, samples, "Ovarian_GPL96_mTOR")
OvarianSelectedA_Jak_STAT_mean_unique_CN = select(OvarianDataA, gene_Jak_STAT, samples, "Ovarian_GPL96_Jak-STAT")
OvarianSelectedA_p53_mean_unique_CN = select(OvarianDataA, gene_p53, samples, "Ovarian_GPL96_p53")
OvarianSelectedA_PI3K_Akt_mean_unique_CN = select(OvarianDataA, gene_PI3K_Akt, samples, "Ovarian_GPL96_PI3K-Akt")

#Normal Cancer, Unique mean gene
OvarianSelectedA_mTOR_mean_unique_NC = select(OvarianDataA, gene_mTOR, samples, "Ovarian_GPL96_mTOR", switchnormal = TRUE)
OvarianSelectedA_Jak_STAT_mean_unique_NC = select(OvarianDataA, gene_Jak_STAT, samples, "Ovarian_GPL96_Jak-STAT", switchnormal = TRUE)
OvarianSelectedA_p53_mean_unique_NC = select(OvarianDataA, gene_p53, samples, "Ovarian_GPL96_p53", switchnormal = TRUE)
OvarianSelectedA_PI3K_Akt_mean_unique_NC = select(OvarianDataA, gene_PI3K_Akt, samples, "Ovarian_GPL96_PI3K-Akt", switchnormal = TRUE)


#Ovarian cancer data management finish, but not include orther dataset.
#Remark: Other Cancer Data is at Data Management.R.
#=============================================================================================================================================================================================================================

