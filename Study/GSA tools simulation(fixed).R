#Hotelling's T^2
#install.packages("Hotelling")
#shrinkage: Shaefer and Strimmer’s James-Stein shrinkage estimator (An efficient way to compute the correlation matrix, even when n < p is also suitable.)
#progBar: progress bar
#Usage: row: sample, column: variable, x: data1, y: data2
library(Hotelling)

#The global test
#source("https://bioconductor.org/biocLite.R")
#biocLite("globaltest")
#biocLite("golubEsets")
library(globaltest)
#response: Y(e.g. tumor vs. normal)
#alternative: design matrix(gene expression data with column: variable, row: sample)

#The N-statistic(Energy Test)
source("https://jeff665547.github.io/Study/Cramer%20test.R")

#Multivariate Kolmogorov-Smirnov test
#biocLite("GSAR")
library(GSAR)
#object: gene expression value.(matrix form, rows: genes(variables), columns: samples)
#group: numeric vector(e.g. normal(1) vs. tumor(2))
#nperm: permutation times

#GSEA-Category package(GSEA like methods)
#source("https://bioconductor.org/biocLite.R")
#biocLite("Category")
#biocLite("Biobase")
library(Category)
library(Biobase)
#create ExpressionSet object
#Assay data: gene expression data with row: gene name, column: sample
#Format: Matrix
# ExpressionSet(assayData=exprs)

#Phenotypic data: information about the samples (e.g., sex, age, and treatment status; referred to as ‘covariates’)
#row: sample, column: Phenotype(covariate)
#The number of rows of phenotypic data should match the number of columns of expression data, 
#and the rowname of phenotypic data and the column name of Assay data should be the same.
#Format: A data.frame with factor and numeric values (Catgorical should be transformed into factor form.)

#GSEA(Subramanian et al 2005) 
source("https://jeff665547.github.io/Study/GSEA.1.0.R")
#input.ds can be a data.frame in R object
#row: gene
#column: sample
#the first column: the first sample gene expression value
#the first row: first gene expression value
#columnname: sample name, rowname: gene name
#
#input.cls can be a list in R object
#the first element name is "phen" which records two phenotypes with two character strings
#the second element name is "class.v" which records values(0, 1) correspond to the two different phenotypes respectively
#0: the first phenotype, 1: the second phenotype
#both two elements in the list are vector.
#e.g.: CLS = list(phen = c("A", "B"), class.v = c(rep(c(0,1), c(15, 20))))
#
#gs.db must be a gmt format
#just like the simulation geneset.gmt file in study folder.
#
#gseadata = as.data.frame(t(obs))
#response = list(phen = c("A", "B"), class.v = c(rep(c(0,1), c(m/2, m/2))))
#resultGSEA = GSEA(input.ds = gseadata, 
#                  input.cls = response, 
#                  gs.db = "C:/Users/jeff/Desktop/Study/simualtion geneset.gmt", 
#                  use.fast.enrichment.routine = FALSE,
#                  nperm = 5000)

#globalANCOVA
#source("https://bioconductor.org/biocLite.R")
#biocLite("GlobalAncova")
library(GlobalAncova)


#SAFE
#biocLite("safe")
library(safe)
#signset = matrix(1, nrow = 1, ncol = probenumber)
#rownames(signset) = "Geneset"
#results <- safe(t(obs), rep(c(0, 1), c(m/2, m/2)), C.mat = t(signset), 
#                Pi.mat = 1000, print.it = FALSE)
#safe.toptable(results, number = 1)




#########################################################################################
#Power analysis
library(Hotelling)
library(globaltest)
library(Category)
library(GlobalAncova)
library(MVN)
library(mvtnorm)
library(energy)

#Random sample from multivariate normal distribution(MVN) with customized parameters settings
mvnormal = function(n, probenumber, rho, MEAN = rep(0, probenumber), rhostructure = "CS"){
  #n: Sample size. (樣本數)
  #probenumber: The number of dimensions(genes). (維度，基因數目)
  #rho: Correlation between any two dimensions(genes). (任意兩個維度，基因之間的相關係數)
  #rhostructure: Structure of correlation matrix. (關係係數矩陣的結構樣式), e.g.:"CS": Compound Symmertry, "UN": Unstructured
  
  #Program start:
  #Create correlation matrix
  if(rhostructure == "CS"){
    SIGMA = matrix(rho, probenumber, probenumber)
    diag(SIGMA) = 1
  }else if(rhostructure == "UN"){
    rhovector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    SIGMA = matrix(c(sample(rhovector, size = probenumber*probenumber, replace = TRUE)), 
                   probenumber, probenumber)
    diag(SIGMA) = 1
  }
  
  #Generate observation from MVN from package "mvtnorm".
  obs = rmvnorm(n, mean = MEAN, sigma = SIGMA)
  return(obs)
}


#Random sample from multivariate t distribution(MVT) with customized parameters settings
mvt = function(n, probenumber, rho, rhostructure = "CS", DF){
  #n: Sample size. (樣本數)
  #probenumber: The number of dimensions(genes). (維度，基因數目)
  #rho: Correlation between any two dimensions(genes). (任意兩個維度，基因之間的相關係數)
  #rhostructure: Structure of correlation matrix. (關係係數矩陣的結構樣式), e.g.:"CS": Compound Symmertry, "UN": Unstructured
  #DF: Degrees of freedom of multivariate t distribution, DF > 2. (多為t分配的自由度,DF 需要大於 2)
  
  #Program start:
  #Create correlation matrtix
  if(rhostructure == "CS"){
    rhomatrix = matrix(rho, probenumber, probenumber)
    diag(rhomatrix) = 1
  }else if(rhostructure == "UN"){
    rhovector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    rhomatrix = matrix(c(sample(rhovector, size = probenumber*probenumber, replace = TRUE)), 
                       probenumber, probenumber)
    diag(rhomatrix) = 1
  }
  
  #Transform correlation matrix into scale matrix. (將相關係數矩陣轉換成多維度t分配中的scale matrix)
  if(DF > 2){
    SIGMA <<- rhomatrix*((DF - 2)/DF) 
  }else{
    cat("error: df must > 2 \n")
  }
  
  #Generate observation from MVT from package "mvtnorm".
  obs = rmvt(n, sigma = SIGMA, df = DF)
  return(obs)
}

#Use the concept of probability integral transformation(PIT) to
#create various scenarios: Mixture of MVN, MVT(df=3), 
#                          MVN(high correlation: 0.9, low correlation: 0.1)
generateobs = function(n, probenumber, alternative = c("normal", "normal"), alternumber = 2, 
                       proportion = rep(1/alternumber, alternumber), MEAN = list(rep(0, probenumber),rep(0, probenumber)),
                       rho, rhostructure = rep("CS", alternumber), DF = Inf){
  #n: Sample size of the (mixture) distribution. 
  #   (從(混合)分配生出來的樣本數)
  #probenumber: The number of dimensions(genes) for the (mixture) distribution. 
  #             ((混合)分配的維度(變數、基因)數目)
  #alternative: The distribution types vector used for creating the (mixture) distribution. 
  #             (創造(混合)分配所需的分配類型)
  #alternumber: The number of components in the mixture distribution. 
  #             (混合分配所要用到的分配個數)
  #proportion: The proportion for each distribution used for creating mixture distribution. 
  #            (各個分配在混合分配中所要的混合比例)
  #MEAN: The mean vectors for all distributions used for creating mixture distribution.
  #      (混合分配中每個分配的期望值向量)
  #rho: Correlation between any two dimensions(genes) for all distributions used for creating mixture distribution. 
  #     (混合分配中每個分配裡任意兩個維度(變數、基因)之間的相關係數)
  #rhostructure: Structure of correlation matrix for all distributions used for creating mixture distribution. 
  #              (混合分配中的每個分配所需要用到的關係係數矩陣的結構樣式), 
  #              e.g.:"CS": Compound Symmertry, "UN": Unstructured
  #DF: Degrees of freedom of multivariate t distribution(if it is needed). 
  #    (多維t分配的自由度(如果需要用到的話), DF 需要大於 2)
  #    DF > 2
  
  #Program start:
  #seed: All data seed
  #obs: All generated observation.
  #j: The value of indicator function: 1 = Distribution 1 of "alternative" vector, 
  #                                    2 = Distribution 2 of "alternative" vector.
  #p: The cumulative probability of the proportion for each distribution in the mixture distribution.
  #per: Each generated observation.
  seed = runif(n) 
  obs = c()
  for(i in seed){
    j = 1 
    p = proportion[1]
    while(TRUE){
      if(j > alternumber){
        break
      }else if(i <= p){
        if(alternative[j] == "normal"){
          per = mvnormal(n = 1, probenumber = probenumber, MEAN = MEAN[[j]],
                         rho = rho[j], rhostructure = rhostructure[j])
        }else if(alternative[j] == "t"){
          per = mvt(n = 1, probenumber = probenumber, rho = rho[j], DF = DF, 
                    rhostructure = rhostructure[j])
        }else{
          per = NA
        }
        obs = rbind(obs, per)
        break
      }else{
        j = j + 1
        p = p + proportion[j]
      }
    }
  }
  return(obs)
}

##########################################################################################
#Mixture of Normal with same mean setting using fixed Mean Difference
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.9)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2)) {
          print(difference)
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, proportion = prop, 
                               MEAN = list(rep(0, probenumber), 
                                           rep(0, probenumber)),
                               rho = rho)
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, proportion = prop,
                               MEAN = list(rep(0 + difference, probenumber), 
                                           rep(0 + difference, probenumber)),
                               rho = rho)
            obs = rbind(obs1, obs2)
            colnames(obs) = paste0("V",(1:ncol(obs)))
            rownames(obs) = paste0("N",(1:nrow(obs)))
            
            sig = 0.05
            print("Hotelling")
            #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
            resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
            #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
            #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
            
            print("Global")
            resultG = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs)
            resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs ,permutations = 1000)
            
            print("GlobalAncova")
            resultGlobalAncova = GlobalAncova(xx = t(obs), group = rep(c(0, 1), c(m/2, m/2)), 
                                              covars=NULL, perm=1000)
            
            print("N-statistic")
            resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
            
            print("GSEA(Category)")
            gseadata = ExpressionSet(t(obs))
            response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
            signset = matrix(1, nrow = 1, ncol = probenumber)
            colnames(signset) = rownames(gseadata)
            rownames(signset) = "geneset"
            resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
            

            if(resultH_Pseudo_Perm$pval < sig){
              countH_Pseudo_Perm = countH_Pseudo_Perm + 1
            }
            
            if(resultG@result[1] < sig){
              countG = countG + 1
            }
            if(resultG_Perm@result[1] < sig){
              countG_Perm = countG_Perm + 1
            }
            
            if(resultGlobalAncova$test.result[2,] < sig){
              countGlobalAncova = countGlobalAncova + 1
            }
            
            if(resultN_Perm$p.value < sig){
              countN_Perm = countN_Perm + 1
            }
            
            if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
              countGSEA_Category = countGSEA_Category + 1
            }
            
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, ",difference:", difference)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(MON with same Mean).csv", row.names = TRUE)

######################################################################################
#Mixture of Normal with different mean setting(0 and 1) using fixed Mean Difference
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.9)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(50, 70)){
        print(m)
        for (difference in c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2)) {
          print(difference)
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, 
                               MEAN = list(rep(0, probenumber),
                                           rep(1, probenumber)),
                               proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, 
                               MEAN = list(rep(0 + difference, probenumber),
                                           rep(1 + difference, probenumber)),
                               proportion = prop, rho = rho)
            obs = rbind(obs1, obs2)
            colnames(obs) = paste0("V",(1:ncol(obs)))
            rownames(obs) = paste0("N",(1:nrow(obs)))
            
            sig = 0.05
            print("Hotelling")
            #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
            resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
            #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
            #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
            
            print("Global")
            resultG = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs)
            resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs ,permutations = 1000)
            
            print("GlobalAncova")
            resultGlobalAncova = GlobalAncova(xx = t(obs), group = rep(c(0, 1), c(m/2, m/2)), 
                                              covars=NULL, perm=1000)
            
            print("N-statistic")
            resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
            
            print("GSEA(Category)")
            gseadata = ExpressionSet(t(obs))
            response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
            signset = matrix(1, nrow = 1, ncol = probenumber)
            colnames(signset) = rownames(gseadata)
            rownames(signset) = "geneset"
            resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
            
            
            if(resultH_Pseudo_Perm$pval < sig){
              countH_Pseudo_Perm = countH_Pseudo_Perm + 1
            }
            
            if(resultG@result[1] < sig){
              countG = countG + 1
            }
            if(resultG_Perm@result[1] < sig){
              countG_Perm = countG_Perm + 1
            }
            
            if(resultGlobalAncova$test.result[2,] < sig){
              countGlobalAncova = countGlobalAncova + 1
            }
            
            if(resultN_Perm$p.value < sig){
              countN_Perm = countN_Perm + 1
            }
            
            if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
              countGSEA_Category = countGSEA_Category + 1
            }
            
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, ",difference:", difference)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(MON with Different Mean(0, 1)).csv", row.names = TRUE)
######################################################################################
#Mixture of Normal with different mean setting(0 and 2) using fixed Mean Difference
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.9)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(50, 70)){
        print(m)
        for (difference in c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2)) {
          print(difference)
          #GSA Tools Start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, 
                               MEAN = list(rep(0, probenumber),
                                           rep(2, probenumber)),
                               proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, 
                               MEAN = list(rep(0 + difference, probenumber),
                                           rep(2 + difference, probenumber)),
                               proportion = prop, rho = rho)
            obs = rbind(obs1, obs2)
            colnames(obs) = paste0("V",(1:ncol(obs)))
            rownames(obs) = paste0("N",(1:nrow(obs)))
            
            sig = 0.05
            print("Hotelling")
            #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
            resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
            #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
            #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
            
            print("Global")
            resultG = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs)
            resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs ,permutations = 1000)
            
            print("GlobalAncova")
            resultGlobalAncova = GlobalAncova(xx = t(obs), group = rep(c(0, 1), c(m/2, m/2)), 
                                              covars=NULL, perm=1000)
            
            print("N-statistic")
            resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
            
            print("GSEA(Category)")
            gseadata = ExpressionSet(t(obs))
            response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
            signset = matrix(1, nrow = 1, ncol = probenumber)
            colnames(signset) = rownames(gseadata)
            rownames(signset) = "geneset"
            resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
            
            
            if(resultH_Pseudo_Perm$pval < sig){
              countH_Pseudo_Perm = countH_Pseudo_Perm + 1
            }
            
            if(resultG@result[1] < sig){
              countG = countG + 1
            }
            if(resultG_Perm@result[1] < sig){
              countG_Perm = countG_Perm + 1
            }
            
            if(resultGlobalAncova$test.result[2,] < sig){
              countGlobalAncova = countGlobalAncova + 1
            }
            
            if(resultN_Perm$p.value < sig){
              countN_Perm = countN_Perm + 1
            }
            
            if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
              countGSEA_Category = countGSEA_Category + 1
            }
            
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, ",difference:", difference)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(MON with Different Mean(0, 2)).csv", row.names = TRUE)

#########################################################################################################################################################
#Multivariate t-distribution using fixed Mean Difference
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.9)



proportionlist = list()
proportionlist[[1]] = c(1)



times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(paste0("rho = ", rho))
    for(DF in 3){
      print(paste0("DF = ", DF))
      for(probenumber in 30){
        print(paste0("probenumber = ", probenumber))
        for(m in c(70)){
          print(paste0("sample = ", m))
          for (difference in c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2)){
            print(difference)
            #GSA Tools start
            names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
            for(i in 1: length(names)) assign(names[i], 0)
            for(s in 1: times){
              print(s)
              set.seed(s)
              obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                                 alternumber = 1, proportion = prop, rho = rho, DF = DF)
              set.seed(s + 2*times)
              obs2 = difference + generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                                              alternumber = 1, proportion = prop, rho = rho, DF = DF)
              obs = rbind(obs1, obs2)
              colnames(obs) = paste0("V",(1:ncol(obs)))
              rownames(obs) = paste0("N",(1:nrow(obs)))
              
              sig = 0.05
              print("Hotelling")
              #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
              resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
              #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
              #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
              
              print("Global")
              resultG = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs)
              resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs ,permutations = 1000)
              
              print("GlobalAncova")
              resultGlobalAncova = GlobalAncova(xx = t(obs), group = rep(c(0, 1), c(m/2, m/2)), 
                                                covars=NULL, perm=1000)
              
              print("N-statistic")
              resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
              
              print("GSEA(Category)")
              gseadata = ExpressionSet(t(obs))
              response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
              signset = matrix(1, nrow = 1, ncol = probenumber)
              colnames(signset) = rownames(gseadata)
              rownames(signset) = "geneset"
              resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
              
              
              if(resultH_Pseudo_Perm$pval < sig){
                countH_Pseudo_Perm = countH_Pseudo_Perm + 1
              }
              
              if(resultG@result[1] < sig){
                countG = countG + 1
              }
              if(resultG_Perm@result[1] < sig){
                countG_Perm = countG_Perm + 1
              }
              
              if(resultGlobalAncova$test.result[2,] < sig){
                countGlobalAncova = countGlobalAncova + 1
              }
              
              if(resultN_Perm$p.value < sig){
                countN_Perm = countN_Perm + 1
              }
              
              if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
                countGSEA_Category = countGSEA_Category + 1
              }
              
              
            }
            count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
            #GSA Tools end
            colnames(count) = names
            rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                     ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                     ",df:", DF,
                                     ",probes:", probenumber, 
                                     ",size:", m, 
                                     ",difference:", difference)
            simulationresult <<- rbind(simulationresult, (count/times))
          }
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(multivariate t).csv", row.names = TRUE)

############################################################################################################################################
#Multivariate Normal Power Simulation(rho=0.1 and 0.9))
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1)
rholist[[2]] = c(0.9)

proportionlist = list()
proportionlist[[1]] = c(1)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(paste0("rho = ", rho))
    for(probenumber in 30){
      print(paste0("probenumber = ", probenumber))
      for(m in c(70)){
        print(paste0("sample = ", m))
        for (difference in c(0, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 2)){
          print(difference)
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, 
                               MEAN = list(rep(0, probenumber)), 
                               proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, 
                               MEAN = list(rep(0 + difference, probenumber)), 
                               proportion = prop, rho = rho)
            obs = rbind(obs1, obs2)
            colnames(obs) = paste0("V",(1:ncol(obs)))
            rownames(obs) = paste0("N",(1:nrow(obs)))
            
            sig = 0.05
            print("Hotelling")
            #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
            resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
            #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
            #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
            
            print("Global")
            resultG = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs)
            resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(obs1),nrow(obs2))), obs ,permutations = 1000)
            
            print("GlobalAncova")
            resultGlobalAncova = GlobalAncova(xx = t(obs), group = rep(c(0, 1), c(m/2, m/2)), 
                                              covars=NULL, perm=1000)
            
            print("N-statistic")
            resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
            
            print("GSEA(Category)")
            gseadata = ExpressionSet(t(obs))
            response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
            signset = matrix(1, nrow = 1, ncol = probenumber)
            colnames(signset) = rownames(gseadata)
            rownames(signset) = "geneset"
            resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
            
            
            if(resultH_Pseudo_Perm$pval < sig){
              countH_Pseudo_Perm = countH_Pseudo_Perm + 1
            }
            
            if(resultG@result[1] < sig){
              countG = countG + 1
            }
            if(resultG_Perm@result[1] < sig){
              countG_Perm = countG_Perm + 1
            }
            
            if(resultGlobalAncova$test.result[2,] < sig){
              countGlobalAncova = countGlobalAncova + 1
            }
            
            if(resultN_Perm$p.value < sig){
              countN_Perm = countN_Perm + 1
            }
            
            if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
              countGSEA_Category = countGSEA_Category + 1
            }
            
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, 
                                   ",size:", m, 
                                   ",difference:", difference)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(multivariate normal distribution(rho=0.1 and 0.9)).csv", row.names = TRUE)

