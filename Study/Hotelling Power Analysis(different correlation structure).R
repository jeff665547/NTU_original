#Power analysis(different correlation structure)
library(Hotelling)
library(globaltest)
library(Category)
library(GlobalAncova)
library(MVN)
library(mvtnorm)
library(energy)

#Random sample from multivariate normal distribution(MVN) with customized parameters settings
mvnormal = function(n, probenumber, covariance, MEAN = rep(0, probenumber), variance = 1, covstructure = "CS"){
  #n: Sample size. (樣本數)
  #probenumber: The number of dimensions(genes). (維度，基因數目)
  #covariance: Covariance between any two dimensions(genes). (任意兩個維度，基因之間的共變異係數)
  #covstructure: Structure of covariance matrix. (變異數矩陣的結構樣式), e.g.:"CS": Compound Symmertry, "UN": Unstructured
  
  #Program start:
  #Create correlation matrix
  if(covstructure == "CS"){
    SIGMA = matrix(covariance, probenumber, probenumber)
    diag(SIGMA) = variance
    AA <<- SIGMA
  }else if(covstructure == "UN"){
    covvector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    SIGMA = matrix(c(sample(covvector, size = probenumber*probenumber, replace = TRUE)), 
                   probenumber, probenumber)
    diag(SIGMA) = variance
    AA <<- SIGMA
  }
  
  #Generate observation from MVN from package "mvtnorm".
  obs = rmvnorm(n, mean = MEAN, sigma = SIGMA)
  return(obs)
}


#Random sample from multivariate t distribution(MVT) with customized parameters settings
mvt = function(n, probenumber, covariance, covstructure = "CS", DF){
  #n: Sample size. (樣本數)
  #probenumber: The number of dimensions(genes). (維度，基因數目)
  #covariance: Covariance between any two dimensions(genes). (任意兩個維度，基因之間的共變異係數)
  #covstructure: Structure of covariance matrix. (變異數矩陣的結構樣式), e.g.:"CS": Compound Symmertry, "UN": Unstructured
  #DF: Degrees of freedom of multivariate t distribution, DF > 2. (多為t分配的自由度,DF 需要大於 2)
  
  #Program start:
  #Create correlation matrtix
  if(covstructure == "CS"){
    covmatrix = matrix(covariance, probenumber, probenumber)
    diag(covmatrix) = variance
  }else if(covstructure == "UN"){
    covvector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    covmatrix = matrix(c(sample(covvector, size = probenumber*probenumber, replace = TRUE)), 
                       probenumber, probenumber)
    diag(covmatrix) = 1
  }
  
  #Transform correlation matrix into scale matrix. (將相關係數矩陣轉換成多維度t分配中的scale matrix)
  if(DF > 2){
    SIGMA = (cov2cor(covmatrix))*((DF - 2)/DF) 
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
                       variance = 1, covariance, covstructure = rep("CS", alternumber), 
                       DF = Inf){
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
  #variance: Variance for each gene under the normal distribution
  #      (常態分佈以及t分佈裡每個基因的variance)
  #covariance: Covariance between any two dimensions(genes) for all distributions used for creating mixture distribution. 
  #     (混合分配中每個分配裡任意兩個維度(變數、基因)之間的共變異數)
  #covstructure: Structure of covariance matrix for all distributions used for creating mixture distribution. 
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
                         covariance = covariance[j], variance = variance[j], covstructure = covstructure[j])
        }else if(alternative[j] == "t"){
          per = mvt(n = 1, probenumber = probenumber, covariance = covariance[j], variance = variance[j], 
                    DF = DF, covstructure = covstructure[j])
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

##########################################################################
#Power
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1)
rholist[[2]] = c(0.5)
rholist[[3]] = c(0.9)

varlist = list()
varlist[[1]] = c(1)
varlist[[2]] = c(2)
varlist[[3]] = c(3)
varlist[[4]] = c(4)
varlist[[5]] = c(5)

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
        for (variance in varlist){
          print(paste0("variance for each gene = ", variance))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, MEAN = list(rep(0, probenumber)), 
                               proportion = prop, 
                               variance = 1, covariance = rho*1, covstructure = rep("CS", 1))
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, MEAN = list(rep(0, probenumber)), 
                               proportion = prop, 
                               variance = variance, covariance = rho*variance, covstructure = rep("CS", 1))
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
            

            
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, 
                                   ",size:", m, 
                                   ",variance:", variance)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "Hotelling power simulation(multivariate normal distribution(different rho structure)).csv", row.names = TRUE)


##########################################################################################################################################################
#Type I error rate
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1)
rholist[[2]] = c(0.5)
rholist[[3]] = c(0.9)

varlist = list()
varlist[[1]] = c(1)
varlist[[2]] = c(2)
varlist[[3]] = c(3)
varlist[[4]] = c(4)
varlist[[5]] = c(5)

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
        for (variance in varlist){
          print(paste0("variance for each gene = ", variance))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, MEAN = list(rep(0, probenumber)), 
                               proportion = prop, 
                               variance = variance, covariance = rho*variance, covstructure = rep("CS", 1))
            set.seed(s + 2*times)
            obs2 = generateobs(n = m/2, probenumber = probenumber, alternative = "normal", 
                               alternumber = 1, MEAN = list(rep(0, probenumber)), 
                               proportion = prop, 
                               variance = variance, covariance = rho*variance, covstructure = rep("CS", 1))
            obs = rbind(obs1, obs2)
            colnames(obs) = paste0("V",(1:ncol(obs)))
            rownames(obs) = paste0("N",(1:nrow(obs)))
            
            sig = 0.05
            print("Hotelling")
            #resultH_Pseudo = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = FALSE, B = 1000, progBar =  FALSE)
            resultH_Pseudo_Perm = hotelling.test(obs1, obs2, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
            #resultH = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = FALSE, B = 1000, progBar =  FALSE)
            #resultH_Perm = hotelling.test(obs1, obs2, shrinkage = FALSE, perm = TRUE, B = 1000, progBar =  FALSE)
            
            print("N-statistic")
            resultN_Perm = eqdist.etest(obs, sizes=c(m/2, m/2) , R = 1000)
            
            
            if(resultH_Pseudo_Perm$pval < sig){
              countH_Pseudo_Perm = countH_Pseudo_Perm + 1
            }
            
            if(resultN_Perm$p.value < sig){
              countN_Perm = countN_Perm + 1
            }
            
          }
          count <<- matrix(c(countH_Pseudo_Perm,  countG, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, 
                                   ",size:", m, 
                                   ",variance:", variance)
          simulationresult <<- rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "Hotelling type I error simulation(multivariate normal distribution).csv", row.names = TRUE)

