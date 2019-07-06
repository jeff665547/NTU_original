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

################################################################################################
#All Test 
Test = function(Control, Case, all.sample.size, probenumber, 
                AA.test, BB.test, CC.test, DD.test, EE.test, sig = 0.05){
  
  ALL = rbind(Control, Case)
  colnames(ALL) = paste0("V",(1:ncol(ALL)))
  rownames(ALL) = paste0("N",(1:nrow(ALL)))
  m = all.sample.size
  
  print("Hotelling")
  resultH_Pseudo_Perm = hotelling.test(Control, Case, shrinkage = TRUE, perm = TRUE, B = 1000, progBar =  FALSE)
  
  print("Global")
  resultG_Perm = gt(rep(c(TRUE, FALSE), c(nrow(Control),nrow(Case))), ALL ,permutations = 1000)
  
  print("GlobalAncova")
  resultGlobalAncova = GlobalAncova(xx = t(ALL), group = rep(c(0, 1), c(m/2, m/2)), 
                                    covars=NULL, perm=1000)
  
  print("N-statistic")
  resultN_Perm = eqdist.etest(ALL, sizes=c(m/2, m/2) , R = 1000)
  
  print("GSEA(Category)")
  gseadata = ExpressionSet(t(ALL))
  response = as.factor(rep(c("A", "B"), c(m/2, m/2)))
  signset = matrix(1, nrow = 1, ncol = probenumber)
  colnames(signset) = rownames(gseadata)
  rownames(signset) = "geneset"
  resultGSEA_Category = gseattperm(gseadata, response, signset, 1000)
  
  if(resultH_Pseudo_Perm$pval < sig){
    AA.test = AA.test + 1
  }
  
  if(resultG_Perm@result[1] < sig){
    BB.test = BB.test + 1
  }
  
  if(resultGlobalAncova$test.result[2,] < sig){
    CC.test = CC.test + 1
  }
  
  if(resultN_Perm$p.value < sig){
    DD.test = DD.test + 1
  }
  
  if((resultGSEA_Category[1] < sig/2) || (resultGSEA_Category[2] < sig/2)){
    EE.test = EE.test + 1
  }
  
  return(c(AA.test, BB.test, CC.test, DD.test, EE.test))
}
##########################################################################################
#Mixture of Normal with same mean using Dynamic Mean Difference 1(Mix 0 0)
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.9)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

set.seed(2)
difference1 = sample(c(rep(0.1, 15), rep(0.3, 7), rep(0.5, 8)))
set.seed(2)
difference2 = sample(c(rep(0.1, 15), rep(0.3, 8), rep(0.5, 7)))
set.seed(2)
difference3 = sample(c(rep(0.1, 7), rep(0.3, 15), rep(0.5, 8)))
set.seed(2)
difference4 = sample(c(rep(0.1, 8), rep(0.3, 15), rep(0.5, 7)))
set.seed(2)
difference5 = sample(c(rep(0.1, 7), rep(0.3, 8), rep(0.5, 15)))
set.seed(2)
difference6 = sample(c(rep(0.1, 8), rep(0.3, 7), rep(0.5, 15)))
set.seed(2)
difference7 = sample(c(rep(0.1, 10), rep(0.3, 10), rep(0.5, 10)))
set.seed(2)
difference8 = sample(c(rep(0, 12), rep(0.5, 18)))
set.seed(2)
difference9 = sample(c(rep(0, 27), rep(3, 3)))
difference10 = sample(c(rep(0, 0), rep(0.3, 30)))





differencelist = list()
differencelist[[1]] = difference1
differencelist[[2]] = difference2
differencelist[[3]] = difference3
differencelist[[4]] = difference4
differencelist[[5]] = difference5
differencelist[[6]] = difference6
differencelist[[7]] = difference7
differencelist[[8]] = difference8
differencelist[[9]] = difference9
differencelist[[10]] = difference10

lapply(X = differencelist, FUN = table)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in differencelist) {
          print(table(difference))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = t(difference + t(generateobs(n = m/2, probenumber = probenumber, proportion = prop, rho = rho)))
            result = Test(Control = obs1, Case = obs2, 
                          all.sample.size = m,
                          probenumber = probenumber, 
                          AA.test = countH_Pseudo_Perm,
                          BB.test = countG_Perm, 
                          CC.test = countGlobalAncova,
                          DD.test = countN_Perm,
                          EE.test = countGSEA_Category)
            
            for(i in 1: length(names)) assign(names[i], result[i])
          }
          count = matrix(c(countH_Pseudo_Perm, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, paste0(",difference(", paste(names(table(difference)), collapse = " "), ")"), 
                                   paste0(as.character((table(difference))), collapse = "", sep = " "))
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(Dynamic Mean Difference(Mix 0, 0, rho 0.9)).csv", row.names = TRUE)
#######################################################################################################################################################
#Mixture of Normal with same mean using Dynamic Mean Difference 2(Mix 0 2)
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.9)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

set.seed(2)
difference1 = sample(c(rep(0.1, 15), rep(0.3, 7), rep(0.5, 8)))
set.seed(2)
difference2 = sample(c(rep(0.1, 15), rep(0.3, 8), rep(0.5, 7)))
set.seed(2)
difference3 = sample(c(rep(0.1, 7), rep(0.3, 15), rep(0.5, 8)))
set.seed(2)
difference4 = sample(c(rep(0.1, 8), rep(0.3, 15), rep(0.5, 7)))
set.seed(2)
difference5 = sample(c(rep(0.1, 7), rep(0.3, 8), rep(0.5, 15)))
set.seed(2)
difference6 = sample(c(rep(0.1, 8), rep(0.3, 7), rep(0.5, 15)))
set.seed(2)
difference7 = sample(c(rep(0.1, 10), rep(0.3, 10), rep(0.5, 10)))
set.seed(2)
difference8 = sample(c(rep(0, 12), rep(0.5, 18)))
set.seed(2)
difference9 = sample(c(rep(0, 27), rep(3, 3)))
difference10 = sample(c(rep(0, 0), rep(0.3, 30)))





differencelist = list()
differencelist[[1]] = difference1
differencelist[[2]] = difference2
differencelist[[3]] = difference3
differencelist[[4]] = difference4
differencelist[[5]] = difference5
differencelist[[6]] = difference6
differencelist[[7]] = difference7
differencelist[[8]] = difference8
differencelist[[9]] = difference9
differencelist[[10]] = difference10

lapply(X = differencelist, FUN = table)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in differencelist) {
          print(table(difference))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, 
                               MEAN = list(rep(0, probenumber),
                                           rep(2, probenumber)),
                               proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = t(difference + t(generateobs(n = m/2, probenumber = probenumber, 
                                                MEAN = list(rep(0, probenumber),
                                                            rep(2, probenumber)),
                                                proportion = prop, rho = rho)))
            result = Test(Control = obs1, Case = obs2, 
                          all.sample.size = m,
                          probenumber = probenumber, 
                          AA.test = countH_Pseudo_Perm,
                          BB.test = countG_Perm, 
                          CC.test = countGlobalAncova,
                          DD.test = countN_Perm,
                          EE.test = countGSEA_Category)
            
            for(i in 1: length(names)) assign(names[i], result[i])
          }
          count = matrix(c(countH_Pseudo_Perm, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, paste0(",difference(", paste(names(table(difference)), collapse = " "), ")"), 
                                   paste0(as.character((table(difference))), collapse = "", sep = " "))
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(Dynamic Mean Difference(Mix 0, 2)).csv", row.names = TRUE)
##########################################################################################################################################################
#MVT using Dynamic Mean Difference 3(MVT 0.1)
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1)

proportionlist = list()
proportionlist[[1]] = c(1)

DF = 3

set.seed(2)
difference1 = sample(c(rep(0.1, 15), rep(0.3, 7), rep(0.5, 8)))
set.seed(2)
difference2 = sample(c(rep(0.1, 15), rep(0.3, 8), rep(0.5, 7)))
set.seed(2)
difference3 = sample(c(rep(0.1, 7), rep(0.3, 15), rep(0.5, 8)))
set.seed(2)
difference4 = sample(c(rep(0.1, 8), rep(0.3, 15), rep(0.5, 7)))
set.seed(2)
difference5 = sample(c(rep(0.1, 7), rep(0.3, 8), rep(0.5, 15)))
set.seed(2)
difference6 = sample(c(rep(0.1, 8), rep(0.3, 7), rep(0.5, 15)))
set.seed(2)
difference7 = sample(c(rep(0.1, 10), rep(0.3, 10), rep(0.5, 10)))
set.seed(2)
difference8 = sample(c(rep(0, 12), rep(0.5, 18)))
set.seed(2)
difference9 = sample(c(rep(0, 27), rep(3, 3)))
difference10 = sample(c(rep(0, 0), rep(0.3, 30)))





differencelist = list()
differencelist[[1]] = difference1
differencelist[[2]] = difference2
differencelist[[3]] = difference3
differencelist[[4]] = difference4
differencelist[[5]] = difference5
differencelist[[6]] = difference6
differencelist[[7]] = difference7
differencelist[[8]] = difference8
differencelist[[9]] = difference9
differencelist[[10]] = difference10

lapply(X = differencelist, FUN = table)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in differencelist) {
          print(table(difference))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                               alternumber = 1, proportion = prop, rho = rho, DF = DF)
            set.seed(s + 2*times)
            obs2 = t(difference + t(generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                                                alternumber = 1, proportion = prop, rho = rho, DF = DF)))
            result = Test(Control = obs1, Case = obs2, 
                          all.sample.size = m,
                          probenumber = probenumber, 
                          AA.test = countH_Pseudo_Perm,
                          BB.test = countG_Perm, 
                          CC.test = countGlobalAncova,
                          DD.test = countN_Perm,
                          EE.test = countGSEA_Category)
            
            for(i in 1: length(names)) assign(names[i], result[i])
          }
          count = matrix(c(countH_Pseudo_Perm, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, paste0(",difference(", paste(names(table(difference)), collapse = " "), ")"), 
                                   paste0(as.character((table(difference))), collapse = "", sep = " "))
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(Dynamic Mean Difference(MVT rho 0.1)).csv", row.names = TRUE)

##########################################################################################################################################################
#MVT using Dynamic Mean Difference 4(MVT rho 0.9)
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.9)

proportionlist = list()
proportionlist[[1]] = c(1)

DF = 3

set.seed(2)
difference1 = sample(c(rep(0.1, 15), rep(0.3, 7), rep(0.5, 8)))
set.seed(2)
difference2 = sample(c(rep(0.1, 15), rep(0.3, 8), rep(0.5, 7)))
set.seed(2)
difference3 = sample(c(rep(0.1, 7), rep(0.3, 15), rep(0.5, 8)))
set.seed(2)
difference4 = sample(c(rep(0.1, 8), rep(0.3, 15), rep(0.5, 7)))
set.seed(2)
difference5 = sample(c(rep(0.1, 7), rep(0.3, 8), rep(0.5, 15)))
set.seed(2)
difference6 = sample(c(rep(0.1, 8), rep(0.3, 7), rep(0.5, 15)))
set.seed(2)
difference7 = sample(c(rep(0.1, 10), rep(0.3, 10), rep(0.5, 10)))
set.seed(2)
difference8 = sample(c(rep(0, 12), rep(0.5, 18)))
set.seed(2)
difference9 = sample(c(rep(0, 27), rep(3, 3)))
difference10 = sample(c(rep(0, 0), rep(0.3, 30)))





differencelist = list()
differencelist[[1]] = difference1
differencelist[[2]] = difference2
differencelist[[3]] = difference3
differencelist[[4]] = difference4
differencelist[[5]] = difference5
differencelist[[6]] = difference6
differencelist[[7]] = difference7
differencelist[[8]] = difference8
differencelist[[9]] = difference9
differencelist[[10]] = difference10

lapply(X = differencelist, FUN = table)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in differencelist) {
          print(table(difference))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                               alternumber = 1, proportion = prop, rho = rho, DF = DF)
            set.seed(s + 2*times)
            obs2 = t(difference + t(generateobs(n = m/2, probenumber = probenumber, alternative = "t", 
                                                alternumber = 1, proportion = prop, rho = rho, DF = DF)))
            result = Test(Control = obs1, Case = obs2, 
                          all.sample.size = m,
                          probenumber = probenumber, 
                          AA.test = countH_Pseudo_Perm,
                          BB.test = countG_Perm, 
                          CC.test = countGlobalAncova,
                          DD.test = countN_Perm,
                          EE.test = countGSEA_Category)
            
            for(i in 1: length(names)) assign(names[i], result[i])
          }
          count = matrix(c(countH_Pseudo_Perm, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, paste0(",difference(", paste(names(table(difference)), collapse = " "), ")"), 
                                   paste0(as.character((table(difference))), collapse = "", sep = " "))
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(Dynamic Mean Difference(MVT rho 0.9)).csv", row.names = TRUE)

##########################################################################################################################################################
#Mixture of Normal with same mean using Dynamic Mean Difference 1(Mix 0 0)
simulationresult = c()

rholist = list()
rholist[[1]] = c(0.1, 0.1)

proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)

set.seed(2)
difference1 = sample(c(rep(0.1, 15), rep(0.3, 7), rep(0.5, 8)))
set.seed(2)
difference2 = sample(c(rep(0.1, 15), rep(0.3, 8), rep(0.5, 7)))
set.seed(2)
difference3 = sample(c(rep(0.1, 7), rep(0.3, 15), rep(0.5, 8)))
set.seed(2)
difference4 = sample(c(rep(0.1, 8), rep(0.3, 15), rep(0.5, 7)))
set.seed(2)
difference5 = sample(c(rep(0.1, 7), rep(0.3, 8), rep(0.5, 15)))
set.seed(2)
difference6 = sample(c(rep(0.1, 8), rep(0.3, 7), rep(0.5, 15)))
set.seed(2)
difference7 = sample(c(rep(0.1, 10), rep(0.3, 10), rep(0.5, 10)))
set.seed(2)
difference8 = sample(c(rep(0, 12), rep(0.5, 18)))
set.seed(2)
difference9 = sample(c(rep(0, 27), rep(3, 3)))
difference10 = sample(c(rep(0, 0), rep(0.3, 30)))





differencelist = list()
differencelist[[1]] = difference1
differencelist[[2]] = difference2
differencelist[[3]] = difference3
differencelist[[4]] = difference4
differencelist[[5]] = difference5
differencelist[[6]] = difference6
differencelist[[7]] = difference7
differencelist[[8]] = difference8
differencelist[[9]] = difference9
differencelist[[10]] = difference10

lapply(X = differencelist, FUN = table)

times = 1000

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(30)){
      print(probenumber)
      for(m in c(70)){
        print(m)
        for (difference in differencelist) {
          print(table(difference))
          #GSA Tools start
          names = c("countH_Pseudo_Perm", "countG_Perm", "countGlobalAncova", "countN_Perm", "countGSEA_Category")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs1 = generateobs(n = m/2, probenumber = probenumber, proportion = prop, rho = rho)
            set.seed(s + 2*times)
            obs2 = t(difference + t(generateobs(n = m/2, probenumber = probenumber, proportion = prop, rho = rho)))
            result = Test(Control = obs1, Case = obs2, 
                          all.sample.size = m,
                          probenumber = probenumber, 
                          AA.test = countH_Pseudo_Perm,
                          BB.test = countG_Perm, 
                          CC.test = countGlobalAncova,
                          DD.test = countN_Perm,
                          EE.test = countGSEA_Category)
            
            for(i in 1: length(names)) assign(names[i], result[i])
          }
          count = matrix(c(countH_Pseudo_Perm, countG_Perm, countGlobalAncova, countN_Perm, countGSEA_Category), 1, length(names))
          #GSA Tools end
          colnames(count) = names
          rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                   ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                   ",probes:", probenumber, ",size:", m, paste0(",difference(", paste(names(table(difference)), collapse = " "), ")"), 
                                   paste0(as.character((table(difference))), collapse = "", sep = " "))
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/Study/Results")
write.csv(simulationresult, file = "GSA tools power simulation(Dynamic Mean Difference(Mix 0, 0, rho 0.1)).csv", row.names = TRUE)
