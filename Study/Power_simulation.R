#Power analysis
library(MVN)
library(mvtnorm)
library(energy)

#multivariate normal distribution
mvnormal = function(n, probenumber, rho, rhostructure = "CS", setseed = FALSE){
  #rho is the correlation among data
  if(rhostructure == "CS"){
    SIGMA = matrix(rho, probenumber, probenumber)
    diag(SIGMA) = 1
  }
  if(rhostructure == "UN"){
    rhovector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    SIGMA = matrix(c(sample(rhovector, size = probenumber*probenumber, replace = TRUE)), 
                   probenumber, probenumber)
    diag(SIGMA) = 1
  }
  if(setseed != FALSE){
    set.seed(setseed)
  }
  obs = rmvnorm(n, mean = rep(0, probenumber), sigma = SIGMA)
  return(obs)
}
mvnormal(1, 2, 0.9)

#multivariate t distribution
mvt = function(n, probenumber, rho, DF, rhostructure = "CS", setseed = FALSE){
  #rho is the correlation among data
  if(rhostructure == "CS"){
    rhomatrix = matrix(rho, probenumber, probenumber)
    diag(rhomatrix) = 1
  }
  if(rhostructure == "UN"){
    rhovector = c(0.1, 0.3, 0.5, 0.7, 0.9)
    rhomatrix = matrix(c(sample(rhovector, size = probenumber*probenumber, replace = TRUE)), 
                   probenumber, probenumber)
    diag(rhomatrix) = 1
  }
  
  if(setseed != FALSE){
    set.seed(setseed)
  }
  
  if(DF > 2){
    SIGMA <<- rhomatrix*((DF - 2)/DF) #transform into scale matrix 
  }else{
    cat("error: df must > 2 \n")
  }
  
  obs = rmvt(n, sigma = SIGMA, df = DF)
  return(obs)
}

#Create alternative: Mixture of multivariate normal, multivariate t distribution
generateobs = function(n, probenumber, alternative = c("normal", "normal"), alternumber = 2, 
                       proportion = rep(1/alternumber, alternumber), rho, 
                       rhostructure = rep("CS", alternumber), DF = Inf){
  seed = runif(n)
  obs = c()
  for(i in seed){
    j = 1
    p = proportion[1]
    while(TRUE){
      if(j > alternumber){
        break
      }else if(i < p){
        if(alternative[j] == "normal"){
          per = mvnormal(n = 1, probenumber = probenumber, rho = rho[j], 
                         rhostructure = rhostructure[j])
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

set.seed(1)
generateobs(n = 4, probenumber = 3, rho = c(0.1, 0.9))

##########################################################################################
#normal mixture(2, 4, 8, 9)

simulationresult = c()

rholist = list()
rholist[[1]] = c(0, 0.9)
rholist[[2]] = c(0, 0.5)


proportionlist = list()
proportionlist[[1]] = c(0.5, 0.5)
proportionlist[[2]] = c(0.1, 0.9)


times = 500

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(2, 4, 8, 9)){
      print(probenumber)
      for(m in c(10, 30, 50, 100)){
        print(m)
        names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
        for(i in 1: length(names)) assign(names[i], 0)
        for(s in 1: times){
          print(s)
          set.seed(s)
          obs = generateobs(n = m, probenumber = probenumber, proportion = prop, rho = rho)
          
          sig = 0.05
          resultM <<- mvn(obs, mvnTest = "mardia")$multivariateNormality
          resultH <<- mvn(obs, mvnTest = "hz")$multivariateNormality
          resultR <<- mvn(obs, mvnTest = "royston")$multivariateNormality
          resultN <<- mvnorm.etest(obs, R = 999)
          
          if(sum(c(resultH$`p value` < sig, resultR$`p value` < sig, resultN$p.value < sig)) >= 2){
            count2 = count2 + 1
          }
          if(as.numeric(as.character(resultM$`p value`[which.min(as.numeric(as.character(resultM$`p value`[-3])))])) < sig){
            countM = countM + 1
          }
          if(resultH$`p value` < sig){
            countH = countH + 1
            if(resultR$`p value` < sig){
              countHR = countHR + 1
              if(resultN$p.value < sig){
                countHRN = countHRN + 1
              }
            }
            if(resultN$p.value < sig){
              countHN = countHN + 1
            }
          }
          if(resultR$`p value` < sig){
            countR = countR + 1
            if(resultN$p.value < sig){
              countRN = countRN + 1
            }
          }
          if(resultN$p.value < sig){
            countN = countN + 1
          }
        }
        count <<- matrix(c(countM, countH, countR, countN, countHR, countHN, countRN, count2, countHRN), 1, length(names))
        colnames(count) = names
        rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                 ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                 ",probes:", probenumber, ",size:", m)
        simulationresult <<- rbind(simulationresult, (count/times))
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_Power(2,4,8,9).csv", row.names = TRUE)

##############################################################################################
#normal mixture(25, 49)

simulationresult = c()

times = 500

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(rho)
    for(probenumber in c(25, 49)){
      print(probenumber)
      for(m in c(50, 100)){
        print(m)
        names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
        for(i in 1: length(names)) assign(names[i], 0)
        for(s in 1: times){
          print(s)
          set.seed(s)
          obs = generateobs(n = m, probenumber = probenumber, proportion = prop, rho = rho)
          
          sig = 0.05
          resultM <<- mvn(obs, mvnTest = "mardia")$multivariateNormality
          resultH <<- mvn(obs, mvnTest = "hz")$multivariateNormality
          resultR <<- mvn(obs, mvnTest = "royston")$multivariateNormality
          resultN <<- mvnorm.etest(obs, R = 999)
          
          if(sum(c(resultH$`p value` < sig, resultR$`p value` < sig, resultN$p.value < sig)) >= 2){
            count2 = count2 + 1
          }
          if(as.numeric(as.character(resultM$`p value`[which.min(as.numeric(as.character(resultM$`p value`[-3])))])) < sig){
            countM = countM + 1
          }
          if(resultH$`p value` < sig){
            countH = countH + 1
            if(resultR$`p value` < sig){
              countHR = countHR + 1
              if(resultN$p.value < sig){
                countHRN = countHRN + 1
              }
            }
            if(resultN$p.value < sig){
              countHN = countHN + 1
            }
          }
          if(resultR$`p value` < sig){
            countR = countR + 1
            if(resultN$p.value < sig){
              countRN = countRN + 1
            }
          }
          if(resultN$p.value < sig){
            countN = countN + 1
          }
        }
        count = matrix(c(countM, countH, countR, countN, countHR, countHN, countRN, count2, countHRN), 1, length(names))
        colnames(count) = names
        rownames(count) = paste0("composition:", paste0(as.character(prop), collapse = "", sep = " "),
                                 ",rho:", paste0(as.character(rho), collapse = "", sep = " "), 
                                 ",probes:", probenumber, ",size:", m)
        simulationresult = rbind(simulationresult, (count/times))
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_Power(25, 49).csv", row.names = TRUE)



#######################################################################################################################
#multivariate t(2, 4, 8, 9)

simulationresult = c()

rholist = list()
rholist[[1]] = c(0.9)
rholist[[2]] = c(0)



proportionlist = list()
proportionlist[[1]] = c(1)



times = 500

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(paste0("rho = ", rho))
    for(DF in c(3, 9, 15)){
      print(paste0("DF = ", DF))
      for(probenumber in c(2, 4, 8, 9)){
        print(paste0("probenumber = ", probenumber))
        for(m in c(10, 30, 50, 100)){
          print(paste0("sample = ", m))
          names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            print(s)
            set.seed(s)
            obs = generateobs(n = m, probenumber = probenumber, alternative = "t", 
                              alternumber = 1, proportion = prop, rho = rho, DF = DF)
            sig = 0.05
            resultM <<- mvn(obs, mvnTest = "mardia")$multivariateNormality
            resultH <<- mvn(obs, mvnTest = "hz")$multivariateNormality
            resultR <<- mvn(obs, mvnTest = "royston")$multivariateNormality
            resultN <<- mvnorm.etest(obs, R = 999)
            
            if(sum(c(resultH$`p value` < sig, resultR$`p value` < sig, resultN$p.value < sig)) >= 2){
              count2 = count2 + 1
            }
            if(as.numeric(as.character(resultM$`p value`[which.min(as.numeric(as.character(resultM$`p value`[-3])))])) < sig){
              countM = countM + 1
            }
            if(resultH$`p value` < sig){
              countH = countH + 1
              if(resultR$`p value` < sig){
                countHR = countHR + 1
                if(resultN$p.value < sig){
                  countHRN = countHRN + 1
                }
              }
              if(resultN$p.value < sig){
                countHN = countHN + 1
              }
            }
            if(resultR$`p value` < sig){
              countR = countR + 1
              if(resultN$p.value < sig){
                countRN = countRN + 1
              }
            }
            if(resultN$p.value < sig){
              countN = countN + 1
            }
          }
          count = matrix(c(countM, countH, countR, countN, countHR, countHN, countRN, count2, countHRN), 1, length(names))
          colnames(count) = names
          rownames(count) = paste0("composition:", prop ,",rho:", rho, ",df:", DF, ",probes:", probenumber, ",size:", m)
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}

setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_tPower().csv", row.names = TRUE)


#########################################################################################################################
#multivariate t(25, 49)

simulationresult = c()


rholist = list()
rholist[[1]] = c(0.9)
rholist[[2]] = c(0)


proportionlist = list()
proportionlist[[1]] = c(1)

times = 500

for(prop in proportionlist){
  print(prop)
  for(rho in rholist){
    print(paste0("rho = ", rho))
    for(DF in c(3, 9, 15)){
      print(paste0("DF = ", DF))
      for(probenumber in c(25, 49)){
        print(paste0("probenumber = ", probenumber))
        for(m in c(50, 100)){
          print(paste0("sample = ", m))
          names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
          for(i in 1: length(names)) assign(names[i], 0)
          for(s in 1: times){
            if(s == 323){
              s = times + 1
            }
            set.seed(s)
            print(s)
            
            obs = generateobs(n = m, probenumber = probenumber, alternative = "t", 
                              alternumber = 1, proportion = prop, rho = rho, DF = DF)
            sig = 0.05
            resultM <<- mvn(obs, mvnTest = "mardia")$multivariateNormality
            resultH <<- mvn(obs, mvnTest = "hz")$multivariateNormality
            resultR <<- mvn(obs, mvnTest = "royston")$multivariateNormality
            resultN <<- mvnorm.etest(obs, R = 999)
            
            if(sum(c(resultH$`p value` < sig, resultR$`p value` < sig, resultN$p.value < sig)) >= 2){
              count2 = count2 + 1
            }
            if(as.numeric(as.character(resultM$`p value`[which.min(as.numeric(as.character(resultM$`p value`[-3])))])) < sig){
              countM = countM + 1
            }
            if(resultH$`p value` < sig){
              countH = countH + 1
              if(resultR$`p value` < sig){
                countHR = countHR + 1
                if(resultN$p.value < sig){
                  countHRN = countHRN + 1
                }
              }
              if(resultN$p.value < sig){
                countHN = countHN + 1
              }
            }
            if(resultR$`p value` < sig){
              countR = countR + 1
              if(resultN$p.value < sig){
                countRN = countRN + 1
              }
            }
            if(resultN$p.value < sig){
              countN = countN + 1
            }
          }
          count = matrix(c(countM, countH, countR, countN, countHR, countHN, countRN, count2, countHRN), 1, length(names))
          colnames(count) = names
          rownames(count) = paste0("composition:", prop ,",rho:", rho, ",df:", DF, ",probes:", probenumber, ",size:", m)
          simulationresult = rbind(simulationresult, (count/times))
        }
      }
    }
  }
}


setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_tPower(25,49).csv", row.names = TRUE)
