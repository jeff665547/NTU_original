#Type I error rate analysis
library(MVN)
library(mvtnorm)
library(energy)
generateobs = function(n, probenumber, rho, rhostructure = "CS", setseed = FALSE){
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

simulationresult = c()

times = 1000
for(rho in c(0.1, 0.9)){
  print(rho)
  for(probenumber in c(2, 4)){
    print(probenumber)
    for(m in c(10, 20, 30, 40, 50, 100)){
      print(m)
      names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
      for(i in 1:length(names)) assign(names[i], 0)
      for(s in 1: times){
        set.seed(1000 + s)
        obs = generateobs(n = m, probenumber = probenumber, rho = rho)
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
      rownames(count) = paste0("size:", m, ",probes:", probenumber, ",rho:", rho)
      simulationresult = rbind(simulationresult, (count/times))
    }
  }
}

setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_1000times.csv", row.names = TRUE)


###########################################################################################
#8, 9 probes for n = 10, 30
simulationresult = c()

times = 1000
for(rho in c(0.1, 0.9)){
  print(rho)
  for(probenumber in c(8, 9)){
    print(probenumber)
    for(m in c(10, 30)){
      print(m)
      names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
      for(i in 1:length(names)) assign(names[i], 0)
      for(s in 1: times){
        set.seed(500 + s)
        obs = generateobs(n = m, probenumber = probenumber, rho = rho)
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
      rownames(count) = paste0("size:", m, ",probes:", probenumber, ",rho:", rho)
      simulationresult = rbind(simulationresult, (count/times))
    }
  }
}

setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_mix(8(10)probes).csv", row.names = TRUE)


#######################################################################################################################
#25, 49 probes for n = 50, 100
simulationresult = c()
times = 1000

for(rho in c(0.1, 0.9)){
  print(rho)
  for(probenumber in c(25, 49)){
    print(probenumber)
    for(m in c(50, 100)){
      print(m)
      names = c("countM", "countH", "countR", "countN", "countHR", "countHN", "countRN", "count2", "countHRN")
      for(i in 1:length(names)) assign(names[i], 0)
      for(s in 1: times){
        set.seed(500 + s)
        print(s)
        obs = generateobs(n = m, probenumber = probenumber, rho = rho)
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
      rownames(count) = paste0("size:", m, ",probes:", probenumber, ",rho:", rho)
      simulationresult = rbind(simulationresult, (count/times))
    }
  }
}


setwd("C:/Users/jeff/Desktop/")
write.csv(simulationresult, file = "Simulation result_mix(25(50)probes).csv", row.names = TRUE)


