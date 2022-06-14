# mvnt with ce
library(copula)
library(copent)
library(MVN)

# statistic
tce<-function(x){
  - 0.5 * log( det(cov(x)) ) - copent(x)
}

tce1 = rep(0,10) # copula entropy
mardia1 = rep(0,10) #mardia
hz1 = rep(0,10) # hz
royston1 = rep(0,10) # royston
dh1 = rep(0,10) # dh
energy1 = rep(0,10) # energy distance
k1 = 60
for(k in 1:k1){
  for(i in 1:10){
    # simulation 1
    # mv.NE <- mvdc(normalCopula(0.8), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = i)))
    # simulation 2
    mv.NE <- mvdc(gumbelCopula(i), c("norm", "norm"), list(list(mean = 0, sd =2), list(mean = 0, sd =2)))
    data1 <- rMvdc(800, mv.NE)
    
    tce1[i] = tce1[i] + tce(data1) / k1
    mardia1[i] = mardia1[i] + as.numeric(as.character(mvn(data1,mvnTest = "mardia")$multivariateNormality$Statistic[1])) / k1
    hz1[i] = hz1[i] + mvn(data1,mvnTest = "hz")$multivariateNormality$HZ / k1
    royston1[i] = royston1[i] + mvn(data1,mvnTest = "royston")$multivariateNormality$H / k1
    dh1[i] = dh1[i] + mvn(data1,mvnTest = "dh")$multivariateNormality$E / k1
    energy1[i] = energy1[i] + mvn(data1,mvnTest = "energy")$multivariateNormality$Statistic / k1
  }
}

# xlab1 = "rho" # simulation1
xlab1 = "alpha" # simulation2
x11();plot(tce1,xlab = xlab1,ylab="Copula Entropy");lines(tce1)
x11();plot(mardia1,xlab = xlab1,ylab="Mardia");lines(mardia1)
x11();plot(hz1,xlab = xlab1,ylab="Royston");lines(hz1)
x11();plot(royston1,xlab = xlab1,ylab="Henze-Zirkler");lines(royston1)
x11();plot(dh1,xlab = xlab1,ylab="Dornik-Haansen");lines(dh1)
x11();plot(energy1,xlab = xlab1,ylab="Energy Distance");lines(energy1)
