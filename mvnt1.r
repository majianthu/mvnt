# mvnt with ce
library(copula)
library(copent)
library(MVN)
library(mvnTest)

tce1 = mardia1 = hz1 = royston1 = dh1 = energy1 = ad1 = cm1 = s2 = y2 = u2 = rep(0,10)
k1 = 10
for(k in 1:k1){
  for(i in 1:10){
    # simulation 1
    mv.NE <- mvdc(normalCopula(0.8), c("norm", "exp"), list(list(mean = 0, sd =2), list(rate = i)))
    # simulation 2
    # mv.NE <- mvdc(gumbelCopula(i), c("norm", "norm"), list(list(mean = 0, sd =2), list(mean = 0, sd =2)))
    data1 <- rMvdc(800, mv.NE)
    
    tce1[i] = tce1[i] + mvnt(data1) / k1
    mardia1[i] = mardia1[i] + as.numeric(as.character(mvn(data1,mvnTest = "mardia")$multivariateNormality$Statistic[1])) / k1
    hz1[i] = hz1[i] + mvn(data1,mvnTest = "hz")$multivariateNormality$HZ / k1
    royston1[i] = royston1[i] + mvn(data1,mvnTest = "royston")$multivariateNormality$H / k1
    dh1[i] = dh1[i] + mvn(data1,mvnTest = "dh")$multivariateNormality$E / k1
    energy1[i] = energy1[i] + mvn(data1,mvnTest = "energy")$multivariateNormality$Statistic / k1
    ad1[i] = ad1[i] + AD.test(data1)@AD / k1
    cm1[i] = cm1[i] + CM.test(data1)@CM / k1
    s2s = S2.test(data1)
    s2[i] = s2[i] + s2s@s2 / k1
    y2[i] = y2[i] + s2s@y2 / k1
    u2[i] = u2[i] + s2s@u2 / k1
  }
}

xlab1 = "rate" # simulation1
# xlab1 = "alpha" # simulation2
x11(width = 12, height = 9);
par(mfrow = c(3,4))
plot(tce1,xlab = xlab1,ylab="statistic", main = "Copula Entropy");lines(tce1)
plot(mardia1,xlab = xlab1,ylab="statistic", main = "Mardia");lines(mardia1)
plot(hz1,xlab = xlab1,ylab="statistic", main = "Royston");lines(hz1)
plot(royston1,xlab = xlab1,ylab="statistic", main = "Henze-Zirkler");lines(royston1)
plot(dh1,xlab = xlab1,ylab="statistic", main = "Dornik-Haansen");lines(dh1)
plot(energy1,xlab = xlab1,ylab="statistic", main = "Energy Distance");lines(energy1)
plot(ad1,xlab = xlab1,ylab="statistic", main = "Anderson-Darling");lines(ad1)
plot(cm1,xlab = xlab1,ylab="statistic", main = "Cramer-von Mises");lines(cm1)
plot(s2,xlab = xlab1,ylab="statistic", main = "McCulloch");lines(s2)
plot(y2,xlab = xlab1,ylab="statistic", main = "Nikulin-Rao-Robson");lines(y2)
plot(u2,xlab = xlab1,ylab="statistic", main = "Dzhaparidze-Nikulin");lines(u2)
