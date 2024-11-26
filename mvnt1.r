# Evaluation of mvnt with ce
library(copula)
library(mvtnorm)
library(copent)
library(MVN)
library(mvnTest)
library(mnt)
library(mvnormtest)
library(latex2exp)

tce1 = rep(0,10) #copent
mardia1 = hz1 = royston1 = dh1 = energy1 = rep(0,10) #MVN
ad1 = cm1 = s2 = y2 = u2 = rep(0,10) #mvnTest
#mnt
cs1 = bhep1 = deht1 = dehu1 = ehs1 = hjg1 = hv1 = hz1 = kkurt1 = makurt1 = maskew1 = mkurt1 = mq1 = mq2 = mrsskew1 = mskew1 = pu1 = sr1 = rep(0,10)
# hjm1 = rep(0,10)
mshapiro1 = rep(0,10) #mvnormtest

# for simulation 3
mu1 = c(0,0); sigma1 = matrix(c(1,0.5,0.5,1),2,2)
mu2 = c(3,3); sigma2 = matrix(c(1,0.8,0.8,1),2,2)
d1 = rmvnorm(800,mu1,sigma1)  
d2 = rmvnorm(800,mu2,sigma2)

for(i in 1:10) {
  # simulation 1
  # mv.NE <- mvdc(normalCopula(0.8), c("norm", "exp"), list(list(mean = 0, sd = 2), list(rate = i)))
  data1 <- rMvdc(800, mv.NE)
  # simulation 2
  ## mv.NE <- mvdc(gumbelCopula(i), c("norm", "norm"), list(list(mean = 0, sd =2), list(mean = 0, sd =2)))
  ## data1 <- rMvdc(800, mv.NE)
  data1 <- rmvt(800,sigma = matrix(c(1,0.5,0.5,1),2,2),df = i)
  # simulation 3
  data1 = d1 * (i-1) / 9 + d2 * (10 -i) / 9
  
  tce1[i] = mvnt(data1)
  mardia1[i] = as.numeric(as.character(
    mvn(data1, mvnTest = "mardia")$multivariateNormality$Statistic[1]
  ))
  hz1[i] = mvn(data1, mvnTest = "hz")$multivariateNormality$HZ
  royston1[i] = mvn(data1, mvnTest = "royston")$multivariateNormality$H
  dh1[i] = mvn(data1, mvnTest = "dh")$multivariateNormality$E
  energy1[i] = mvn(data1, mvnTest = "energy")$multivariateNormality$Statistic
  ad1[i] = AD.test(data1)@AD
  cm1[i] = CM.test(data1)@CM
  s2s = S2.test(data1)
  s2[i] = s2s@s2
  y2[i] = s2s@y2
  u2[i] = s2s@u2
  cs1[i] = CS(data1)
  bhep1[i] = BHEP(data1)
  deht1[i] = DEHT(data1)
  dehu1[i] = DEHU(data1, a = 1)
  ehs1[i] = EHS(data1)
  hjg1[i] = HJG(data1)
  #hjm1[i] = HJM(data1,a=1)
  hv1[i] = HV(data1)
  hz1[i] = HZ(data1)
  kkurt1[i] = KKurt(data1)
  makurt1[i] = MAKurt(data1)
  maskew1[i] = MASkew(data1)
  mkurt1[i] = MKurt(data1)
  mq1[i] = MQ1(data1)
  mq2[i] = MQ2(data1)
  mrsskew1[i] = MRSSkew(data1)
  mskew1[i] = MSkew(data1)
  pu1[i] = PU(data1)
  sr1[i] = SR(data1) # comment this line in simulation 2
  mshapiro1[i] = mshapiro.test(t(data1))$statistic
}#i

# xlab1 = TeX(r'($\lambda$)') # simulation 1
xlab1 = TeX(r'($\nu$)') # simulation 2
# xlab1 = TeX(r'($\beta$)') # simulation 3
x11(width = 8, height = 12);
par(mfrow = c(6,5))
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
plot(bhep1,xlab = xlab1,ylab="statistic", main = "BHEP");lines(bhep1)
plot(cs1,xlab = xlab1,ylab="statistic", main = "Cox-Small");lines(cs1)
plot(deht1,xlab = xlab1,ylab="statistic", main = "DEHT");lines(deht1)
plot(dehu1,xlab = xlab1,ylab="statistic", main = "DEHU");lines(dehu1)
plot(ehs1,xlab = xlab1,ylab="statistic", main = "EHS");lines(ehs1)
plot(hjg1,xlab = xlab1,ylab="statistic", main = "HJG");lines(hjg1)
# plot(hjm1,xlab = xlab1,ylab="statistic", main = "HJM");lines(hjm1)
plot(hv1,xlab = xlab1,ylab="statistic", main = "HV");lines(hv1)
plot(hz1,xlab = xlab1,ylab="statistic", main = "HZ");lines(hz1)
plot(kkurt1,xlab = xlab1,ylab="statistic", main = "KKurt");lines(kkurt1)
plot(makurt1,xlab = xlab1,ylab="statistic", main = "MAKurt");lines(makurt1)
plot(maskew1,xlab = xlab1,ylab="statistic", main = "MASkew");lines(maskew1)
plot(mkurt1,xlab = xlab1,ylab="statistic", main = "MKurt");lines(mkurt1)
plot(mq1,xlab = xlab1,ylab="statistic", main = "MQ1");lines(mq1)
plot(mq2,xlab = xlab1,ylab="statistic", main = "MQ2");lines(mq2)
plot(mrsskew1,xlab = xlab1,ylab="statistic", main = "MRSSkew");lines(mrsskew1)
plot(mskew1,xlab = xlab1,ylab="statistic", main = "MSkew");lines(mskew1)
plot(pu1,xlab = xlab1,ylab="statistic", main = "PU");lines(pu1)
plot(sr1,xlab = xlab1,ylab="statistic", main = "SR");lines(sr1)
plot(mshapiro1,xlab = xlab1,ylab="statistic", main = "Shapiro-Wilk");lines(mshapiro1)
