library(CARBayesST)
library(sp)
library(sf)
library(dplyr)
library(spdep)
library(lme4)
library(lmtest)

NYC <- st_read("MODZCTA_2010.shx")
zipdata <- read.csv('COVIDdataNYC.csv', header = T)
scaled.dat <- zipdata
scaled.dat[,c(10:21)] <- data.frame(scale(zipdata[,c(10:21)]));#standardize data

NYC <- NYC[-178,]
NYC1 <- as(NYC, 'Spatial')

Trans.av <- summarise(group_by(zipdata,zip), Trans.mean = mean(localtrans))
#get adjacency matrix
W.nb <- poly2nb(NYC1, row.names = Trans.av$zip)
W.list <- nb2listw(W.nb, style = "B", zero.policy = TRUE)
W <- nb2mat(W.nb, style = "B", zero.policy = TRUE)

#CARar
#fix roosevelt island
W[38,19] <- 1
W[38,83] <- 1
W[38,88] <- 1
W[19,38] <- 1
W[83,38] <- 1
W[88,38] <- 1

################################### No lag for within-zip code transmission

formula <- localtrans ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                   data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                   thin = 20, AR = 2)
#replace ST.CARar with ST.CARanova if run using the ST.CARanova model

model$samples <- NULL # if you want to save all samples for posterior parameters, comment out this line; otherwise, samples will be dropped to save space
modellist <- list(model)

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtests)

################################### 1-week lag for within-zip code transmission

formula <- localtrans_1 ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                  data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                  thin = 20, AR = 2)
model$samples <- NULL
modellist <- append(modellist,list(model))

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(moranspvalue, morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtestpvalue,dwtests)

################################### 2-week lag for within-zip code transmission

formula <- localtrans_2 ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                  data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                  thin = 20, AR = 2)

model$samples <- NULL
modellist <- append(modellist,list(model))

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(moranspvalue, morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtestpvalue,dwtests)

################################### no lag for cross-zip code transmission

formula <- crosstrans ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                  data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                  thin = 20, AR = 2)
# model$samples <- NULL
modellist <- append(modellist,list(model))

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(moranspvalue, morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtestpvalue,dwtests)

################################### 1-week lag for cross-zip code transmission

formula <- crosstrans_1 ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                  data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                  thin = 20, AR = 2)
model$samples <- NULL
modellist <- append(modellist,list(model))

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(moranspvalue, morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtestpvalue,dwtests)

################################### 2-week lag for cross-zip code transmission

formula <- crosstrans_2 ~ offset(logpop) + logpopdensity + logcase + logtest + logcumucaserate + black + hispanic + age65plus + householdincome + bachelor + householdsize + vaccoverage  + totalvisitor

model <- ST.CARar(formula = formula, family = "poisson",
                  data = scaled.dat, W = W, burnin = 20000, n.sample = 420000,
                  thin = 20, AR = 2)

model$samples <- NULL
modellist <- append(modellist,list(model))

sink(file = "modelsummary.txt", append = TRUE, type = "output")
print(model)
sink()

resid.CAR <- residuals(model)
#test spatial autocorrelation
morans <- rep(0,31)
for (t in 1:31) {
  moran <- moran.mc(x = resid.CAR[((t-1)*177+1):(t*177)], listw = W.list, nsim = 10000, zero.policy = TRUE)
  morans[t] <- moran$p.value
}
moranspvalue <- rbind(moranspvalue, morans)

#test temporal autocorrelation
dwtests <- rep(0,177)
for (t in 1:177) {
  y1 <- resid.CAR[seq(t,5487,177)]
  dw <- dwtest(y1 ~ 1)
  dwtests[t] <- dw$p.value
}

dwtestpvalue <- rbind(dwtestpvalue,dwtests)

#####################output
save(modellist,file="modellist.RData")
save(moranspvalue,file="moranspvalue.RData")
save(dwtestpvalue,file="dwtestpvalue.RData")
