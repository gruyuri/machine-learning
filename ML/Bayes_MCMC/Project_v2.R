setwd("D:\\Data\\Bayes_MCMC\\Data\\Marathon")

dat <- read.csv(file="napa_marathon_fm2015.csv", header=TRUE)
str(dat)
head(dat)

df <- dat[, c(1:3)]
df$isM <- as.logical(df$Gender == "M")
df$is_35 <- as.logical(df$Age <= 35)
df$is36_45 <- as.logical(df$Age > 35 & df$Age <= 45)
df$is46_55 <- as.logical(df$Age > 45 & df$Age <= 55)
df$is_60 <- as.logical(df$Age > 55)

df$ageGroup <- as.numeric(0)
df$ageGroup[df$is_35 == T] <- 1
df$ageGroup[df$is36_45 == T] <- 2
df$ageGroup[df$is46_55 == T] <- 3
df$ageGroup[df$is_60 == T] <- 4

df$sexGroup <- as.numeric(0)
df$sexGroup[df$isM == T] <- 1


str(df)
head(df)
tail(df)

idx1 <- df$isM == T & df$is_35 == T
dfM_1 <- df[idx1,]

idx2 <- df$isM == T & df$is36_45 == T
dfM_2 <- df[idx2,]

idx3 <- df$isM == T & df$is46_55 == T
dfM_3 <- df[idx3,]

idx4 <- df$isM == T & df$is_60 == T
dfM_4 <- df[idx4,]

idxM <- df$isM == T;
dfM <- df[idxM,]

hist(dfM_1$Hours)
hist(dfM_2$Hours)
hist(dfM_3$Hours)
hist(dfM_4$Hours)

# find dependency for linear regression
library(ggplot2)

attach(dfM_1)
plot(x = (df$Age - mean(df$Age))/sqrt(var(df$Age)) + 3, y = df$Hours)
plot(x = (df$Age - mean(df$Age))/sqrt(var(df$Age)) + 3, y = log(df$Hours)^.5)
detach(dfM_1)

plot(jitter(y) ~ jitter(normAge), data = dfM_1, xlab = "Norm Age", ylab="sqrt(log(Hours))")


df$normAge <- (df$Age - mean(df$Age))/sqrt(var(df$Age)) + 3
df$y <- log(df$Hours)^.5

lin.model <- lm(y ~ normAge, data = dfM_1)

summary(lin.model)
plot(lin.model, add = T)

df <- df[, c(-1)]

hist(df$Hours[df$isM == F])

####################################################

df_model <- df[, c("Hours", "normAge", "y", "ageGroup", "sexGroup", "Age")]

library(rjags)
mod_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[i], prec)
      mu[i] = b0 + b1*normAge[i]
    }
    
  b0 ~ dnorm(0.0, 1.0/1.0e6)
  b1 ~ dnorm(0.0, 1.0/1.0e6)

  prec ~ dgamma(1.0/2.0, 1.0*300/2.0)
  sig2 = 1.0 / prec
  sig = sqrt(sig2)
} "

set.seed(75);
data_jags <- as.list(dfM);
params <- c("b0", "b1");
mod <- jags.model(textConnection(mod_string), data = data_jags, n.chains = 3);
update(mod, 1000);
mod_sim <- coda.samples(model = mod, variable.names = params, n.iter = 5000);
mod_csim <- as.mcmc(do.call(rbind, mod_sim))

plot(mod_sim);
gelman.diag(mod_sim);
# autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

summary(mod_sim);

dic1 <- dic.samples(mod, n.iter = 1e4) ## 686

################################################
mod2_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[i], prec)
      mu[i] = b0 + b1*normAge[i] +b2*ageGroup[i] +b3*sexGroup[i]
    }
    
  b0 ~ dnorm(0.0, 1.0/1.0e6)
  b1 ~ dnorm(0.0, 1.0/1.0e6)
  b2 ~ dnorm(0.0, 1.0/1.0e6)
  b3 ~ dnorm(0.0, 1.0/1.0e6)

  prec ~ dgamma(1.0/2.0, 1.0*300.0/2.0)
  sig2 = 1.0 / prec
  sig = sqrt(sig2)
} "

set.seed(75);
data2_jags <- as.list(df_model);
params2 <- c("b0", "b1", "b2", "b3");
mod2 <- jags.model(textConnection(mod2_string), data = data2_jags, n.chains = 3);
update(mod2, 1000);
mod2_sim <- coda.samples(model = mod2, variable.names = params2, n.iter = 1e4);
mod2_csim <- as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim);
gelman.diag(mod2_sim);
# autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

summary(mod2_sim);

dic2 <- dic.samples(mod2, n.iter = 1e4) ## 132

################################################
mod3_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[ageGroup[i]], prec)
    }
    
    for (j in 1:4) {
      mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }

    prec ~ dgamma(1.0/2.0, 1.0*300.0/2.0)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

set.seed(75);
data3_jags <- as.list(df_model);
params3 <- c("mu", "sig");
mod3 <- jags.model(textConnection(mod3_string), data = data3_jags, n.chains = 3);
update(mod3, 1000);
mod3_sim <- coda.samples(model = mod3, variable.names = params3, n.iter = 1e4);
mod3_csim <- as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim);
gelman.diag(mod3_sim);
# autocorr.diag(mod2_sim)
autocorr.plot(mod3_sim)
effectiveSize(mod3_sim)

summary(mod3_sim);

dic3 <- dic.samples(mod3, n.iter = 1e4) ## 139

################################################
mod4_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[ageGroup[i]], prec[ageGroup[i]])
    }

    for (j in 1:4) {
      mu[j] ~ dnorm(0.0, 1.0/1.0e6)
      prec[j] ~ dgamma(1/2.0, 1.0*300.0/2.0)
    }

    sig = sqrt(1.0 / prec)
} "

set.seed(75);
data4_jags <- as.list(df_model);
params4 <- c("mu", "sig");
mod4 <- jags.model(textConnection(mod4_string), data = data4_jags, n.chains = 3);
update(mod4, 1000);
mod4_sim <- coda.samples(model = mod4, variable.names = params4, n.iter = 1e4);
mod4_csim <- as.mcmc(do.call(rbind, mod4_sim))

plot(mod4_sim);
gelman.diag(mod4_sim);
# autocorr.diag(mod4_sim)
autocorr.plot(mod4_sim)
effectiveSize(mod4_sim)

summary(mod4_sim);

dic4 <- dic.samples(mod4, n.iter = 1e4) ## 2533

################################################
mod5_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dnorm(mu[4*sexGroup[i]+ageGroup[i]], prec[4*sexGroup[i]+ageGroup[i]])
    }

    for (j in 1:8) {
      mu[j] ~ dnorm(0.0, 1.0/1.0e6)
      prec[j] ~ dgamma(1/2.0, 1.0*300.0/2.0)
    }

    sig = sqrt(1.0 / prec)
} "

set.seed(75);
data5_jags <- as.list(df_model);
params5 <- c("mu", "sig");
mod5 <- jags.model(textConnection(mod5_string), data = data5_jags, n.chains = 3);
update(mod5, 1000);
mod5_sim <- coda.samples(model = mod5, variable.names = params5, n.iter = 1e4);
mod5_csim <- as.mcmc(do.call(rbind, mod5_sim))

plot(mod5_sim);
gelman.diag(mod5_sim);
# autocorr.diag(mod5_sim)
autocorr.plot(mod5_sim)
effectiveSize(mod5_sim)

summary(mod5_sim);

dic5 <- dic.samples(mod5, n.iter = 1e4) ## 6814

################################################
mod6_string = " model {
    for (i in 1:length(y)) {
      y[i] ~ dnorm(mu[ageGroup[i]], prec[sexGroup[i] + 1])
    }

    for (j in 1:4) {
      mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }

    for (j in 1:2) {
      prec[j] ~ dgamma(1/2.0, 1.0*300.0/2.0)
    }

    sig = sqrt(1.0 / prec)
} "

set.seed(75);
data6_jags <- as.list(df_model);
params6 <- c("mu", "sig");
mod6 <- jags.model(textConnection(mod6_string), data = data6_jags, n.chains = 3);
update(mod6, 1000);
mod6_sim <- coda.samples(model = mod6, variable.names = params6, n.iter = 1e4);
mod6_csim <- as.mcmc(do.call(rbind, mod6_sim))

plot(mod6_sim);
gelman.diag(mod6_sim);
# autocorr.diag(mod6_sim)
autocorr.plot(mod6_sim)
effectiveSize(mod6_sim)

summary(mod6_sim);

dic6 <- dic.samples(mod6, n.iter = 1e4) ## 1382

################################################
mod7_string = " model {
    for (i in 1:length(Hours)) {
      Hours[i] ~ dnorm(mu[Age[i]], prec)
    }

    for (j in 1:100) {
      mu[j] ~ dnorm(0.0, 1.0/1.0e6)
    }

    prec ~ dgamma(1.0/2.0, 1.0*300.0/2.0)
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
} "

set.seed(75);
data7_jags <- as.list(df_model);
params7 <- c("mu", "sig");
mod7 <- jags.model(textConnection(mod7_string), data = data7_jags, n.chains = 3);
update(mod7, 1000);
mod7_sim <- coda.samples(model = mod7, variable.names = params7, n.iter = 1e4);
mod7_csim <- as.mcmc(do.call(rbind, mod7_sim))

plot(mod7_sim);
gelman.diag(mod7_sim);
# autocorr.diag(mod7_sim)
autocorr.plot(mod7_sim)
effectiveSize(mod7_sim)

summary(mod7_sim);

dic7 <- dic.samples(mod7, n.iter = 1e4) ## 4389

