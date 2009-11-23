library(R2WinBUGS)
library(rjags)
library(foreign)
library(reshape)
library(ggplot2)
source("utils.R")


## read data and subset
dnow <- read.dta("dataforanalysis2.dta")
dnow <- subset(dnow, country!=4, select=c(caseno, party, nptys, seats, portshare, ptyprop, ppunwtd))

dnow$partyR <- with(dnow, ave(caseno, caseno, FUN=function(x) 1:length(x)))

## split it by number of parties
## we need to do this in order to account for the structural zeroes
dnow.split <- split(dnow, dnow$nptys)

## now reshape, creating matrices with portfolio share and seat share
caseno <- lapply(dnow.split, function(x) as.matrix(recast(x, caseno ~ partyR, measure.var="ppunwtd")[,1]))
pshare <- lapply(dnow.split, function(x) as.matrix(recast(x, caseno ~ partyR, measure.var="ppunwtd")[,-1]))
cshare <- lapply(dnow.split, function(x) log(as.matrix(recast(x, caseno ~ partyR, measure.var="ptyprop")[,-1])))
names(pshare) <- paste("p", names(pshare), sep='')
names(cshare) <- paste("z", names(cshare), sep='')
## number of obs in each split
N <- as.numeric(sapply(pshare, nrow))

## data for the bugs model
jags.data <- c(pshare,cshare, N2=N[1], N3=N[2], N4=N[3], N5=N[4], N6=N[5],
               n.parties2=2,n.parties3=3,n.parties4=4, n.parties5=5, n.parties6=6)





## bugs model
model.dirichlet <- function() {
    ## note that we code each split separately
    ## but the coefficients are shared
    for (i in 1:N2) {
        for (j in 1:n.parties2) {
            alpha.hat2[i,j] <- b[1]+b[2]*z2[i, j]
        }
        p2[i,  ] ~ ddirch(exp(alpha.hat2[i,]))
    }
    for (i in 1:N3) {
        for (j in 1:n.parties3) {
            alpha.hat3[i,j] <- b[1]+b[2]*z3[i, j]
        }
        p3[i,  ] ~ ddirch(exp(alpha.hat3[i,]))
    }
    for (i in 1:N4) {
        for (j in 1:n.parties4) {
            alpha.hat4[i,j] <- b[1]+b[2]*z4[i, j]
        }
        p4[i,  ] ~ ddirch(exp(alpha.hat4[i,]))
    }
    for (i in 1:N5) {
        for (j in 1:n.parties5) {
            alpha.hat5[i,j] <- b[1]+b[2]*z5[i, j]
        }
        p5[i,  ] ~ ddirch(exp(alpha.hat5[i,]))
    }
    for (i in 1:N6) {
        for (j in 1:n.parties6) {
            alpha.hat6[i,j] <- b[1]+b[2]*z6[i, j]
        }
        p6[i,  ] ~ ddirch(exp(alpha.hat6[i,]))
    }
    ## priors
    for (i in 1:2) {
        b[i] ~ dnorm(0, .0001)
    }
}

write.model(model.dirichlet, con="model.bug")
parameters.to.save <- c("b", paste(c("alpha.hat"), 2:6, sep=''))


n.chains <- 2
n.iter <- 1000
n.burnin <- 1000
n.thin <- 1

## compile the model and run burnin
jags1 <- jags.model(file="model.bug", data=jags.data,
                    n.adapt=n.burnin,n.chains=n.chains)

## get the samples
jags1.s <- coda.samples(jags1,variable.names=parameters.to.save, n.iter=n.iter,thin=n.thin)

## coda2bugs is a function to translate mcmc objects to
## (nicer) gelman's bugs objects
bugs <- coda2bugs(jags1.s)

## function to merge the predicted values together
getalpha <- function() {
    res <- NULL
    for (x in 2:6)  {
        alpha <- exp(bugs$mean[[paste("alpha.hat",x,sep='')]])
        df <- data.frame(caseno=caseno[[paste(x)]],
                         alpha/apply(alpha,1,sum))
        df <- melt(df, id.var="caseno")
        df$partyR <- as.numeric(gsub("X", "", df$variable))
        df$variable <- NULL
        if (is.null(res)) {
            res <- df
        } else {
            res <- merge(res, df, all=TRUE)
        }
    }
}
res <- getalpha()
res <- merge(res, dnow)

## the dirichlet model is slightly better
with(subset(res), cor(cbind(value, ppunwtd, ptyprop)))^2

with(res, qplot(value, ppunwtd, colour=nptys))+geom_abline()

with(res, qplot(ptyprop, value, colour=nptys))+geom_abline()


















tmp <- lapply(c("rjagsCluster.R","utils.R","mcmc-cluster.R"),
              function(x) try(source(x)))








## samples from a gamsons-like dist
sim <- function() {
    K <- 10 ## number of countries
    P <- 5 ## number of parties per country
    T <- 5 ## number of periods
    df <- expand.grid(group=1:K, period=1:T, party=1:P)
    N <- nrow(df)    
    df$z <- rnorm(N)
    df$e <- rnorm(N,0,1)
    df$delta <- with(df, z+e)
    df$p <- with(df, ave(exp(delta), group, period, FUN=function(x) x/sum(x)))
    df
}

s1 <- sim()
##s1$party.f <- with(s1, factor(paste(group, party, sep=";")))
##s1$party.i <- as.numeric(s1$party.f)
##s1party <- recast(s1, group+period ~ party + variable, measure.var=c("party.i"), id.var=c("group", "party", "period"))
s1p <- recast(s1, group+period ~ party + variable, measure.var=c("p"), id.var=c("group", "party", "period"))
s1z <- recast(s1, group+period ~ party + variable, measure.var=c("z"), id.var=c("group", "party", "period"))
n.b <- 6 ## number of covariates
jags.data <- list(p=as.matrix(s1p[,grepl("_p$", colnames(s1p))]),
                  z=as.matrix(s1z[,grepl("_z$", colnames(s1z))]),
                  ##party.i=as.matrix(s1party[,grepl("_party.i$", colnames(s1party))]),
                  ##n.countries=max(s1$group),
                  ##group=s1$group,
                  N=nrow(s1z),
                  n.parties=max(s1$party),
                  n.b=n.b                  
                  )



write.model(model.dirichlet, con="model.bug")
parameters.to.save <- c("b", "alpha.hat" )
n.chains <- 2
n.iter <- 1000
n.burnin <- 1000
n.thin <- 1

jags1 <- jags.model(file="model.bug", data=jags.data,
                    inits=function() {
                        list(b=rnorm(n.b))
                    },
                    n.adapt=n.burnin,n.chains=n.chains)

jags1.s <- coda.samples(jags1,variable.names=parameters.to.save, n.iter=n.iter,thin=n.thin)

tmp <- lapply(c("rjagsCluster.R","utils.R","mcmc-cluster.R"),
              function(x) try(source(x)))











## samples from a dirichlet
library(VGAM)
p <- rdiric(10, c(1, 2, 3), dimension = NULL)


N <- nrow(p)
n.parties <- ncol(p)
n.b <- 2
jags.data <- list(p=p, n.b=n.b, N=N, n.parties=n.parties, z=matrix(rnorm(n.parties*N), ncol=n.parties))

write.model(model.dirichlet, con="model.bug")
parameters.to.save <- c("b", "b.out")

jags1 <- jags.model(file="model.bug", data=jags.data,
                    inits=function() {
                        list(b=rnorm(n.b))
                    },
                    n.adapt=n.burnin,n.chains=n.chains)

jags1.s <- coda.samples(jags1,variable.names=parameters.to.save, n.iter=n.iter,thin=n.thin)




library(R2WinBUGS)
library(rjags)
n.chains <- 2
n.burnin <- 1000
n.iter <- 1000
n.thin <- 1
write.model(model.dirichlet, con="model.bug")
parameters.to.save <- c("b", "b.out")
jags1 <- jags.model(file="model.bug", data=jags.data,
                    inits=function() {
                        list(b=rnorm(n.b))
                    },
                    n.adapt=n.burnin,n.chains=n.chains)
jags1.s <- coda.samples(jags1,variable.names=parameters.to.save, n.iter=n.iter,thin=n.thin)














## dirichlet - gamma parameterization
model.gamma2 <- function() {
    for (i in 1:N) {
        for (j in 1:n.parties) {
            p[i, j] <- delta[i, j] / sum(delta[i,])
            delta[i, j] ~ dgamma(exp(b[1]), 1)
        }
    }
    ## priors
    for (i in 1:2) {
        b[i] ~ dnorm(0, .001)
    }
    ## for (j in 1:n.parties) {
    ##     alpha.party[j] <- b[1]
    ## }
}

## normal
model.normal <- function() {
    for (i in 1:N) {
        for (j in 1:n.parties) {
            p[i, j] ~ dnorm(b[1]+b[2]*z[i, j], tau.p)
        }
    }
    ## priors
    for (i in 1:2) {
        b[i] ~ dnorm(0, .001)
    }
    tau.p ~ dgamma(.01, .01)
}

sim <- function() {
    K <- 20 ## number of countries
    P <- 3 ## number of parties per country
    T <- 10 ## number of periods
    df <- expand.grid(group=1:K, period=1:T, party=1:P)
    N <- nrow(df)    
    df$z <- rnorm(N)
    df$e <- rnorm(N,0,1)
    ##df$delta <- with(df, z+e)
    df$delta <- with(df, 1+e)
    df$p <- with(df, ave(exp(delta), group, period, FUN=function(x) x/sum(x)))
    df
}

s1 <- sim()
s1$party.f <- with(s1, factor(paste(group, party, sep=";")))
s1$party.i <- as.numeric(s1$party.f)


library(reshape)
s1p <- recast(s1, group+period ~ party + variable, measure.var=c("p"), id.var=c("group", "party", "period"))
s1z <- recast(s1, group+period ~ party + variable, measure.var=c("z"), id.var=c("group", "party", "period"))
s1party <- recast(s1, group+period ~ party + variable, measure.var=c("party.i"), id.var=c("group", "party", "period"))


jags.data <- list(p=as.matrix(s1p[,grepl("_p$", colnames(s1p))]),
                  z=as.matrix(s1z[,grepl("_z$", colnames(s1z))]),
                  party.i=as.matrix(s1party[,grepl("_party.i$", colnames(s1party))]),
                  N=nrow(s1z),
                  n.parties=max(s1$party),
                  n.countries=max(s1$group),
                  group=s1$group                  
                  )


library(R2WinBUGS)
library(rjags)
n.chains <- 2
n.burnin <- 1000
n.iter <- 10000
n.thin <- 10
##write.model(model.normal, con="model.bug")
#write.model(model.gamma, con="model.bug")
write.model(model.gamma2, con="model.bug")
parameters.to.save <- c("b")

system.time(jags1 <- jags.model(file="model.bug", data=jags.data,
                                inits=function() {
                                    list(b=rnorm(2))
                                },
                                n.adapt=n.burnin,n.chains=n.chains))

system.time(jags1.s <- coda.samples(jags1,variable.names=parameters.to.save, n.iter=n.iter,thin=n.thin))


                  










sim <- function() {
    N <- 100
    K <- 10 ## number of countries
    P <- 5 ## number of parties per country
    T <- 3 ## number of periods
    group <- rep(1:K, length=N)
    party <- rep(1:P, length=N)
    period <- rep(1:T, length=N)
    z <- rnorm(N)
    e <- rnorm(N,0,1)
    delta <- z+e
    p <- ave(exp(delta), group, period, FUN=function(x) x/sum(x))
    summary(lm(log(p) ~ z))
}
res <- replicate(1000, coef(sim())[2,1])
s1 <- summary(res)

sim2 <- function() {
    N <- 100
    K <- 10 ## number of countries
    P <- 5 ## number of parties per country
    T <- 3 ## number of periods
    group <- rep(1:K, length=N)
    party <- rep(1:P, length=N)
    period <- rep(1:T, length=N)
    z <- rnorm(N)
    e <- rnorm(N,0,1)
    delta <- z+e
    p <- ave(exp(delta), group, period, FUN=function(x) x/sum(x))
    summary(lm(p ~ z))
}
res2 <- replicate(1000, coef(sim2())[2,1])
s2 <- summary(res2)

sim3 <- function() {
    N <- 100
    K <- 10 ## number of countries
    P <- 5 ## number of parties per country
    T <- 3 ## number of periods
    group <- rep(1:K, length=N)
    party <- rep(1:P, length=N)
    period <- rep(1:T, length=N)
    z <- rnorm(N)
    e <- rnorm(N,0,1)
    delta <- z+e
    p <- ave(exp(delta), group, period, FUN=function(x) x/sum(x))
    summary(lm(log(p[group!=1]) ~ z[group!=1]))
}
res3 <- replicate(1000, coef(sim3())[2,1])
s3 <- summary(res3)
















### logs 

sim <- function(N=1000, ## n. obs
                K=100 ## n. clusters
                ) {
    require(reshape)
    group <- sample(1:K, N, replace=TRUE)
    log.seats <- rnorm(N)
    log.shares <- log.seats+rnorm(N,0,1)
    seats <- ave(log.seats, group, FUN=function(x) exp(x)/(1+sum(exp(x))))
    shares <- ave(log.shares, group, FUN=function(x) exp(x)/(1+sum(exp(x))))
    df <- data.frame(seats,shares,log.seats, log.shares, group)
    df0 <- recast(df, group~variable, measure.var=c("seats", "shares"), fun.aggregate=function(x) 1-sum(x))
    df <- rbind(df, data.frame(df0,
                               log.seats=NA,
                               log.shares=NA))
    df$nparties <- with(df, ave(group, group, FUN=length))
    df <- df[order(df$group, df$shares),]
    df
}

df <- sim(N=15000, K=5000)

j <- 2
df$log.seats2 <- with(df, ave(seats, group, FUN=function(x) log(x/x[j])))
df$log.shares2 <- with(df, ave(shares, group, FUN=function(x) log(x/x[j])))

df0 <- subset(df, !is.na(log.seats))
## DGP
m1 <- summary(lm(log.shares~log.seats, data=df0))
coef(m1)[,1:3]
## ARBITRARY LOG SCALE
m2 <- summary(lm(log.shares2~log.seats2, data=subset(df, log.seats2!=0)))
coef(m2)[,1:3]
## LM, DROP base
m3 <- summary(lm(shares~seats, data=df0))
coef(m3)[,1:3]
## LM, DROP LARGEST
df$largest <- with(df, seats==ave(seats, group, FUN=max))
m4 <- summary(lm(shares~seats, data=subset(df, !largest)))
coef(m4)[,1:3]
## LM, ALL
m5 <- summary(lm(shares~seats, data=df))
coef(m5)[,1:3]

m4










## leg seats
n <- 100
## portfolios
m <- 10
## prop of seats
s <- runif(k)
s <- s/sum(s)




## prop of porfolio
p <- s+runif(k,-.1,.1)
p <- ifelse(p>0, p, 0)
p <- p/sum(p)
## draw votes
seats <- sample(1:k, n, replace=TRUE, prob=s)
## draw port shares
ports <- sample(1:k, m, replace=TRUE, prob=p)



## summaries
sseats <- as.numeric(table(factor(seats, levels=1:k)))
pports <- as.numeric(table(factor(ports, levels=1:k)))
