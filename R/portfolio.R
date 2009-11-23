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


