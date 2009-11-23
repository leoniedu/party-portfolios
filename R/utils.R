clean<-function(x){
    #CLEANS ACCENTS AND OTHER MARKS FROM FIELD
    y<-toupper(x)
    y<-gsub("Â","A", y) 
    y<-gsub("Á","A", y)
    y<-gsub("Ã","A", y)
    y<-gsub("É","E", y)
    y<-gsub("Ê","E", y)
    y<-gsub("Í","I", y)
    y<-gsub("Ó","O", y)
    y<-gsub("Ô","O", y)
    y<-gsub("Õ","O", y)
    y<-gsub("Ú","U", y)
    y<-gsub("Ü","U", y)
    y<-gsub("Ç","C", y)
    y<-gsub("'"," ", y)
    y<-gsub("."," ", y, fixed=TRUE)  
    y<-gsub("-"," ", y, fixed=TRUE)    
    return(y)
  }

clean.df <- function(df) {
    foo <- function(x) {
        if ("character"%in%class(x)) {
            clean(x)
        } else if ("factor"%in%class(x)) {
            levels(x) <- clean(levels(x))
            x
        } else {
            x
        }
    }
    data.frame(lapply(df,foo))
}


## twitter poster for R
## allows delays
##based on http://www.cerebralmastication.com/?p=274
tweet.now <- function(status,userpwd="leoniedu:e321109",shorten=NULL) {
  shorturl <- function(x) {
    require(rjson)
    require(RCurl)
    res <- res <- getURL(paste("http://api.bit.ly/shorten?version=2.0.1&longUrl=",x,"&login=leoniedu&apiKey=R_75ef8862cba341b64f0d7ad92a02ad66&history=1&format=json",sep=''))
    fromJSON(res)$results[[1]]$shortUrl
  }
  if (!is.null(shorten)) {
    print(status)    
    status <- gsub("#1",shorturl(shorten),status)
    print(status)
  }
  if (nchar(status)>140) stop("status is too long!")
  require(tcltk)
  require("RCurl")
  opts <- curlOptions(header = FALSE, userpwd = userpwd, netrc = FALSE)
  method <- "http://twitter.com/statuses/update.xml?status="
  encoded_status <- URLencode(status)
  request <- paste(method,encoded_status,sep = "")
  postForm(request,.opts = opts)
  print(paste("You just twitted: ",status))
}

tweet <- function(status,wait=1,shorten=NULL,...) {
  require(tcltk)
  ## only one tweet permitted
  if (exists(".id")&&(!is.null(.id))) stop("another tweet is in the queue!")
  ##wait is in minutes
  ## convert to ms
  waitS <- wait*60*1000
  z <- function(...) {
    tweet.now(status,shorten=shorten,...)
    tweet.cancel()
  }
  .id <<- tcl("after", waitS, z)
  attr(.id,"status") <- status
}


tweet.cancel <- function() {
  tcl("after","cancel",.id)
  .id <<- NULL
}
tweet.info <- function()  {
  tcl("after","info",.id)
  print(attr(.id,"status"))
}


## summarise <- function(.data, ...) { 
##   eval(substitute(data.frame(...)), .data, parent.frame()) 
## }


##check if file is update
## file.updated <- function(x,y) {
##   ix <- file.info(x)
##   iy <- file.info(y)
##   !((ix$size==iy$size)&(ix$mtime==iy$mtime))
## }


wordwrap<-function(x,len,collapse="\n",...) paste(strwrap(x,width=len,...),collapse=collapse)

ls.size <- function(n=10) {
  x <- ls(name=".GlobalEnv")
  df <- NULL
  j <- 1
  for (i in x) {
    df <- rbind(df,data.frame(i,c(object.size(get(i)))))
    j <- j+1
  }
  names(df) <- c("object","size")
  df <- df[order(df$size,decreasing=TRUE),]
  df[1:n,]
}
##source("http://bioconductor.org/biocLite.R")
##biocLite()

ls.size <- function(n=10) {
  x <- ls(name=".GlobalEnv")
  df <- NULL
  j <- 1
  for (i in x) {
    df <- rbind(df,data.frame(i,c(object.size(get(i)))))
    j <- j+1
  }
  names(df) <- c("object","size")
  df <- df[order(df$size,decreasing=TRUE),]
  df[1:n,]
}
##source("http://bioconductor.org/biocLite.R")
##biocLite()

fastestmirrors <- function(country="USA") {
    ##http://cran.cnr.berkeley.edu/bin/windows/base/NEWS.R-2.7.2
    ##http://cran.stat.ucla.edu/bin/windows/base/rw-FAQ.html
    tmp <- url("http://cran.r-project.org/mirrors.html")
    mirrors <- readLines(tmp)
    close(tmp)
    ## mirror list between "here" (hyperlinked) and "Many of these"
    l1 <- grep("<a href=\"mirmon_report.html\">here</a>",mirrors,value=FALSE)
    l2 <- grep("Many of these",mirrors,value=FALSE)
    mirrors <- mirrors[(l1+1):(l2-1)]    
    ##country list (added the last row, so that the last country is properly matched)
    clist <- c(grep("<dt>",mirrors,value=FALSE),length(mirrors))
    cnow <- grep(country,mirrors,value=FALSE)
    cnext <- match(cnow,clist)+1
    mirrors <- mirrors[cnow:clist[cnext]]
    mirrors <- grep("http",mirrors[grep("href",mirrors,value=FALSE)+1],
                    value=TRUE)
    ##urls <- paste(mirrors,"src/contrib/",package,sep="")
    urls <- paste(mirrors,"bin/windows/base/rw-FAQ.html",sep="")
    times <- sapply(urls,function(x) {
        tm <- try(system.time(download.file(x,"tmp"))[3])
        ifelse(class(tm)=="numeric",tm,NA)
    }
                    )
    tmp <- data.frame(mirror=mirrors,time=round(times,1))[order(times),]
    rownames(tmp) <- NULL
    tmp
}

t.rollcall <- function(x) {
    ## "transpose" roll call object: legs to votes; votes to legs
    newx <- x
    newx$votes <- t(x$votes)
    newx$n <- x$m
    newx$m <- x$n
    newx$legis.data <- x$vote.data
    newx$vote.data <- x$legis.data
    newx
}


tracex2 <- function (object, legis = NULL, d = 1, conf.int = 0.95, showAll = FALSE, 
                     burnin = NULL, span = 0.25, legendLoc = "topright",...) 
{
    warnOption <- options()$warn
    options(warn = -1)
    if (class(object) != "ideal") 
        stop("object passed to function tracex must be of class ideal")
    if (is.null(d)) 
        stop("default value must be supplied for dimension to trace\n")
    if (length(d) > 2) 
        stop("tracex only works with up to 2 dimensions\n")
    if (!is.character(legis)) 
        stop("legis must be character (names of legislators)\n")
    Rv <- as.numeric(version$major) + 0.1 * as.numeric(version$minor)
    if (Rv >= 2.4) 
        old.par <- par(no.readonly = TRUE)
    else old.par <- par()
    legis.names <- as.vector(dimnames(object$x)[[2]])
    nLegis <- length(legis)
    p <- list()
    for (j in 1:nLegis) {
        p[[j]] <- grep(pattern = paste("^", legis[j], sep = ""), 
            x = legis.names)
        if (length(p[[j]]) == 0) {
            cat(paste("could not find legislator", legis[j], 
                "\n"))
            p[[j]] <- NA
        }
        if (!is.null(p[[j]]) & length(p[[j]]) > object$d) {
            cat(paste("no unique match for legislator", legis[j], 
                "\n"))
            cat("try providing more of the unique identifier for the legislator\n")
            p[[j]] <- NA
        }
    }
    plotName <- rep(NA, length(p))
    for (j in 1:nLegis) {
        if (length(d) == 1) {
            if (!is.na(p[[j]])) {
                p[[j]] <- p[[j]][1]
            }
        }
        if (!is.na(p[[j]])) {
            foo <- pmatch(x = legis[j], table = as.vector(dimnames(object$xbar)[[1]]))
            if (!is.na(foo)) 
                plotName[j] <- dimnames(object$xbar)[[1]][foo]
            cat(paste("matching", legis[j], "with", plotName[j], 
                "\n"))
        }
    }
    p <- p[!is.na(p)]
    plotName <- plotName[!is.na(plotName)]
    names(p) <- plotName
    nLegis <- length(p)
    if (is.null(burnin)) 
        keep <- pscl:::checkBurnIn(object, eval(object$call$burnin))
    else keep <- pscl:::checkBurnIn(object, burnin)
    start <- object$x[keep, 1][1]
    if (length(d) == 1) {
        options(warn = 0)
        checkD(object, d)
        if (span <= 0 | span >= 1) 
            stop("span must be between 0 and 1")
        if ((conf.int <= 0) || (conf.int >= 1)) 
            stop("conf.int must be between 0 and 1")
        rw <- 3
        if (nLegis < 3) 
            rw <- nLegis
        count <- 0
        if (rw < 4 & dev.interactive()) 
            par(mfrow = c(rw, 1))
        if (length(legendLoc) == 1) 
            legendLoc <- rep(legendLoc, nLegis)
        for (i in 1:nLegis) {
            meat <- object$x[keep, p[[i]] + d - 1]
            iter <- object$x[keep, 1]
            par(mar = c(4, 4, 4, 2) + 0.1)
            mainText <- plotName[i]
            if (object$d > 1) 
                mainText <- paste(mainText, ", Dimension ", d, 
                  sep = "")
            plot(y = meat, x = iter, type = "l", xlab = "Iteration", 
                ylab = "", main = mainText,...)
            runmean <- cumsum(meat)/1:length(iter)
            lines(iter, runmean, col = "red", lwd = 3)
            lf <- loess(meat ~ iter, span = span)
            lines(iter, predict(lf), col = "blue", lwd = 3)
            xbar <- mean(meat)
            q <- c((1 - conf.int)/2, 1 - ((1 - conf.int)/2))
            q <- quantile(meat, q)
            abline(h = xbar, lwd = 3, col = "grey")
            abline(h = q[1], lty = 2, col = "grey")
            abline(h = q[2], lty = 2, col = "grey")
            count <- count + 1
            if (!is.null(legendLoc[i])) 
                legend(x = legendLoc[i], bg = "white", ncol = 1, 
                  legend = c("Trace", "Cumulative Mean", paste("Moving Average (loess, span=", 
                    round(span, 2), ")", sep = ""), "Posterior Mean", 
                    paste(round(100 * conf.int), "% Confidence Interval", 
                      sep = "")), lty = c(1, 1, 1, 1, 2), lwd = c(2, 
                    3, 3, 3, 2), col = c("black", "red", "blue", 
                    "grey", "grey"), yjust = 0, cex = 0.65)
            if ((count == 3) & (nLegis > 3) & (dev.interactive())) {
                count <- 0
                readline("Press return/enter to see next set of plots: ")
            }
        }
    }
    if (length(d) == 2) {
        goodD <- d %in% (1:object$d)
        if (!all(goodD)) 
            stop("invalid dimensions requested in tracex")
        col <- rainbow(nLegis)
        meat <- list()
        for (i in 1:nLegis) {
            xTraces <- object$x[keep, p[[i]][1]]
            yTraces <- object$x[keep, p[[i]][2]]
            meat[[i]] <- list(x = xTraces, y = yTraces, col = col[i])
        }
        if (showAll) {
            xRange <- range(unlist(lapply(meat, function(x) x$x)), 
                na.rm = TRUE)
            yRange <- range(unlist(lapply(meat, function(x) x$y)), 
                na.rm = TRUE)
            require(graphics)
            layout(mat = matrix(c(1, 2), 1, 2, byrow = TRUE), 
                width = c(0.7, 0.3))
            par(mar = c(4, 4, 1, 1))
            plot(x = xRange, y = yRange, type = "n", axes = FALSE, 
                xlab = paste("Dimension", d[1]), ylab = paste("Dimensions", 
                  d[2]),...)
            axis(1)
            axis(2)
            lineFunc <- function(obj) {
                points(obj$x[1], obj$y[1], pch = 16, col = obj$col)
                lines(obj$x, obj$y, col = obj$col)
            }
            lapply(meat, lineFunc)
            mtext(side = 3, outer = FALSE, line = -0.5, cex = 0.75, 
                paste("Two-dimensional trace plots, MCMC iterations,\n", 
                  eval(object$call$object)$desc, ", Iterations ", 
                  start, " to ", object$call$maxiter, " thinned by ", 
                  object$call$thin, sep = ""))
            par(mar = c(3, 0, 1, 0))
            plot(x = c(0, 1), y = c(0.5, nLegis + 0.5), xlab = "", 
                ylab = "", xaxs = "i", yaxs = "i", axes = F, 
                type = "n",...)
            for (i in 1:nLegis) {
                lines(x = c(0, 0.15), y = rep(i, 2), lwd = 2, 
                  col = col[i])
                text(x = 0.25, y = i, cex = 0.75, plotName[i], 
                  adj = 0)
            }
        }
        if (!showAll) {
            par(mfrow = c(2, 2))
            count <- 0
            for (i in 1:nLegis) {
                plot(x = meat[[i]]$x, y = meat[[i]]$y, type = "l", 
                  xlab = paste("Ideal Point, Dimension ", d[1], 
                    sep = ""), ylab = paste("Ideal Point, Dimension ", 
                    d[2], sep = ""),...)
                title(plotName[i])
                count <- count + 1
            }
            if (count == 3 & dev.interactive()) {
                count <- 0
                readline("Press any key to see next set of plots:  ")
            }
        }
    }
    par(old.par)
    options(warn = warnOption)
    invisible(NULL)
}

jos.ls <- function(envir = as.environment(-1)) {
    ##http://jaeweb.cantr.net/2007/04/29/a-nicer-ls-for-r/
    df <- NULL
    names <- .Internal(ls(envir, all.names=T))
    for (item in names) {
        l1 <- length(get(item))
        l2 <- 0
        if (!is.null(dim(get(item)))) {
            l1 <- dim(get(item))[1]
            l2 <- dim(get(item))[2]
        }
        tm <- data.frame(item,
                         paste(class(get(item)),collapse=","),
                         mode(get(item)),
                         l1, l2
                         ,check.names = TRUE,
                         stringsAsFactors = FALSE)
        ltm <- length(tm) 
        ##print(tm)
        if (ltm>5) stop(print(item))
        df <- rbind(df,tm)
        ##         cat(sprintf("%-30s  %-10s  %-10s  %10d %10sn",
        ##                     item, class(get(item)), mode(get(item)), l1, l2))
    }
    rownames(df) <- NULL
    df <- data.frame(df)
    names(df) <- c("item","class","mode","dim1","dim2")
    df
}

rename.var <- function(x,...) {
    ## from memisc
    subst <- c(...)
    for (i in 1:length(subst)) {
        names(x)[names(x) == names(subst[i])] <- subst[i]
    }
    return(x)
}


ssource <- function(x,echo=TRUE) {
    Stangle(x)
    source(gsub("Rnw$","R",x),echo=echo)
}
toZ <- function(x) (x-mean(x))/sd(x)

dupall <- function(x) {
    duplicated(x) | duplicated(x,fromLast=TRUE)
}

pad0 <- function(x,mx=NULL,fill=0) {
  lx <- nchar(as.character(x))
  mx.calc <- max(lx,na.rm=TRUE)
  if (!is.null(mx)) {
    if (mx<mx.calc) {
      stop("number of maxchar is too small")
    }
  } else {
    mx <- mx.calc
  }
  px <- mx-lx
  paste(sapply(px,function(x) paste(rep(fill,x),collapse="")),x,sep="")
}


var.abbreviate <- function(x) {
    x <- trim(tolower(x))
    o <- c("/","\\(|\\)","percent","district","native american","attendance","composite","exceeds","region","enrollment","graduation"," +","-","%","&","\\.+")
    a <- c(" or ","_","pct","dist","ntv american","attend","comp","excd","reg","enrl","grad",".",".","pct","and",".")
    for (i in 1:length(o)) {
        x <- gsub(o[i],a[i],x)
}
    x
}

parMCMCirt1d <- function (..., n.chains=2, inits.fun=NULL) {
  post <- as.list(1:n.chains)
  for (i in 1:n.chains) {
    if (is.null(inits.fun)) {
      post[[i]] <- MCMCirt1d(...,seed=list(NA,i))
    } else {
      post[[i]] <- MCMCirt1d(...,seed=list(NA,i),
                             alpha.start=inits[[i]]$alpha,
                             beta.start=inits[[i]]$beta,
                             theta.start=inits[[i]]$theta
                             )
    }
    ## FIX: for some reason gamma has different names depending on the run.
    ##post[[i]] <- post[[i]][,-grep("^gamma",varnames(post[[i]]))]
  }
  mcmc.list(post)
}


"as.bugs.array2" <- function (sims.array, pts, DIC=FALSE,engine="bugs")
{
  ## 'sims.array' is supposed to be a 3-way array with
  # n.sims*n.chains*n.parameters simulations, and
  # the 3rd component of dimnames(x) should have the parameter names.

  ## From Andrew Gelman's bugs.r function
  ## a couple of lines commented out by Eduardo Leoni (see comment below)
  require("R2WinBUGS")
  d <- dim(sims.array)
  n.burnin     <- 0
  n.keep       <- d[1]
  n.chains     <- d[2]
  n.parameters <- d[3]
  n.sims       <- n.keep*n.chains
  n.iter       <- n.keep
  n.thin       <- 1
  #
  parameter.names <- dimnames(sims.array)[[3]]
  if (is.null(parameter.names)) {
    parameter.names <- paste("P", 1:n.parameters, sep="")
    dimnames(sims.array)[[3]] <- parameter.names
  }
  parameters.to.save <- unique(sapply(strsplit(parameter.names, "\\["), "[", 1))
  #
  sims <- matrix(NA, n.sims, n.parameters)
  root.long <- character(n.parameters)
  indexes.long <- vector(n.parameters, mode = "list")
  for (i in 1:n.parameters) {
    temp <- R2WinBUGS:::decode.parameter.name(parameter.names[i])
    root.long[i] <- temp$root
    indexes.long[[i]] <- temp$indexes
  }
  n.roots <- length(parameters.to.save)
  left.bracket.short <- as.vector(regexpr("[[]", parameters.to.save))
  right.bracket.short <- as.vector(regexpr("[]]", parameters.to.save))
  root.short <- ifelse(left.bracket.short == -1, parameters.to.save,
      substring(parameters.to.save, 1, left.bracket.short -
          1))
  dimension.short <- rep(0, n.roots)
  indexes.short <- vector(n.roots, mode = "list")
  n.indexes.short <- vector(n.roots, mode = "list")
  long.short <- vector(n.roots, mode = "list")
  length.short <- numeric(n.roots)
  for (j in 1:n.roots) {
      long.short[[j]] <- (1:n.parameters)[root.long == root.short[j]]
      length.short[j] <- length(long.short[[j]])
      if (length.short[j] == 0)
          stop(paste("parameter", root.short[[j]], "is not in the model"))
      else if (length.short[j] > 1) {
          dimension.short[j] <- length(indexes.long[[long.short[[j]][1]]])
          n.indexes.short[[j]] <- numeric(dimension.short[j])
          for (k in 1:dimension.short[j]) n.indexes.short[[j]][k] <- length(unique(unlist(lapply(indexes.long[long.short[[j]]],
              .subset, k))))
          length.short[j] <- prod(n.indexes.short[[j]])
          ## Modified by Eduardo Leoni
          ## this check fails if you take out a part of the simulations
          ## (for example, you don't want the array to have some of the
          ## parameters) so I took them out.

          ## if (length(long.short[[j]]) != length.short[j])
          ##   stop(paste("error in parameter", root.short[[j]],
          ##   "in parameters.to.save"))
          indexes.short[[j]] <- as.list(numeric(length.short[j]))
          for (k in 1:length.short[j]) indexes.short[[j]][[k]] <- indexes.long[[long.short[[j]][k]]]
      }
  }
  rank.long <- unlist(long.short)
  # -----
  # yes, it's inefficient to do this, but for now I'm just letting this be as it is:
  for (k in 1:n.parameters) {
    sims[,k] <- as.vector(sims.array[,,k])
  }
  # ----
  dimnames(sims) <- list(NULL, parameter.names)
  summary <- R2WinBUGS:::monitor(sims.array, n.chains, keep.all = TRUE)
  last.values <- as.list(numeric(n.chains))
  for (i in 1:n.chains) {
    n.roots.0 <- if (DIC)
      n.roots - 1
    else n.roots
    last.values[[i]] <- as.list(numeric(n.roots.0))
    names(last.values[[i]]) <- root.short[1:n.roots.0]
    for (j in 1:n.roots.0) {
      if (dimension.short[j] <= 1) {
        last.values[[i]][[j]] <- sims.array[n.keep, i,
                                            long.short[[j]]]
        names(last.values[[i]][[j]]) <- NULL
      }
      else if (engine=="jags") last.values[[i]][[j]] <- array(sims.array[n.keep,
                 i, long.short[[j]]], n.indexes.short[[j]])
      ## only winbugs have to permute the array.
      else last.values[[i]][[j]] <- aperm(array(sims.array[n.keep,
                                                           i, long.short[[j]]], rev(n.indexes.short[[j]])),
                                          dimension.short[j]:1)
    }
  }
  sims <- sims[sample(n.sims), ]
  sims.list <- summary.mean <- summary.sd <- summary.median <- summary.025 <-  summary.975 <- vector(n.roots,
                                                                      mode = "list")
  names(sims.list) <- names(summary.mean) <- names(summary.sd) <- names(summary.median) <- names(summary.025) <- names(summary.975) <- root.short
  for (j in 1:n.roots) {
    if (length.short[j] == 1) {
        sims.list[[j]] <- sims[, long.short[[j]]]
        summary.mean[[j]] <- summary[long.short[[j]], "mean"]
        summary.sd[[j]] <- summary[long.short[[j]], "sd"]
        summary.median[[j]] <- summary[long.short[[j]], "50%"]
        ##ell: added 025 and 975
        summary.025[[j]] <- summary[long.short[[j]], "2.5%"]
        summary.975[[j]] <- summary[long.short[[j]], "97.5%"]
    }
    else if (engine=="bugs") {
      temp2 <- dimension.short[j]:1
      sims.list[[j]] <- aperm(array(sims[, long.short[[j]]],
                                    c(n.sims, rev(n.indexes.short[[j]]))), c(1, (dimension.short[j] +
                                                                                 1):2))
      summary.mean[[j]] <- aperm(array(summary[long.short[[j]],
                                               "mean"], rev(n.indexes.short[[j]])), temp2)
      summary.sd[[j]] <- aperm(array(summary[long.short[[j]],
                                             "sd"], rev(n.indexes.short[[j]])), temp2)
      summary.median[[j]] <- aperm(array(summary[long.short[[j]],
                                                 "50%"], rev(n.indexes.short[[j]])), temp2)
      ##ell: added 025 and 975
      summary.025[[j]] <- aperm(array(summary[long.short[[j]],
                                              "2.5%"], rev(n.indexes.short[[j]])), temp2)
      summary.975[[j]] <- aperm(array(summary[long.short[[j]],
                                              "97.5%"], rev(n.indexes.short[[j]])), temp2)
    } else if (engine=="jags") {
      ##fix this list
      sims.list[[j]] <- sims[, long.short[[j]]]
      summary.mean[[j]] <- array(summary[long.short[[j]],"mean"],n.indexes.short[[j]])
      summary.sd[[j]] <- array(summary[long.short[[j]],"sd"],n.indexes.short[[j]])
      summary.median[[j]] <- array(summary[long.short[[j]],"50%"],n.indexes.short[[j]])
      ##ell: added 025 and 975
      summary.025[[j]] <- array(summary[long.short[[j]],"2.5%"],n.indexes.short[[j]])
      summary.975[[j]] <- array(summary[long.short[[j]],"97.5%"],n.indexes.short[[j]])
    }
  }
  summary <- summary[rank.long, ]
  all <- list(n.chains = n.chains, n.iter = n.iter, n.burnin = n.burnin,
        n.thin = n.thin, n.keep = n.keep, n.sims = n.sims, sims.array = sims.array[,
            , rank.long, drop = FALSE], sims.list = sims.list,
        sims.matrix = sims[, rank.long], summary = summary, mean = summary.mean,
        sd = summary.sd, median = summary.median, root.short = root.short,
        long.short = long.short, dimension.short = dimension.short,
        indexes.short = indexes.short, last.values = last.values, is.DIC=DIC,p02.5=summary.025,p97.5=summary.975)
    if (DIC) {
        deviance <- all$sims.array[, , dim(sims.array)[3], drop = FALSE]
        dim(deviance) <- dim(deviance)[1:2]
        pD <- numeric(n.chains)
        DIC <- numeric(n.chains)
        for (i in 1:n.chains) {
            pD[i] <- var(deviance[, i])/2
            DIC[i] <- mean(deviance[, i]) + pD[i]
        }
        all <- c(all, list(pD = mean(pD), DIC = mean(DIC)))
    }
  class(all) <- "bugs"
  return(all)
}


codadrop <- function(object.list,bugs=TRUE,parms.todrop=NA,parms.tokeep=NA,engine=NULL) {    
  n.chains <- length(object.list)
  if (n.chains<2) stop ("n.chains must be at least 2")
  parnames <- colnames(object.list[[1]])
  if (!is.na(parms.tokeep[1])) {
      c.tokeep <- NULL
      for (i in 1:length(parms.tokeep)) {
          c.tokeep <- c(c.tokeep,grep(parms.tokeep[i],parnames))
      }
      c.tokeep <- unique(c.tokeep)
  } else {
      c.tokeep <- 1:length(parnames)
  }
  res <- mcmc.list(lapply(object.list,function(x) x[,c.tokeep]))
  parnames <- colnames(object.list[[1]])[c.tokeep]
  c.todrop <- NULL
  if (!is.na(parms.todrop[1])) {
      for (i in 1:length(parms.todrop)) {
          c.todrop <- c(c.todrop,grep(parms.todrop[i],parnames))
      }
      c.todrop <- unique(c.todrop)
      res <- mcmc.list(lapply(res,function(x) x[,-c.todrop]))
  }
  res
}

"coda2bugs" <- function(object.list,bugs=TRUE,parms.todrop=NA,parms.tokeep=NA,engine=NULL) {
  ## try to get engine from object.list attribute
  if (is.null(engine)) engine <- attr(object.list,"engine")
  ## if empty assume it is jags
  if (is.null(engine)) engine <- "jags"
  ## One more piece borrowed from bugs.r
  ## This time from the mcsamp function
  ## this function gets a list of mcmc objects
  ## and transforms it to a bugs object (so we can get
  ## the nice bugs.r summary and plots)
  n.chains <- length(object.list)
  parnames <- colnames(object.list[[1]])
  if (!is.na(parms.tokeep[1])) {
    c.tokeep <- NULL
    for (i in 1:length(parms.tokeep)) {
      c.tokeep <- c(c.tokeep,grep(parms.tokeep[i],parnames))
    }
    c.tokeep <- unique(c.tokeep)
  } else {
    c.tokeep <- 1:length(parnames)
  }
  if (n.chains<2) stop ("n.chains must be at least 2")
  first.chain <- object.list[[1]][,c.tokeep]
  n.parameters <- ncol(first.chain)
  n.iter <- nrow(first.chain)
  sims <- array (NA, c(n.iter, n.chains, n.parameters))
  par.names <- dimnames(first.chain)[[2]]
  sims[,1,] <- first.chain
  for (k in 2:n.chains){
    sims[,k,] <- object.list[[k]][,c.tokeep]
  }
  dimnames(sims) <- list (NULL, NULL, par.names)
  if (!is.na(parms.todrop[1])) {
    c.todrop <- NULL
    for (i in 1:length(parms.todrop)) {
      c.todrop <- c(c.todrop,grep(parms.todrop[i],par.names))
    }
    c.todrop <- unique(c.todrop)
    sims <- sims[,,-c.todrop]
  }
  if (bugs){
    r <- as.bugs.array2 (sims,engine=engine)
    r$isDIC <- FALSE
    return (r)
  }
  else {
    return (sims)
  }
}




##- var y[N],p[N],theta[n.legs],beta[n.rolls,2],alpha.party[n.parties],alpha.state[n.states],leg[

mcarray2coda <- function(model.samples,thin=1,start=0) {
  require(rjags)
  nc <- dim(model.samples[[1]])["chain"]
  niter <- dim(model.samples[[1]])["iteration"]
  ans <- vector("list", nc)
  for (ch in 1:nc) {
    ans.ch <- vector("list", length(model.samples))
    vnames.ch <- NULL
    for (i in seq(along = model.samples)) {
      varname <- names(model.samples)[[i]]
      d <- dim(model.samples[[i]])
      if (length(d) < 3) {
        stop("Invalid dimensions for sampled output")
      }
      vardim <- d[1:(length(d) - 2)]
      nvar <- prod(vardim)
      niter <- d[length(d) - 1]
      nchain <- d[length(d)]
      values <- as.vector(model.samples[[i]])
      var.i <- matrix(NA, nrow = niter, ncol = nvar)
      for (j in 1:nvar) {
        var.i[, j] <- values[j + (0:(niter - 1)) * nvar +
                             (ch - 1) * niter * nvar]
      }
      vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
      ans.ch[[i]] <- var.i
    }
    ans.ch <- as.matrix(data.frame(ans.ch))
    colnames(ans.ch) <- vnames.ch
    ans[[ch]] <- mcmc(ans.ch, start = start, thin = thin)
  }
  mcmc.list(ans)
}

mcarray2bugs <- function(post,...) {
  coda2bugs(mcarray2coda(post),...)
}


array2mat <- function(x) {
  nc <- dim(x)["chain"]
  newx <- NULL
  for (i in 1:nc) {
    if (dim(x)[1]>1) {
      newx <- cbind(newx,x[,,i])
    } else {
      ## single parameter posterior
      newx <- c(newx,x[1,,i])
    }
  }
  newx
}

f2n <- function(f) {
    if (is.factor(f)) {
        as.numeric(levels(f))[f]
    } else {
        as.numeric(f)
    }    
}

melt.ideal <- function(x,sims.to.keep=NULL) {
    require(reshape)
    if (!is.null(sims.to.keep)) {
        x <- x[x[,"Iteration"]%in%sims.to.keep,-1]
    } else {
        x <- x[,-1]
    }
    x.samples <- melt(x)
    nc <- nchar(as.character(x.samples$X2))
    x.samples$j1 <- as.character(substr(x.samples$X2,nc-1,nc))
    x.samples$i <- gsub("Vote |Difficul","",as.character(substr(x.samples$X2,1,nc-2)))    
    ## if is difficulty, let it be dim 0
    x.samples$j1[x.samples$j1=="ty"] <- "d0"    
    x.samples <- subset(x.samples,select=-X2)    
    names(x.samples)[1] <- "sim"
    ##x.samples$i <- f2n(x.samples$i)
    x.samples$sim <- f2n(x.samples$sim)
    x.samples$j1 <- factor(x.samples$j1)
    x.samples <- data.frame(x.samples)
    x.samples <- recast(x.samples,sim+i~j1,id.var=c("i","j1","sim"))
    x.samples
}

melt.mcarray <- function(x) {
  require(reshape)
  dimx <- dim(x)
  x.samples <- melt(array(x,dim=dimx))
  jnames <- NULL
  if(length(dimx)>3) {
    jnames <- paste ("j", 1:(dimx[2]-1),sep="")
  }
  names(x.samples) <- c("i",jnames,"sim","chain","value")
  x.samples
}




sc <- function(x,r=c(0,1)) {
    r <- sort(r)
    x <- x-min(x,na.rm=TRUE)
    x <- x/max(x,na.rm=TRUE)
    x*diff(r,na.rm=TRUE)+r[1]
}

pties <- cbind(party=c("pcb", "pcdob", "pdt", "pfl", "phdbs", "phs", "pl", "pmdb",
               "pmn", "pp", "pps", "prb", "prona", "prp", "prs", "psb", "psc",
               "psdb", "psdc", "psl", "psol", "pst", "pstu", "pt", "ptb", "ptc",
               "ptn", "ptr", "pv", "s.part."),
               color=rainbow(30,alpha=.5))



session2year <- function(x) {
    x*4+1794
}

## read ipea
read.ipea <- function(file,comma=FALSE) {
    zz <- paste("unaccent ", file," | sed s/\\\\.//g",sep="")
    if (comma) {
        d <- read.csv(pipe(zz),dec=",",header=FALSE,skip=1)
    } else {
        d <- read.csv2(pipe(zz),dec=",",header=FALSE,skip=1)
    }
    newname <- paste("V",1:ncol(d),sep="")
    noname <- d[1,]==""|is.na(d[1,])
    names(d)[noname] <- newname[noname]
    names(d)[!noname] <- d[1,!noname]
    d <- d[-1,]
    d
}


read.ipea.series <- function(dir="/Users/eduardo/projects/BrazilianPolitics/trunk/electoral/ipea/",
                             pattern="pop*.csv",
                             var="population",gu="state",comma=FALSE) {
    require(reshape)
    zz <- pipe(paste("ls ",dir,pattern,sep=""))
    files <- readLines(zz)
    close(zz)
    d <- NULL
    for (i in files) {
        dnow <- read.ipea(i,comma=comma)
        if (is.null(d)) {
            d <- dnow
        } else {
            d <- cbind(d,dnow[,-c(1:3)])
        }
    }
    ##return(d)
    lgu <- length(gu)
    ##d <- d[,-c((lgu+1):(lgu+2))]
    d <- d[,-c((lgu+1))]
    ##if (gu=="state") {
    names(d)[1:lgu] <- gu
    d <- melt(d,id.var=gu)
    d$variable <- gsub("V","",d$variable)
    d$state <- tolower(d$state)
    d$variable <- as.numeric(d$variable)
    d <- memisc::rename(d,variable="year",value=var)
    ##}
    data.frame(d)
}


##brazilian parties recode
recode.sit <- function(x) {
    factor(car::recode(x,"c(1,5)='titular'; 2='suplente'; c(3,7)='renuncia/cassacao'; 4='nao eleito'; 6='segundo turno'; c(8,9)='registro negado'; 10='renuncia/falecimento com substituicao'"))
}

recode.parties <- function (x) {
  x <- as.character(x)
  x <- gsub(" +","",x)
  p1 <- c("pdc","pds","ppr","pp","ppb","pmr","pan","psd","ptc")
  p2 <- c("pp" ,"pp" ,"pp" ,"pp","pp" ,"prb","ptb","ptb","prn")
  for (i in 1:length(p1)) x[x==p1[i]] <- p2[i]
  x <- factor(x)
  x
}


write.dta2 <- function(dataframe,version=8,...) {
  nms <- names(dataframe)
  nms <- gsub(" +|\\.","_",nms)
  n1 <- grep("^[0-9]+",nms)
  if (length(n1)>0) nms[n1] <- paste("v",nms[n1],sep="")
  nms <- iconv(nms,from="UTF-8",to="ASCII//TRANSLIT")
  names(dataframe) <- nms
  for (i in nms[sapply(dataframe,is.factor)]) {
    levels(dataframe[,i]) <- iconv(levels(dataframe[,i]),
                                   from="UTF-8",to="ASCII//TRANSLIT")
  }
  write.dta(dataframe=dataframe,version=version,...)
}


write.csv2 <- function(dataframe,...) {
  nms <- names(dataframe)
  nms <- iconv(nms,from="UTF-8",to="macroman")
  names(dataframe) <- nms
  for (i in nms[sapply(dataframe,function(x) is.factor(x))]) {
    levels(dataframe[,i]) <- iconv(levels(dataframe[,i]),
                                   from="UTF-8",to="macroman")
  }
  for (i in nms[sapply(dataframe,function(x) is.character(x))]) {
    dataframe[,i] <- iconv(dataframe[,i],
                           from="UTF-8",to="macroman")
  }
  write.csv(dataframe,...)
}


library(foreign)
library(grDevices)
library(lattice)
##try(rcompgen::rc.settings(file = TRUE))

## ltheme <- canonical.theme(color = FALSE)     ## in-built B&W theme
## ltheme$strip.background$col <- "transparent" ## change strip bg
## lattice.options(default.theme = ltheme)      ## set as default

sumC <- function(X,...) apply(X,2,sum,...)




pdf4 <- function(...) pdf(...,version="1.4")

invlogit <- function(x) plogis(x)

sort.df <- function (x,by,...) {
  x.o <- x[,by]
  if(is.factor(x[,by]))  {
    x[,by] <- as.character(x[,by])
  }
  ##sorts data.frame x by by
  ## ... options passed to "order", (so you can use something like decreasing=TRUE)
  x <- x[order(x[,by],...),]
  x[,by] <- x.o
  x
}

sort.s <- function (x,by,...) {
  ##sorts x and assigns new value to x
  ## ... options passed to "order", (so you can use something like decreasing=TRUE)
  subx <- substitute(x)
  if (is.name(subx))
    subx <- deparse(subx)
  if (!is.character(subx) || length(subx) != 1)
    stop("'sort.s' requires a name")
  parent <- parent.frame()
  x <- get(subx, envir = parent)
  x <- x[order(x[,by],...),]
  assign(subx, x, env = .GlobalEnv)
}


normalize.factor <- function(x,exclude=NULL) {
  x <- factor(gsub(" +","_",trim(tolower(as.character(x)))))
  if (!missing(exclude)) {
    for (i in 1:length(exclude)) {
      x[x==exclude[i]] <- NA
    }
  }
  factor(x)
}


##personal finance

###process wesabe csv download



get.finance <- function(file.loc="/Users/eduardo/Desktop/downloads/accounts.csv") {
  tmp <- read.csv(file.loc)
  tmp <- tmp[!tmp$Tags=="transfer",]
  tmp$up <- with(tmp,ifelse(Amount>0,Amount,0))
  tmp$down <- with(tmp,ifelse(Amount<0,Amount,0))
  tmp$year <- substr(as.character(tmp$Transaction.Date),27,30)
  tmp$month <- tolower(substr(as.character(tmp$Transaction.Date),5,7))
  tmp$date <- with(tmp,paste(year,month,28,sep=""))
  ##tmp

  tmp0 <- tmp

  tmp <- aggregate(list(tmp$Amount,tmp$up,tmp$down),by=list(tmp$date),FUN=sum)
  names(tmp) <- c("date","balance","up","down")
  tmp$date <- as.Date(tmp$date,format="%Y%b%d")
  tmp <- tmp[order(tmp$date),]
  tmp <- tmp[tmp$date>as.Date("20070101",format="%Y%m%d"),]

  tmp1 <- tmp

  tmp <- tmp0
  tmp <- aggregate(list(tmp$Amount,tmp$up,tmp$down),by=list(tmp$date,tmp$Tags),FUN=sum)
  names(tmp) <- c("date","tag","balance","up","down")
  tmp$date <- as.Date(tmp$date,format="%Y%b%d")
  tmp <- tmp[order(tmp$tag,tmp$date),]
  tmp <- tmp[tmp$date>as.Date("20070101",format="%Y%m%d"),]

  list(tmp1,tmp)

}



coef.lmer <- function(object,digits=NULL) {
  fcoef <- .Call("mer_fixef", object, PACKAGE = "lme4")
  useScale <- object@devComp[8]
  corF <- vcov(object)@factors$correlation
  coefs <- cbind(fcoef, corF@sd)
  if (length (fcoef) > 0){
    dimnames(coefs) <- list(names(fcoef), c("coef.est", "coef.se"))
  }
  if (is.null(digits)) {
    coefs
  } else {
    round(coefs,digits)
  }
}


ns <- function(x,y) {
  gsub(" ","_",paste(x,y,sep="_"))
}

fix.parties <- function(x) {
  require(car)
  x <- toupper(x)
  x <- car:::recode(x,"'PPR'='PP'; 'PPB'='PP'")
  parties <- c("PT", "PCDOB", "PSB", "PDT", "PV", "PPS", "PL", "PMDB", "PP",
               "PTB", "PSDB", "PFL", "PRESIDENT")
  x <- factor(parties[(match(x,parties))],levels=parties)
  levels(x) <- toupper(levels(x))
  x
}


.eval <- function(evaltext,envir=sys.frame()) {
  eval(parse(text=evaltext), envir=envir)
}


##viewColors by Kevin Wright
viewColors = function (ind = 1:length(colors()), identify = FALSE, sort = TRUE)
{
  ind = setdiff(ind, grep("grey", colors()[ind]))
  if (sort) {
    indNames = colors()[ind]
    indRGB = col2rgb(indNames)
    indHSV = rgb2hsv(indRGB)
    hueOrder = order(indHSV[1, ], indHSV[2, ], indHSV[3,
      ])
    ind = ind[hueOrder]
  }
  indNames = colors()[ind]
  num = ceiling(sqrt(length(ind)))
  m = matrix(1:num^2, nrow = num)
  x = y = 1:num
  image(x, y, m, axes = F, col = c(colors()[ind], rep(colors()[1],
                             num^2 - length(ind))), xlab = "", ylab = "")
  for (i in 1:length(ind) - 1) text(i%%num + 1, i%/%num + 1,
                                    ind[i + 1], cex = 0.6)
  if (identify) {
    title("Left-click to choose colors.  Right-click to exit")
    cols = locator()
    x = round(cols$x)
    y = round(cols$y)
    ind = (y - 1) * num + x
    return(indNames[ind])
  }
}

substr.e <- function(x,start,stop) {
  if (start<1) {
    start <- nchar(as.character(x))-start
  }
  substr(x,start,stop)
}


sweatex <- function(file,bibtex=FALSE,open=FALSE,cache=FALSE) {
  file.root <- gsub("\\.Rnw$","",file)
  if (cache) {
      require("cacheSweave")
      Sweave(file,driver = cacheSweaveDriver)
  } else {
      Sweave(file)
  }
  system(paste("pdflatex ",file.root,".tex",sep=""))
  if (bibtex) {
    system(paste("bibtex ",file.root,sep=""))
    system(paste("pdflatex ",file.root,".tex",sep=""))
    system(paste("pdflatex ",file.root,".tex",sep=""))
  }
  if(open) system(paste("open -a Skim ",file.root,".pdf",sep=""))
}


states <- c("AC","AL","AP","AM","BA","CE","DF","ES","GO","MA","MT","MS","MG",      
            "PA","PB","PR","PE","PI","RJ","RN","RS","RO","RR","SC","SP","SE","TO") 



demean <- function(x,index) {
  for (i in unique(index)) {
    for (j in 1:ncol(x)) {
      x.now <- x[index==i,j]
      x[index==i,j] <- x.now-mean(x.now,na.rm=TRUE)
    }
  }
  x
}

mylag <- function(x) {
  l <- length(x)
  x <- c(NA,x)[1:l]
  x
}

mydif <- function(x) {
  l <- length(x)
  x <- x-c(NA,x)[1:l]
  x
}

decode <- function(x) {
  ##if (is.factor(x)) as.numeric(levels(x))[x]
  if (is.factor(x)) {
    as.numeric(as.character(x))
  }  else x
}

corrections <- function (x) {
  x <- as.character(x)
  x <- trim(x)
  x <- tolower(x)
  x <- gsub("'"," ",x)
  x <- gsub("-"," ",x)
  x <- gsub("sao vicente do serido","serido",x)
  x <- gsub("de goias"," ",x)
  x <- gsub("do tocantins"," ",x)
  x <- gsub("do maranhao"," ",x)
  x <- gsub("de umbuzeiro"," ",x)
  x <- gsub("do anaua"," ",x)
  x <- gsub("senador"," ",x)
  x <- gsub("stssma","santissima",x)
  x <- gsub("tejucuoca","tejussuoca",x)
  x <- gsub("acu","assu",x)
  x <- gsub("brodosqui","brodowski",x)
  x <- gsub("seb\\.","sebastiao",x)
  trim(x)
}

states.l <- c('acre','alagoas','amazonas','amapa','bahia','ceara','distrito federal','espirito santo','goias','maranhao','minas gerais','mato grosso do sul','mato grosso','para','paraiba','pernambuco','piaui','parana','rio de janeiro','rio grande do norte','rondonia','roraima','rio grande do sul','santa catarina','sergipe','sao paulo','tocantins')


trim <-  function (s)
{
  s <- sub("^\t+","", s)
  s <- sub("^ +", "", s)
  s <- sub(" +$", "", s)
  s
}

normalize <-  function (s)
{
  s <- sub("^ +", "", s)
  s <- sub(" +$", "", s)
  s <- gsub(" +"," ", s)
  s <- gsub("%","pct_",s)
  s <- gsub("\\.+","_",s)
  tolower(s)
}


state.a2l <- function(object) {
  object <- as.character(tolower(object))
  require(car)
  car:::recode(object,"'ac' ='acre';
                       'al' ='alagoas';
                       'am' ='amazonas';
                       'ap' ='amapa';
                       'ba' ='bahia';
                       'ce' ='ceara';
                       'df' ='distrito federal';
                       'es' ='espirito santo';
                       'go' ='goias';
                       'ma' ='maranhao';
                       'mg' ='minas gerais';
                       'ms' ='mato grosso do sul';
                       'mt' ='mato grosso';
                       'pa' ='para';
                       'pb' ='paraiba';
                       'pe' ='pernambuco';
                       'pi' ='piaui';
                       'pr' ='parana';
                       'rj' ='rio de janeiro';
                       'rn' ='rio grande do norte';
                       'ro' ='rondonia';
                       'rr' ='roraima';
                       'rs' ='rio grande do sul';
                       'sc' ='santa catarina';
                       'se' ='sergipe';
                       'sp' ='sao paulo';
                       'to' ='tocantins'")
}

state.a2L <- function(object) {
  object <- as.character(tolower(object))
  require(car)
  car:::recode(object,"'ac' ='Acre';
                       'al' ='Alagoas';
                       'am' ='Amazonas';
                       'ap' ='Amapá';
                       'ba' ='Bahia';
                       'ce' ='Ceará';
                       'df' ='Distrito Federal';
                       'es' ='Espírito Santo';
                       'go' ='Goiás';
                       'ma' ='Maranhão';
                       'mg' ='Minas Gerais';
                       'ms' ='Mato Grosso do Sul';
                       'mt' ='Mato Grosso';
                       'pa' ='Pará';
                       'pb' ='Paraíba';
                       'pe' ='Pernambuco';
                       'pi' ='Piauí';
                       'pr' ='Paraná';
                       'rj' ='Rio de Janeiro';
                       'rn' ='Rio Grande do Norte';
                       'ro' ='Rondônia';
                       'rr' ='Roraima';
                       'rs' ='Rio Grande do Sul';
                       'sc' ='Santa Catarina';
                       'se' ='Sergipe';
                       'sp' ='São Paulo';
                       'to' ='Tocantins'")
}

state.a2ln <- function(object) {
  object <- as.character(tolower(object))
  require(car)
  car:::recode(object,"'ac' = 1;
                       'al' = 2;
                       'am' = 3;
                       'ap' = 4;
                       'ba' = 5;
                       'ce' = 6;
                       'df' = 7;
                       'es' = 8;
                       'go' = 9;
                       'ma' = 10;
                       'mg' = 11;
                       'ms' = 12;
                       'mt' = 13;
                       'pa' = 14;
                       'pb' = 15;
                       'pe' = 16;
                       'pi' = 17;
                       'pr' = 18;
                       'rj' = 19;
                       'rn' = 20;
                       'ro' = 21;
                       'rr' = 22;
                       'rs' = 23;
                       'sc' = 24;
                       'se' = 25;
                       'sp' = 26;
                       'to' = 27")
}


## function to recode long state names to short
state.l2a <- function(object) {
  ##require(car)
  object <- tolower(as.character(object))
  car:::recode(object,"
                       'acre'                  = 'ac';
                       'alagoas'               = 'al';
                       'amazonas'              = 'am';
                       c('amapa','amapá')      = 'ap';
                       'bahia'                 = 'ba';
                       c('ceara','ceará')      = 'ce';
                       'distrito federal'      = 'df';
                       'espirito santo'        = 'es';
                       'espitito santo'        = 'es';
                       c('goias','goiás')      = 'go';
                       c('maranhao','maranhão')= 'ma';
                       'minas gerais'          = 'mg';
                       'mato grosso do sul'    = 'ms';
                       'mato g sul'            = 'ms';
                       'm. g. do sul'          = 'ms';
                       'mato grosso'           = 'mt';
                       'para'                  = 'pa';
                       'pará'                  = 'pa';
                       'paraiba'               = 'pb';
                       'paraíba'               = 'pb';
                       'pernambuco'            = 'pe';
                       'piaui'                 = 'pi';
                       'piauí'                 = 'pi';
                       'parana'                = 'pr';
                       'paraná'                = 'pr';
                       'rio de janeiro'        = 'rj';
                       'r. g. do norte'        = 'rn';
                       'rio grande do norte'   = 'rn';
                       'rondonia'              = 'ro';
                       'rondônia'              = 'ro';
                       'roraima'               = 'rr';
                       'rio grande do sul'     = 'rs';
                       'r. g. do sul'          = 'rs';
                       'santa catarina'        = 'sc';
                       'sergipe'               = 'se';
                       'sao paulo'             = 'sp';
                       'são paulo'             = 'sp';
                       'tocantins'              = 'to'")
}





## amatch <- function(names1,names2,max) {
##   ##convert to char
##   names1 <- as.character(names1)
##   names2 <- as.character(names2)
##   ##hold unmatched names1
##   ##hold matches
##   matched <- NULL
##   for (j in max) {
##     names1u <- NULL
##     for (i in names1) {
##       match.now <- agrep(i,names2,value=FALSE,max=j)[1]
##       if (length(match.now)!=0&(!is.na(match.now))) {
##         matched <- rbind(matched,cbind(i,names2[match.now],j))
##         names2 <- names2[-match.now]
##       } else {
##         names1u <- c(names1u,i)
##       }
##     }
##     names1 <- unique(names1u)
##     if (length(names1u)==0) break
##   }
##   list(
##        matched=data.frame(names1=matched[,1],names2=matched[,2],max=matched[,3]),
##        unmatched1=unique(names1),unmatched2=unique(names2))
## }

## amatch2 <- function(names1,names2,max,max2=1) {
##   ## if column 1 matches, check column 2
##   ##convert to char
##   names1a <- as.character(names1[,1])
##   names1b <- as.character(names1[,2])
##   names2a <- as.character(names2[,1])
##   names2b <- as.character(names2[,2])
##   index1 <- 1:nrow(names1)
##   index2 <- 1:nrow(names2)
##   ##hold matches
##   matched <- NULL
##   for (j in max) {
##     cat("using max=",j,"\n")
##     index1u <- NULL     ## holds unmatched index1
##     index2u <- NULL     ## holds unmatched index2
##     for (i in index1) {
##       ##names2a available to match
##       names2tomatch <- names2a[index2]
##       ##match  column 1
##       match.now <- agrep(names1a[i],names2tomatch,value=FALSE,max=j)
##       ##match.now is the index of the names2match vector (which has length index2)
##       ##the index2 matches are therefore
##       index2matched <- index2[match.now]
##       ##match column 2 within matches
##       match2 <- agrep(names1b[i],names2b[index2matched],max=max2)[1]
##       ## and hold 1st value only
##       ## is there a match in match2?
##       match2T <- length(match2)>0&(!is.na(match2))
##       ##if match2T get match2 from match.now
##       match.now <- ifelse(match2T,
##                           match.now[match2],
##                           match.now[1])##else get 1st value of match1
##       ##update index2matched
##       index2matched <- index2[match.now]
##       match1T <- length(match.now)>0&(!is.na(match.now))
##       if ((match1T)        ## if existe match.now and (match2 exists
##           & (match2T|is.na(names1b[i]))) { ## or names1b is NA
##         matched <- rbind(matched,c(i,index2matched,1-j))
##         index2 <- index2[!(index2==index2matched)]
##       }
##       else { ## there was no match
##         index1u <- c(index1u,i)
##       }
##     }
##     ##update index1
##     index1 <- index1u
##     if (length(index1)==0) break
##   }
##   list(matched1=matched[,1],matched2=matched[,2],confidence=matched[,3],unmatched1=index1,unmatched2=index2)
## }



## agrepwrap <- function(x,y,max) {
##   matchedy <- rep(NA,length(y))
##   matchedx <- rep(NA,length(x))
##   indexx <- 1:length(x)
##   indexy <- 1:length(y)

##   for (i in 1:length(unmatchedx)) {
##     x2match <- indexx[is.na(matchedx)][1]
##     matched.now <- agrep(x[x2match],y[is.na(matchedy)],max=max)
##     if (length(matched.now)==1) {
##       matchedx[x2match] <- 1
##     } else {
##       matchedx[x2match] <- 0
##     }
##   }
## }


## amatch2b <- function(x,y,by1,by2) {
## ### match datasets exactly on by1 and approximate match on by2
##   ##
##   x <- as.data.frame(x)
##   y <- as.data.frame(y)
##   x[,by2] <- as.character(x[,by2])
##   y[,by2] <- as.character(y[,by2])

##   ##insert indices in x and y
##   x <- x$index <- 1:nrow(x)
##   y <- y$index <- 1:nrow(y)

##   ##available on by1
##   by1u <- unique(as.character(x[,by1]))
##   for (i in by1u) {
##     ##exact match on by1
##     xnow <- x[x[,by1]==i,]
##     ynow <- y[y[,by1]==i,]
##     ##approximate match on by2
##     for (j in 1:nrow(xnow)) {
##       matchnow <- c(xnow[j,],grep(xnow[j,by2],ynow[,by2],max=max))
##     }

##   }
## }

user.prompt <- function (message="\nPress <return> to continue: ")
  silent <- readline(message)
##from Zelig


index <- function(value,object) {
  index.i <- 1:length(object)
  index.i[object==value]
}

states <- (c("ac",
             "al",
             "am",
             "ap",
             "ba",
             "ce",
             "df",
             "es",
             "go",
             "ma",
             "mg",
             "ms",
             "mt",
             "pa",
             "pb",
             "pe",
             "pi",
             "pr",
             "rj",
             "rn",
             "ro",
             "rr",
             "rs",
             "sc",
             "se",
             "sp",
             "to"))


parse2 <- function(x) eval(parse(text=x))

capwords <- function(s, strict = TRUE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}


coefplot.polr.e <- function (object, varnames = NULL, ...) {
  coefs <- summary(object)$coef[, 1]
  sds <- summary(object)$coef[, 2]
  ifelse(is.null(varnames), varnames <- names(coefs), varnames <- varnames)
  coefplot.default(coefs, sds, varnames = varnames, CI = 2,
                   vertical = TRUE, cex.var = 0.8, cex.pts = 0.9, col.pts = 1,
                   var.las = 2,...)
}





pcor <- function (x) {
  require(polycor)
  if (is.data.frame(x))     x <- as.matrix(x)
  ncy <- ncx <- ncol(x)
  r <- matrix(0, nrow = ncx, ncol = ncy)
  for (i in seq(2, length = ncx - 1)) {
    for (j in seq_len(i - 1)) {
      x2 <- x[, i]
      y2 <- x[, j]
      ok <- complete.cases(x2, y2)
      r[i, j] <- polychor(x2, y2)
    }
  }
  r <- r + t(r)
  diag(r) <- 1
  rownames(r) <- colnames(x)
  colnames(r) <- colnames(x)
  r
}




legend2 <- function (x, y = NULL, legend, fill = NULL, col = par("col"),
                     lty, lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"),
                     box.lwd = par("lwd"), box.lty = par("lty"), pt.bg = NA, cex = 1,
                     pt.cex = cex, pt.lwd = lwd, xjust = 0, yjust = 1, x.intersp = 1,
                     y.intersp = 1, adj = c(0, 0.5), text.width = NULL, text.col = par("col"),
                     merge = do.lines && has.pch, trace = FALSE, plot = TRUE,
                     ncol = 1, horiz = FALSE, title = NULL, inset = 0, fill.lwd=1,bcol.fill="transparent")
{
  if (missing(legend) && !missing(y) && (is.character(y) ||
                                         is.expression(y))) {
    legend <- y
    y <- NULL
  }
  mfill <- !missing(fill) || !missing(density)
  title <- as.graphicsAnnot(title)
  if (length(title) > 1)
    stop("invalid title")
  legend <- as.graphicsAnnot(legend)
  n.leg <- if (is.call(legend))
    1
  else length(legend)
  if (n.leg == 0)
    stop("'legend' is of length 0")
  auto <- if (is.character(x))
    match.arg(x, c("bottomright", "bottom", "bottomleft",
                   "left", "topleft", "top", "topright", "right", "center"))
  else NA
  if (is.na(auto)) {
    xy <- xy.coords(x, y)
    x <- xy$x
    y <- xy$y
    nx <- length(x)
    if (nx < 1 || nx > 2)
      stop("invalid coordinate lengths")
  }
  else nx <- 0
  xlog <- par("xlog")
  ylog <- par("ylog")
  rect2 <- function(left, top, dx, dy, density = NULL, angle,
                    ...) {
    r <- left + dx
    if (xlog) {
      left <- 10^left
      r <- 10^r
    }
    b <- top - dy
    if (ylog) {
      top <- 10^top
      b <- 10^b
    }
    rect(left, top, r, b, angle = angle, density = density,
         ...)
  }
  segments2 <- function(x1, y1, dx, dy, ...) {
    x2 <- x1 + dx
    if (xlog) {
      x1 <- 10^x1
      x2 <- 10^x2
    }
    y2 <- y1 + dy
    if (ylog) {
      y1 <- 10^y1
      y2 <- 10^y2
    }
    segments(x1, y1, x2, y2, ...)
  }
  points2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    points(x, y, ...)
  }
  text2 <- function(x, y, ...) {
    if (xlog)
      x <- 10^x
    if (ylog)
      y <- 10^y
    text(x, y, ...)
  }
  if (trace)
    catn <- function(...) do.call("cat", c(lapply(list(...),
                                                  formatC), list("\n")))
  cin <- par("cin")
  Cex <- cex * par("cex")
  if (is.null(text.width))
    text.width <- max(abs(strwidth(legend, units = "user",
                                   cex = cex)))
  else if (!is.numeric(text.width) || text.width < 0)
    stop("'text.width' must be numeric, >= 0")
  xc <- Cex * xinch(cin[1], warn.log = FALSE)
  yc <- Cex * yinch(cin[2], warn.log = FALSE)
  if (xc < 0)
    text.width <- -text.width
  xchar <- xc
  xextra <- 0
  yextra <- yc * (y.intersp - 1)
  ymax <- yc * max(1, strheight(legend, units = "user", cex = cex)/yc)
  ychar <- yextra + ymax
  if (trace)
    catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra,
                                                   ychar))
  if (mfill) {
    xbox <- xc * 0.8
    ybox <- yc * 0.5
    dx.fill <- xbox
  }
  do.lines <- (!missing(lty) && (is.character(lty) || any(lty >
                                                          0))) || !missing(lwd)
  n.legpercol <- if (horiz) {
    if (ncol != 1)
      warning("horizontal specification overrides: Number of columns := ",
              n.leg)
    ncol <- n.leg
    1
  }
  else ceiling(n.leg/ncol)
  if (has.pch <- !missing(pch) && length(pch) > 0) {
    if (is.character(pch) && !is.na(pch[1]) && nchar(pch[1],
                                                     type = "c") > 1) {
      if (length(pch) > 1)
        warning("not using pch[2..] since pch[1] has multiple chars")
      np <- nchar(pch[1], type = "c")
      pch <- substr(rep.int(pch[1], np), 1:np, 1:np)
    }
    if (!merge)
      dx.pch <- x.intersp/2 * xchar
  }
  x.off <- if (merge)
    -0.7
  else 0
  if (is.na(auto)) {
    if (xlog)
      x <- log10(x)
    if (ylog)
      y <- log10(y)
  }
  if (nx == 2) {
    x <- sort(x)
    y <- sort(y)
    left <- x[1]
    top <- y[2]
    w <- diff(x)
    h <- diff(y)
    w0 <- w/ncol
    x <- mean(x)
    y <- mean(y)
    if (missing(xjust))
      xjust <- 0.5
    if (missing(yjust))
      yjust <- 0.5
  }
  else {
    h <- (n.legpercol + (!is.null(title))) * ychar + yc
    w0 <- text.width + (x.intersp + 1) * xchar
    if (mfill)
      w0 <- w0 + dx.fill
    if (has.pch && !merge)
      w0 <- w0 + dx.pch
    if (do.lines)
      w0 <- w0 + (2 + x.off) * xchar
    w <- ncol * w0 + 0.5 * xchar
    if (!is.null(title) && (tw <- strwidth(title, units = "user",
                                           cex = cex) + 0.5 * xchar) > w) {
      xextra <- (tw - w)/2
      w <- tw
    }
    if (is.na(auto)) {
      left <- x - xjust * w
      top <- y + (1 - yjust) * h
    }
    else {
      usr <- par("usr")
      inset <- rep(inset, length.out = 2)
      insetx <- inset[1] * (usr[2] - usr[1])
      left <- switch(auto, bottomright = , topright = ,
                     right = usr[2] - w - insetx, bottomleft = , left = ,
                     topleft = usr[1] + insetx, bottom = , top = ,
                     center = (usr[1] + usr[2] - w)/2)
      insety <- inset[2] * (usr[4] - usr[3])
      top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3] +
                    h + insety, topleft = , top = , topright = usr[4] -
                    insety, left = , right = , center = (usr[3] +
                                                         usr[4] + h)/2)
    }
  }
  if (plot && bty != "n") {
    if (trace)
      catn("  rect2(", left, ",", top, ", w=", w, ", h=",
           h, ", ...)", sep = "")
    rect2(left, top, dx = w, dy = h, col = bg, density = NULL,
          lwd = box.lwd, lty = box.lty)
  }
  xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1),
                                              rep.int(n.legpercol, ncol)))[1:n.leg]
  yt <- top - 0.5 * yextra - ymax - (rep.int(1:n.legpercol,
                                             ncol)[1:n.leg] - 1 + (!is.null(title))) * ychar
  if (mfill) {
    if (plot) {
      fill <- rep(fill, length.out = n.leg)
      rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox,
            col = fill, density = density, angle = angle,lwd=fill.lwd,
            border = bcol.fill)
    }
    xt <- xt + dx.fill
  }
  if (plot && (has.pch || do.lines))
    col <- rep(col, length.out = n.leg)
  if (missing(lwd))
    lwd <- par("lwd")
  if (do.lines) {
    seg.len <- 2
    if (missing(lty))
      lty <- 1
    lty <- rep(lty, length.out = n.leg)
    lwd <- rep(lwd, length.out = n.leg)
    ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
    if (trace)
      catn("  segments2(", xt[ok.l] + x.off * xchar, ",",
           yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
    if (plot)
      segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len *
                xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l],
                col = col[ok.l])
    xt <- xt + (seg.len + x.off) * xchar
  }
  if (has.pch) {
    pch <- rep(pch, length.out = n.leg)
    pt.bg <- rep(pt.bg, length.out = n.leg)
    pt.cex <- rep(pt.cex, length.out = n.leg)
    pt.lwd <- rep(pt.lwd, length.out = n.leg)
    ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
    x1 <- (if (merge)
           xt - (seg.len/2) * xchar
    else xt)[ok]
    y1 <- yt[ok]
    if (trace)
      catn("  points2(", x1, ",", y1, ", pch=", pch[ok],
           ", ...)")
    if (plot)
      points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok],
              bg = pt.bg[ok], lwd = pt.lwd[ok])
    if (!merge)
      xt <- xt + dx.pch
  }
  xt <- xt + x.intersp * xchar
  if (plot) {
    if (!is.null(title))
      text2(left + w/2, top - ymax, labels = title, adj = c(0.5,
                                                      0), cex = cex, col = text.col)
    text2(xt, yt, labels = legend, adj = adj, cex = cex,
          col = text.col)
  }
  invisible(list(rect = list(w = w, h = h, left = left, top = top),
                 text = list(x = xt, y = yt)))
}


edvreg2 <- function(y,x,G) {
  require(MASS)
  tr <- function(x) sum(diag(x))
  X <- cbind(1,x)
  N <- nrow(X)
  K <- ncol(X)
  sumv <-  sum(residuals(lm(y~x))^2)
  ssq <- sumv/(N-K)
  omegasq <- diag(G)
  W <- ssq*(diag(N)-G)
  sigmahatsq <- (sumv-sum(omegasq)-
                 tr(ginv(t(X)%*%X)%*%t(X)%*%W%*%X
                    ))/(N-K)
  if (sigmahatsq<0) sigmahatsq <- 0
  S <- G+diag(N)*sigmahatsq
  ## now the final regression
  Sinv <- ginv(S)
  xSx <- ginv(t(X)%*%Sinv%*%X)
  beta <- xSx%*%t(X)%*%Sinv%*%y
  yhat <- X%*%beta
  res <- y-yhat
  s2 <- sum(t(res)%*%Sinv%*%(res))/(N-K)
  vc <- s2*xSx
  cat("Estimated Sigma =",sqrt(sigmahatsq),"\n")
  ##model.Fgls
  list(beta=beta,vbeta=vc)
}




## jackknife se estimator for Design::lrm
lrm.all <- function(formula,data,cluster) {
  library(Design)
  u.g <- unique(cluster)
  N <- length(u.g)
  pool <- lrm(formula,data=data,x=TRUE,y=TRUE)
  pool.b <- coef(pool)
  J <- length(pool.b)
  jack.gamma <- sapply(u.g,function(x) lrm(formula,data=data[cluster!=x,])$coefficients)
  jack.gamma <- N*pool.b-(N-1)*jack.gamma
  get.se <- function(x) sqrt(diag(x))
  r.vcov <- robcov(pool,cluster)
  jack.se <- sqrt(apply(jack.gamma,1,var)/N)
  res <- cbind(coef(pool),get.se(pool$var),get.se(r.vcov$var),jack.se)
  colnames(res) <- c("coef","se","robust.se","jackknife.se")
  rownames(res) <- names(pool.b)
  print(round(res,4))
  list(pool,r.vcov,jack.se,res)
}



latex.table <- function(coef.mat,se.mat,digits=3,table.command=TRUE) {
  nc <- ncol(coef.mat)
  coef.mat <- round(coef.mat,3)
  ##se.mat <- round(se.mat,3)
  se.mat <- format(se.mat,digits=1)
  text.now <- NULL
  if (table.command) {
    text.now <- c(text.now,"\\begin{table}\n")
    text.now <- c(text.now,"\\centering\n")
  }
  text.now <- c(text.now,"\\begin{tabular}[R]{",rep("c",nc+1),"}\n")
  text.now <- c(text.now,"\\hline\n")
  for (j in 1:ncol(coef.mat)) {
    text.now <- c(text.now," & ", colnames(coef.mat)[j])
  }
  text.now <- c(text.now,"\\\\\n\\hline\n")
  for (i in 1:nrow(coef.mat)) {
    ##print coef estimates
    text.now <- c(text.now,rownames(coef.mat)[i])
    for (j in 1:ncol(coef.mat)) {
      if (is.na(coef.mat[i,j])) {
        text.now <- c(text.now," & ")
      } else {
        text.now <- c(text.now," & ", coef.mat[i, j])
      }
    }
    text.now <- c(text.now,"\\\\\n")
    ##  print SEs
    for (j in 1:ncol(coef.mat)) {
      if (trim(se.mat[i,j])=="NA") {
        text.now <- c(text.now," & ")
      } else {
        text.now <- c(text.now," & ",
                      ##format(
                      sprintf("(%s)",
                      se.mat[i,j])
                      ##,digits=1
                      )
      }
    }
    text.now <- c(text.now,"\\\\[1mm]\n")
  }
  text.now <- c(text.now,"\\\\\n")
  text.now <- c(text.now,"\\hline")
  text.now <- c(text.now,"\n")
  text.now <- c(text.now,"\\end{tabular}\n")
  if (table.command) text.now <- c(text.now,"\\end{table}\n")
  paste(text.now,collapse="")
}



plot.se <- function(x,
                    y,
                    se,
                    lo,
                    hi,
                    identify=FALSE,
                    label.points=y,lowess=FALSE,fgls=FALSE,pvalue=TRUE,padj=0.6,
                    lm=FALSE,ylim,digits=2,f=function(x) x,
                    ...) {
  if (missing(ylim)) ylim <- range(c(lo,hi))
  model.lm <- lm(y~x)
  model.h <- edvreg(y,x,se)
  s.h <- summary(model.h)$coefficients
  p <- s.h[2,4]
  b <- s.h[2,1]
  o.x <- order(x)
  ##xlab.text=paste("hat(gamma)=",round(b,2),"\n p value=",round(p,2)),
  plot(x=f(x),y=y,type="p",xlab="",
       ##xlab=xlab.text,
       ylim=ylim,...)
  if (pvalue) {
    xlab.text <- substitute(list(hat(gamma)==b,pvalue==p),list(b=round(b,digits),p=round(p,digits)))
    mtext(xlab.text,side=1,cex=.6,line=2)
  } else {
    xlab.text <- substitute(list(hat(gamma)==b),list(b=round(b,digits)))
    xlab.text2 <- substitute(list((se)),list(se=round(s.h[2,2],digits)))
    mtext(xlab.text,side=1,cex=.7,line=2,adj=0.5)
    mtext(xlab.text2,side=1,cex=.7,line=3,adj=padj)
  }
  ##points(x,y)
  for (i in 1:length(x)) {
    lines(y=c(lo[o.x[i]],hi[o.x[i]]),x=f(c(x[o.x[i]],x[o.x[i]])),lwd=0.35,col="gray30")
  }
  ##mtext(side=1,text=paste("p-value",round(p,2)),line=-1,cex=0.75,adj=0)
  if (lowess) {
    lines(lowess(x,y), col = 2)
  }
  PtoPlot <- o.x[c(1,length(o.x))]
  if (lm) {
    lines(x=f(x[PtoPlot]),
          y=model.lm$fitted.values[PtoPlot],lty=8)
  }
  if (fgls) {
    lines(x=f(x[PtoPlot]),
          y=model.h$fitted.values[PtoPlot],lty=1)
  }
  if (identify) {
    identities <- identify(f(x),y,label=label.points)
    identities
  }
}

edvreg <- function(y,x,se) {
  require(MASS)
  X <- cbind(1,x)
  N <- nrow(X)
  sumv <-      sum(residuals(lm(y~x))^2)
  omegasq <- se^2
  G <- diag(N)*omegasq
  sigmahatsq <-
    (
     sumv-sum(omegasq)
     +sum(diag(
         ginv(t(X)%*%X)%*%t(X)%*%G%*%X
         ))
     )/
       (N-ncol(X))
  if (sigmahatsq<0) sigmahatsq <- 0
  model.Fgls <- lm(y~x,weights=sqrt(1/(omegasq+sigmahatsq)))
  ##cat("Estimated Sigma =",sqrt(sigmahatsq),"\n")
  model.Fgls
}


## circle
makecirclepoly <- function(center,radius,numpts=1000) {
  require(gpclib)
  theta <- seq(0,2*pi,length=numpts)
  verts <- cbind(center[1]+radius*cos(theta),center[2]+radius*sin(theta))
  return( as(verts,"gpc.poly") )
}


## circle
circle <- function(center,radius,numpts=1000) {
  theta <- seq(0,2*pi,length=numpts)
  verts <- cbind(center[1]+radius*cos(theta),center[2]+radius*sin(theta))
  verts
}










## ####### wnominate plots
## plot.cutlines2 <- function (x, main.title = "Cutting Lines", d1.title = "First Dimension",
##     d2.title = "Second Dimension", lines = 50, dims = c(1, 2),
##     lwd = 2, sset=NULL, ...)
## {
##     if (!class(x) == "nomObject")
##         stop("Input is not of class 'nomObject'.")
##     if (x$dimensions == 1)
##         stop("All angles in 1D NOMINATE are 90 degrees.")
##     if (length(dims) != 2)
##         stop("'dims' must be an integer vector of length 2.")
##     if (lines < 1)
##         stop("'Lines' must be less than 1.")
##     if (!is.null(sset)) {
##         x$rollcalls <- x$rollcalls[sset,]
##         if (length(sset)==1) {
##             x$rollcalls <- t(x$rollcalls)
##         }
##         ##print(sset)
##     }
##     constrained <- ((abs(x$rollcalls[, "spread1D"]) > 0 | abs(x$rollcalls[,
##         "spread2D"]) > 0) & (x$rollcalls[, "midpoint1D"]^2 +
##         x$rollcalls[, "midpoint2D"]^2) < 0.95)
##     cutlineData <- cbind(x$rollcalls[constrained, paste("midpoint",
##         dims[1], "D", sep = "")], x$rollcalls[constrained, paste("spread",
##         dims[1], "D", sep = "")], x$rollcalls[constrained, paste("midpoint",
##         dims[2], "D", sep = "")], x$rollcalls[constrained, paste("spread",
##         dims[2], "D", sep = "")])
##     cutlineData <- na.omit(cutlineData)
##     suppressWarnings(symbols(x = 0, y = 0, circles = 1, inches = FALSE,
##         asp = 1, main = main.title, xlab = d1.title, ylab = d2.title,
##         xlim = c(-1, 1), ylim = c(-1, 1), cex.main = 1.2, cex.lab = 1.2,
##         font.main = 2, lwd = 2, fg = "grey", frame.plot = FALSE,
##         ...))
##     if (lines < dim(cutlineData)[1])
##         cutlineData <- cutlineData[sample(1:dim(cutlineData)[1],
##             lines), ]
##     suppressWarnings(apply(cutlineData, 1, wnominate:::add.cutline, weight = x$weights[dims[2]]/x$weights[dims[1]],
##         lwd = lwd))
## }

dropRollCall2 <- function (object, dropList = NULL)
{
    if (class(object) != "rollcall") {
        stop("dropRollCall only works for objects of class rollcall.")
    }
    tmpRollCall <- object
    if (!is.list(dropList)) {
        cat("dropList must be a non-null list or alist.\nNo subsetting will occur.\n")
        return(object)
    }
    cat("Dropping elements of rollcall matrix\nusing the following dropList:\n")
    print(dropList)
    flag <- TRUE
    flag2 <- TRUE
    counter <- 1
    while (flag) {
        cat(paste("pass number", counter, "over roll call object\n"))
        v <- tmpRollCall$votes
        dimOld <- dim(v)
        if (!is.null(dropList$codes) & length(dropList$codes) >
            0) {
            cat("Processing dropList voting codes...\n")
            dc <- dropList$codes
            dCodes <- NULL
            if (all(is.character(dc))) {
                dropCodes <- match(dc, names(tmpRollCall$codes))
                dropCodes <- dropCodes[!is.na(dropCodes)]
                if (length(dropCodes) > 0) {
                  for (j in dropCodes) dCodes <- c(dCodes, tmpRollCall$codes[j])
                  keepCodes <- !(names(tmpRollCall$codes) %in%
                    dc)
                  keepCodes <- tmpRollCall$codes[keepCodes]
                  tmpRollCall$codes <- keepCodes
                }
            }
            if (is.numeric(dc)) {
                dCodes <- dc[dc %in% unique(as.vector(v))]
            }
            bad <- v %in% dCodes
            cat(paste("Will set", sum(bad), "voting decisions to NA.\n"))
            tmpRollCall$votes[bad] <- NA
            rm(bad)
        }
        dropLegis <- rep(FALSE, dim(v)[1])
        dropVotes <- rep(FALSE, dim(v)[2])
        if (!is.null(dropList$legisMin)) {
            legisMin <- dropList$legisMin
            if (length(legisMin) != 1 | is.na(legisMin) | !is.numeric(legisMin) |
                legisMin >= tmpRollCall$m)
                stop("bad value for legisMin in drop list.")
            vtmp <- convertCodes(tmpRollCall)
            goodCount <- apply(vtmp, 1, function(x) sum(!is.na(x)))
            dropLegis <- dropLegis | goodCount < legisMin
        }
        if (!is.null(dropList$dropLegis)) {
            r <- pscl:::dropRollCallViaData(dropList$dropLegis, object = tmpRollCall,
                d = expression(legis.data))
            if (!is.null(r))
                dropLegis <- dropLegis | r
        }
        if (!is.null(dropList$lop)) {
            cat("Processing lop-sided restrictions...\n")
            lop <- dropList$lop
            if (length(lop) != 1 | is.na(lop) | !is.numeric(lop) |
                lop < 0 | lop >= tmpRollCall$n) {
                print(paste("lop=",lop,", n=",tmpRollCall$n,sep=""))
                stop("Invalid value for lop")
            }
            ## FIX THIS (ELL: Shouldn't we always compute the margins? this is for
            ## estimation, not description
            ##if (is.null(tmpRollCall$voteMargins)) {
            cat("Computing vote margins...\n")
            tmpRollCall <- computeMargins(tmpRollCall, dropList = NULL)
            ##}
            r <- tmpRollCall$voteMargins[, "Min"] <= lop
            cat(paste("Will drop", sum(r), "roll calls.\n"))
            dropVotes <- dropVotes | r
            cat("Finished processing lop-sided restrictions.\n")
        }
        if (!is.null(dropList$dropVotes)) {
            r <- dropRollCallViaData(dropList$dropVotes, object = tmpRollCall,
                d = expression(vote.data))
            if (!is.null(r))
                dropVotes <- dropVotes | r
        }
        cat(paste("pass number ", counter, ": dropRollCall will drop ",
            sum(dropLegis), " legislators & ", sum(dropVotes),
            " rollcalls.\n", sep = ""))
        if (sum(dropLegis) > 0) {
            ##cat("Dropped Legislators:\n")
            ##print(dimnames(tmpRollCall$votes)[[1]][dropLegis])
        }
        if (sum(dropVotes) > 0) {
            ##cat("Dropped Votes:\n")
            ##print(dimnames(tmpRollCall$votes)[[2]][dropVotes])
        }
        if (sum(dropLegis > 0) | sum(dropVotes) > 0)
            tmpRollCall$votes <- tmpRollCall$votes[!dropLegis,
                !dropVotes]
        if (!is.null(tmpRollCall$legis.data) & sum(dropLegis) >
            0) {
            tmpRollCall$legis.data <- tmpRollCall$legis.data[!dropLegis,
                ]
            if (class(object$legis.data) == "data.frame") {
                class(tmpRollCall$legis.data) <- "data.frame"
                names(tmpRollCall$legis.data) <- names(object$legis.data)
            }
        }
        if (!is.null(tmpRollCall$vote.data) & sum(dropVotes) >
            0) {
            tmpRollCall$vote.data <- tmpRollCall$vote.data[!dropVotes,
                ]
            if (class(object$vote.data) == "data.frame") {
                class(tmpRollCall$vote.data) <- "data.frame"
                names(tmpRollCall$vote.data) <- names(object$vote.data)
            }
        }
        if (!is.null(tmpRollCall$voteMargins) & sum(dropVotes) >
            0) {
            tmpRollCall$voteMargins <- tmpRollCall$voteMargins[!dropVotes, ]
        }
        dimNew <- dim(tmpRollCall$votes)
        tmpRollCall$n <- dimNew[1]
        tmpRollCall$m <- dimNew[2]
        if (all(dimNew == dimOld)) {
            cat("dropRollCall has finished processing rollcall object\n")
            flag <- FALSE
        }
        counter <- counter + 1
    }
    tmpRollCall
}


plot.coords2 <- function (x, main.title = "W-NOMINATE Coordinates", d1.title = "First Dimension",
    d2.title = "Second Dimension", dims = c(1, 2), plotBy = "party",
    color = TRUE, shape = TRUE, cutline = NULL, Legend = TRUE,
    legend.x = 0.8, legend.y = 1, rsize=1,
                              colorlist = c("darkblue", "firebrick", "darkcyan", "darkgreen",
        "darkmagenta", "darkolivegreen", "darkorange", "darkorchid",
        "darkred", "darksalmon", "darkseagreen", "darkslateblue",
        "darkslategray", "darkturquoise", "darkviolet", "deeppink",
        "deepskyblue", "dodgerblue"),...) {
    if (!class(x) == "nomObject")
        stop("Input is not of class 'nomObject'.")
    if (!any(colnames(x$legislators) == plotBy)) {
        warning("Variable '", plotBy, "' does not exist in your W-NOMINATE object.")
        types <- rep("Leg", dim(x$legislators)[1])
    }
    else {
        types <- x$legislators[, plotBy]
    }
    if (length(dims) != 2 & x$dimensions != 1)
        stop("'dims' must be an integer vector of length 2.")
    nparties <- length(unique(types))

    shapes <- rep(c(16, 15, 17, 18, 19, 3, 4, 8), 3)
    if (color == FALSE)
        colorlist <- sample(colors()[160:220], 50)
    if (shape == FALSE)
        shapes <- rep(16, 50)
    if (x$dimensions == 1) {
        coord1D <- x$legislators[, "coord1D"]
        ranking <- rank(x$legislators[, "coord1D"])
        plot(seq(-1, 1, length = length(coord1D)), 1:length(coord1D),
            type = "n", cex.main = 1.2, cex.lab = 1.2, font.main = 2,
            xlab = "First Dimension Nominate", ylab = "Rank",
            main = "1D W-NOMINATE Plot")
        if (Legend)
            legend(0.67, 0.7 * length(coord1D), unique(types),
                pch = shapes[1:nparties], col = colorlist[1:nparties],
                cex = 0.7)
        for (i in 1:nparties) suppressWarnings(points(coord1D[types ==
            unique(types)[i]], ranking[types == unique(types)[i]],
            pch = shapes[i], col = colorlist[i], cex = rsize, lwd = 2))
    }
    else {
        coord1D <- x$legislators[, paste("coord", dims[1], "D",
            sep = "")]
        coord2D <- x$legislators[, paste("coord", dims[2], "D",
            sep = "")]
        suppressWarnings(symbols(x = 0, y = 0, circles = 1, inches = FALSE,
            asp = 1, main = main.title, xlab = d1.title, ylab = d2.title,
            xlim = c(-1, 1), ylim = c(-1, 1), cex.main = 1.2,
            cex.lab = 1.2, font.main = 2, lwd = 2, fg = "grey",
            frame.plot = FALSE, ...))
        if (!is.null(cutline)) {
            for (i in 1:length(cutline)) {
                if (all(is.na(x$rollcalls[cutline[i], ])))
                  stop("Roll call for cutline did not meet minimum lopsidedness requirements.")
                add.cutline(c(x$rollcalls[cutline[i], paste("midpoint",
                  dims[1], "D", sep = "")], x$rollcalls[cutline[i],
                  paste("spread", dims[1], "D", sep = "")], x$rollcalls[cutline[i],
                  paste("midpoint", dims[2], "D", sep = "")],
                  x$rollcalls[cutline[i], paste("spread", dims[2],
                    "D", sep = "")]), weight = x$weights[dims[2]]/x$weights[dims[1]],
                  lwd = 2)
            }
        }
        if (Legend)
            legend(legend.x, legend.y, unique(types), pch = shapes[1:nparties],
                bty = "n", col = colorlist[1:nparties], cex = 0.7)
        for (i in 1:nparties) suppressWarnings(points(coord1D[types ==
            unique(types)[i]], coord2D[types == unique(types)[i]],
            pch = shapes[i], col = colorlist, cex = rsize, lwd = 2))
    }
}


plot.cutlines2 <- function (x, main.title = "Cutting Lines", d1.title = "First Dimension",
                            d2.title = "Second Dimension", lines = 850, dims = c(1, 2), col=1,
                            col.vec,lwd = 2, sset=NULL, ...)
{
    if (!class(x) == "nomObject")
        stop("Input is not of class 'nomObject'.")
    if (x$dimensions == 1)
      stop("All angles in 1D NOMINATE are 90 degrees.")
    if (length(dims) != 2)
        stop("'dims' must be an integer vector of length 2.")
    if (lines < 1)
        stop("'Lines' must be less than 1.")
    x$rollcalls <- cbind(x$rollcalls,col)
    colnames(x$rollcalls)[ncol(x$rollcalls)] <- "col"
    if (!is.null(sset)) {
        x$rollcalls <- x$rollcalls[sset,]
        if (length(sset)==1) {
            x$rollcalls <- t(x$rollcalls)
        }
        ##print(sset)
    }
    constrained <- ((abs(x$rollcalls[, "spread1D"]) > 0 | abs(x$rollcalls[,
        "spread2D"]) > 0) & (x$rollcalls[, "midpoint1D"]^2 +
        x$rollcalls[, "midpoint2D"]^2) < 0.95)
    cutlineData <- cbind(x$rollcalls[constrained, paste("midpoint",
        dims[1], "D", sep = "")], x$rollcalls[constrained, paste("spread",
        dims[1], "D", sep = "")], x$rollcalls[constrained, paste("midpoint",
        dims[2], "D", sep = "")], x$rollcalls[constrained, paste("spread",
        dims[2], "D", sep = "")],x$rollcalls[constrained,"col"])
    ##cutlineData <- cbind(cutlineData,col)
    cutlineData <- na.omit(cutlineData)
    suppressWarnings(symbols(x = 0, y = 0, circles = 1, inches = FALSE,
        asp = 1, main = main.title, xlab = d1.title, ylab = d2.title,
        xlim = c(-1, 1), ylim = c(-1, 1), cex.main = 1.2, cex.lab = 1.2,
        font.main = 2, lwd = 2, fg = "grey", frame.plot = FALSE,
        ...))
    if (lines < dim(cutlineData)[1])
        cutlineData <- cutlineData[sample(1:dim(cutlineData)[1],
            lines), ]
    suppressWarnings(apply(cutlineData, 1, add.cutline2, weight = x$weights[dims[2]]/x$weights[dims[1]],lwd = lwd,col.vec=col.vec))
  }



add.cutline2 <- function (cutData, weight, lwd = 2,col,col.vec) {
  slope <- -cutData[2]/(cutData[4] * weight)
    if (is.na(slope)) {
        x <- c(cutData[1], cutData[1])
        y <- c(sqrt(1 - cutData[1]^2), -sqrt(1 - cutData[1]^2))
        slope <- NA
        intercept <- NA
    }
    else {
        intercept <- -slope * cutData[1] + cutData[3]
        x <- c((-slope * intercept + sqrt((slope * intercept)^2 -
            (1 + slope * slope) * (intercept * intercept - 1)))/(1 +
            slope * slope), (-slope * intercept - sqrt((slope *
            intercept)^2 - (1 + slope * slope) * (intercept *
            intercept - 1)))/(1 + slope * slope))
        if (is.na(x[1])) {
            warning("Couldn't solve for points on the unit circle!\n")
            x <- NA
            y <- NA
            slope <- NA
            intercept <- NA
        }
        else {
            y <- intercept + slope * x
            y[y < -1] <- -sqrt(1 - x[y < 1]^2)
            y[y > 1] <- sqrt(1 - x[y > 1]^2)
          }
      }
    lines(x, y, lwd = lwd, col=col.vec[cutData[length(cutData)]])
  }


plot.wnom2 <- function(x,col,col.vec) {
  with(data.frame(x$legislators),{
    plot(coord1D,coord2D,type="n",xlim=c(-1.2,1.2),ylim=c(-1.2,1.2))
    ##plot.cutlines2(portugat.w,add=TRUE,col="red")
    plot.cutlines2(x,add=TRUE,
                   ##lines=800,
                   col.vec=col.vec,col=col)
                   ###col=data[id.data$id,j])
    ##text(coord1D,coord2D,labels=abbreviate(sport,10),cex=1.2)
  })
}


meansd <- function(x) {
  mx <- mean(c(x),na.rm=TRUE)
  sx <- sd(c(x),na.rm=TRUE)
  if (is.na(sx)) sx <- 0
  c(mean=mx,sd=sx)
}




ptable <- function(x) prop.table(table(x))




packages <- c("Amelia","AnnotationDbi","Biobase","Biostrings","BradleyTerry","CDNmoney","Cairo","CarbonEL","DAAG","DBI","Design","DynDoc","Ecdat","FAiR","FactoMineR","FinTS","GO.db","GPArotation","Hmisc","KEGG.db","KernSmooth","MASS","MBESS","MCMCpack","MEMSS","MLDS","MNP","MSBVAR","Matching","Matrix","Mcomp","MiscPsycho","PBSmapping","PTAk","R.matlab","R.methodsS3","R.oo","R.utils","R2WinBUGS","R2jags","RArcInfo","RBGL","RColorBrewer","RGtk2","RMySQL","ROC","ROCR","RODBC","RSQLite","RSvgDevice","RTisean","RUnit","RWeka","RaschSampler","Rdbi","RdbiPgSQL","Rmpi","SensoMineR","Snowball","SparseM","VGAM","XML","Zelig","abind","accuracy","acepack","ada","ade4","ade4TkGUI","affy","affyPLM","affyQCReport","affydata","affyio","akima","amap","anacor","annaffy","annotate","ape","apsrtable","arm","aroma.light","arules","bayesm","betareg","bitops","bivpois","boot","bootstrap","brglm","ca","caTools","cacheSweave","cairoDevice","car","cba","cfa","chron","class","cluster","cocorresp","coda","codetools","coin","colorspace","combinat","concord","corpcor","crq","ctv","degreenet","deldir","digest","distr","distrEx","dlm","doBy","dse1","dse2","dyn","dynlm","e1071","eRm","eba","ecodist","ellipse","energy","ergm","evd","expsmooth","fBasics","fCalendar","fEcofin","fImport","fSeries","fUnitRoots","fUtilities","fda","ff","fields","filehash","fma","forecast","foreign","fpc","gWidgets","gWidgetsRGtk2","gam","gcrma","gdata","gee","geepack","genefilter","geneplotter","ggplot","ggplot2","gnm","gpclib","gplots","graph","gsl","gstat","gtools","hgu95av2.db","homals","hopach","ifa","igraph","ineq","irr","irtoys","its","jit","kernlab","kinship","labdsv","latentnet","lattice","leaps","limma","lme4","lmtest","locfit","ltm","mFilter","mapproj","maps","maptools","marray","matchprobes","mcgibbsit","mclust","mclust02","meboot","memisc","mgcv","micEcon","mice","misc3d","mitools","mlbench","mlica","mlmRev","mmlcr","modeltools","mokken","moments","mprobit","msm","multcomp","multcompView","multilevel","multinomRob","multtest","mvnmle","mvnormtest","mvtnorm","nFactors","network","networksis","nlme","nnet","np","numDeriv","nws","odesolve","odfWeave","oz","pan","partitions","party","pastecs","pcaMethods","pcaPP","pgirmess","pixmap","plRasch","plink","plm","pls","plyr","pmml","poLCA","polycor","polynom","prabclus","preprocessCore","profileModel","proto","proxy","pscl","psy","psych","psychometric","psyphy","pwt","quadprog","quantreg","qvcalc","rJava","randomForest","randomLCA","randomSurvivalForest","rattle","reldist","relimp","reshape","rgdal","rgenoud","rggobi","rgl","rjags","rlecuyer","robustbase","rpanel","rpart","rpvm","rrcov","rv","sandwich","scapeMCMC","scatterplot3d","sde","segmented","sem","setRNG","sfsmisc","shapes","simpleaffy","simpleboot","sm","smacof","sna","snow","snowFT","sp","spam","spatial","spatstat","spdep","splancs","startupmsg","stashR","statmod","statnet","strucchange","survey","survival","systemfit","tensor","tensorA","tframe","tkrplot","tm","trimcluster","tripack","tsDyn","tseries","tseriesChaos","tsfa","urca","uroot","vars","vcd","vegan","vsn","waveslim","wnominate","xgobi","xtable","zicounts","zoo","KernSmooth","MASS","R2WinBUGS","Rmpi","base","boot","class","cluster","coda","codetools","datasets","foreign","grDevices","graphics","grid","lattice","methods","mgcv","nlme","nnet","rpart","snow","spatial","splines","stats","stats4","survival","tcltk","tools","utils")

jags2coda <- function(out) {
    nc <- dim(out[[1]])
    nc <- nc[length(nc)]
    ans <- vector("list", nc)
    for (ch in 1:nc) {
        ans.ch <- vector("list", length(out))
        vnames.ch <- NULL
        for (i in seq(along = out)) {
            varname <- names(out)[[i]]
            d <- dim(out[[i]])
            if (length(d) < 3) {
                stop("Invalid dimensions for sampled output")
            }
            vardim <- d[1:(length(d) - 2)]
            nvar <- prod(vardim)
            niter <- d[length(d) - 1]
            nchain <- d[length(d)]
            values <- as.vector(out[[i]])
            var.i <- matrix(NA, nrow = niter, ncol = nvar)
            for (j in 1:nvar) {
                var.i[, j] <- values[j + (0:(niter - 1)) * nvar + 
                  (ch - 1) * niter * nvar]
            }
            vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
            ans.ch[[i]] <- var.i
        }
        ans.ch <- as.matrix(data.frame(ans.ch))
        colnames(ans.ch) <- vnames.ch
        ans[[ch]] <- mcmc(ans.ch, start = 1, thin = 1)
    }
    mcmc.list(ans)
}


whos <- function() sort( sapply(ls(),function(x){object.size(get(x))}))
