library(ggplot2)
library(colorblindr)
library(grid)
library(gridExtra)
library(coda)

# clear workspace
rm(list = ls())
# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


doCompare <- function(dsFileName, mcmcFileName) {
    dfs <- read.table(paste0("simulator/", dsFileName), header=T)
    df <- read.table(paste0("operators/", mcmcFileName), header=T)

    # Remove 10% for burnin
    N <- dim(df)[1]
    df <- df[-(1:ceiling(0.1*N)),]

    # Tree age/length plot
    pdf(gsub(".log", ".pdf", mcmcFileName), width = 10, height = 10) 
    

    par(mfcol=c(3,1))

    maxLength <- max(quantile(df$network.totalLength, probs=0.99),
                     quantile(dfs$network.totalLength, probs=0.99),
                     quantile(df$network.height, probs=0.99),
                     quantile(dfs$network.height, probs=0.99))

    maxDensity <- max(density(dfs$network.height)$y,
                      density(df$network.height)$y,
                      density(dfs$network.totalLength)$y,
                      density(df$network.totalLength)$y)

    plot(density(df$network.height), 'l', col='red', lwd=2, lty=2,
         xlim=c(0,maxLength), ylim=c(0,maxDensity),
         xlab="Statistic", ylab="Density",
         main="")
    lines(density(dfs$network.height), col='blue', lwd=2, lty=2)
    lines(density(df$network.totalLength), col='red', lwd=2)
    lines(density(dfs$network.totalLength), col='blue', lwd=2)
    legend("topright", inset=0.05,
       c("MCMC", "Direct simulation","Network length", "Network height"),
       lty=c(1,1,1,2), lwd=2, col=c("red","blue","black","black"))

    # Network node count plot

    maxCount <- max(df$network.reassortmentNodeCount, dfs$network.reassortmentNodeCount)+1
    h <- hist(df$network.reassortmentNodeCount, plot=F,
              breaks=seq(-0.5,maxCount+0.5,by=1))
    hs <- hist(dfs$network.reassortmentNodeCount, plot=F,
               breaks=seq(-0.5,maxCount+0.5,by=1))

    maxDensity <- max(h$density, hs$density)

    plot(h$mids, h$density, 'o', col='red', lwd=4,
        xlab="Reassortment Count",
        ylab="Posterior probability",
        ylim=c(0,maxDensity))
    lines(hs$mids, hs$density, 'o', col='blue', lwd=2)

    legend("topright", inset=0.04,
          c("MCMC", "Direct simulation"), lty=1, pch=1, lwd=2,
          col=c("red","blue"))

  # Kolmogorov-Smirnov test for network statistics

  H <- c()
  L <- c()
  R <- c()
  l <- length(df$network.height)

  # For comparison to run faster, calculations need to be done at fewer ponts.
  # It can be achieved by taking a fraction smaller than 0.01 below.
  nPoints <- round(l*0.01)
  by <- as.integer(l/nPoints)

  for (i in seq(l, 0, by=-by)){
      H <- append(H, log(ks.test(df$network.height[1:(l-i)], dfs$network.height)$statistic))
      L <- append(L, log(ks.test(df$network.totalLength[1:(l-i)], dfs$network.totalLength)$statistic))
      R <- append(L, log(ks.test(df$network.reassortmentNodeCount[1:(l-i)], dfs$network.reassortmentNodeCount)$statistic))
  }
  plot(seq(1, length(L)), L, type='l', col='red',  lwd=2,
       xlab=paste("Number of logged iterations (minus burnin) / ", by),
       ylab="Kolmogorov-Smirnov stat. (log scale)", 
       main="Kolmogorov-Smirnov test for network statistics",
       ylim=c(min(L, H, R), max(L, H, R)))

  lines(seq(1, length(H)), H, col='blue', lwd=2, lty=2)
  lines(seq(1, length(R)), R, lwd=2, lty=3)
  legend("topright", inset=0.05,
         c("Network length", "Network height", "Reassortment count"),
         lty=c(2,2, 3), lwd=2, col=c("red","blue","black"))
  dev.off()
}

doCompare("simulate_serial5taxon8seg.log",
          "testAddRemove_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+Unif_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+DS_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+NS_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+NarrowExchange_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+WideExchange_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+Slide_serial5taxon8seg.log")

doCompare("simulate_serial5taxon8seg.log",
          "testAR+Gibbs_serial5taxon8seg.log")


