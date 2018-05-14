library(ggplot2)
library(reshape2)

df <- read.table("coalescent_density.log", header=T)


ggplot(df, aes(x=popSize, y=reassortmentRate)) +
    geom_raster(aes(fill=density)) +
    geom_contour(aes(z=density), colour="white", bins=20)
