library(adegenet)
library(ade4)
library(ggplot2)
library(poppr)
library(tidyverse)


#save nona object
l1.nona <- read_rds("data/genind/l1.nona.RDS")
l2.nona <- read_rds("data/genind/l2.nona.RDS")

#find clusters L1
l1.nona_Clust <- find.clusters(l1.nona, max.n.clust=40,n.pca=40) #create BIC values versus nb of clusters plot




#find clusters L2
l2.nona_Clust <- find.clusters(l2.nona, max.n.clust=40, n.pca = 60) #create BIC values versus nb of clusters plot


