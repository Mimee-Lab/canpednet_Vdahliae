library(adegenet)
library(ade4)
library(ggplot2)
library(poppr)
library(tidyverse)

setwd("N:/Actif-Active/J-002391_CanPEDnet/")

#load genind objects
L1 <- readRDS("Data/processed_data/verti_wgs/genind/l1.genind.rds")
L2 <- readRDS("Data/processed_data/verti_wgs/genind/l2.genind.rds")


#loading meta file
popmap <- read.table("Data/processed_data/verti_wgs/genind/pop.data2.tsv", header = T) %>% filter(complete.cases(.)) %>% 
                                                                        mutate(across(c(pop,province,lineage),factor))


#removing loci containing NA. Takes a while
l1.nona <- missingno(L1, cutoff = 0, type = "loci")
#Found 51522 missing values.
#1819 loci contained missing values greater than 0%

#save nona object
saveRDS(l1.nona, "l1.nona.RDS")

#find clusters
l1.nona_Clust <- find.clusters(l1.nona, max.n.clust=40)#PC = 40, clusters = 3
head(l1.nona_Clust$Kstat, 20)
head(l1.nona_Clust$grp, 74)#see in which group each samples are.
head(l1.nona_Clust$size, 20)#groups sizes



d.dapcL1 <- dapc(l1.nona, l1.nona_Clust$grp)#PC = 40, discr. func = 3
scatter(d.dapcL1, grp=l1.nona_Clust$grp,clab = 0.85,
        posi.da="bottomright", scree.pca = F, clabel = 0) #posi.pca = "bottomleft"


#MSN plot by genetic clusters
library(igraph)
pop(l1.nona) <- l1.nona_Clust$grp
popNames(l1.nona) <- c("group_3", "group_1", "group_2")

l1.nona_dist_clust <- bitwise.dist(l1.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)
l1_min_span_net_clust <- poppr.msn(l1.nona, l1.nona_dist_clust, showplot = FALSE, include.ties = TRUE)


set.seed(48)
#jpeg("l1_msn_norep_clusters.jpeg", width = 480, height = 620, quality = 99)
plot_poppr_msn(l1.nona,
               l1_min_span_net_clust,
               inds = "NONE",
               gadj = 800,
               gweight = 1,
               glim = c(0,0.8),
               nodescale = 60,
               palette = pal,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = T,
               layfun = igraph::layout_with_fr)
#dev.off()


#MSN plot by Locations
pop(l1.nona) <- popmap %>% filter(.,lineage=="L1") %>% pull(province)
popNames(l1.nona) <- c("AB", "MB", "NB", "ON", "PEI", "QC")

l1.nona_dist <- bitwise.dist(l1.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)
l1_min_span_net <- poppr.msn(l1.nona, l1.nona_dist, showplot = FALSE, include.ties = TRUE)
#
library(RColorBrewer)
#jpeg("l1_msn_norep_location.jpeg", width = 480, height = 620, quality = 99)
set.seed(48)
plot_poppr_msn(l1.nona,
               l1_min_span_net,
               inds = "NONE",
               gadj = 800,
               gweight = 1,
               glim = c(0,0.8),
               nodescale = 60,
               palette = pal,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = T,
               layfun = igraph::layout_with_fr)
#dev.off()

#Joel conversion list to ggplot object for multipanel option

theme_set(theme_classic())

layout <- layout_with_kk(l1_min_span_net_clust$graph) %>%   
  as_tibble() %>%   
  mutate(name = names(l1.nona_Clust$grp), grp = l1.nona_Clust$grp) %>%
  rename(x = V1, y = V2) %>%   
  mutate(province = str_remove(name, "[0-9]+.*"))

g <- get.data.frame(l1_min_span_net_clust$graph) %>%   
  as_tibble() %>%   
  left_join(.,layout,by = c("from"="name")) %>%   
  left_join(.,layout, by = c("to" = "name"),suffix = c("_from","_to"))

#graph for groups
a <- ggplot(layout, aes(x = x, y = y, color = grp)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to), inherit.aes = FALSE,size = 2)+  
  geom_point( size = 3) +  scale_color_manual(values = c("blue","red","darkgrey"))+
  labs(x="", y="", color="Group")

#graph for provinces
b <- ggplot(layout, aes(x = x, y = y, color = province)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to), inherit.aes = FALSE,size = 2)+  
  geom_point( size = 3) +  scale_color_manual(values = pal)+
  labs(x="", y="", color="Province")


#######################lineage 2##########################

#removing loci containing NA
l2.nona <- missingno(L2, cutoff = 0, type = "loci")
#Found 75465 missing values.
#1949 loci contained missing values greater than 0%

#save nona object
saveRDS(l2.nona, "l2.nona.RDS")


#find clusters
l2.nona_Clust <- find.clusters(l2.nona, max.n.clust=40)#PC = 60, k = 6
head(l2.nona_Clust$Kstat, 20)
head(l2.nona_Clust$grp, 118)#see in what group each samples are.
head(l2.nona_Clust$size, 20)#4 11 29 63  2  9. nb in each of the 6 groups



d.dapcL2 <- dapc(l2.nona, l2.nona_Clust$grp)#PC = 40, discr. func = 3
scatter(d.dapcL2, grp=l2.nona_Clust$grp,clab = 0.85,
        posi.da="bottomleft", scree.pca = F, clabel = 0) #posi.pca = "bottomleft"



#MSN plot by genetic clusters
library(igraph)
pop(l2.nona) <- l2.nona_Clust$grp
popNames(l2.nona) <- c("group_3", "group_4", "group_2","group_6","group_5","group_1")

l2.nona_dist_clust <- bitwise.dist(l2.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)
l2_min_span_net_clust <- poppr.msn(l2.nona, l2.nona_dist_clust, showplot = FALSE, include.ties = TRUE)
#
#pdf("l2_msn_norep_clusters.pdf", width = 8, height = 6)
set.seed(48)
plot_poppr_msn(l2.nona,
               l2_min_span_net_clust,
               inds = "NONE",
               gadj = 800,
               gweight = 1,
               glim = c(0,0.8),
               nodescale = 60,
               palette = pal,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = T,
               layfun = igraph::layout_with_fr)
#dev.off()


#MSN plot by Locations
pop(l2.nona) <- popmap %>% filter(.,lineage=="L2") %>% pull(province)
popNames(l2.nona) <- c("QC", "PEI", "AB", "NB", "ON", "MB")

l2.nona_dist <- bitwise.dist(l2.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)
l2_min_span_net <- poppr.msn(l2.nona, l2.nona_dist, showplot = FALSE, include.ties = TRUE)
#
library(RColorBrewer)
pal <- brewer.pal(6, "Set1")
#pdf("l2_msn_norep_location.pdf", width = 8, height = 6)
set.seed(48)
plot_poppr_msn(l2.nona,
               l2_min_span_net,
               inds = "NONE",
               gadj = 800,
               gweight = 1,
               glim = c(0,0.8),
               nodescale = 60,
               palette = pal,
               pop.leg = TRUE,
               size.leg = FALSE,
               scale.leg = T,
               layfun = igraph::layout_with_fr)
#dev.off()

#Joel conversion list to ggplot object

theme_set(theme_classic())

layout <- layout_with_kk(l2_min_span_net_clust$graph) %>%   
  as_tibble() %>%   
  mutate(name = names(l2.nona_Clust$grp), grp = l2.nona_Clust$grp) %>%
  rename(x = V1, y = V2) %>%   
  mutate(province = str_remove(name, "[0-9]+.*"))

g <- get.data.frame(l2_min_span_net_clust$graph) %>%   
  as_tibble() %>%   
  left_join(.,layout,by = c("from"="name")) %>%   
  left_join(.,layout, by = c("to" = "name"),suffix = c("_from","_to"))

#graph for groups
c <- ggplot(layout, aes(x = x, y = y, color = grp)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to), inherit.aes = FALSE,size = 2)+  
  geom_point( size = 3) +  scale_color_manual(values = c("blue","red","darkgrey", "darkgreen", "orange","grey25"))+
  labs(x="", y="", color="Group")

#graph for provinces
d <- ggplot(layout, aes(x = x, y = y, color = province)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to), inherit.aes = FALSE,size = 2)+  
  geom_point( size = 3) +  scale_color_manual(values = pal)+
  labs(x="", y="", color="Province")


#multiplot arrangment
library(ggpubr)

ggarrange(a,b,c,d,ncol = 2,nrow = 2, labels = c("A","B","C","D"))
ggsave("multipanel.jpeg", width = 14, height = 8,dpi = 300)

#####################################################################################################


#adding pop slot
pop(L1) <- popmap %>% filter(.,lineage=="L1") %>% pull(province)
pop(L2) <- popmap %>% filter(.,lineage=="L2") %>% pull(province)



#Na
sum(is.na(L1$tab))
#51522

sum(is.na(L2$tab))
#75465

#remove NA
LL1 <- scaleGen(L1, NA.method="mean")
class(LL1)#"matrix" "array"

LL1_gi <- df2genind(LL1, ploidy = 1, type = "codom", sep = "\t")

#find clusters
LL1_Clust <- find.clusters(LL1, max.n.clust=40)#PC = 20, k = 4
head(LL1_Clust$Kstat, 20)
head(LL1_Clust$grp, 74)#see in what group each samples are.
head(LL1_Clust$size, 20)#1  1 63 9. nb in each of the 4 groups


#remove NA
LL2 <- scaleGen(L2, NA.method="mean")
class(LL2)##"matrix" "array"

#find clusters
LL2_Clust <- find.clusters(LL2, max.n.clust=40)#PC = 20, k = 2
head(LL2_Clust$Kstat, 20)
head(LL2_Clust$grp, 119)#see in what group each samples are.
head(LL2_Clust$size, 21)#1 118. nb in each of the 21 groups. doit oter outgroup



d.dapcL1 <- dapc(LL1, LL1_Clust$grp)#PC = 40, discr. func = 3
scatter(d.dapcL1, grp=LL1_Clust$grp,clab = 0.85,
        posi.da="bottomleft", scree.pca = F, clabel = 0) #posi.pca = "bottomleft"


d.dapcL2 <- dapc(LL2, LL2_Clust$grp)#PC = 40, discr. func = 3
scatter(d.dapcL2, grp=LL2_Clust$grp, clab = 0.85,
        posi.da="topright", scree.pca = F, clabel = 0) #posi.pca = "bottomleft"




pcaL1 <- dudi.pca(LL1,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pcaL1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

s.label(pcaL1$li)
title("L1")
add.scatter.eig(pcaL1$eig[1:20], 3,1,2,posi = "topright")

#group names avec dernier de chaque cat?gorie. Modif .gen original
s.class(pcaL1$li,pop(L1)) #x = genind object used for pop names clustering
add.scatter.eig(pcaL1$eig[1:20], 3,1,2, posi = "topright")






######OnlyQc######
x <- read.genepop("populations.snps.genepop.gen")


table(pop(x))


###essai PCA 17 nov 20###
sum(is.na(x$tab))
#416402 NA's

#diff results si use tab() function...
X <- scaleGen(x, NA.method="mean")
class(X)

#find clusters
#NumClust <- find.clusters(X, max.n.clust=20)
#head(NumClust$Kstat, 20)
#plot(head(NumClust$Kstat, 20))


pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

s.label(pca1$li)
title("onlyQc")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

#group names avec dernier de chaque cat?gorie. Modif .gen original
s.class(pca1$li, pop(x)) #x = genind object used for pop names clustering
add.scatter.eig(pca1$eig[1:20], 3,1,2, posi = "topleft")


#s.class(pca1$li,pop(x),xax=1,yax=3,sub="PCA 1-3",csub=2)
#title("OnlyQc axes 1-3")
#add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)


col <- deepseasun(12)
s.class(pca1$li, pop(x),xax=1,yax=2, col=transp(col,0.7), axesell=FALSE,
        cstar=0, cpoint=3, grid=TRUE, label = NA, sub = pop(x))


colorplot(pca1$li, pca1$li, transp=TRUE, cex=2, xlab="PC 1", ylab="PC 2")
title("GBS2")
abline(v=0,h=0,col="grey", lty=2)
legend(pop(x))


#colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
#title("OnlyQc axes 1-3")
#abline(v=0,h=0,col="grey", lty=2)



ggplot(pca1$li, aes(x=Axis1, y=Axis2,color=pop(x)))+
  geom_point(aes(color=pop(x)), size=3)+
  #geom_text(aes(label=pop(x)))+
  geom_hline(yintercept=0,linetype=3) + 
  geom_vline(xintercept=0,linetype=3) +
  scale_color_discrete(name="Groups")+
  ggtitle("                                                   merge_2pop")



d.dapc <- dapc(x, n.pca = 20, n.da = 2)
scatter(d.dapc, clab = 0.85, col = funky(18),
        posi.da="topright", scree.pca = F) #posi.pca = "bottomleft"



#13192_138.01 13192_138.04 46139_146.01 46139_146.02 57849_154.01 57849_154.03


temp    <- seploc(x)     
snp13192 <- tab(temp[["13192_138"]])

freq1041 <- apply(snp13192, 2, function(e) tapply(e, pop(x), mean, na.rm = TRUE))

