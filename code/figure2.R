library(adegenet)
library(ade4)
library(ggplot2)
library(poppr)
library(tidyverse)
library(igraph)
library(RColorBrewer)
library(ggpubr)

#set colors
pal <- brewer.pal(6, "Set1")

#load genind objects
l1.nona <- readRDS("data/genind/l1.nona.RDS")
l2.nona <- readRDS("data/genind/l2.nona.RDS")

#loading meta file
popmap <- read.table("data/pop.data.tsv", header = T) %>% filter(complete.cases(.)) %>% 
  mutate(across(c(pop,province,lineage),factor))


#### Lineage 1 ####
#find clusters
l1.nona_Clust <- find.clusters(l1.nona, max.n.clust=40, n.pca=40, n.clust =3) #create BIC values versus nb of clusters plot if n.pca and n.clust parameters omitted

head(l1.nona_Clust$grp, 74)#see in which group each samples are.
head(l1.nona_Clust$size, 20)#groups sizes


#MSN plot by genetic clusters

pop(l1.nona) <- l1.nona_Clust$grp
popNames(l1.nona) <- c("group_3", "group_1", "group_2")

#dissimilarity and Euclidean distances calculation
l1.nona_dist_clust <- bitwise.dist(l1.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)

#Create a minimum spanning network of selected populations using a distance matrix
l1_min_span_net_clust <- poppr.msn(l1.nona, l1.nona_dist_clust, showplot = FALSE, include.ties = TRUE)


#Converting list to ggplot objects for easier manipulations
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
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to),
               inherit.aes = FALSE,linewidth = 1)+  
  geom_point( size = 1.5) +  
  scale_color_manual(values = c("darkgrey","red","blue"))+
  labs(x="", y="", color="Group")

#graph for provinces
b <- ggplot(layout, aes(x = x, y = y, color = province)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to),
               inherit.aes = FALSE,linewidth = 1)+  
  geom_point( size = 1.5) +  
  scale_color_manual(values = pal)+
  labs(x="", y="", color="Province")



#### Lineage 2 ####
#find clusters
l2.nona_Clust <- find.clusters(l2.nona, max.n.clust=40, n.pca=60, n.clust =6) #create BIC values versus nb of clusters plot if n.pca and n.clust parameters omitted

head(l2.nona_Clust$grp, 118)#see in which group each samples are.
head(l2.nona_Clust$size, 20)#groups sizes


#MSN plot by genetic clusters

pop(l2.nona) <- l2.nona_Clust$grp
popNames(l2.nona) <- c("group_3", "group_4", "group_2","group_6","group_5","group_1")

#dissimilarity and Euclidean distances calculation
l2.nona_dist_clust <- bitwise.dist(l2.nona, percent = TRUE, mat = FALSE, missing_match = TRUE, scale_missing = FALSE, euclidean = FALSE, differences_only = FALSE, threads = 0)

#Create a minimum spanning network of selected populations using a distance matrix
l2_min_span_net_clust <- poppr.msn(l2.nona, l2.nona_dist_clust, showplot = FALSE, include.ties = TRUE)


#Converting list to ggplot objects for easier manipulations
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
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to),
               inherit.aes = FALSE,linewidth = 1)+  
  geom_point( size = 1.5) +  
  scale_color_manual(values = c("blue","red","darkgrey", "darkgreen", "orange","grey25"))+
  labs(x="", y="", color="Group")

#graph for provinces
d <- ggplot(layout, aes(x = x, y = y, color = province)) +   
  geom_segment(data=g,aes(x=x_from,xend = x_to, y=y_from,yend = y_to),
               inherit.aes = FALSE,linewidth = 1)+  
  geom_point( size = 1.5) +  
  scale_color_manual(values = pal)+
  labs(x="", y="", color="Province")


#multiplot arrangment
windowsFonts("Arial" = windowsFont("Arial"))
ggarrange(a,b,c,d,ncol = 2,nrow = 2, labels = c("A","B","C","D"),
          font.label = list(family = "Arial",size = 14,face="bold"))
ggsave("fig_2.tiff",path = "results",compression = "lzw", width = 7, height = 5,
       units = "in",dpi = 600)
