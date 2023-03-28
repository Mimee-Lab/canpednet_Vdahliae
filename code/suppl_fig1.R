library(tidyverse)
library(ggtree)
library(tidytree)#load after tidyverse

#vcf allpop
tree <- read.tree("data/tree/bautista2020_merge_tree.txt")

#bar indicator per node
#if pop.data not loaded
pop.path <- "data/pop.data.tsv"
pop.data <- read.table(pop.path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#create sample_origin column
tidy.tree <- tree %>%
as_tibble() %>%
mutate(sample_origin = if_else(label %in% pop.data$pop,"canpednet","Bautista")) %>%
as.treedata()


#tree plot
pf <- ggtree(tidy.tree) + 
  theme_tree() +
  geom_tippoint(aes(color = sample_origin), size = 1, shape =15)+
  #geom_tiplab() + 
  geom_cladelabel(node=378, label="4B", 
                  color="black", offset=0.01,offset.text = 0.01,hjust = 0.2, align = T) + #\nAB=2\nMB=1\nON=8\nQC=18\nPEI=48\nNB=41
  geom_cladelabel(node=364, label="2A", 
                  color="black", offset=0.01,offset.text = 0.01,hjust = 0.2) +
  geom_cladelabel(node=599, label="4A", 
                  color="black", offset=0.01,offset.text = 0.01,hjust = 0.2) + #\nAB=24\nMB=10\nON=14\nQC=1\nPEI=15\nNB=10
  geom_cladelabel(node=695, label="2B/1", 
                  color="black", offset=0.01,offset.text = 0.01,hjust = 0.2) +
  xlim(0, 0.58) +
  scale_color_manual(labels = c("Bautista-Jalón*et al*. (2021)","This study"),
                       values = c("red", "blue"))+
  theme(legend.text = ggtext::element_markdown(),
        legend.title = element_blank())+
  guides(color = guide_legend(override.aes = list(size = 3)))

pf
