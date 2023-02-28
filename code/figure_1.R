library(ggpubr)
library(tidyverse)
library(adegenet)
library(ape)
library(ggtree)
library(tidytree)
library(treeio)
library(magick)



theme_set(theme_classic())

map.raw <- image_read("data/map/map_Vdahliae_600dpi.tif")

tree <- read.tree("data/tree/tree_aboot_gi_fullvcf2.tree")
pop.data <- read_tsv(file = "data/pop.data.tsv")



out.node <- 1
vdahliae.node <- 195
lin1.node <- 313
lin2.node <- 196



ggt <- tree %>%
  as_tibble() %>%
  left_join(., pop.data, by = c("label" = "pop")) %>%
  mutate(province = factor(province, levels = c("AB","MB","ON","QC","NB","PEI"))) %>% 
  as.treedata() 

lin.df <- ggt %>% 
  as_tibble() %>% 
  drop_na(province) %>% 
  mutate(lineage == if_else(node >= lin1.node, "L1","L2")) %>% 
  dplyr::select(label, lineage)

ggt <- ggt %>% 
  as_tibble() %>% 
  left_join(., lin.df, by = "label") %>% 
  as.treedata()
  

# print node label
ggtree(ggt) +
  geom_text(aes(label = node), size = 3.5)





### Only V.dahliae
tree.vdah <- tree_subset(ggt, node = vdahliae.node, levels_back = 0)

p.vdah <- ggtree(tree.vdah)+
  geom_label2(aes(subset=node == lin1.node-1), fill = "lightgrey", label = "Lineage 1", hjust = 1.5) + 
  geom_label2(aes(subset=node == lin2.node-1), fill = "lightgrey", label = "Lineage 2", hjust = 1.5) 

p.vdah



# Lineage ratio
p.lin.ratio <- ggt %>% 
  as_tibble() %>%
  drop_na(province) %>% 
  mutate(lineage = if_else(parent >= lin1.node, "L1", "L2")) %>% 
  group_by(lineage, province) %>%
  summarize(n = n()) %>% 
  ungroup %>% 
  group_by(province) %>% 
  mutate(lin_ratio = n / sum(n), prov_n = sum(n)) %>% 
  ggplot(aes(x = province,fill = lineage,y = lin_ratio))+
  geom_label(aes(label = paste0("n=",prov_n)), y = 1, vjust = -0.1, color = "black", fill = "grey")+
  geom_col() + 
  ylim(c(0,1.1))+
  scale_fill_manual(values = c("#c8c8c8","#646464"))+
  labs(y = "Sample ratio")

p.lin.ratio


image.map = image_scale(map.raw, "1100")
print(image.map)


out <- image_composite(image = fig, composite_image=image.map,offset = "+50-50")
print(out)









