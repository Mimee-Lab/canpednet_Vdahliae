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

tree <- read.tree("data/tree/tree_aboot_gi_vcf.tree")
pop.data <- read_tsv(file = "data/pop.data.tsv")



### Split lineage
vdahliae.node <- 193
lin1.node <- 311
lin2.node <- 194



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


### Full tree

p.vdah <- ggtree(ggt,size = 0.2) +
  geom_label2(aes(subset=node == lin1.node), fill = "lightgrey", label = "Lineage 1",
              hjust = 1.5, vjust = -0.2) + 
  geom_label2(aes(subset=node == lin2.node), fill = "lightgrey", label = "Lineage 2",
              hjust = 1.5, vjust = -0.2) 

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
  geom_label(aes(label = paste0("n=",prov_n)), y = 1, vjust = -0.1,
             color = "black", fill = "grey", size = 2.5)+
  geom_col() + 
  scale_fill_manual(values = c("#c8c8c8","#646464"))+
  scale_y_continuous(breaks=seq(0,1.2,0.25),limits=c(0, 1.1))+
  labs(y = "Sample ratio")

p.lin.ratio


image.map = image_scale(map.raw,geometry = "4000x1425")
print(image.map)

windowsFonts("Arial" = windowsFont("Arial"))


fig <- image_graph(width = 4200, height = 3500, res = 600)
ggarrange(ggplot()+theme_void(),
          ggpubr::ggarrange(p.vdah, p.lin.ratio, ncol = 2,labels = c("B","C"),font.label = list(family = "Arial",font = "Arial",size = 14)),
          nrow = 2, labels = c("A",""),font.label = list(family = "Arial",size = 14,face="bold"))

dev.off()
out <- image_composite(image = fig, composite_image=image.map,offset = "+100-200")


image_ggplot(out) %>% 
  ggsave(filename = "fig_1.tiff", path = "results",width = 7, height = 5.833,
         units = "in",device = "tiff",dpi = 600,compression = "lzw")






