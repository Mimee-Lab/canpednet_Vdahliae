library(vcfR)
library(adegenet)
library(factoextra)
library(poppr)
library(ade4)
library(tidyverse)
library(ggtree)
library(tidytree)
library(treeio)

# set ggplot theme
theme_set(theme_bw())

# path to vcf file
vcf.path <- "../canpednet_verti/data/full_vcf/verti_wgs_freebayes.snpeff.vcf"

#creating vcfR object
vcf.full <- read.vcfR(vcf.path)

# extract province from pop name
pop.data <- colnames(vcf.full@gt)[-1]
pop.data <- pop.data %>%
  as_tibble %>% 
  rename(pop = value) %>% 
  mutate(province = str_remove_all(pop,"[0-9]+.*")) %>% 
  mutate(province = factor(province, level = c("AB","MB","ON","QC","NB","PEI","V_alfalfae"))) %>%
  dplyr::select(pop, province)

# samples by province
pop.data %>% 
  count(province)

#### vcf filtering ####


# remove all loci with variant QUAL of 0
var.id.qual <- vcfR2tidy(vcf.full,info_only = TRUE)$fix %>% 
  mutate(variantid = paste0(CHROM, "_", POS)) %>% 
  select(variantid, QUAL) %>% 
  filter(QUAL > 20) %>% 
  pull(variantid)



var.id.full <- extract.gt(vcf.full, "GT") %>%
  as_tibble(rownames = "variantid") %>%
  pull(variantid)

vcf <- vcf.full[which(var.id.full %in% var.id.qual),]


#### vcf analysis ####


# build genind object
gi <- vcfR2genind(vcf)
#gl <- vcfR2genlight(vcf) # light version of genind



### summary vcf
print(vcf)

### summary genind
print(gi)


### Tree 
set.seed(2023)
tree.gi <- aboot(gi)

#Write tree
treeio::write.tree(tree.gi, file = "data/new_data/tree_aboot_gi_vcf.tree")

# read tree
tree <- treeio::read.tree("data/new_data/tree_aboot_gi_vcf.tree")

ggt <- tree %>%
  as_tibble() %>%
  left_join(., pop.data, by = c("label" = "pop")) %>%
  mutate(province = factor(province, levels = c("AB","MB","ON","QC","NB","PEI"))) %>% 
  as.treedata() 

# print node label
ggtree(ggt) +
  geom_text(aes(label = node), size = 3)


### Split lineage
vdahliae.node <- 193
lin1.node <- 311
lin2.node <- 194






lin.df <- ggt %>% 
  as_tibble() %>% 
  drop_na(province) %>% 
  mutate(lineage = if_else(parent >= lin1.node, "L1","L2")) %>% 
  dplyr::select(label, lineage, province)


### Build lineage specific genind
l1 <- lin.df %>% filter(lineage == "L1") %>% pull(label) 
l2 <- lin.df %>% filter(lineage == "L2") %>% pull(label)

l1.gi <- gi[rownames(gi@tab) %in% l1 ,]
l2.gi <- gi[rownames(gi@tab) %in% l2,]



### Write lineage specific genind
write_rds(l1.gi, "data/new_data/l1.genind.rds")
write_rds(l2.gi, "data/new_data/l2.genind.rds")



#removing loci containing NA. Takes a while
l1.nona <- missingno(l1.gi, cutoff = 0, type = "loci")
l2.nona <- missingno(l2.gi, cutoff = 0, type = "loci")

#save no NA object
write_rds(l1.nona, "data/new_data/l1.nona.RDS")
write_rds(l2.nona, "data/new_data/l2.nona.RDS")
