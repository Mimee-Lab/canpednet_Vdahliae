library(vcfR)
library(adegenet)
library(factoextra)
library(poppr)
library(ade4)
library(tidyverse)


# set ggplot theme
theme_set(theme_bw())

# path to vcf file
vcf.path <- "../canpednet_verti/data/full_vcf/verti_wgs2_freebayes.snpeff.vcf"

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

# remove all pop with same genotype
var.id <- extract.gt(vcf.full, "GT") %>%
  as_tibble(rownames = "variantid") %>%
  filter(!if_all(-variantid, ~ .x == 0)) %>%  
  filter(!if_all(-variantid, ~.x ==1)) %>%
  pull(variantid)

var.id.full <- extract.gt(vcf.full, "GT") %>%
  as_tibble(rownames = "variantid") %>%
  pull(variantid)

vcf <- vcf.full[which(var.id.full %in% var.id),]


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

write.tree(tree.gi, file = "analysis/wgs2/tree_aboot_gi_fullvcf2.tree")


### Split lineage

ggt %>% 
  as_tibble() %>%
  filter(node < 194) %>% 
  mutate(lineage = if_else(parent >= lin1.node, "L1", "L2")) %>% 
  mutate(lineage = if_else(label == "V_alfalfae", NA_character_, lineage)) %>% 
  select(pop = label , province, lineage) %>% 
  write_tsv("analysis/wgs2/pop.data2.tsv")

pop.data <- read_tsv("analysis/wgs2/pop.data2.tsv")

l1 <- pop.data %>% filter(lineage == "L1") %>% pull(pop) 
l2 <- pop.data %>% filter(lineage == "L2") %>% pull(pop)

l1.gi <- gi[rownames(gi@tab) %in% l1 ,]
l2.gi <- l2.gi[rownames(l2.gi@tab) %in% l2,]




write_rds(l1.gi, "l1.genind.rds")
write_rds(l2.gi, "l2.genind.rds")



#removing loci containing NA. Takes a while
l1.nona <- missingno(L1, cutoff = 0, type = "loci")


#save nona object
saveRDS(l1.nona, "l1.nona.RDS")