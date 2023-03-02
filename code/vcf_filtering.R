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

df <- vcfR2tidy(vcf.full,info_only = TRUE,
                info_fields = c("DP","MQM","TYPE"))

df$fix

df$fix %>% 
  janitor::clean_names() %>% 
  mutate(type = str_remove_all(type, ",.*")) %>% 
  group_by(type) %>% 
  summarize(mqm_40 = sum(mqm > 40),
            qual_100 = sum(qual>100),
            total = n_distinct(paste0(chrom,pos)))

df$fix %>% 
  janitor::clean_names() %>% 
  mutate(type = str_remove_all(type, ",.*")) %>% 
  count(type)
