### Script explaining how and why to make reproducible figures

### loading R packages ## also includes installation information

library(tidyverse) #install.packages("tidyverse") # this one takes a long time to install...
library(ggpubr) #install.packages("ggpubr")
#library(devtools) #install.packages("devtools")
#library(patchwork) #devtools::install_github("thomasp85/patchwork")
#library(limma)  #BiocManager::install("limma")

#setwd("~/iclouddrive/Documents/labmeetings/2020_tutorials/reproducible-figs")
#setwd("...") ## set to your working directory (this repository)

## Read in data
## Cheating here...this is already analyzed data
## We went from peptide intensities to the following 4 files
# tidied peptide level GSVA scores
pep_saha <- read.table("./data/saha_gsva_peptide.csv", header=T, row.names=1, sep="\t") %>%
  as.data.frame()
# tidied protein level GSVA scores
prot_saha <- read.table("./data/saha_gsva_protein.csv", header=T, row.names=1, sep="\t") %>% 
  as.data.frame()

## Significantly altered KEGG pathways 
## peptide level data
sig_pep_saha <- read.table("./data/saha_sig_peps.txt", header=T, sep="\t") %>% 
  as.data.frame()
## protein level data
sig_prot_saha <- read.table("./data/saha_sig_prot.txt", header=T, sep="\t") %>% 
  as.data.frame()


### shaping data into format for plotting 
## What are the significantly altered KEGG pathways identified at both levels?
sig_saha <- c(sig_pep_saha$x %>% as.character(), sig_prot_saha$x %>% as.character()) %>% 
  unique()

## Subsetting/reoriganizing peptide GSVA scores to match sig_saha
## wide -> long format for ggplot
pep_saha <- pep_saha[pep_saha$Pathways %in% sig_saha,] %>% 
  reshape2::melt(id = "Pathways", variable.name = "Sample", value.name = "GSVA")

## Subsetting/reoriganizing protein GSVA scores to match sig_saha
## wide -> long format for ggplot
prot_saha <- prot_saha[prot_saha$Pathways %in% sig_saha,] %>% 
  reshape2::melt(id = "Pathways", variable.name = "Sample", value.name = "GSVA")

## Merging peptide and protein data
saha <- merge(pep_saha, prot_saha, by=c("Pathways", "Sample"))
# Changing long name to short
saha <- saha %>% mutate(Sample =ifelse(
  Sample == "Xu_20181212_SAHA_rapidAIM_Low_4_181214210722", "Low_4", as.character(Sample)))
# Adding condition information
saha$Condition <- factor(substr(saha$Sample, 1, nchar(saha$Sample %>% as.character())-2))

# Finding median values per conditon
saha <- saha %>% dplyr::group_by(Pathways, Condition) %>% 
  dplyr::summarise(med_pep = median(GSVA.x), med_prot = median(GSVA.y))

# Completing the correlation test
saha_cor <- cor.test(saha$med_pep, saha$med_prot %>% as.numeric(), method=c("pearson"))


## Plot options
pointcolour <- "#363635"
#pointcolour <- "purple"
linecolour <- "#FFBA49"
#linecolour <- "green"

(cor1 <- ggplot(saha, aes(x=med_pep, y = med_prot)) +
    geom_point(shape=16, alpha = 0.9, size = 5, color = pointcolour) + 
    stat_cor(method = "pearson", size = 8) + 
    geom_smooth(method=lm, se=T, color = linecolour, size = 2) + theme_bw(base_size=26) +
    xlab("GSVA enrichment scores (peptides)") + ylab("GSVA enrichment scores (proteins)")) 
ggsave("./figs/saha_cor_example.pdf", p, width = 5, height = 5, units = "in")

(cor2 <- ggplot(saha, aes(x=med_pep, y = med_prot)) +
    geom_point(shape=16, alpha = 0.9, size = 10, color = pointcolour) + 
    stat_cor(method = "pearson", size = 8) + 
    geom_smooth(method=lm, se=T, color = linecolour, size = 2) + theme_dark(base_size=26) +
    xlab("PEPTIDE GSVA") + ylab("PROTEIN GSVA")) 


##### Combining figures...
## Putting heatmaps together
## Analysis in ./scripts/gsvaProteinLevel.R and ./scripts/gsvaPeptideLevel.R
sahapep_data <- read.table("saha_pepplotdata.txt", header=T, row.names = 1, sep = "\t") %>% as.data.frame()
sahaprot_data <- read.table("saha-prot_plotdata.txt", header=T, row.names = 1, sep = "\t") %>% as.data.frame()
sahapep_dataclust<- read.table("saha-pep_dataclust.txt", sep="\t", skip=1)
sahaprot_dataclust<- read.table("saha-prot_dataclust.txt", sep="\t", skip=1)

intersect(sahapep_data$Pathway %>% unique(), sahaprot_data$Pathway %>% unique()) %>% length()

saha_data <- rbind(sahapep_data, sahaprot_data) %>% 
  mutate(Type =  c(rep("Peptide", nrow(sahapep_data)), rep("Protein", nrow(sahaprot_data)))) #%>% na.omit()

saha_data$Pathway <- factor(saha_data$Pathway %>% as.character(), levels = c(sahapep_dataclust$V1 %>% as.character(), 
                                                                             sahaprot_dataclust$V1 %>% as.character()) %>% 
                              unique()) #%>% na.omit()





saha_data <- saha_data %>% mutate(GSVA = ifelse(significance == 0 & Conditions != "DMSO", NA, GSVA))
#saha_data <- saha_data[complete.cases(saha_data), ]
