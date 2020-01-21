## want to compare the peptide-based gsva to the protein based....
## make this a general script please
library(tidyverse)
library(data.table)
library(plyr)
library(DESeq2)
library(GSVA)
library(limma)
library(reshape2)
library(qvalue)

## Xu's data
prot_data <- "~/projects/annotation/SAHA_rapidAIM/proteinGroups.txt" 
out_file <- "../figs/gsva_prot_SAHA.pdf"
control_cond <- "DMSO"
condition1 <- "Low"
condition2 <- "High"
condition_opts <- c(control_cond, condition1, condition2)


### Leyuan's data
prot_data <- "~/iclouddrive/Documents/Manuscripts/figeys/peptide_centric/data/20190819_MG_MP/proteinGroups.txt"
out_file <- "../figs/gsva_prot_mtfm.pdf"
control_cond <- "CTRL"
condition1 <- "MTFM"
condition_opts <- c(control_cond, condition1)

### kegg
gene_kegg <- fread("~/iclouddrive/Documents/Manuscripts/figeys/peptide_centric/annotation/prot_kegg_annotated.txt.gz")
## Data
# Kegg full database information (L4 kegg with L3 pathways for gene set enrichment...)
kegg_L3 <- read.delim("~/iclouddrive/Documents/Manuscripts/figeys/peptide_centric/annotation/kegg_L3.txt",
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc'),
                      colClasses=c('character','character','character','character')) %>% as.data.frame()
## Creating gene sets
pathways <- kegg_L3$L4_desc %>% unique() 
pathways <- pathways[-c(208:229)] #removing KO that are not in brite/pathway
pathway_kegg <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))
### functions
match_pathway <- function(df){
  subset <- protidInfo[protidInfo$V5 %in% df$L3,]
  subset <- subset$newIDname
  return(subset)
} 


exp_data <- read.delim(prot_data, sep = '\t') %>%
  as.data.frame() %>%  dplyr::select(., Protein.IDs, Number.of.proteins, starts_with('LFQ'))
#write.table(exp_data, file="../testing/test_proteins.txt", sep="\t", row.names = T, col.names = T, quote=F)

## MTFR median
techrep <- as.matrix(exp_data[,10:12]) 
techrep_median <- rowMedians(techrep)
exp_data <- exp_data[,-c(10:12)]
exp_data <- data.frame(exp_data, LFQ.intensity.MTFM_3median=techrep_median)
#


## Remove peptides that are only sometimes identified...
## Data filtering function
## https://datascienceplus.com/proteomics-data-analysis-2-3-data-filtering-and-missing-value-imputation/
filter_valids = function(df, conditions, min_count, at_least_one = TRUE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  df <- df %>% select(-Number.of.proteins) %>% column_to_rownames(., var = "Protein.IDs")
  df[df==0] <- NA
  all_names <- colnames(df) 
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, all_names, value = TRUE, perl = TRUE))
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  #return(df) # No rows are omitted, filter rules are listed in the KEEP column
  df[is.na(df)] <- 0 
  return(df %>%  rownames_to_column(., var='Protein.IDs') %>% filter(KEEP) %>% dplyr::select(-KEEP) %>%
           column_to_rownames(., var='Protein.IDs'))  # only keeping rows that meet the criteria!
}


## Apply filtering
exp_data = filter_valids(exp_data,
                         conditions = condition_opts,
                         #min_count = c(4,5,4), #want peptide to have been identified in at least half of the samples 
                          min_count = c(5,5),
                         at_least_one = TRUE)  #but if it is consistently identified in a condition, keep it


#######
# PCA #
#######

log_data <- data.frame(exp_data) %>% dplyr::select(starts_with('LFQ')) %>%
  mutate_all(., funs(log2(1 + .)))


pca<- prcomp(t(log_data), center=T, scale=F)
sampleVals<-data.frame(pca$x)
exprVals<-data.frame(pca$rotation)
PoV <- (pca$sdev^2/sum(pca$sdev^2))*100
conditions <- c(rep("CTRL", 5), rep("MTFM", 5))
conditions <- colnames(log_data) %>% substr(., 15, nchar(.)-2)
#conditions[10] <- "MTFM"
conditions[13] <- "Low"
coords<-data.frame(sampleVals, Condition = conditions,
                   samplename = rownames(sampleVals))
numPCs <- 1:length(PoV)

for (i in 1:length(PoV)) {
  percent <- paste0("(", round(PoV[i],2), "%)")
  name <- paste0("PC", i, "per")
  assign(name, percent)
}

(mtfm_pcaplot <- ggplot(coords, aes(x = PC1, y = PC2)) + #accept selectInput to choose axes!
    geom_point(size=8, aes(fill=Condition, shape=Condition)) + 
    stat_ellipse(geom = "polygon", alpha=.2, aes(color=Condition, fill=Condition)) +
    #scale_color_manual(values=c(input$control_col, input$cond1_col, input$cond2_col)) + #pick colours for colour picker
    scale_color_manual(values=c("#FFBA49", "#20A39E")) +
    scale_fill_manual(values=c("#FFBA49", "#20A39E")) +
    scale_shape_manual(values=c(22, 21)) + #24 after
    #scale_shape_manual(values=shapes2use) +
    scale_x_continuous(name= paste0("PC1", " ", PC1per))+ # labels depend on selected PCs
    scale_y_continuous(name= paste0("PC2", " ", PC2per))+ theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=24)) )
#ggsave("../figs/mtfmprotpca.pdf", mtfm_pcaplot,  width = 9.5, height = 8.7, units = "in")

(saha_pcaplot <- ggplot(coords, aes(x = PC1, y = PC2)) + #accept selectInput to choose axes!
    geom_point(size=8, aes(fill=Condition, shape=Condition)) + 
    stat_ellipse(geom = "polygon", alpha=.2, aes(color=Condition, fill=Condition)) +
    #scale_color_manual(values=c(input$control_col, input$cond1_col, input$cond2_col)) + #pick colours for colour picker
    scale_color_manual(values=c("#FFBA49", "#20A39E", "#363635")) +
    scale_fill_manual(values=c("#FFBA49", "#20A39E", "#363635")) +
    scale_shape_manual(values=c(22, 21, 24) )+ #24 after
    #scale_shape_manual(values=shapes2use) +
    scale_x_continuous(name= paste0("PC1", " ", PC1per))+ # labels depend on selected PCs
    scale_y_continuous(name= paste0("PC2", " ", PC2per))+ theme_bw() +
    theme(legend.position = "bottom", legend.title = element_blank(), text=element_text(size=24)) )
#ggsave("../figs/sahaprotpca.pdf", saha_pcaplot, width = 9.5, height = 8.7, units = "in")


exp_data <- exp_data %>% rownames_to_column("Protein.IDs")
exp_data$protgroupID <- seq.int(nrow(exp_data)) %>% as.character()

protidInfo <- data.frame(IDnum = exp_data$protgroupID, protID = exp_data$Protein.IDs)

protidInfo <- separate_rows(protidInfo, c('protID'), sep = ';')
protidInfo <- merge(protidInfo, gene_kegg, by.x='protID', by.y="V1", all.x=T) 
#test
x <- protidInfo %>%
  dplyr::group_by(IDnum, V5) %>% dplyr::select(IDnum, V5)  %>%
  dplyr::count() %>% drop_na() 

protidInfo <-  x %>% dplyr::group_by(IDnum) %>%
  dplyr::mutate(total = sum(n), prop = n/total)
protidInfo$newIDname <- make.names(protidInfo$IDnum,unique=T) ##unque protID names so we can multiply the intensities...
 
log_exp <- exp_data %>% 
  merge(., protidInfo, by.x = 'protgroupID', by.y = 'IDnum', all.x = T) %>%
  mutate(prop=replace(prop, is.na(prop), 1)) %>%
  mutate_each(funs(.*prop), starts_with('LFQ')) %>% #multiplies the intensities by the proportion 
  mutate(correct_ID = case_when(is.na(newIDname) ~ protgroupID,
                                 !is.na(newIDname) ~ newIDname)) %>%
  dplyr::select(starts_with('LFQ'), correct_ID) %>%
  mutate_each(funs(log2(.+1)), starts_with('LFQ')) %>%
  column_to_rownames('correct_ID')
 
 colnames(log_exp) <- substr(colnames(log_exp), 15, nchar(log_exp))
 ## must be changed by hand
 #colnames(log_exp)[13] <- "Low_4"
 
 ## normalize by library size
 #norm_pep <- estimateSizeFactorsForMatrix(log_exp)
 #exp_data <- sweep(log_exp, 2, norm_pep, "/")
 
 
 ## data is now all set up for GSVA
 
 ## GSVA gene sets must be created
 pathway_genesets <- lapply(pathway_kegg, match_pathway) 
 
 ## running GSVA
 gsva_kegg <- gsva(as.matrix(log_exp),pathway_genesets, min.sz=10,
                   kcdf='Gaussian')
 
# gsva_tofile <- gsva_kegg %>% as.data.frame %>% rownames_to_column("Pathways")
# write.table(gsva_tofile, file="mtfm_gsva_protein.csv", sep="\t", row.names = T, col.names = T, quote=F)
conditions <- colnames(log_exp) %>% substr(., 1, nchar(.)-2)
conditions[13] <- "Low"
conditions[10] <- "MTFM"
#conditions <- c(rep(control_cond, 5), rep(condition1, 7))
 cond <- factor(conditions) %>% relevel(control_cond)
 #print(cond)
 design <- model.matrix(~  cond) # we are comparing all to DMSO which is our control
#
colnames(design)[2:ncol(design)] <- substr(colnames(design)[2:ncol(design)], 5, 
                                            nchar(colnames(design)[2:ncol(design)])) #just removing "cond"
colnames(design)[1] <- control_cond
 
 #design <- model.matrix(~  cond)
 #colnames(design) <- c(control_cond, condition1)
 fit <- lmFit(gsva_kegg, design)
 fit <- eBayes(fit, trend=T)
 allGeneSets <- topTable(fit, number=Inf)
 DEgeneSets <- topTable(fit,  number=Inf,
                        p.value=0.05, adjust="BH")
 res <- decideTests(fit, p.value=0.05, adjust = "BH")
 summary(res)
 sigdrugs <- res[,abs(res) %>% colSums(.) > 0]
 
 res <- res %>% as.data.frame()
 sig_tests <- res[abs(res[,2:3]) %>% rowSums(.) > 0,]
 sig_tests <- res[abs(res[,2])  > 0,]
 sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(sig_tests),]
 gsva_tofile <- sig_gsva %>% as.data.frame %>% rownames_to_column("Pathways")
 #write.table(gsva_tofile$Pathways, file = "saha_sig_prot.txt", sep = "\t", row.names = F, col.names = T, quote = F)
 write.table(gsva_tofile$Pathways, file = "mtfm_sig_prot.txt", sep = "\t", row.names = F, col.names = T, quote = F)
 
 gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
   reshape2::melt(., id='Pathway', variable.name = "Samples", value.name="GSVA") 
 
 gsvaplot_data <- gsvaplot_data %>% mutate(Conditions = case_when(
   str_detect(Samples, control_cond) ~ control_cond,
   str_detect(Samples, condition1) ~ condition1)) #, #) # , 
   str_detect(Samples, condition2) ~ condition2))
 clusterdata <- rownames(gsva_kegg)[hclust(dist(gsva_kegg))$order]
# write.table(clusterdata, file="mtfm-prot_dataclust.txt", sep = " ", append = FALSE, dec = ".",
#            row.names = F, col.names = F, quote=F)
 
#write.table(clusterdata, file="saha-prot_dataclust.txt", sep = " ", append = FALSE, dec = ".",
#             row.names = F, col.names = F, quote=F) 
gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
 
 

clustersamples <- colnames(gsva_kegg)[hclust(dist(gsva_kegg %>% t()))$order]
 #write.csv(clustersamples, file="mtfm-prot_sampleclust.txt", row.names = T, col.names = T, quote=F)
 #write.csv(clustersamples, file="saha-prot_sampleclust.txt", row.names = T, col.names = T, quote=F) 
 gsvaplot_data$Samples<- factor(gsvaplot_data$Samples, levels = clustersamples)
 
 sig_tests <- sig_tests %>% rownames_to_column(., var = "Pathway") %>%  
   reshape2::melt(. ,id="Pathway", variable.name= "Conditions", value.name = "significance")
 gsvaplot_data <- gsvaplot_data %>% merge(., sig_tests, by=c("Pathway", "Conditions")) %>%
   dplyr::mutate(GSVA = ifelse(significance == 0 & Conditions != control_cond, NA, GSVA))

# write.table(gsvaplot_data, file="mtfm-prot_plotdata.txt", sep="\t", row.names = T, col.names = T, quote=F)
# write.table(gsvaplot_data, file="saha-prot_plotdata.txt", sep="\t", row.names = T, col.names = T, quote=F)
 
 
# sigpathways <- as.data.frame(sigdrugs %>% abs())   
# sigpathways <- sigpathways[sigpathways > 0,, drop=F] %>% as.data.frame() %>% 
#   drop_na() %>% rownames_to_column(., var='Pathway')     
# 
# sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% sigpathways$Pathway,]
# gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
#   melt(., id='Pathway') %>%
#   mutate(Condition = case_when(
#     str_detect(variable, control_cond) ~ control_cond,
#     str_detect(variable, condition1) ~ condition1))
# 
# clusterdata <- rownames(gsva_kegg)[hclust(dist(gsva_kegg))$order]
# gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
#
# clustersamples <- colnames(gsva_kegg)[hclust(dist(gsva_kegg %>% t()))$order]
# gsvaplot_data$Samples<- factor(gsvaplot_data$Samples, levels = clustersamples)
 
 (proteinhm <- ggplot(data = gsvaplot_data, mapping = aes(x = Samples, y = Pathway, fill = GSVA)) + 
   facet_grid(~ Conditions, switch='x', scales = "free") +
   scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
   #scale_fill_gradientn(colours=c(input$low_col, "white", input$high_col),
                        space = "Lab", name="GSVA enrichment score", na.value="#64686F") + 
   geom_tile(na.rm = TRUE) +
   xlab(label = "Sample") +
   ylab(label="") +
   theme(axis.text.x = element_text(angle = 45, hjust = 1)) )
 
 #ggsave(out_file, width=7.78, height=2.5, units = 'in')
                   