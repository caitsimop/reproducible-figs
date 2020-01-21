## Testing peptide-centric workflow!
## peptide gene set variation analysis (GSVA)
## enrichment in KEGG terms. KEGG intensity is summed for each sample
##    - weighted to the proportion of kegg assignment
##    - this is so we consider all possibilities of the function (and they are distributed 
##       by our confidence in the assigned function. )
library(tidyverse)
library(plyr)

# Set working directory 
setwd("~/projects/annotation/test/")

## Data
# Kegg full database information (L4 kegg with L3 pathways for gene set enrichment...)
kegg_L3 <- read.delim("~/projects/annotation/kegg/kegg_L3.txt", sep='\t', header=F, 
                      col.names=c('L3', 'L3_desc', 'L4', 'L4_desc'),
                      colClasses=c('character','character','character','character')) %>% as.data.frame()

# Cog fill database
cog_database <- read.delim("~/projects/annotation/cog/cognames2003-2014.tab", sep='\t', col.names = c('cog', 'func', 'desc')) %>%
  separate_rows(., 'func', sep='(?<=.)(?=.)')

# Peptides matched with core peptides 
core_test <- read.delim("core_peptides-kegg.txt", row.names = 1, na.strings=c(""," ","NA")) %>% as.data.frame()
as.character(core_test$kegg) %>% na.omit() %>% length() #8078
## only expression data
exp_data <- core_test %>% dplyr::select(c(starts_with("DMSO"), starts_with("High"), starts_with("Low")))
exp_annot <- core_test %>% dplyr::select(-c(starts_with("DMSO"), starts_with("High"), starts_with("Low")))
exp_data[exp_data==1] <-NA
sapply(exp_data, is.numeric) #checking that all the intensity columns are numeric
cond <- colnames(exp_data) %>% substr(., 1, nchar(.)-1)

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
  return(df %>%  rownames_to_column(., var='peptides') %>% filter(KEEP) %>% dplyr::select(-KEEP) %>%
         column_to_rownames(., var='peptides'))  # only keeping rows that meet the criteria!
}


## Apply filtering
exp_data = filter_valids(exp_data,
                     conditions = c('DMSO', 'High', 'Low'),
                     min_count = c(2, 3, 2), #want peptide to have been identified in at least half of the samples 
                     at_least_one = TRUE)  #but if it is consistently identified in a condition, keep it



## Creating gene sets
# first make a list that contains each pathway and corresponding 
pathways <- kegg_L3$L4_desc %>% unique() 
pathways <- pathways[-c(208:229)] #removing KO that are not in brite/pathway
pathway_kegg <- dlply(kegg_L3 %>% dplyr::select(L4_desc, L3), .(L4_desc))

pathway_cog <- dlply(cog_database %>% dplyr::select(func, cog), .(func))

##########
## KEGG ##
##########
# weigh intensity of functional annotation by the identified proportions
# all core peptides!
core_pep_kegg <- read.delim2("~/projects/annotation/db/core_annot/core_pep_kegg.csv", 
                             sep=",", header=F, col.names = c("pep", "kegg", "count", "eval"))
## get information on the number of proteins fora peptide...
core_pep_kegg_extrainfo <- core_pep_kegg %>% dplyr::group_by(pep) %>% 
  dplyr::summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  dplyr::mutate(prop = count/total) #%>% dplyr::select(pep, kegg, prop)

core_pep_kegg_extrainfo[grep("DGNKLDLYGK", core_pep_kegg_extrainfo$pep),]


core_pep_kegg <- core_pep_kegg %>% dplyr::group_by(pep) %>% 
  dplyr::summarize(total = sum(count))  %>%
  merge(., core_pep_kegg, by='pep', all.y=T) %>%
  dplyr::mutate(prop = count/total) %>% dplyr::select(pep, kegg, prop)



## update pep if it is a duplicate...
core_pep_kegg$newpep_name <- make.names(core_pep_kegg$pep,unique=T)
## do the same thing to the test data

#########
## COG ##
#########
core_pep_cog <- read.delim2("~/projects/annotation/db/core_annot/core_pep_cog.csv", 
                             sep=",", header=F, col.names = c("pep", "cog", "count"))
core_pep_cog <- core_pep_cog %>% dplyr::group_by(pep) %>% 
  dplyr::summarize(total = sum(count)) %>% merge(., core_pep_cog, by='pep', all.y=T) %>%
  dplyr::mutate(prop = count/total) %>% dplyr::select(pep, cog, prop)

core_pep_cog$newpep_name <- make.names(core_pep_cog$pep,unique=T)

## create a list that has the peptide name instead of kegg (just like "pathway_genesets)
# function for matching peptides with kegg
#match_pathway <- function(df){
#  subset <- core_test[core_test$kegg %in% df$L3,]
#  subset <- data.frame(pep = rownames(subset), kegg = subset$kegg)
#  return(subset)
#} 


## intensities normalized by the proportion of functional annotation
## if there are more than one functional annotation, the peptide will have a suffix added to the end (i.e. .1, .2, .3...etc)
## $newpep_name
test_core_kegg <- exp_data %>% 
  rownames_to_column(., var='pep') %>%
  merge(., core_pep_kegg, by='pep') %>% # this is how I know I'm not 
  mutate(prop=replace(prop, is.na(prop), 1)) %>%
  mutate_each(funs(.*prop), c(starts_with("DMSO"), starts_with("Low"), 
                              starts_with("High"))) %>% #multiplies the intensities by the proportion 
  column_to_rownames(., var='newpep_name') %>%
  dplyr::select(c(starts_with("DMSO"), starts_with("Low"), starts_with("High")))

test_core_cog <- exp_data  %>% 
  rownames_to_column(., var='pep') %>%
  merge(., core_pep_cog, by='pep') %>% # this is how I know I'm not 
  mutate(prop=replace(prop, is.na(prop), 1)) %>%
  mutate_each(funs(.*prop), c(starts_with("DMSO"), starts_with("Low"), 
                              starts_with("High"))) %>% #multiplies the intensities by the proportion 
  column_to_rownames(., var='newpep_name') %>%
  dplyr::select(c(starts_with("DMSO"), starts_with("Low"), starts_with("High")))


expression_data <- function(annotationtype){
  if (annotationtype=='COG'){
   exp_data <- test_core_cog %>% dplyr::select(starts_with("DMSO"), starts_with("High"), starts_with("Low")) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), c(starts_with("DMSO"), starts_with("High"), starts_with("Low"))) %>% ##should be log10 data...
      column_to_rownames(., var='peptide') 
  return(exp_data)
    } else if (annotationtype=='kegg'){
    exp_data <- test_core_kegg %>% dplyr::select(starts_with("DMSO"), starts_with("High"), starts_with("Low")) %>%
      rownames_to_column(., var='peptide') %>% 
      mutate_each(., funs(log(1 + .)), c(starts_with("DMSO"), starts_with("High"), starts_with("Low"))) %>% ##should be log10 data...
      column_to_rownames(., var='peptide')
  return(exp_data)    
  } else {
    stop("Not an accepted functional annotation type.")
    }}


exp_data <- expression_data('kegg')

## write a function that will do this for you!  
norm_pep <- estimateSizeFactorsForMatrix(exp_data)
exp_data <- sweep(exp_data, 2, norm_prot, "/") ##peptides are normalized
match_pathway <- function(df){
  subset <- core_pep_kegg[core_pep_kegg$kegg %in% df$L3,]
  subset <- subset$newpep_name
  return(subset)
} 
match_pathway_cog <- function(df){
  subset <- core_pep_cog[core_pep_cog$cog %in% df$cog,]
  subset <- subset$newpep_name
  return(subset)
} 

# applying function over our pathway list
pathway_genesets <- lapply(pathway_kegg, match_pathway) 
COG_genesets <- lapply(pathway_cog, match_pathway_cog) 
## filtering our data...did not actually do this
#nsFilter(leukemia_eset, require.entrez=TRUE, remove.dupEntrez=TRUE, var.func=IQR, var.filter=TRUE, var.cutoff=0.5, filterByQuantile=TRUE,
#         feature.exclude="^AFFX")


## this is testing only
gsva_kegg <- gsva(as.matrix(exp_data),pathway_genesets, min.sz=10,
                  kcdf='Gaussian') ## rnaseq=F because we have continuous data
cond<- factor(cond) %>% relevel("DMSO") # DMSO is the control
design <- model.matrix(~  cond)
colnames(design) <- c("DMSO", "HighVsDMSO", "LowVsDMSO")
fit <- lmFit(gsva_kegg, design)
fit <- eBayes(fit, trend=T)
allGeneSets <- topTable(fit, coef=c("HighVsDMSO", "LowVsDMSO"), number=Inf)
DEgeneSets <- topTable(fit, coef=c("HighVsDMSO", "LowVsDMSO"), number=Inf,
              p.value=0.05, adjust="BH")
res <- decideTests(fit, p.value=0.05)
summary(res)

## How to plot this??
## plot the GSVA statistics (KS ranking statistics) for each sample... (the significant ones!)
## only significant differences in High vs DMSO

sig_gsva <- gsva_kegg[rownames(gsva_kegg) %in% rownames(DEgeneSets),]
gsvaplot_data <- data.frame(sig_gsva) %>% rownames_to_column(., var="Pathway") %>%
  melt(., id='Pathway') 
gsvaplot_data$condition <- substr(gsvaplot_data$variable, 1, nchar(as.character(gsvaplot_data$variable))-1)

clusterdata <- rownames(sig_gsva)[hclust(dist(sig_gsva))$order]
gsvaplot_data$Pathway<- factor(gsvaplot_data$Pathway, levels = clusterdata)
(highplot <- ggplot(data = gsvaplot_data, mapping = aes(x = variable, y = Pathway, fill = value)) + 
    facet_grid(~ condition, switch = "x", scales = "free_x", space = "free_x") +
    scale_fill_gradientn(colours=c("#67A7C1","white","#FF6F59"),
                         space = "Lab", name="GSVA enrichment score") + 
    geom_tile(na.rm = TRUE) +
    xlab(label = "Sample") +
    ylab(label="") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave("gsva_heatmap.pdf", highplot)


## results are unexpected...why? Let's look at the expression of the unexpected results!
bactmot <- pathway_genesets[["Bacterial motility proteins [BR:ko02035]"]]
bactmotpep <- exp_data[rownames(exp_data) %in% bactmot,]
