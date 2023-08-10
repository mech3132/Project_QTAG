##### WORKSPACE SETUP #####
library(phyloseq)
library(tidyverse)
library(vegan)
library(ranacapa)
library(stats)

#for merging
library(speedyseq)
library(janitor)

## tell R where to get the files
setwd("~/1. FilteringMergingAddingDescription")

## read in files
df1 = readRDS("hu_salinity_16S_unfiltered.RDS")
df2 = readRDS("wang_salinity_unfiltered_phyloseq.RDS")
df3 = readRDS("fraser_salinity_16S_unfiltered.RDS")
df4 = readRDS("campbell_salinity_16S_unfiltered.RDS")


######FIXING WANG ############
#rename columns to match
renamed2 <- df2 %>%
  rename_sample_data(salinity = water.salinity) %>%
  rename_sample_data(Lat = lat) %>%
  rename_sample_data(Long = long) %>%
  rename_sample_data(temperature = water.temperature) %>%
  rename_sample_data(Year = ena_year) %>%
  rename_sample_data(Season = ena_season) 

#view(renamed2@sam_data)
#To see columns for wang
renamed2 %>% sample_variables
#For comparison
df3 %>% sample_variables

#Add missing columns
## extract sample data dataframe from phyloseq object
whywang <- as.data.frame(as.matrix(sample_data(renamed2)))

#add columns
whywang$study = "wang"
whywang$Region = "Chesapeake"

## get your metadata 
samples = as.data.frame(as.matrix(whywang))
samples = sample_data(samples)

## you need to check if the rows are the ASV sequence (TRUE) or not (FALSE)
otu_mat <- as.data.frame(t(as.matrix(otu_table(renamed2))))
#View(otu_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat = renamed2@tax_table
TAX = tax_table(tax_mat)  

# make new phyloseq object
wangwang = phyloseq(OTU, TAX, samples)

#Merge data
df = merge_phyloseq(df1, wangwang, df3, df4)
#view(df@sam_data)

#to view
df %>% sample_variables




########## ADD SALINITY RANGES####################
#Add fresh, breackish, marine
##remove null salinity values
ww1 = subset_samples(df1, as.double(salinity) >= 0)
ww2 = subset_samples(wangwang, as.double(salinity) >= 0)
ww3 = subset_samples(df3, as.double(salinity) >= 0)
ww4 = subset_samples(df4, as.double(salinity) >= 0)
#combined data
ww = subset_samples(df, as.double(salinity) >= 0)

## extract sample data dataframe from phyloseq object
##### NEED TO MANUALLY CHANGE LINE 84 and 176 #####
study <- ww1
#transposes (use for 1, 3, 4)
{
#extract dataframe
whywater <- as.data.frame(as.matrix(sample_data(study)))
whywater$Range = 0


#view(whywater)
{
getRangeResult <- function(salinity) {
  if (salinity < 5) {
    result <- "fresh"
  }
  else if (salinity > 18) {
    result <- "marine"
  }
  else {
    result <- "brackish"
  }
  return(result)
}

maxRows = nrow(whywater)

for(i in 1:maxRows) {    
  rangeResult <- getRangeResult(as.double(whywater$salinity[i]))
  whywater$Range[i]=rangeResult
}
}



## get your metadata 
samples = as.data.frame(as.matrix(whywater))
samples = sample_data(samples)

## you need to check if the rows are the ASV sequence (TRUE) or not (FALSE)
otu_mat <- as.data.frame(t(as.matrix(otu_table(study))))
#View(otu_mat)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat = study@tax_table
TAX = tax_table(tax_mat)  
}

#doesnt transpose (use for 2)
{
  #extract dataframe
  whywater <- as.data.frame(as.matrix(sample_data(study)))
  whywater$Range = 0
  
  
  #view(whywater)
  {
    getRangeResult <- function(salinity) {
      if (salinity < 5) {
        result <- "fresh"
      }
      else if (salinity > 18) {
        result <- "marine"
      }
      else {
        result <- "brackish"
      }
      return(result)
    }
    
    maxRows = nrow(whywater)
    
    for(i in 1:maxRows) {    
      rangeResult <- getRangeResult(as.double(whywater$salinity[i]))
      whywater$Range[i]=rangeResult
    }
  }
  
  
  
  ## get your metadata 
  samples = as.data.frame(as.matrix(whywater))
  samples = sample_data(samples)
  
  ## you need to check if the rows are the ASV sequence (TRUE) or not (FALSE)
  otu_mat <- as.data.frame(as.matrix(otu_table(study)))
  #View(otu_mat)
  OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
  
  tax_mat = study@tax_table
  TAX = tax_table(tax_mat)  
}

# make new phyloseq object
phyloseq1 = phyloseq(OTU, TAX, samples)

#Merge data
phyloseq_all = merge_phyloseq(phyloseq1, phyloseq2, phyloseq3, phyloseq4)

#to view
#view(all_phyloseq@sam_data)
#all_phyloseq %>% sample_variables


##### FILTERING - REMOVE OFF-TAGET TAXA AND MAJOR CONTAMINATION FROM TAXONOMY FILE #####

#paste study of interest here, phyloseq1, 2, 3, 4, or all
df = phyloseq1

df = subset_taxa(df, domain == "Bacteria")
df = subset_taxa(df, order != "Chloroplast")
df = subset_taxa(df, family != "Mitochondria")
#view(df@tax_table)

{
##### FILTERING - REMOVE SAMPLES WITH LESS THAN N READS #####
# my N is usually >1000
## Remove samples with less than N reads 
sample_sums(df)
plot(sort(sample_sums(df)))

## add number of total sequences in a sample (Read depth) to the metadata 
df@sam_data$read_depth_noofftargets = sample_sums(df) 
sample_variables(df)
## check which samples have less than N reads
which(df@sam_data$read_depth_noofftargets < 1000) 

## remove the samples by your threshold 
df.pruned <- prune_samples(sample_sums(df) >= 1000, df)

## write file to know which samples were lost here. This is important for the methods section. 
#df.below1000 <- prune_samples(sample_sums(df) < 1000, df)
#df.below1000 = as.matrix(df.below1000@sam_data)
#write_rds(df.below1000, "Hu_salinity_pruned_SF_samples_less_than_1000.csv")

##### FILTERING - REMOVE INDIVIDUAL ASVS WITH LESS THAN N READS #####
## N is probably around 100 
## extract OTU dataframe from phyloseq object
otu.pruned <- as.data.frame(t(as.matrix(otu_table(df.pruned))))

## remove ASVs (rows) with less than N reads accross whole dataset but keep all samples
#!# make sure asv sequence is rownames and sample id is column name
otu.pruned$rowsum = rowSums(otu.pruned)

## remove low frequency asvs
otu.pruned = subset(otu.pruned, otu.pruned$rowsum>9)
#view(otu.pruned)

## remove rowsum column from your OTU table 
otu = subset(otu.pruned, select=-c(rowsum))


##### FILTERING - REMOVE ASVs FOUND IN N SAMPLES OR LESS #####
## your N is probably between 2 and 5

# has sample ID as column name and asv is as row name. Needs to be this way to use richness function 
# function to calculate richness, sums along a row (OTU)
richness = function(x){return(sum(x>0))}

## calculate richness on entire dataframe
otu$richness = apply(otu,1,richness) # use all columns of otu dataframe
summary(otu$richness)

## remove OTU (rows) with richness N or less (found in two samples or less) but keep all samples (columns)
otu = subset(otu, otu$richness>2)
## check that it worked
summary(otu$richness)
## remove richness column
otu = subset(otu, select=-c(richness))


##### DENOISING - MAKE ALL THE CELLS IN THE OTU TABLE WITH VALUES N OR LESS 0 #####
## your N is probably between 2 and 10

otu <- mutate_all(otu, funs(ifelse(. < 5, 0, .)))
#summary(otu)


##### CREATE AND READ BACK IN FILTERED BUT NOT RAREFIED PHYLOSEQ OBJECT ######

## format for phyloseq
otu_mat = t(as.matrix(otu))
## you need to check if the rows are the ASV sequence (TRUE) or not (FALSE)
OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
view(OTU)

## use the taxonomy table you had in your original phyloseq object
tax_mat = df.pruned@tax_table
  #as.data.frame(t(as.matrix(df.pruned@tax_table)))
TAX = tax_table(tax_mat)  

## get your metadata 
samples = as.data.frame(as.matrix(df.pruned@sam_data))

## get metadata ready for phyloseq
samples = sample_data(samples)

# make new phyloseq object
df.prerarefaction = phyloseq(OTU, TAX, samples)

#view(TAX)
# check that the phyloseq object was made correctly
df.prerarefaction
rank_names(df.prerarefaction)
sample_variables(df.prerarefaction)




##### FINAL DENOISING AND SAVE FILTERED DATA #####

## remove taxa from taxonomy table with empty otus (from denoising)
df.prerarefaction = prune_taxa(taxa_sums(df.prerarefaction)>0, df.prerarefaction)

## get final read depth
df.prerarefaction@sam_data$read_depth_filtered = sample_sums(df.prerarefaction)

## remove samples with very highest sample counts (might not be required for your dataset)
histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)
df.prerarefaction = subset_samples(df.prerarefaction, read_depth_filtered<300000)
histogram(df.prerarefaction@sam_data$read_depth_filtered, breaks=100)

## save your filtered but not rarefied phyloseq object 
#view(df.prerarefaction@sam_data)

}

view(df.prerarefaction@otu_table)
#write_rds(df.prerarefaction, "Hu_20230320_salinity_ranges_filtered_not_rareified.RDS", path = deprecated())