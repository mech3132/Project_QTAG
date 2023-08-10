### SET UP ###
{
  ## load packages
  library(phyloseq)
  library(tidyverse)
  library(plyr)
  library(ggplot2); theme_set(theme_bw()+
                                theme(panel.grid = element_blank(),
                                      strip.background = element_rect(fill="white"),
                                      axis.text.y = element_text(size = 12, colour = "black"),
                                      axis.title = element_text(size=15, face="bold"),
                                      strip.text = element_text(color="black", size=10),
                                      legend.text=element_text(size=10),
                                      axis.line = element_line(colour = "black"),
                                      axis.text.x = element_blank(),))
  
  ## load function
  dephyloseq = function(phylo_obj){
    
    ## get the metadata
    meta = as.data.frame(as.matrix(phylo_obj@sam_data))
    
    ## how many metadata columns you have 
    metacols = ncol(meta)+1
    
    ## get out the otu table 
    ## if your metadta is empty after running this, you need to use otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
    otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
    
    ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
    mo = merge(meta, otu, by=0)
    
    ## get out the taxonomy file 
    tax = as.data.frame(phylo_obj@tax_table)
    
    ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
    tax = tax %>% rownames_to_column(var="asv_id")
    
    ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
    mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance")
    
    ## Join the metadata and otu table with the taoxnomy table 
    mot = full_join(mo, tax)
    
    ## Specify the output for the dephyloseq funciton 
    output = mot
  }
  

}

setwd("~/2. GetTopTaxa")

## read in the data
#change read_rds for other datasets
phyloseq = read_rds("fraser_20230320_salinity_ranges_filtered_not_rareified.RDS")
phyloseq = subset_taxa(phyloseq, genus != "NA")


## get the number of reads in each sample in your phyloseq object
phyloseq@sam_data$number_of_reads = sample_sums(phyloseq)

## subset the phyloseq object to only keep the samples you want
sample_names(phyloseq)
phyloseq = subset_samples(phyloseq)

## get the hu data out of phyloseq
## I always call this dataframe mot becayse it's the Metadata, OTU table, and Taxonomy, in that order
mot = dephyloseq(phyloseq)

## get relative abundance of otu
#removes the Na 
mot = subset(mot, (mot$asv_abundance > 0) > 0)


#calculate relative abundance
mot$relative_abundance = as.numeric(mot$asv_abundance)/as.numeric(mot$number_of_reads)

#round salinity values
mot$Sal = 0.0

{
  maxRows = nrow(mot)
  for(i in 1:maxRows) {
    rounded <- round(as.double(mot$salinity[i]), digits = 0)
    rounded <- format(rounded, nsmall = 0)
    mot$Sal[i]=rounded }
}

print(mot$Sal)


## make taxa names for plot column
mot$plot_names = paste0(mot$order, "; ", mot$genus)



##### LOOP TO SELECT TOP 15 GENERA IN SAMPLES ######
## make variable to track sample types and salinities
mot$group = paste0(mot$Range)#, "-", mot$Region)

## summarize mot
mot.sum = ddply(mot, c("group", "plot_names"),
                   summarise,
                   N=length(relative_abundance),
                   sum = sum(relative_abundance))

## get list of all samples
samplegroups = unique(mot.sum$group)
samplegroups

## sort data by relative abundance
sorted = mot.sum[order(-mot.sum$sum),]

mot.top = mot

## make empty dataframe to store taxa
top.df = NULL

## start loop
for(i in samplegroups) {
  for(j in i) {
    # subset dataframe by samples
    sample = subset(sorted, sorted$group %in% c(j))

    ## get top 15 genera
    top = sample[c(1:10),]

    # save list of top  abundance taxa
    t.tmp <- top
    top.df <- rbind.fill(top.df, t.tmp)

    ## close loop
  }
}

## add identifier for top and bottom taxa
top.df$place = "top_10"


##### SET UP TAXA PLOT #####
## join the top taxa and existing mot dataframe
mot.top = full_join(mot, top.df)
mot.top = subset(mot.top, place != "NA")

top_asv <- as.list(mot.top$asv_id)
OTU <- subset(otu_table(phyloseq), rownames(otu_table(phyloseq))%in% c(top_asv))
rephyloseqed <- phyloseq(OTU, tax_table(phyloseq), sample_data(phyloseq))


print(unique(top.df$plot_names))

#save
#change write_rds for other datasets
write_rds(rephyloseqed, "Fraser_20230328_TOP_OTU.RDS")