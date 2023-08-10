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

setwd("~/3.Input")
filepath ="~/TaxaplotOutput"

## read in the data
phyloseq1 = read_rds("Campbell_20230328_TOP_OTU.RDS")
phyloseq2 = read_rds("Fraser_20230328_TOP_OTU.RDS")
phyloseq3 = read_rds("Hu_20230328_TOP_OTU.RDS")
phyloseq4 = read_rds("Wang_20230328_TOP_OTU.RDS")
phyloseq = merge_phyloseq(phyloseq1, phyloseq2, phyloseq3, phyloseq4)

phyloseq = subset_taxa(phyloseq, genus != "NA")

## get the number of reads in each sample in your phyloseq object
phyloseq@sam_data$number_of_reads = sample_sums(phyloseq)

## get the hu data out of phyloseq
## I always call this dataframe mot because it's the Metadata, OTU table, and Taxonomy, in that order
mot = dephyloseq(phyloseq)
#sample_variables(mot)

## get relative abundance of otu
#remove NA
mot = subset(mot, asv_abundance >= 0)
#print(sort(unique(mot$relative_abundance)))

#calculate relative abundance
mot$relative_abundance = as.numeric(mot$asv_abundance)/as.numeric(mot$number_of_reads)

#round salinity values
#change mot$salinity to temperature for taxaplots with temperature
mot$Sal = 0.0
maxRows = nrow(mot)

{
  for(i in 1:maxRows) {    
    rounded <- round(as.double(mot$salinity[i]), digits = 0)
    mot$Sal[i]=rounded }
}

#print(mot$Sal)

## make taxa names for plot column
mot$plot_names = paste0(mot$order, "; ", mot$genus)


##rename study for nice axis titles
mot$renamed = ""
{
  renamestudy <- function(study) {
    if (study == "wang") {
      result <- "Chesapeake Bay"
    }
    else if (study == "hu2016") {
      result <- "Baltic Sea"
    }
    else if (study == "chen") {
      result <- "Fraser River"
    }
    else if (study == "campbell2013") {
      result <- "Delaware Bay"
    }
    else {
      result <- "NA"
    }
    return(result)
  }
  
  maxRows = nrow(mot)
  
  for(i in 1:maxRows) {    
    renames <- renamestudy(as.character(mot$study[i]))
    mot$renamed[i]=renames
  }
}


##### LOOP TO SELECT TOP 15 GENERA IN SAMPLES ######
## make variable to track sample types and salinities
mot$group = paste0(mot$renamed)#, "-", mot$Region)

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
top.df$place = "top_15"


##### SET UP TAXA PLOT #####
## join the top taxa and existing mot dataframe
mot.top = full_join(mot, (top.df))

## make the empty "place" cells say bottom
mot.top$place = replace(mot.top$place, is.na(mot.top$place), "bottom")

## replace plot_names for bottom taxa with Other
mot.top[mot.top$place == "bottom",]$plot_names <- "Others"

## order mot.top by decreasing relative abundance
mot.top = mot.top[order(-mot.top$relative_abundance),]#(mot.top$salinity)]

## get list of factors in order
natural.genus.order = as.list(c(unique(mot.top$plot_names)))

## remove others from list #!#
no.others=natural.genus.order[!natural.genus.order == 'Others']

## add Others to end of list
plot.order = append(no.others, "Others")

## set plot_names levels
plot.order = unlist(plot.order)

## order dataframe by relative abundance
mot.top$plot_names = factor((mot.top$plot_names), levels=c(plot.order))

alltaxa <- unique(mot.top$plot_names)


### NEED TO LOAD ALL DATASETS FIRST ###
##### MAKE COLOR LIST #####
# write a csv with the names
#write.csv(alltaxa, "allcolors.csv")
allcolors = read.csv("allcolors.csv")
colors <- filter(allcolors, x %in% plot.order)

##### MAKE TAXAPLOT #####
## get list to cycle through
mot.top$taxaplotgroup = paste0(mot.top$renamed)#, "-", mot.top$salinity)
taxaplot.groups = unique(mot.top$taxaplotgroup)

#changing order
mot.top <- mot.top[order(as.double(mot.top$Sal)),]

## field taxaplot loop
for (i in taxaplot.groups){
  for (j in i){
    
    sub.df = subset(mot.top, mot.top$taxaplotgroup== c(j)) 
      a <- (unique(sub.df$plot_names))
      colorsa <- filter(allcolors, x %in% a)
      colorsa$x <- factor(colorsa$x, levels=plot.order)
      colorsa <- colorsa[order(colorsa$x),]
      scolorsa <- colorsa$plot_colors
    
    myplot=ggplot(sub.df, aes(x=as.character(Row.names), 
                              y=as.numeric(relative_abundance), 
                              fill=as.factor(plot_names)))+
      geom_bar(stat = "identity")+
      theme(text = element_text(size=10),
            legend.title = element_text(size=10, face = "bold"),
            legend.text = element_text(size=5),
#            legend.position = "none",
            strip.text.x = element_text(size=5),
            axis.text.y = element_text(size=10),
            axis.title.y = element_text(size=14),
            title = element_text(size=10, face = "bold"))+
      scale_fill_manual(values=scolorsa)+
      facet_grid(.~format(as.double(Sal), nsmall = 0), scales="free", space="free")+
      labs(x= "", 
           y="Relative abundance", 
#           title= j,
           fill="Order; Genus")+
      guides(fill=guide_legend(ncol=3))
    
    
    myplot
    
    ## save plot
    ggsave(myplot, filename=paste0(j,"_taxaplot_fix",".png",sep=""), width=15, heigh=5)
  }
}



