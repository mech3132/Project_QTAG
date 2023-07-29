#!bin/bash
library(ape)
# library(phytools)
# library(castor)
# library(VennDiagram)
library(MASS)
library(vegan)
library(tidyverse)
# VEN DIAGRAM


dir.create("07e_MLclustering")

###### Loading data ##########
tax <- read.delim("06_processing_QTAG/all_taxonomy.txt")
tre16 <- read.tree("06_processing_QTAG/tree16_filt.tre")
dat16 <- read.delim("06_processing_QTAG/QTAG_output_with_taxa.txt")
# dat16season <- read.delim("06_processing_QTAG_season/QTAG_output_with_taxa.txt") %>%
  # filter(TaxaID %in% tre16$tip.label) 
# tax16 <- tax %>% filter(TaxaID %in% tre16$tip.label)

# tre16.filt <- keep.tip(tre16, tip=dat16$TaxaID)
# remove(tre16)
# NH16 <- node.depth.edgelength(tre16.filt) # Get depths of all nodes
# treeH16 = max(NH16) # Get max height of tree

############ Cluster ASVs based on distribution #########
allStudies <- unique(dat16$Study)

if (file.exists("07e_MLclustering/asv_clustering.RData")) {
  load("asv_clustering.RData")
} else {
  asv_clustering <- data.frame(TaxaID=NA)
  # asv_clustering_LOW <- data.frame(TaxaID=NA)
  for ( s in allStudies) {
    mat16_salinity <- dat16 %>%
      filter(Study==paste0(s), Classification !="Unclassified", !is.na(Genus), Genus !="g__uncultured") %>%  
      select(starts_with("Salinity_")) %>%
      as.matrix()
    # Filter to include only abundant things, because otherwise there are too many
    # And also, low abundance things have too many zeros to work sometimes
    maxRow <- apply(mat16_salinity, 1, FUN=function(x) max(x, na.rm=TRUE))
    preva <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE)/length(x))
    nobs <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE))
    
    # Tried not using preva here
    toKeepRows <- which(maxRow>0.001 & nobs>2)
    # toKeepRowsLOW <- which(maxRow<=0.001 & preva<=0.1)
    newNames <- dat16 %>%   
      filter(Study==paste0(s), Classification !="Unclassified", !is.na(Genus), Genus !="g__uncultured") %>%  
      select(TaxaID, Classification, Study) %>%
      unite(., col="newNames") %>% pull(newNames)
    rownames(mat16_salinity) <- newNames
    # make it so rows sum to 1 (each ASV has 100% total)
    mat16_salinity <- mat16_salinity/rowSums(mat16_salinity, na.rm=TRUE)
    mat16_salinity_filt <- mat16_salinity[toKeepRows,which(!is.na(colSums(mat16_salinity)))]
    # mat16_salinity_filtLOW <- mat16_salinity[toKeepRowsLOW,which(!is.na(colSums(mat16_salinity)))]
    # Distance matrix
    dist_mat16 <- vegdist(mat16_salinity_filt, method="bray", na.rm=TRUE)
    # dist_mat16LOW <- vegdist(mat16_salinity_filtLOW, method="bray", na.rm=TRUE)
    nmds_dist_mat16 <- isoMDS(dist_mat16)
    # nmds_dist_mat16_LOW <- isoMDS(dist_mat16LOW)
    nmds_dat <- nmds_dist_mat16$points %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study")) 
    # nmds_dat_LOW <- nmds_dist_mat16_LOW$points %>%
    # as.data.frame() %>% rownames_to_column(var="temp") %>%
    # separate(temp, into=c("TaxaID","Classification","Study")) 
    pca_dist_mat16 <- prcomp(dist_mat16)
    # pca_dist_mat16_LOW <- prcomp(dist_mat16_LOW)
    # pca_dat_LOW <-pca_dist_mat16_LOW$x %>%
    # as.data.frame() %>% rownames_to_column(var="temp") %>%
    # separate(temp, into=c("TaxaID","Classification","Study"))
    pca_dat <-pca_dist_mat16$x %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study"))
    asv_clustering <- full_join(asv_clustering, full_join(nmds_dat, pca_dat))
    # asv_clustering_LOW <- full_join(asv_clustering_LOW, full_join(nmds_dat_LOW, pca_dat_LOW))
  }
  asv_clustering <- asv_clustering %>%filter(!is.na(TaxaID))
  
  save(asv_clustering, file="07e_MLclustering/asv_clustering.RData")
}

asv_clustering %>%
  ggplot() + 
  geom_point(aes(x=V1, y=V2, col=Classification)) +
  scale_color_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  facet_wrap(.~Study,scales="free")
asv_clustering %>%
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, col=Classification)) +
  scale_color_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  facet_wrap(.~Study,scales="free")


asv_clustering %>%
  left_join(dat16 %>% select(TaxaID,Study, Class)) %>%
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, col=Class), show.legend=FALSE) +
  # scale_color_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  facet_wrap(.~Study,scales="free")

###
##### Use kmeans clustering to group into similar clusters of appearance #######
if ( file.exists()) {
  load("07e_MLclustering/kmeans_clustering.RData")
  
} else {
  kmeans_clustering <- data.frame(TaxaID=NA)
  for ( s in allStudies) {
    mat16_salinity <- dat16 %>%
      filter(Study==paste0(s), Classification !="Unclassified", !is.na(Genus), Genus !="g__uncultured") %>%  
      select(starts_with("Salinity_")) %>%
      as.matrix()
    # Filter to include only abundant things, because otherwise there are too many
    # And also, low abundance things have too many zeros to work sometimes
    maxRow <- apply(mat16_salinity, 1, FUN=function(x) max(x, na.rm=TRUE))
    preva <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE)/length(x))
    nobs <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE))
    
    # Tried not using preva here
    toKeepRows <- which(maxRow>0.001 & nobs>2)
    # toKeepRowsLOW <- which(maxRow<=0.001 & preva<=0.1)
    newNames <- dat16 %>%   
      filter(Study==paste0(s), Classification !="Unclassified", !is.na(Genus), Genus !="g__uncultured") %>%  
      select(TaxaID, Classification, Study) %>%
      unite(., col="newNames") %>% pull(newNames)
    rownames(mat16_salinity) <- newNames
    # make it so rows sum to 1 (each ASV has 100% total)
    mat16_salinity <- mat16_salinity/rowSums(mat16_salinity, na.rm=TRUE)
    # mat16_salinity[mat16_salinity>0] <-1
    # get rid of column NAs
    # Filter to include only abundant things, because otherwise there are too many
    # And also, low abundance things have too many zeros to work sometimes
    maxRow <- apply(mat16_salinity, 1, FUN=function(x) max(x, na.rm=TRUE))
    preva <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE)/length(x))
    nobs <- apply(mat16_salinity, 1, FUN=function(x) sum(x>0, na.rm = TRUE))
    
    # Tried not using preva here
    toKeepRows <- which(maxRow>0.001 & nobs>2)
    toKeepCol <- !is.na(colSums(mat16_salinity))
    mat16_salinity_filt <- mat16_salinity[toKeepRows, toKeepCol]
    
    tempkmeans2<- kmeans(mat16_salinity_filt, centers =2)
    tempkmeans3<- kmeans(mat16_salinity_filt, centers =3)
    tempkmeans4<- kmeans(mat16_salinity_filt, centers =4)
    tempkmeans5<- kmeans(mat16_salinity_filt, centers =5)
    
    ## Find average salinity 
    maxLoading2 <- apply( tempkmeans2$centers , MARGIN=1, function(x) which.max(x))
    sal2 <- data.frame(cluster=names(maxLoading2), Salinity=colnames(tempkmeans2$centers)[maxLoading2]) %>%
      separate(Salinity, into=c("msc","Sal")) %>% mutate(Sal = as.numeric(Sal), cluster=as.numeric(cluster)) %>%
      select(cluster, Sal)
    maxLoading3 <- apply( tempkmeans3$centers , MARGIN=1, function(x) which.max(x))
    sal3 <- data.frame(cluster=names(maxLoading3), Salinity=colnames(tempkmeans3$centers)[maxLoading3]) %>%
      separate(Salinity, into=c("msc","Sal")) %>% mutate(Sal = as.numeric(Sal), cluster=as.numeric(cluster)) %>%
      select(cluster, Sal)
    maxLoading4 <- apply( tempkmeans4$centers , MARGIN=1, function(x) which.max(x))
    sal4 <- data.frame(cluster=names(maxLoading4), Salinity=colnames(tempkmeans4$centers)[maxLoading4]) %>%
      separate(Salinity, into=c("msc","Sal")) %>% mutate(Sal = as.numeric(Sal), cluster=as.numeric(cluster)) %>%
      select(cluster, Sal)
    maxLoading5 <- apply( tempkmeans5$centers , MARGIN=1, function(x) which.max(x))
    sal5 <- data.frame(cluster=names(maxLoading5), Salinity=colnames(tempkmeans5$centers)[maxLoading5]) %>%
      separate(Salinity, into=c("msc","Sal")) %>% mutate(Sal = as.numeric(Sal), cluster=as.numeric(cluster)) %>%
      select(cluster, Sal)

    kmeans2 <- tempkmeans2$cluster %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study"), sep="_") %>%
      rename(cluster=".") %>%
      mutate(withinSS = mean(tempkmeans2$withinss), betweenSS = tempkmeans2$betweenss, k=2) %>%
      left_join(sal2)
    kmeans3 <- tempkmeans3$cluster %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study"), sep="_") %>%
      rename(cluster=".") %>%
      mutate(withinSS = mean(tempkmeans3$withinss), betweenSS = tempkmeans3$betweenss, k=3) %>%
      left_join(sal3)
    kmeans4 <- tempkmeans4$cluster %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study"), sep="_") %>%
      rename(cluster=".") %>%
      mutate(withinSS = mean(tempkmeans4$withinss), betweenSS = tempkmeans4$betweenss, k=4) %>%
      left_join(sal4)
    kmeans5 <- tempkmeans5$cluster %>%
      as.data.frame() %>% rownames_to_column(var="temp") %>%
      separate(temp, into=c("TaxaID","Classification","Study"), sep="_") %>%
      rename(cluster=".") %>%
      mutate(withinSS = mean(tempkmeans5$withinss), betweenSS = tempkmeans5$betweenss, k=5) %>%
      left_join(sal5)
    allkmeanstemp <- rbind(kmeans2, kmeans3, kmeans4, kmeans5)
    
    kmeans_clustering <- full_join(kmeans_clustering, allkmeanstemp )
  }
  kmeans_clustering <- kmeans_clustering %>% filter(!is.na(TaxaID))
  
  save(kmeans_clustering, file="07e_MLclustering/kmeans_clustering.RData")
}

allClustering <- asv_clustering %>% full_join(kmeans_clustering) %>%
  group_by(Study, k) %>%
  mutate(cluster_ranked = rank(Sal)) %>%
  ungroup() 

### Need to find categories for each cluster
allClustering %>%
  filter(k==5) %>%
  ggplot() + 
  geom_point(aes(x=PC1, y=PC2, col=factor(cluster), fill=Classification), pch=21, stroke=2) +
  # geom_point(aes(x=V1, y=V2, col=(cluster_named), fill=Classification), pch=21, stroke=3) +
  # scale_color_manual(values=c(`1`="blue", `2`="purple", `3`="red")) +
  # scale_color_gradient2(low="blue",mid="purple",high="red") +
  scale_fill_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  # scale_color_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  facet_wrap(.~Study,scales="free")



## Legend




kmeans_clustering %>%
  select(Study, withinSS, betweenSS, k) %>% distinct() %>%
  mutate(rSS = withinSS/betweenSS) %>%
  ggplot() + 
  geom_line(aes(x=k, y=rSS, col=Study))

#### Combine kmeans and pca #########
allClustering <- asv_clustering %>% full_join(kmeans_clustering) %>%
  group_by(Study, k) %>%
  mutate(cluster_ranked = rank(Sal)) %>%
  ungroup() 

allClustering %>%
  filter(k==3) %>%
  ggplot() + 
  # geom_point(aes(x=V1, y=V2, col=factor(cluster), fill=Classification), pch=21, stroke=3) +
  geom_point(aes(x=V1, y=V2, col=(cluster_named), fill=Classification), pch=21, stroke=3) +
  # scale_color_manual(values=c(`1`="blue", `2`="purple", `3`="red")) +
  # scale_color_gradient2(low="blue",mid="purple",high="red") +
  scale_fill_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  scale_color_manual(values=c(Fresh="blue", Brackish="purple", Marine="red")) +
  facet_wrap(.~Study,scales="free")



#### First, let's find out which ASVs are consistent in classification #####

### Venn Diagram to see what overlap is ###
allFresh <- dat16 %>% filter(Classification=="Fresh") %>% pull(TaxaID) %>% unique()
allBrackish <- dat16 %>% filter(Classification=="Brackish") %>% pull(TaxaID) %>% unique()
allMarine <- dat16 %>% filter(Classification=="Marine") %>% pull(TaxaID) %>% unique()
listTypes <- list(Fresh=allFresh, Brackish=allBrackish, Marine= allMarine)

length(allFresh)
length(intersect(allFresh, allBrackish))
length(allBrackish)
length(intersect(allBrackish, allMarine))
length(allMarine)
length(intersect(allFresh, allMarine))

