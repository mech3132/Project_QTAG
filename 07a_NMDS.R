#!bin/bash
# library(phytools)
# library(GUniFrac)
# library(ape)
library(vegan)
library(MASS)
library(ggbiplot)
library(tidyverse)
dir.create("07a_NMDS")

######## Load ##########
mf <- read.delim("04_adjust_files_for_QTAG/all_meta.txt")
otucount <- read.delim("04_adjust_files_for_QTAG/otu.txt")
taxa <- read.delim("04_adjust_files_for_QTAG/taxonomy.txt")
dmWU16 <- read.delim("06_diversity/dm_WU_16S/distance-matrix.tsv", row.names = 1) %>% as.matrix()
dmUWU16 <- read.delim("06_diversity/dm_UWU_16S/distance-matrix.tsv", row.names = 1) %>% as.matrix()

### Filter dm by samples in mf
# mf$SampleID %in% colnames(dmWU16)
toKeepSamps <- intersect(mf$SampleID, colnames(dmWU16))

dmWU16_filt <- dmWU16[toKeepSamps,toKeepSamps]
dmUWU16_filt <- dmUWU16[toKeepSamps,toKeepSamps]
# Sanity check; the below should be redundant anyway
mf <- mf %>% filter(SampleID %in% toKeepSamps)
mf16 <- mf %>% filter(target_region=="16S")

###### NMDS 16S #########
nmdsWU16 <- isoMDS(dmWU16_filt, k=2)
nmdsUWU16 <- isoMDS(dmUWU16_filt, k=2)

gg_WU16 <- nmdsWU16$points %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
  left_join(mf) %>%
  ggplot() +
  geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity), shape=Study), cex=0.5) +
  ylab("NMDS2") + xlab("NMDS1") +
  scale_color_gradient(low="blue",high="darkred") +
  labs(title=paste0("Stress: ",round(nmdsUWU16$stress,2))) +
  labs(col="Salinity (ppt)")

gg_WU16
ggsave(filename = "07a_NMDS/nmds_wu_16.png",
       gg_WU16, height=4, width=6)

# dev.off()
ggUWU16 <- nmdsUWU16$points %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
  left_join(mf) %>%
  ggplot() +
  geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity), shape=Study)) +
  ylab("NMDS2") + xlab("NMDS1") +
  scale_color_gradient(low="blue",high="darkred") +
  labs(title=paste0("Stress: ",round(nmdsUWU16$stress,2))) +
  labs(col="Salinity (ppt)")
ggUWU16
ggsave(filename = "07a_NMDS/nmds_uwu_16.png",
       gg_WU16, height=4, width=6)


#### PCA ####
pcaWU16 <- prcomp(dmWU16_filt)
pcaUWU16 <- prcomp(dmUWU16_filt)

### WORKING HERE ##

var_expl_WU16 <- (pcaWU16$sdev^2/sum(pcaWU16$sdev^2))
var_expl_UWU16 <- (pcaUWU16$sdev^2/sum(pcaUWU16$sdev^2))

gg_WU16_pca <- pcaWU16$x %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
  left_join(mf) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, col=as.numeric(Salinity), shape=Study), cex=0.5) +
  ylab(paste0("PC2 (",round(var_expl_WU16[2]*100,1),"% var. expl)")) + 
  xlab(paste0("PC1 (",round(var_expl_WU16[1]*100,1),"% var. expl)")) +
  scale_color_gradient(low="blue",high="darkred") +
  labs(col="Salinity (ppt)")

gg_WU16_pca
ggsave(filename = "07a_NMDS/pca_wu_16.png",
       gg_WU16_pca, height=4, width=6)

# dev.off()
ggUWU16_pca <- pcaUWU16$x %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
  left_join(mf) %>%
  ggplot() +
  geom_point(aes(x=PC1, y=PC2, col=as.numeric(Salinity), shape=Study)) +
  ylab(paste0("PC2 (",round(var_expl_UWU16[2]*100,1),"% var. expl)")) + 
  xlab(paste0("PC1 (",round(var_expl_UWU16[1]*100,1),"% var. expl)")) +  scale_color_gradient(low="blue",high="darkred") +
  labs(col="Salinity (ppt)")
ggUWU16_pca
ggsave(filename = "07a_NMDS/pca_uwu_16.png",
       ggUWU16_pca, height=4, width=6)


#### Each season/study separately ####
allStudies <- unique(mf$Study)
allStudySeason <- unique(mf$StudySeason)

for ( d in c("dmUWU16_filt", "dmWU16_filt")) {
  for (s in allStudies ) {
    tempmf <- mf16 %>% filter(Study==s) 
    tempDM <- get(d)[tempmf$SampleID, tempmf$SampleID]
    tempnmds <- isoMDS(tempDM)
    temppca <- prcomp(tempDM)
    
    temp_gg <- tempnmds$points %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
      left_join(tempmf) %>%
      ggplot() +
      geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity), pch=Season) )+
      ylab("NMDS2") + xlab("NMDS1") +
      scale_color_gradient(low="blue",high="darkred") +
      labs(title=paste0("Stress: ",round(tempnmds$stress,2))) +
      labs(col="Salinity (ppt)")
    ggsave(filename = paste0("07a_NMDS/nmds_",s,"_",d,".png"), temp_gg
           , height=4, width=6)
    
    var_expl_temp <- (temppca$sdev^2/sum(temppca$sdev^2))
    temp_gg_pca <- temppca$x %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
      left_join(tempmf) %>%
      ggplot() +
      geom_point(aes(x=PC1, y=PC2, col=as.numeric(Salinity), pch=Season) )+
      ylab(paste0("PC2 (",round(var_expl_UWU16[2]*100,1),"% var. expl)")) + 
      xlab(paste0("PC1 (",round(var_expl_UWU16[1]*100,1),"% var. expl)")) +  
      scale_color_gradient(low="blue",high="darkred") +
      labs(col="Salinity (ppt)")
    ggsave(filename = paste0("07a_NMDS/pca_",s,"_",d,".png"), temp_gg_pca
           , height=4, width=6)

    temp_gg_season <- tempnmds$points %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
      left_join(tempmf) %>%
      ggplot() +
      geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity))) +
      ylab("NMDS2") + xlab("NMDS1") +
      scale_color_gradient(low="blue",high="darkred") +
      labs(title=paste0("Stress: ",round(tempnmds$stress,2))) +
      labs(col="Salinity (ppt)") +
      facet_wrap(.~Season, nrow=1)
    nSeasons <- length(unique(tempmf$Season))
    ggsave(filename = paste0("07a_NMDS/nmds_",s,"_",d,"_season.png"), temp_gg_season
           , height=4, width=5*nSeasons)
    
    temp_gg_season_pca <- temppca$x %>%as.data.frame() %>% rownames_to_column(var="SampleID") %>%
      left_join(tempmf) %>%
      ggplot() +
      geom_point(aes(x=PC1, y=PC2, col=as.numeric(Salinity), pch=Season) )+
      ylab(paste0("PC2 (",round(var_expl_UWU16[2]*100,1),"% var. expl)")) + 
      xlab(paste0("PC1 (",round(var_expl_UWU16[1]*100,1),"% var. expl)")) +  
      scale_color_gradient(low="blue",high="darkred") +
      labs(col="Salinity (ppt)") +
      facet_wrap(.~Season, nrow=1)
    # nSeasons <- length(unique(tempmf$Season))
    ggsave(filename = paste0("07a_NMDS/pca_",s,"_",d,"_season.png"), temp_gg_season_pca
           , height=4, width=5*nSeasons)
    
  }
}



