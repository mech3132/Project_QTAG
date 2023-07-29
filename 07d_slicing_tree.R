#!bin/bash
library(ape)
library(phytools)
library(castor)
library(tidyverse)

##### SLICING TREE ########

# tre <- read.tree("03_phylogenetic_tree/sepp/phylogenetic_tree_16S_2022-08-14/exported_16S_tree/tree.nwk")
# tax2 <- read.delim("02_adjust_files_for_processing/taxa_adj.txt")

tax <- read.delim("06_processing_QTAG/all_taxonomy.txt")
tre16 <- read.tree("06_processing_QTAG/tree16_filt.tre")
dat16 <- read.delim("06_processing_QTAG/QTAG_output_with_taxa.txt") %>%
  filter(TaxaID %in% tre16$tip.label) 
tax16 <- tax %>% filter(TaxaID %in% tre16$tip.label)

tre16.filt <- keep.tip(tre16, tip=dat16$TaxaID)
NH16 <- node.depth.edgelength(tre16.filt) # Get depths of all nodes
treeH16 = max(NH16) # Get max height of tree

### Can't just correlate distances overall because clade-specific patterns make this messy.

allStudies16 <- dat16 %>% select(Study) %>% pull() %>% unique()
allLevels <- c("Phylum","Class","Order","Family","Genus","Species")
slices = seq(0.1,treeH16, length.out=20)

# iterate through studies
slice_collapsed_dat_filt <- data.frame()
nn_mindist <- data.frame()
n5n_mindist <- data.frame()
# nnratio_mindist <- data.frame()
tax_collapsed_data_filt <- data.frame()
for ( s in allStudies16) {
  print(paste0("doing ",s))
  
  dat_filt <- dat16 %>% filter(Study ==s)
  tre.filt_temp <- keep.tip(tre16.filt, tip=dat_filt$TaxaID)
  ######## Try slicing tree ##########
  all_slice_collapsed <- data.frame()
  for (slice in slices) {
    collapsed <- collapse_tree_at_resolution(tre.filt_temp, resolution = slice
                                             , rename_collapsed_nodes = TRUE
                                             , shorten=FALSE
                                             , criterion = 'max_tip_pair_dist')
    
    collapsedtips <- get_tips_for_mrcas(tre.filt_temp, mrca_nodes = collapsed$collapsed_nodes)
    for ( t in collapsedtips) {
      node_collapsed <- get_mrca_of_set(tre.filt_temp, t) 
      tiptemp <- getDescendants(tre.filt_temp,get_mrca_of_set(tre.filt_temp, t),node_collapsed)
      all_slice_collapsed <- rbind(all_slice_collapsed
                                   , cbind(slice=slice, mrca = node_collapsed, tips = tiptemp, TaxaID=tre.filt_temp$tip.label[tiptemp]))
    }
  }
  all_slice_collapsed_dat_filt <- all_slice_collapsed %>%
    mutate(Study=s) %>%
    filter(!is.na(TaxaID)) %>%
    left_join(dat_filt) %>%
    select(Study, slice, mrca, TaxaID, Classification, aveAbund, SalMax, lwrObs, lwrTol, uprObs, uprTol, Domain_NA, Phylum_NA, Class_NA, Order_NA, Family_NA, Genus_NA, Species_NA, ScientificName)
  slice_collapsed_dat_filt <- rbind(slice_collapsed_dat_filt, all_slice_collapsed_dat_filt)
  
  ######### nearest neighbour approach ###########
  disttips <- cophenetic.phylo(tre.filt_temp)
  disttips[disttips==0] <- NA
  
  # Filter frames to be fresh, marine, or brackish only
  allFresh <- dat_filt%>% filter(Classification=="Fresh") %>% pull(TaxaID)
  allMarine <- dat_filt%>% filter(Classification=="Marine") %>% pull(TaxaID)
  allBrackish <- dat_filt%>% filter(Classification=="Brackish") %>% pull(TaxaID)
  distFresh <- disttips[,allFresh]
  distMarine <- disttips[,allMarine]
  distBrackish <- disttips[,allBrackish]
  
  ## Nearest neightbour ##
  allMinDist <- data.frame(Study=s
                           , TaxaID = rownames(disttips)
                           , minDistFresh=as.numeric(apply(distFresh,1, function(x) min(x,na.rm=TRUE)))
                           , minDistMarine=as.numeric(apply(distMarine,1, function(x) min(x,na.rm=TRUE)))
                           , minDistBrackish=as.numeric(apply(distBrackish,1, function(x) min(x,na.rm=TRUE)))
  ) %>%
    rowwise() %>%
    mutate(nn_class = c("Fresh","Marine","Brackish")[which.min(c(minDistFresh, minDistMarine, minDistBrackish))]) %>%
    mutate(FM_ratio = (minDistFresh-minDistMarine)/mean(c(minDistFresh, minDistMarine))) %>%
    ungroup() %>% mutate(nn_class_FM = c("Fresh","Marine")[which.min(c(minDistFresh, minDistMarine))]) %>%
    left_join(dat_filt %>% select(TaxaID, Classification, aveAbund, SalMax, lwrTol, uprTol))

  nn_mindist <- rbind(nn_mindist, allMinDist)
  
  ######### nearest 5 neighbours approach ###########

  allMin5Dist <- data.frame(Study=s
                           , TaxaID = rownames(disttips)
                           , minDistFresh=as.numeric(apply(distFresh,1, function(x) median(sort(x)[1:5])))
                           , minDistMarine=as.numeric(apply(distMarine,1, function(x) median(sort(x)[1:5])))
                           , minDistBrackish=as.numeric(apply(distBrackish,1, function(x) median(sort(x)[1:5])))
                           , sdDistFresh=as.numeric(apply(distFresh,1, function(x) sd(sort(x)[1:5])))
                           , sdDistMarine=as.numeric(apply(distMarine,1, function(x) sd(sort(x)[1:5])))
                           , sdDistBrackish=as.numeric(apply(distBrackish,1, function(x) sd(sort(x)[1:5])))
  ) %>%
    rowwise() %>%
    mutate(nn_class = c("Fresh","Marine","Brackish")[which.min(c(minDistFresh, minDistMarine, minDistBrackish))]) %>%
    mutate(FM_ratio = (minDistFresh-minDistMarine)/mean(c(minDistFresh, minDistMarine))) %>%
    mutate(nn_class_FM = c("Fresh","Marine")[which.min(c(minDistFresh, minDistMarine))]) %>%
    ungroup() %>%
    left_join(dat_filt %>% select(TaxaID, Classification, aveAbund, SalMax, lwrTol, uprTol))
  
  n5n_mindist <- rbind(n5n_mindist, allMin5Dist)
  # 
  # ######### nearest neighbour ratio ###########
  # 
  # allMinratioDist <- data.frame(Study=s
  #                           , TaxaID = rownames(disttips)
  #                           , minDistFresh=as.numeric(apply(distFresh,1, function(x) min(x,na.rm = TRUE)))
  #                           , minDistMarine=as.numeric(apply(distMarine,1, function(x) min(x,na.rm = TRUE)))
  #                           , minDistBrackish=as.numeric(apply(distBrackish,1, function(x) min(x,na.rm = TRUE)))
  # ) %>%
  #   rowwise() %>%
  #   mutate(nn_class = c("Fresh","Marine","Brackish")[which.min(c(minDistFresh, minDistMarine, minDistBrackish))]) %>%
  #   mutate(FM_ratio = (minDistFresh-minDistMarine)/mean(c(minDistFresh, minDistMarine))) %>%
  #   ungroup() %>% mutate(nn_class_FM = c("Fresh","Marine")[which.min(c(minDistFresh, minDistMarine))]) %>%
  #   left_join(dat_filt %>% select(TaxaID, Classification, aveAbund, SalMax, lwrTol, uprTol))
  # 
  # nnratio_mindist <- rbind(nnratio_mindist, allMinratioDist)
  # 
  ######### taxonomic approach ###########
  allleveldat <- data.frame()
  for ( l in allLevels) {
    level_col <- dat_filt %>% group_by_at(l) %>%
      summarise(nFresh=sum(class=="low"), nMarine=sum(class=="high"), nBrackish=sum(class=="intermediate"), nUnassigned=sum(class=="unclassified"), totalCount=n()) %>%
      rename_at(vars(all_of(l)), ~"Group")%>%
      mutate(Level=l, Study=s)
    allleveldat <- rbind(allleveldat, level_col)
  }
  tax_collapsed_data_filt <- rbind(tax_collapsed_data_filt, allleveldat)
}

####### Look at results ##########

slice_collapsed_dat_filt 
nn_mindist
n5n_mindist
# nnratio_mindist
tax_collapsed_data_filt

############ Visualization #######

#### Slicing tree ####
slice_summary <- slice_collapsed_dat_filt %>%
  group_by(slice, mrca, Study) %>%
  summarise(n=n(), nFresh=sum(Classification%in% "Fresh"), nMarine = sum(Classification=="Marine"), nBrackish=sum(Classification=="Brackish")
            , se_average = sd(aveAbund)/sqrt(n())
            , whichFM = ifelse(nFresh>nMarine, "Fresh", 
                               ifelse(nMarine>nFresh,"Marine", 
                                      ifelse(nBrackish>nMarine & nBrackish>nFresh, "Brackish","Equal")))
            , se_salmax = sd(SalMax)/sqrt(n())) %>%
  ungroup() %>%rowwise() %>% mutate(propMax = max(c(nFresh, nMarine, nBrackish))/n) %>% ungroup() %>%
  mutate(slice = as.numeric(slice)) 


ggplot(slice_summary) +
  geom_jitter(aes(x=(slice), y=propMax, col=whichFM), width=0.1)  +
  geom_smooth(aes(x=(slice), y=propMax, col=whichFM), method="loess") +
  ylim(0,1) +
  facet_grid(.~Study)+
  scale_color_manual(values=c(Fresh="blue", Marine="red", Brackish="purple", Equal="grey")) +
  xlab("Slicing depth into phylogenetic tree") +
  ylab("Maximum proportion of clade that is one classificaion") +
  labs(col="Dominant classification\nin each clade")

# ### Proportion ###
# ggplot(slice_summary) +
#   geom_jitter(aes(x=(slice), y=propMax, col=whichFM), width=0.1)  +
#   geom_smooth(aes(x=(slice), y=propMax, col=whichFM)) 

# ggplot(slice_summary) +
#   geom_bar(aes(x=factor(round(slice)), y=propMax, fill=whichFM, group=mrca), stat="identity", position="dodge")  
# slice_summary$se_average
ggplot(slice_summary) +
  geom_jitter(aes(x=(slice), y=se_average, col=whichFM)) +
  geom_smooth(aes(x=(slice), y=se_average, col=whichFM)) +
  ylim(0,15) +
  facet_grid(.~Study)+
  scale_color_manual(values=c(Fresh="blue", Marine="red", Brackish="purple", Equal="grey"))


ggplot(slice_summary) +
  geom_jitter(aes(x=(slice), y=se_salmax, col=whichFM)) +
  geom_smooth(aes(x=(slice), y=se_salmax, col=whichFM)) +
  ylim(0,20) +
  facet_grid(.~Study)+
  scale_color_manual(values=c(Fresh="blue", Marine="red", Brackish="purple", Equal="grey"))


#### Nearest Neighbour ####
# nn_mindist %>% ungroup() %>%
#   # filter(Study=="wang2021") %>%
#   # filter(Classification=="Brackish") %>%
#   # mutate(MF = aveAbund>10) %>%
#   ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=Classification), alpha=0.25) +
#   geom_abline(aes(slope=1, intercept=0), col="red") +
#   scale_color_manual(values=c(Fresh="blue", Marine="red", Brackish="purple", Equal="grey"))



gg_freshmarinedist <- nn_mindist %>%
  filter(Classification != "Unclassified") %>%
  ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=SalMax)) +
  geom_abline(aes(slope=1, intercept=0), col="grey") +
  facet_grid(Classification~Study) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)+
  ylab("Distance to nearest\nmarine neighbour") + xlab("Distance to nearest\nfresh neighbour") +
  labs(col="Salinity where\nabundance peaks\nfor ASV")
gg_freshmarinedist

gg_freshmarinedist2 <- nn_mindist %>%
  filter(Classification != "Unclassified") %>%
  ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=aveAbund)) +
  geom_abline(aes(slope=1, intercept=0), col="grey") +
  facet_grid(Classification~Study) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)+
  ylab("Distance to nearest\nmarine neighbour") + xlab("Distance to nearest\nfresh neighbour") +
  labs(col="Salinity midpoint of\nobserved range")
gg_freshmarinedist2


### Nearest 5 neighbours, filtered ###
n5n_mindist %>% 
  # filter(Classification != "Unclassified", sdDistFresh<0.5, sdDistMarine<0.5) %>%
  filter(Classification != "Unclassified") %>%
  ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=aveAbund)) +
  geom_abline(aes(slope=1, intercept=0), col="grey") +
  facet_grid(Classification~Study, scales="free") +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)+
  ylab("Distance to nearest\nmarine neighbour") + xlab("Distance to nearest\nfresh neighbour") +
  labs(col="Salinity midpoint of\nobserved range")

n5n_mindist %>% 
  # filter(Classification != "Unclassified", sdDistFresh<0.5, sdDistMarine<0.5) %>%
  filter(Classification != "Unclassified") %>%
  ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=SalMax)) +
  geom_abline(aes(slope=1, intercept=0), col="grey") +
  facet_grid(Classification~Study, scales="free") +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)+
  ylab("Distance to nearest\nmarine neighbour") + xlab("Distance to nearest\nfresh neighbour") +
  labs(col="Salinity where\nabundance peaks\nfor ASV")

# What about brackish-brackish relationship?
n5n_mindist %>%
  mutate(tolRange = uprTol-lwrTol) %>%
  select(Study, TaxaID, minDistFresh, minDistMarine, minDistBrackish, Classification, SalMax) %>%
  pivot_longer(c(minDistFresh, minDistMarine, minDistBrackish), names_to="Comparison", values_to="PD") %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
  mutate(Comparison = factor(Comparison, levels=c("minDistFresh","minDistBrackish","minDistMarine"))) %>%
  filter(Classification!="Unclassified") %>%
  ggplot() +
  geom_boxplot(aes(x=Comparison, y=PD))+
  facet_grid(.~Classification)+
  theme(axis.text.x = element_text(angle=90))

nn_mindist %>%
  mutate(tolRange = uprTol-lwrTol) %>%
  select(Study, TaxaID, minDistFresh, minDistMarine, minDistBrackish, Classification, SalMax) %>%
  pivot_longer(c(minDistFresh, minDistMarine, minDistBrackish), names_to="Comparison", values_to="PD") %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
  mutate(Comparison = factor(Comparison, levels=c("minDistFresh","minDistBrackish","minDistMarine"))) %>%
  filter(Classification!="Unclassified") %>%
  ggplot() +
  geom_boxplot(aes(x=Comparison, y=(PD)))+
  facet_grid(Classification~Study)+
  theme(axis.text.x = element_text(angle=90))

#### Filter by abundance #####
maxObs <- apply(dat16 %>% select(starts_with("Salinity_")) ,1,FUN=function(x) max(as.numeric(x), na.rm=TRUE))
dat16_highabund <- dat16[maxObs>0.01,]

dat16_highabund %>% select(TaxaID, Study) %>% left_join(nn_mindist) %>%
  filter(Classification != "Unclassified") %>%
  ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=aveAbund)) +
  geom_abline(aes(slope=1, intercept=0), col="grey") +
  facet_grid(Classification~Study) +
  scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)+
  ylab("Distance to nearest\nmarine neighbour") + xlab("Distance to nearest\nfresh neighbour") +
  labs(col="Salinity midpoint of\nobserved range")

### Histogram peak salinity
dat16 %>%
  ggplot() + 
  geom_histogram(aes(x=SalMax, y=after_stat(density), group=Classification, fill=Classification)
                 , alpha=0.5, position="identity") +
  geom_density(aes(x=SalMax, y=after_stat(density), group=Classification, fill=Classification)
                 , alpha=0.5, position="identity") +
  facet_grid(Classification~Study, scales="free")

# n5n_mindist %>%
#   ggplot() + geom_point(aes(x=minDistFresh, y=minDistMarine, col=SalMax)) +
#   geom_abline(aes(slope=1, intercept=0), col="red") +
#   facet_grid(Classification~Study) +
#   scale_color_gradient2(low="blue", mid="white", high="red", midpoint=15)

#### Nearest Neighbour Ratio ####
nn_mindist %>%
  filter(Classification !="Unclassified") %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine"))) %>%
  ggplot(aes(x=Classification, y=(FM_ratio))) +
  geom_jitter(aes(col=SalMax)) +
  facet_grid(.~Study) 


#### Sal max vs distance to fresh and marine
gg_distfresh_sal <- nn_mindist %>%
  filter(Classification%in%c("Fresh","Marine")) %>%
  # mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine"))) %>%
  ggplot(aes(x=SalMax, y=minDistFresh, col=Classification)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(.~Study) +
  ylab("Distance to nearest\nfreshwater neighbour") +
  xlab("Salinity at which population\nis highest")
gg_distfresh_sal


gg_distmarine_sal <- nn_mindist %>%
  filter(Classification%in%c("Fresh","Marine")) %>%
  # mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine"))) %>%
  ggplot(aes(x=SalMax, y=minDistMarine, col=Classification)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(.~Study) +
  ylab("Distance to nearest\nmarine neighbour")+
  xlab("Salinity at which population\nis highest")
gg_distmarine_sal

gridExtra::grid.arrange(gg_distfresh_sal, gg_distmarine_sal, nrow=2)
gg_distfresh_sal_brackish <- nn_mindist %>%
  filter(Classification =="Brackish") %>%
  ggplot(aes(x=SalMax, y=minDistFresh)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Study~., scales= "free_y") +
  ylab("Distance to nearest\nfreshwater neighbour") +
  xlab("Salinity at which population\nis highest")
gg_distfresh_sal_brackish

gg_distmarine_sal_brackish <- nn_mindist %>%
  filter(Classification =="Brackish") %>%
  ggplot(aes(x=SalMax, y=minDistMarine)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Study~., scales= "free_y") +
  ylab("Distance to nearest\nmarine neighbour")+
  xlab("Salinity at which population\nis highest")
gg_distmarine_sal_brackish


gg_distbrack_sal_brackish <- nn_mindist %>%
  filter(Classification =="Brackish") %>%
  ggplot(aes(x=SalMax, y=minDistBrackish)) + 
  geom_point() +
  geom_smooth(method="lm") +
  facet_grid(Study~., scales= "free_y") +
  ylab("Distance to nearest\nbrackish neighbour")+
  xlab("Salinity at which population\nis highest")
gg_distbrack_sal_brackish

gridExtra::grid.arrange(gg_distfresh_sal_brackish, gg_distbrack_sal_brackish, gg_distmarine_sal_brackish, ncol=3)

#### Try looking at one slice at a time
# all_slice_collapsed_dat_filt %>% filter(slice==slices[3]) %>% View()

#### Taxonomy collapse ####
tax_collapsed_data_filt %>%
  mutate(propFresh=nFresh/(nFresh+nMarine), propMarine=nMarine/(nFresh+nMarine)) %>%
  filter(Level=="Class") %>%
  select(Group, Study, propFresh, propMarine) %>%
  pivot_longer(c(propFresh, propMarine), names_to="Classification", values_to="proportion") %>%
  ggplot() +
  geom_bar(aes(x=Group, y=proportion, fill=Classification), stat="identity")+
  facet_grid(Study~.) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  scale_fill_manual(values=c(propBrackish="purple", propFresh="blue", propMarine = "red"))

### Calculate dominance over taxa
tax_collapsed_adj <- tax_collapsed_data_filt %>%
  mutate(Level=factor(Level, levels=allLevels)) %>%
  mutate(propFresh=nFresh/(nFresh+nMarine), propMarine=nMarine/(nFresh+nMarine)) %>%
  mutate(prFresh = nFresh/totalCount, prMarine = nMarine/totalCount, prBrackish = nBrackish/totalCount) 

# filter(totalCount>1) %>%
  # select(Level, Study, pFresh, pMarine, pBrackish) %>%
  # pivot_longer(c(pFresh, pMarine, pBrackish), names_to="Classification", values_to="Proportion") %>%
  # ggplot(aes(x=Level, y=maxProp, col=maxClass)) +
  # geom_violin() 

  
allTaxaCollapse <- tax16
for ( l in allLevels ) {
  # l = "Phylum"
  tax_temp <- tax_collapsed_adj %>% filter(Level==l) %>%
    select(Group, Study, propFresh, propMarine, prFresh, prMarine, prBrackish) %>%
    rename_at(vars(Group), ~paste0(l)) %>%
    rename_at(vars(starts_with("pr")), ~paste0(l,"_",.))
  allTaxaCollapse <- left_join(allTaxaCollapse, tax_temp)
  }
  
allTaxaCollapse %>%
  select(Study, TaxaID, one_of(allLevels), ends_with("propMarine")) %>%
  pivot_longer(c(paste0(allLevels, "_propMarine")), names_to="Group", values_to="Proportion") %>%
  separate(Group, into=c("Level","ID"), sep="_") %>%
  mutate(Level=factor(Level, levels=allLevels)) %>%
  ggplot(aes(x=Level, y=Proportion, group=TaxaID)) + 
  geom_point() +
  geom_line() +
  facet_grid(Study~.)

#################
allStudies16
tempdat <- dat16 %>% 
  filter(Study=="chen16") %>%
  filter(Class%in% c("c__Gammaproteobacteria", "c__Alphaproteobacteria"))
tre.filt2 <- keep.tip(tre16, tip=tempdat$TaxaID)

classcol <- tempdat[match(tre.filt2$tip.label, tempdat$TaxaID),"Class"]
classcol[classcol=="c__Alphaproteobacteria"] <- "darkred"
classcol[classcol=="c__Gammaproteobacteria"] <- "blue"

plot(tre.filt2, tip.color = classcol)

######### Simple heatmap plots for brackish
onlyBrackish_dat <- dat16 %>% filter(Classification=="Brackish", Study=="herleman")
tre16_brack <- keep.tip(tre16, onlyBrackish_dat$TaxaID)
dm_brack_phylo <- cophenetic.phylo(tre16_brack)
heatmap(dm_brack_phylo)
