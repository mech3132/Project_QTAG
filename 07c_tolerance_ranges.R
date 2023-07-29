#!bin/bash

####### tolerance ranges ####
library(tidyverse)
dir.create("07c_tolerance_ranges")

# load 
dat <- read.delim("06_processing_QTAG/QTAG_output_with_taxa.txt")
meta <- read.delim("04_adjust_files_for_QTAG/all_meta.txt") %>%
  select(StudySeason, Study, Season)

# dat %>% View()
gg_tolrange <- dat %>%
  mutate(TaxaID=make.unique(TaxaID)) %>%
  mutate(midSal = (lwrTol+uprTol)/2) %>%
  mutate(Classification=ifelse(class=="high","Marine"
                               , ifelse(class=="intermediate","Brackish",
                                                              ifelse(class=="low","Fresh","Unclassified")))) %>%
  mutate(Classification=factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
  arrange( Classification, midSal) %>%
  mutate(TaxaID=factor(TaxaID, levels=unique(TaxaID)))  %>%
  ggplot() +
  geom_segment(aes(x=lwrObs, xend=uprObs, y=TaxaID, yend=TaxaID,col=Classification), alpha=0.25) +
  geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID,col=Classification)) +
  facet_wrap(.~Study, drop=TRUE, scales="free_y")+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
  xlab("Salinity") 
gg_tolrange
ggsave(filename = "07c_tolerance_ranges/tolerance_range_plot_all.png"
       ,gg_tolrange, height=8, width=12)


# # dat %>% View()
# gg_tolrange_rare <- dat %>%
#   mutate(sumRare=rowSums(across(starts_with("SalinityRare_")), na.rm=TRUE)) %>%
#   filter(sumRare>0) %>%
#   mutate(TaxaID=make.unique(TaxaID)) %>%
#   mutate(midSal = (lwrTol+uprTol)/2) %>%
#   mutate(Classification=ifelse(class=="high","Marine"
#                                , ifelse(class=="intermediate","Brackish",
#                                         ifelse(class=="low","Fresh","Unclassified")))) %>%
#   mutate(Classification=factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
#   arrange( Classification, midSal) %>%
#   mutate(TaxaID=factor(TaxaID, levels=unique(TaxaID)))  %>%
#   ggplot() +
#   geom_segment(aes(x=lwrObs, xend=uprObs, y=TaxaID, yend=TaxaID,col=Classification), alpha=0.25) +
#   geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID,col=Classification)) +
#   facet_wrap(.~Study, drop=TRUE, scales="free_y")+
#   theme(axis.text.y = element_blank())+
#   scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
#   xlab("Salinity") 
# gg_tolrange_rare
# ggsave(filename = "07c_tolerance_ranges/tolerance_range_plot_rare.png"
#        ,gg_tolrange_rare, height=8, width=12)

gg_tolrange_facet_no18 <-  dat %>%
  filter(Study!="chen18") %>%
  mutate(TaxaID=make.unique(TaxaID)) %>%
  mutate(midSal = (lwrTol+uprTol)/2) %>%
  mutate(Classification=ifelse(class=="high","Marine"
                               , ifelse(class=="intermediate","Brackish",
                                        ifelse(class=="low","Fresh","Unclassified")))) %>%
  mutate(Classification=factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
  arrange( Classification, midSal) %>%
  mutate(TaxaID=factor(TaxaID, levels=unique(TaxaID)))  %>%
  ggplot() +
  geom_segment(aes(x=lwrObs, xend=uprObs, y=TaxaID, yend=TaxaID,col=Classification), alpha=0.25) +
  geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID,col=Classification)) +
  facet_grid(.~Study, drop=TRUE)+
  theme(axis.text.y = element_blank())+
  scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
  xlab("Salinity") 
gg_tolrange_facet_no18
ggsave(filename = "07c_tolerance_ranges/tolerance_range_facet.png"
       ,gg_tolrange_facet_no18, height=6, width=12)

######## Other plots

# 
# dat %>%
#   filter(Study=="chen16") %>%
#   mutate(abund=rowSums(across(starts_with("Salinity_")), na.rm=TRUE)) %>%
#   filter(abund>0) %>%
#   select(TaxaID,Classification, starts_with("Salinity_")) %>%
#   pivot_longer(-c(TaxaID,Classification), names_to="Salinity", values_to="RA")%>%
#   mutate(Salinity=as.numeric(gsub("Salinity_","",Salinity))) %>%
#   filter(!is.na(RA)) %>%
#   group_by(Salinity) %>%
#   mutate(maxRA=max(RA)) %>%
#   ungroup() %>%
#   mutate(relRA = RA/maxRA) %>%
#   ggplot() +
#   geom_smooth(aes(x=Salinity, y=relRA, group=TaxaID, col=Classification), se=FALSE)

## Calculate rank and tolerance range
# rank_and_maxabund <- dat %>%
#   select(TaxaID, Study, Classification,lwrTol, uprTol,SalMax, starts_with("SAMPLE_")) %>%
#   pivot_longer(c(-TaxaID, -Study, -Classification, -lwrTol, -uprTol,-SalMax),names_to="Sample", values_to="RA") %>%
#   drop_na() %>%
#   group_by(Sample) %>%
#   mutate(rank_each = rank(RA,ties.method = "first"), rank_rel = 1-rank(RA, ties.method="first")/n()) %>%
#   ungroup() %>%
#   # group_by(TaxaID, Study, Classification, lwrTol, uprTol, SalMax) %>%
#   group_by(TaxaID, Study, lwrTol, uprTol, SalMax) %>%
#   summarise(maxAbund = max(RA), sumAbund = sum(RA), minRank_each = min(rank_each), maxRank_rel = max(rank_rel), medRank_rel = median(rank_rel)) %>% ungroup() %>%
#   mutate(tolRangeSize = uprTol-lwrTol, aveSal = mean(c(uprTol, lwrTol)))
rank_and_maxabund <- dat %>%
  select(TaxaID, Study,starts_with("SAMPLE_")) %>%
  pivot_longer(c(-TaxaID, -Study),names_to="Sample", values_to="RA") %>%
  drop_na() %>%
  group_by(Sample) %>%
  # filter(RA>0) %>%
  mutate(rank_each = rank(-RA, ties.method = "average"), rank_rel = rank(RA, ties.method="average")/n()) %>%
  ungroup() %>%
  # group_by(TaxaID, Study, Classification, lwrTol, uprTol, SalMax) %>%
  group_by(TaxaID, Study) %>%
  summarise(maxAbund = max(RA), sumAbund = sum(RA), minRank_each = min(rank_each), maxRank_rel = max(rank_rel), medRank_rel = median(rank_rel)) %>% ungroup() 

rank_and_maxabund <- rank_and_maxabund %>%
  left_join(dat %>% select(TaxaID, Study, Classification, uprTol, lwrTol, uprObs, lwrObs, SalMax)) %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unassigned")))%>%
  mutate(tolRangeSize = uprTol-lwrTol, aveSal = mean(c(uprTol, lwrTol)), obsRangeSize = uprObs-lwrObs)

# Correlation between tolerance range size and best rank?
rank_and_maxabund %>%
  filter(Classification !="Unassigned") %>%
  # filter(maxAbund>0.01) %>%
  ggplot() + 
  geom_point(aes(x=(maxRank_rel), y=tolRangeSize, col=SalMax)) +
  # geom_smooth(aes(x=(maxRank_rel), y=tolRangeSize), method="lm")+
  facet_grid(Study~Classification, scales="free")+
  # xlab("Relative rank in community\n(1 = most abundant)") +
  scale_color_gradient2(low="blue", high="red", midpoint=15)


# Check ranks are correlated; sanity check
rank_and_maxabund %>%
  ggplot() +
  geom_point(aes(x=maxRank_rel, y=minRank_each))+
  facet_wrap(.~Study)
# Check distr of maxAbund
rank_and_maxabund %>%
  ggplot() +
  geom_histogram(aes(x=(maxAbund), y=after_stat(density))) +
  facet_wrap(.~Study)+
  scale_x_log10()

## Tolerance range size and salinity
dat %>%
  mutate(tolRangeSize = uprTol-lwrTol) %>%
  ggplot() +
  geom_point(aes(x=SalMax, y=tolRangeSize, col=Classification))+
  geom_smooth(aes(x=SalMax, y=tolRangeSize), method="lm")+
  facet_grid(Classification~Study)
# Correlation between tolerance range and maximum observed abundance?
rank_and_maxabund %>%
  ggplot() + 
  geom_point(aes(x=(maxAbund), y=tolRangeSize, col=Classification)) +
  geom_smooth(aes(x=(maxAbund), y=tolRangeSize), method="lm")+
  facet_grid(Classification~Study)
# Correlation between tolerance range and summed observed abundance?
rank_and_maxabund %>%
  ggplot() + 
  geom_point(aes(x=log(sumAbund), y=tolRangeSize, col=Classification)) +
  geom_smooth(aes(x=log(sumAbund), y=tolRangeSize), method="lm")+
  facet_grid(Classification~Study)

# Correlation between tolerance range size and best rank?
rank_and_maxabund %>%
  filter(Classification !="Unassigned") %>%
  ggplot() + 
  geom_point(aes(x=(maxRank_rel), y=tolRangeSize, col=SalMax)) +
  geom_smooth(aes(x=(maxRank_rel), y=tolRangeSize), method="lm")+
  facet_grid(Study~Classification, scales="free")+
  xlab("Relative rank in community\n(1 = most abundant)") +
  scale_color_gradient2(low="blue", high="red", midpoint=15)

# Correlation between tolerance range size and average rank?
rank_and_maxabund %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unassigned"))) %>%
  filter(Classification !="Unassigned") %>%
  ggplot() + 
  geom_point(aes(x=(medRank_rel), y=tolRangeSize, col=SalMax)) +
  geom_smooth(aes(x=(medRank_rel), y=tolRangeSize), method="lm")+
  facet_grid(Study~Classification, scales="free")+
  xlab("Relative rank in community\n(1 = most abundant)") +
  scale_color_gradient2(low="blue", high="red", midpoint=15)


### The best plot, but remove things that never reach <0.01 abundance
# Correlation between tolerance range size and best rank?
# rank_and_maxabund %>%
#   filter(maxAbund<0.001) %>%
#   mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unassigned"))) %>%
#   filter(Classification !="Unassigned") %>%
#   ggplot() + 
#   geom_point(aes(x=(maxRank_rel), y=tolRangeSize, col=SalMax)) +
#   geom_smooth(aes(x=(maxRank_rel), y=tolRangeSize), method="lm")+
#   facet_grid(Study~Classification, scales="free")+
#   xlab("Relative rank in community\n(1 = most abundant)") +
#   scale_color_gradient2(low="blue", high="red", midpoint=15)


gg_rank_vs_tol_highlowabund <- rank_and_maxabund %>%
  mutate(Abundant=ifelse(maxAbund<0.01, "Low abundance\n(<1% in all samples)", "High abundance\n(>1% in at least one sample)")) %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unassigned"))) %>%
  filter(Classification !="Unassigned") %>%
  ggplot() + 
  geom_point(aes(x=(maxRank_rel), y=tolRangeSize, col=SalMax, pch=Classification)) +
  # geom_smooth(data=rank_and_maxabund %>% filter(Classification%in%c("Marine"))
              # ,aes(x=(maxRank_rel), y=tolRangeSize, group=Classification), method="lm", se = FALSE, col="darkred")+
  # geom_smooth(data=rank_and_maxabund %>% filter(Classification%in%c("Fresh"))
              # ,aes(x=(maxRank_rel), y=tolRangeSize, group=Classification), method="lm", se = FALSE, col="darkblue")+
  facet_grid(Study~Abundant, scales="free")+
  xlab("Maximum observed\nrelative rank in community\n(1 = most abundant)") +
  ylab("Tolerance range size") +
  scale_color_gradient2(low="blue", high="red", midpoint=15)+
  scale_shape_manual(values=c(Fresh=25, Brackish=19, Marine=24))+
  theme_grey()
gg_rank_vs_tol_highlowabund
ggsave(filename = "07c_tolerance_ranges/gg_rank_vs_tol_highlowabund.png", gg_rank_vs_tol_highlowabund, 
       height=8, width=8)

gg_rank_vs_tol <- rank_and_maxabund %>%
  mutate(Abundant=ifelse(maxAbund<0.01, "Low abundance\n(<1% in all samples)", "High abundance\n(>1% in at least one sample)")) %>%
  mutate(Classification = factor(Classification, levels=c("Fresh","Brackish","Marine","Unassigned"))) %>%
  filter(Classification !="Unassigned") %>%
  ggplot() + 
  geom_point(aes(x=(maxRank_rel), y=tolRangeSize, col=SalMax, pch=Classification)) +
  # geom_smooth(data=rank_and_maxabund %>% filter(Classification%in%c("Marine"))
  #             ,aes(x=(maxRank_rel), y=tolRangeSize, group=Classification), method="lm", se = FALSE, col="darkred")+
  # geom_smooth(data=rank_and_maxabund %>% filter(Classification%in%c("Fresh"))
  #             ,aes(x=(maxRank_rel), y=tolRangeSize, group=Classification), method="lm", se = FALSE, col="darkblue")+
  facet_grid(Study~Classification, scales="free")+
  xlab("Maximum observed\nrelative rank in community\n(1 = most abundant)") +
  ylab("Tolerance range size") +
  scale_color_gradient2(low="blue", high="red", midpoint=15)+
  scale_shape_manual(values=c(Fresh=25, Brackish=19, Marine=24))+
  theme_grey()
gg_rank_vs_tol
ggsave(filename = "07c_tolerance_ranges/gg_rank_vs_tol.png", gg_rank_vs_tol, 
       height=8, width=6)
# 
# 
# 
# ## Tolerance plot, but intensity is the max rank within communities
# dat_wranks <- dat %>%
#   left_join(rank_and_maxabund) %>%
#   mutate(TaxaID=make.unique(TaxaID)) %>%
#   mutate(midSal = (lwrTol+uprTol)/2) %>%
#   mutate(Classification=ifelse(class=="high","Marine"
#                                , ifelse(class=="intermediate","Brackish",
#                                         ifelse(class=="low","Fresh","Unclassified")))) %>%
#   mutate(Classification=factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
#   arrange( Classification, midSal) %>%
#   mutate(TaxaID=factor(TaxaID, levels=unique(TaxaID)))  
# 
# dat_wranks%>%
#   # filter(Classification=="Fresh") %>%
#   ggplot() +
#   geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID, col=(maxRank_rel))) +
#   # facet_grid(.~Study, drop=TRUE, scales="free")+
#   theme(axis.text.y = element_blank())+
#   scale_color_gradient(low="white", high="darkgreen") +
#   # scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
#   xlab("Salinity") 
# 
# dat_wranks%>%
#   filter(Classification=="Brackish") %>%
#   ggplot() +
#   geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID, col=(maxRank_rel))) +
#   # facet_grid(.~Study, drop=TRUE, scales="free")+
#   theme(axis.text.y = element_blank())+
#   scale_color_gradient(low="white", high="darkgreen") +
#   # scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
#   xlab("Salinity") 

## Cluster based on observed places

allStudies <- unique(dat$Study)
allClusterDat <- data.frame()
for (s in allStudies) {
  tempmeta <- dat %>% filter(Study==s) %>%
    select(TaxaID, Study, Classification, SalMax)
  tempdat <- dat %>% filter(Study==s) %>%
    select(starts_with("SAMPLE_"))
  tempdat<- tempdat[,which(!is.na(colSums(tempdat)))]
  temp_kmeans <- kmeans((tempdat), centers =4)
  tempmeta$cluster <- temp_kmeans$cluster
  allClusterDat <- rbind(allClusterDat, tempmeta)
}
summaryClust <- allClusterDat %>% 
  group_by(Study, cluster) %>%
  summarise(aveSalMax = mean(SalMax)) %>% ungroup() %>%
  group_by(Study) %>%
  mutate(rankClust = rank(aveSalMax)) %>% ungroup()
# 
# allClusterDat_adj <- allClusterDat %>%
#   # left_join(summaryClust, relationship = "many-to-many") %>%
#   mutate(ClustClass = ifelse(rankClust==1,"Fresh", 
#                              ifelse(rankClust==2, "Lower Brackish",
#                                     ifelse(rankClust==3, "Upper Brackish",
#                                            ifelse(rankClust==4,"Marine", "TIE"))))) %>%
#   ungroup() %>%
#   select(TaxaID, Study, ClustClass, cluster)
  
dat %>%
  left_join(allClusterDat)  %>%
  mutate(TaxaID=make.unique(TaxaID)) %>%
  mutate(midSal = (lwrTol+uprTol)/2) %>%
  # mutate(Classification=ifelse(class=="high","Marine"
                               # , ifelse(class=="intermediate","Brackish",
                                        # ifelse(class=="low","Fresh","Unclassified")))) %>%
  mutate(Classification=factor(Classification, levels=c("Fresh","Brackish","Marine","Unclassified"))) %>%
  arrange( Classification, midSal) %>%
  mutate(TaxaID=factor(TaxaID, levels=unique(TaxaID)))  %>%
  ggplot() +
  geom_segment(aes(x=lwrObs, xend=uprObs, y=TaxaID, yend=TaxaID,col=factor(cluster)), alpha=0.25) +
  geom_segment(aes(x=lwrTol, xend=uprTol, y=TaxaID, yend=TaxaID,col=factor(cluster))) +
  facet_wrap(.~Study, drop=TRUE, scales="free_y")+
  theme(axis.text.y = element_blank())+
  # scale_color_manual(values=c(Fresh="blue",Brackish="purple",Marine="red",Unclassified="grey")) +
  xlab("Salinity") 
  
