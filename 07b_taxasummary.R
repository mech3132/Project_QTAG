#!bin/bash
library(tidyverse)
dir.create("07b_taxasummary")
dir.create("07b_taxasummary/taxasummaries_study")
dir.create("07b_taxasummary/taxasummaries_study_bysal")
dir.create("07b_taxasummary/taxasummaries_study_bysal_rare")
dir.create("07b_taxasummary/taxasummaries_bylevel")

# 
# mf <- read.delim("04_adjust_files_for_QTAG/all_meta.txt")
otucount_initial <- read.delim("04_adjust_files_for_QTAG/otu.txt")
colnames(otucount_initial) <- paste0("SAMPLE_",colnames(otucount_initial))
taxa <- read.delim("04_adjust_files_for_QTAG/taxonomy.txt")

meta <- read.delim("06_processing_QTAG/allMeta.txt")
dat <- read.delim("06_processing_QTAG/QTAG_output_with_taxa.txt")

### Filter out missing things
toKeep <- intersect(paste0("SAMPLE_",meta$SampleID), colnames(dat))
otucount <- otucount_initial %>% select(SAMPLE_TaxaID, one_of(toKeep)) %>% rename(TaxaID=SAMPLE_TaxaID)
meta <- meta %>% filter(paste0("SAMPLE_",SampleID) %in% colnames(otucount))
taxa <- filter(taxa, TaxaID %in% otucount$TaxaID & TaxaID %in% dat$TaxaID)

### Taxonomic summary across salinity
### FIX


set.seed(234)
allColorsRandomClass <- sample(size=length(unique(dat$Class)), colors()[-grep("white|beige|grey|gray",colors())], replace=TRUE)
names(allColorsRandomClass) <- unique(dat$Class)
allStudies <- unique(dat$Study)
for ( s in allStudies) {
  # s="chen16"
  lvl="Class"
  dat_for_bar <- dat %>% 
    filter(Study==s) %>%
    select(paste0(lvl), Study, one_of(meta$SampleID)) %>%
    # group_by(Study) %>% summarise(sum(Salinity_1))
    pivot_longer(-c(paste0(lvl), Study), names_to="SampleID", values_to="RA") %>%
    filter(!is.na(RA)) %>%
    left_join(meta %>% select(SampleID, Salinity, Season)) %>%
    # mutate(Salinity=round(Salinity)) %>%
    arrange(Salinity) %>%
    mutate(SampleID = factor(SampleID, levels=unique(SampleID))) %>%
    mutate(RoundedSal = factor(round(Salinity))) %>%
    full_join(data.frame(RoundedSal=as.character(seq(1,35)), RA = 0, SampleID = "MOCK", Class = "c_Bacteroidia")) %>%
    mutate(RoundedSal = factor(RoundedSal, levels=as.character(seq(0,35)))) %>%
    filter(!is.na(Season))
  # View(dat_for_bar)
  # bar_width_manual <- dat_for_bar %>% select(Study, SampleID, Salinity) %>% distinct() %>%
  #   group_by(Study, Salinity) %>%
  #   mutate(count=n()) %>% mutate(colWidth=1/count) %>% ungroup() %>% select(Study, SampleID, Salinity, count, colWidth)
  gg_indivtemp <- dat_for_bar %>%
    ggplot()+ geom_bar(aes(x=SampleID, y=RA, fill=get(lvl)), stat = "identity"
                       , show.legend = FALSE
    ) +
    scale_fill_manual(values=allColorsRandomClass)+
    theme_dark()+theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
    # scale_x_discrete(name="Salinity", labels=bar_width_manual$Salinity) +
    # theme(axis.text.x = element_text(angle=90,hjust=1, vjust=0.5)) +
    facet_grid(Season~RoundedSal, drop=TRUE, scale="free") +
    xlab("") +
    labs(title=paste0(lvl))
  nSeason <- length(unique(dat_for_bar$Season))
  # gg_indivtemp
  ggsave(filename = paste0("07b_taxasummary/taxasummaries_study/",s,".png"), gg_indivtemp
         , height=4*nSeason, width=15)
  
  ## Keep all salinities
  gg_indivtemp2 <- dat_for_bar %>%
    mutate(RoundedSal = factor(RoundedSal, levels=paste0(seq(0,35)))) %>%
    ggplot()+ geom_bar(aes(x=SampleID, y=RA, fill=get(lvl)), stat = "identity"
                       , show.legend = FALSE
    ) +
    scale_fill_manual(values=allColorsRandomClass)+
    theme_dark()+theme(axis.text.x = element_blank(), axis.ticks = element_blank())+
    # scale_x_discrete(name="Salinity", labels=bar_width_manual$Salinity) +
    # theme(axis.text.x = element_text(angle=90,hjust=1, vjust=0.5)) +
    facet_grid(Season~RoundedSal, drop=FALSE, scale="free") +
    xlab("") +
    labs(title=paste0(lvl))
  nSeason <- length(unique(dat_for_bar$Season))
  # gg_indivtemp
  ggsave(filename = paste0("07b_taxasummary/taxasummaries_study/",s,"_allsal.png"), gg_indivtemp2
         , height=4*nSeason, width=15)
  
}

## for by salinity
for ( s in allStudies) {
  lvl = "Class"
  dat_for_bar <- dat %>% 
    filter(Study==s) %>%
    select(paste0(lvl), Study, starts_with("Salinity_")) %>%
    pivot_longer(-c(paste0(lvl), Study), names_to="Salinity", values_to="RA") %>%
    filter(!is.na(RA)) %>%
    mutate(Salinity=gsub("Salinity_","",Salinity)) %>%
    arrange(as.numeric(Salinity)) 
    # left_join(meta %>% select(SampleID, Salinity, Season) %>% mutate(Salinity=as.numeric(Salinity))) 
  bar_width_manual <- dat_for_bar %>% select(Study, Salinity) %>% distinct() %>%
    group_by(Study, Salinity) %>% 
    mutate(count=n()) %>% mutate(colWidth=1/count) %>% ungroup() %>% select(Study, Salinity, count, colWidth)
  gg_indivtemp <- dat_for_bar %>%
    ggplot()+ geom_bar(aes(x=Salinity, y=RA, fill=get(lvl)), stat="identity"
                       , show.legend = FALSE
    ) +
    scale_fill_manual(values=allColorsRandomClass)+
    scale_x_discrete(name="Salinity", labels=bar_width_manual$Salinity) +
    theme(axis.text.x = element_text(angle=90,, hjust=1, vjust=0.5))
  # gg_indivtemp
  ggsave(filename = paste0("07b_taxasummary/taxasummaries_study_bysal/",s,".png"), gg_indivtemp
         , height=4, width=6)
}
# 
# # RArefied salinity
# for ( s in allStudies) {
#   # s="chen16"
#   dat_for_bar <- dat %>%
#     filter(Study==s) %>%
#     select(Class, Study, starts_with("SalinityRare_")) %>%
#     # group_by(Study) %>% summarise(sum(SalinityRare_0, na.rm=TRUE))
#     pivot_longer(-c(Class, Study), names_to="Salinity", values_to="RA") %>%
#     filter(!is.na(RA)) %>%
#     mutate(Salinity=gsub("SalinityRare_","",Salinity)) %>%
#     # left_join(meta %>% select(SampleID, Salinity)) %>%
#     # mutate(Salinity=round(Salinity)) %>%
#     arrange(as.numeric(Salinity)) 
#   # mutate(SampleID = factor(SampleID, levels=unique(SampleID))) 
#   bar_width_manual <- dat_for_bar %>% select(Study, Salinity) %>% distinct() %>%
#     group_by(Study, Salinity) %>% 
#     mutate(count=n()) %>% mutate(colWidth=1/count) %>% ungroup() %>% select(Study, Salinity, count, colWidth)
#   gg_indivtemp <- dat_for_bar %>%
#     ggplot()+ geom_bar(aes(x=Salinity, y=RA, fill=Class), stat="identity"
#                        , show.legend = FALSE
#     ) +
#     scale_fill_manual(values=allColorsRandomClass)+
#     scale_x_discrete(name="Salinity", labels=bar_width_manual$Salinity) +
#     theme(axis.text.x = element_text(angle=90,, hjust=1, vjust=0.5))
#   ggsave(filename = paste0("07b_taxasummary/taxasummaries_study_bysal_rare/",s,".png"), gg_indivtemp
#          , height=4, width=6)
# }


gg_forlegendClass <- ggplot(data.frame(Taxa=names(allColorsRandomClass))) + 
  geom_bar(aes(x=Taxa, fill=Taxa), col=NA) +
  scale_fill_manual(values=allColorsRandomClass)
ggsave(filename = paste0("07b_taxasummary/taxasummaries_study/forlegend.png"), gg_forlegendClass
       , height=10, width=10)


#### CLASS Grouped, but by sal summaries ####

##### Just class, all studies
dat16_long_bysal <- dat %>% 
  filter(Study!="chen18") %>%
  select(Class, Study, starts_with("Salinity_")) %>%
  # group_by(Study) %>% summarise(sum(Salinity_1))
  pivot_longer(-c(Class, Study), names_to="Salinity", values_to="RA") %>%
  filter(!is.na(RA)) %>%
  mutate(Salinity=as.numeric(gsub("Salinity_","",Salinity)))

highabundclasses <- dat16_long_bysal %>% filter(RA>=0.01) %>% pull(Class) %>% unique()
highabundcol_only <- allColorsRandomClass[highabundclasses]
# allColorsRandomClass_grey <- allColorsRandomClass
# allColorsRandomClass_grey[lowabundclasses] <- NA
gg_bySal_class16 <- dat16_long_bysal %>%
  ggplot()+ geom_bar(aes(x=Salinity, y=RA, fill=Class), stat="identity", show.legend = FALSE) +
  facet_wrap(.~Study) +
  scale_fill_manual(values=highabundcol_only)+theme_dark()+
  ylab("Relative Abundance") + xlab("Salinity (ppt)")
gg_bySal_class16
ggsave(filename = paste0("07b_taxasummary/summary_class16_bysal.png"), gg_bySal_class16
       , height=4, width=8)

# allColorsRandomClass_greyna <- allColorsRandomClass_grey[!is.na(allColorsRandomClass_grey)]
gg_forlegendClass_grey <- ggplot(data.frame(Taxa=names(highabundcol_only))) + 
  geom_bar(aes(x=Taxa, fill=Taxa), col=NA) +
  scale_fill_manual(values=highabundcol_only)
ggsave(filename = paste0("07b_taxasummary/forlegend_classgrey_bysal.png"), gg_forlegendClass_grey
       , height=10, width=10)


#####3 18S

dat18_long_bysal <- dat %>% 
  filter(Study=="chen18") %>%
  select(Class, Study, starts_with("Salinity_")) %>%
  # group_by(Study) %>% summarise(sum(Salinity_1))
  pivot_longer(-c(Class, Study), names_to="Salinity", values_to="RA") %>%
  filter(!is.na(RA)) %>%
  mutate(Salinity=as.numeric(gsub("Salinity_","",Salinity)))
highabundclasses18 <- dat18_long_bysal %>% filter(RA>=0.01) %>% pull(Class) %>% unique()
highabundcol_only18 <- allColorsRandomClass[highabundclasses18]
# allColorsRandomClass_grey <- allColorsRandomClass
# allColorsRandomClass_grey[lowabundclasses] <- NA
gg_bySal_class18 <- dat18_long_bysal %>%
  ggplot()+ geom_bar(aes(x=Salinity, y=RA, fill=Class), stat="identity", show.legend = FALSE) +
  facet_wrap(.~Study) +
  scale_fill_manual(values=highabundcol_only18)+theme_dark()+
  ylab("Relative Abundance") + xlab("Salinity (ppt)")
gg_bySal_class18
ggsave(filename = paste0("07b_taxasummary/summary_class18_bysal.png"), gg_bySal_class18
       , height=4, width=8)

# allColorsRandomClass_greyna <- allColorsRandomClass_grey[!is.na(allColorsRandomClass_grey)]
gg_forlegendClass_grey18 <- ggplot(data.frame(Taxa=names(highabundcol_only18))) + 
  geom_bar(aes(x=Taxa, fill=Taxa), col=NA) +
  scale_fill_manual(values=highabundcol_only18)
ggsave(filename = paste0("07b_taxasummary/forlegend_classgrey18_bysal.png"), gg_forlegendClass_grey18
       , height=10, width=10)

###### Go through all levles #########
set.seed(8019)
colRandTemp <- sample(size=length(unique(dat$Class)), colors()[-grep("white|beige|grey|gray",colors())], replace=TRUE)
for ( l in c("Class","Order","Family","Genus")) {
  dat16_long_bysal <- dat %>% 
    filter(Study!="chen18") %>%
    select(paste0(l), Study, starts_with("Salinity_")) %>%
    # group_by(Study) %>% summarise(sum(Salinity_1))
    pivot_longer(-c(paste0(l), Study), names_to="Salinity", values_to="RA") %>%
    filter(!is.na(RA)) %>%
    rename(Group=paste0(l)) %>%
    mutate(Salinity=as.numeric(gsub("Salinity_","",Salinity)))
  highabundgroups <- dat16_long_bysal %>% group_by(Group, Study) %>% 
    summarise(maxRA = max(RA)) %>% ungroup() %>% 
    group_by(Study) %>% mutate(abndRank = rank(-maxRA, ties.method = "max")) %>%
    ungroup() %>% filter(abndRank<=20) %>% pull(Group) %>% unique()
  highabundcol_only <- colRandTemp[1:length(highabundgroups)]
  names(highabundcol_only) <- highabundgroups
  # allColorsRandomClass_grey <- allColorsRandomClass
  # allColorsRandomClass_grey[lowabundclasses] <- NA
  gg_bySal_temp <- dat16_long_bysal %>%
    ggplot()+ geom_bar(aes(x=factor(Salinity), y=RA, fill=Group), stat="identity", show.legend = FALSE) +
    facet_wrap(.~Study) +
    scale_fill_manual(values=highabundcol_only)+theme_bw()+
    ylab("Relative Abundance") + xlab("Salinity (ppt)")+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
  gg_bySal_temp
  ggsave(filename = paste0("07b_taxasummary/taxasummaries_bylevel/summary_bysal_",l,".png"), gg_bySal_temp
         , height=4, width=8)
  gg_forlegend16 <- ggplot(data.frame(Taxa=names(highabundcol_only))) + 
    geom_bar(aes(x=Taxa, fill=Taxa), col=NA) +
    scale_fill_manual(values=highabundcol_only)
  ggsave(filename = paste0("07b_taxasummary/taxasummaries_bylevel/forlegend_bysal_",l,".png"), gg_forlegend16
         , height=10, width=10)
  
}

# 
# for ( l in c("Class","Order","Family","Genus")) {
#   # l="Genus"
#   # allGroups <- dat %>% select(paste0(l)) %>% pull() %>% unique()
#   # allColsTemp <- sample(size = length(allGroups), allColorsRandomClass, replace = TRUE)
#   # names(allColsTemp) <- NULL
#   # names(allColsTemp) <- allGroups
# gg_bySal <- dat %>% 
#     select(paste0(l), Study, starts_with("Salinity")) %>%
#     # group_by(Study) %>% summarise(sum(Salinity_1))
#     pivot_longer(-c(paste0(l), Study), names_to="Salinity", values_to="RA") %>%
#     filter(!is.na(RA)) %>%
#     mutate(Salinity=as.numeric(gsub("Salinity_","",Salinity))) %>%
#     # group_by(Genus, Study, Salinity) %>%
#     # summarise(RA=mean(RA)) %>% ungroup() %>%
#     ggplot()+ geom_bar(aes(x=Salinity, y=RA, fill=get(l)), stat="identity", show.legend = FALSE) +
#     facet_wrap(.~Study) 
# gg_bySal
# ggsave(filename = paste0("07b_taxasummary/summary_",l,".png"), gg_bySal
#        , height=4, width=8)
# }



# otucount_mat <- as.matrix(otucount[,-1])
# rownames(otucount_mat) <- otucount[,1]
# 
# otu_mat_RA <- t(t(otucount_mat)/colSums(otucount_mat))
# # Collapse by level
# for ( l in c("Class","Order","Family","Genus")) {
#   # l="Family"
#   tax_col <- taxa %>% 
#     filter(TaxaID %in% rownames(otu_mat_RA)) %>%
#     select(TaxaID, paste0(l)) %>% distinct() %>%
#     mutate(pres=1) %>% 
#     pivot_wider(names_from=paste0(l), values_from=pres, values_fill=0)
#   tax_col_mat <- as.matrix(tax_col[,-1])
#   rownames(tax_col_mat) <- pull(tax_col[,1])
#   OTU_collapsed <- t(t(otu_mat_RA)%*% tax_col_mat)
#   save(OTU_collapsed, file = paste0("07b_taxasummary/OTU_collapsed_",l,".RData"))
#   # Plot
#   gg_otucollapsed <- OTU_collapsed %>% as.data.frame() %>%
#     rownames_to_column(var="Group") %>%
#     pivot_longer(-Group, names_to="SampleID", values_to="RA") %>%
#     left_join(meta) %>%
#     mutate(Salinity=as.numeric(Salinity)) %>%
#     arrange(as.numeric(Salinity)) %>%
#     mutate(SampleID = factor(SampleID, levels=unique(SampleID)))%>%
#     mutate(SalGroup=ifelse(Salinity<10, "Low", ifelse(Salinity <20,"Brackish", "Marine"))) %>%
#     ggplot() +
#     geom_bar(aes(x=SampleID, y=RA, fill=Group), stat="identity",  show.legend = FALSE) +
#     facet_wrap(.~Study, drop=TRUE, scales="free") +
#     labs(title=paste0(l))
#   gg_otucollapsed
#   ### WORKING; NEED TO PLTO  
# }


