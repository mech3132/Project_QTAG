#!bin/bash
library(ape)
library(tidyverse)

dir.create("04_adjust_files_for_QTAG")
##### Load  Clustered ###########
# meta$Study %>% unique()
allStudies <- c("chen16","chen18","herleman","wang","hu16","campbell2013")

##### Load ASVs ###########
meta <- read.delim("02b_adjust_files_for_processing/metadata_adj.txt") %>%
  filter(!is.na(Salinity))
otu16 <- read.delim("03a_assign_taxonomy/otu_16S_exported/feature-table.txt", skip=1) 
otu18 <- read.delim("03a_assign_taxonomy/otu_18S_exported/feature-table.txt", skip=1)
tree16 <- read.tree("03b_phylogenetic_tree/filtered-tree_16S/exported_16S_tree/tree.nwk")$tip.label
tree18 <- read.tree("03b_phylogenetic_tree/filtered-tree_18S/exported_18S_tree/tree.nwk")$tip.label

# taxa <- read.delim("02_adjust_files_for_processing/taxa_adj.txt")
taxa16 <- read.delim("03a_assign_taxonomy/vsearch_classify_id_80_taxonomy16/export/taxonomy.tsv") %>%
  rename(TaxaID=Feature.ID) %>%
  filter(TaxaID %in% otu16$X.OTU.ID) %>%
  # rowwise() %>%
  # filter(length(grep("Chloro|Mitochon|Eukary", Taxon))>0) %>% ungroup()
  separate(Taxon, into=c("Domain_NA","Phylum_NA","Class_NA","Order_NA","Family_NA","Genus_NA","Species_NA"), sep ="; ", remove=FALSE) %>%
  mutate(Phylum = ifelse(is.na(Phylum_NA), Domain_NA, Phylum_NA)
         , Class = ifelse(is.na(Class_NA), Phylum, Class_NA)
         , Order = ifelse(is.na(Order_NA), Class, Order_NA)
         , Family = ifelse(is.na(Family_NA), Order, Family_NA)
         , Genus = ifelse(is.na(Genus_NA), Family, Genus_NA)
         , Species = ifelse(is.na(Species_NA), Genus, Species_NA)
  )
taxa18 <- read.delim("03a_assign_taxonomy/vsearch_classify_id_80_taxonomy18/export/taxonomy.tsv") %>%
  rename(TaxaID=Feature.ID) %>%
  filter(TaxaID %in% otu18$X.OTU.ID) %>%
  separate(Taxon, into=c("Domain_NA","Phylum_NA","Class_NA","Order_NA","Family_NA","Genus_NA","Species_NA"), sep ="; ", remove=FALSE) %>%
  mutate(Phylum = ifelse(is.na(Phylum_NA), Domain_NA, Phylum_NA)
         , Class = ifelse(is.na(Class_NA), Phylum, Class_NA)
         , Order = ifelse(is.na(Order_NA), Class, Order_NA)
         , Family = ifelse(is.na(Family_NA), Order, Family_NA)
         , Genus = ifelse(is.na(Genus_NA), Family, Genus_NA)
         , Species = ifelse(is.na(Species_NA), Genus, Species_NA)
  )
# Merge 16 and 18S taxonomies
allTaxa <- rbind(taxa16, taxa18) %>% distinct()
# Readjust metadata for X's
colnames(otu16) <- gsub("[.]","-",colnames(otu16))
colnames(otu16) <- gsub("Baltic[-]","Baltic.", colnames(otu16))
colnames(otu18) <- gsub("[.]","-",colnames(otu18))
colnames(otu18) <- gsub("Baltic[-]","Baltic.", colnames(otu18))

##### ASV filter chloro, mito, euk ###########

# Look for chloro/mito/euk
taxcolmat16 <- taxa16 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse(length(grep("Chloropl|Mitochon", Taxon))>0 | Domain_NA %in% c("d__Archaea","Unassigned","d__Eukaryota"),"Chloro_Mito_Euk_Archaea_Unassigned","Bacteria")) %>%
  ungroup() %>% select(TaxaID, REMOVE) %>%
  mutate(presence=1) %>%
  pivot_wider(names_from=REMOVE, values_from=presence, values_fill = 0)
taxcolmat16_mat <- as.matrix(taxcolmat16[,-1])
rownames(taxcolmat16_mat) <-pull( taxcolmat16[,1])
otu16mat <- as.matrix(otu16[,-1])
rownames(otu16mat) <- otu16[,1]
otu16mat <- otu16mat[rownames(taxcolmat16_mat),]
otu16COLLAPSED <- t(taxcolmat16_mat)%*%otu16mat
gg_chloromitoeuk <- t(t(otu16COLLAPSED)/colSums(otu16COLLAPSED)) %>%
  as.data.frame() %>% rownames_to_column(var="ASVtype") %>%
  pivot_longer(-ASVtype, names_to = "X.SampleID", values_to = "RelativeAbundance") %>%
  left_join(meta) %>%
  arrange(-as.numeric(Salinity)) %>% mutate(X.SampleID=factor(X.SampleID, levels=unique(X.SampleID))) %>%
  ggplot() + 
  geom_bar(aes(x=X.SampleID, y=RelativeAbundance, fill=ASVtype), stat="identity") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(.~Study, scales = "free", drop=TRUE)
gg_chloromitoeuk
ggsave(filename = "04_adjust_files_for_QTAG/gg_plot_eukmitochlor_remove16S.png",gg_chloromitoeuk, height=8, width=8)

# Look for unassigned euks
taxa18$Domain_NA %>% unique()
taxcolmat18 <- taxa18 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse( Domain_NA %in% c("Unassigned","d__Bacteria","d__Archaea"),"Unassigned","Eukaryote")) %>%
  ungroup() %>% select(TaxaID, REMOVE) %>%
  mutate(presence=1) %>%
  pivot_wider(names_from=REMOVE, values_from=presence, values_fill = 0)
taxcolmat18_mat <- as.matrix(taxcolmat18[,-1])
rownames(taxcolmat18_mat) <-pull( taxcolmat18[,1])
otu18mat <- as.matrix(otu18[,-1])
rownames(otu18mat) <- otu18[,1]
otu18mat <- otu18mat[rownames(taxcolmat18_mat),]
otu18COLLAPSED <- t(taxcolmat18_mat)%*%otu18mat
gg_eukunassigned <- t(t(otu18COLLAPSED)/colSums(otu18COLLAPSED)) %>%
  as.data.frame() %>% rownames_to_column(var="ASVtype") %>%
  pivot_longer(-ASVtype, names_to = "X.SampleID", values_to = "RelativeAbundance") %>%
  left_join(meta) %>%
  arrange(-as.numeric(Salinity)) %>% mutate(X.SampleID=factor(X.SampleID, levels=unique(X.SampleID))) %>%
  ggplot() + 
  geom_bar(aes(x=X.SampleID, y=RelativeAbundance, fill=ASVtype), stat="identity") +
  theme(axis.text.x = element_blank()) +
  facet_wrap(.~Study, scales = "free", drop=TRUE)
gg_eukunassigned
ggsave(filename = "04_adjust_files_for_QTAG/gg_plot_unassigned_remove18S.png",gg_chloromitoeuk, height=8, width=8)

### STATS FOR REMOVAL
taxa16_nremoved <- taxa16 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse(length(grep("Chloropl|Mitochon", Taxon))>0 | Domain_NA %in% c("d__Archaea","Unassigned","d__Eukaryota"),"Chloro_Mito_Euk_Archaea_Unassigned","Bacteria")) %>%
  ungroup() %>% select(REMOVE) %>% table()
taxa16_nremoved
taxa18_nremoved <- taxa18 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse( Domain_NA %in% c("Unassigned"),"Unassigned","Eukaryote")) %>%
  ungroup() %>% select(REMOVE) %>% table()
taxa18_nremoved

# PROPORITON OF DATASET
meta16studycol <- meta %>% filter(target_region=="16S") %>% select(X.SampleID, Study) %>%
  mutate(value=1) %>% pivot_wider(names_from=Study, values_from=value, values_fill = 0)
meta16studycol_mat <- as.matrix(meta16studycol[,-1])
rownames(meta16studycol_mat) <- pull(meta16studycol[,1])
# Remove extra sample from otu table
otu16COLLAPSED_filt <- otu16COLLAPSED[,which(colnames(otu16COLLAPSED) %in% rownames(meta16studycol_mat))]
meta16studycol_mat <- meta16studycol_mat[match(colnames(otu16COLLAPSED_filt),rownames(meta16studycol_mat)),]
nseqs_removed_16S <- otu16COLLAPSED_filt %*% meta16studycol_mat
prop_removed_16S <- t(nseqs_removed_16S)/colSums(nseqs_removed_16S)
prop_removed_16S



## for 18S
meta18studycol <- meta %>% filter(target_region=="18S") %>% select(X.SampleID, Study) %>%
  mutate(value=1) %>% pivot_wider(names_from=Study, values_from=value, values_fill = 0)
meta18studycol_mat <- as.matrix(meta18studycol[,-1])
rownames(meta18studycol_mat) <- pull(meta18studycol[,1])
# Remove extra sample from otu table
otu18COLLAPSED_filt <- otu18COLLAPSED[,which(colnames(otu18COLLAPSED) %in% rownames(meta18studycol_mat))]
meta18studycol_mat <- meta18studycol_mat[match(colnames(otu18COLLAPSED_filt),rownames(meta18studycol_mat)),]
nseqs_removed_18S <- otu18COLLAPSED_filt %*% meta18studycol_mat
prop_removed_18S <- t(nseqs_removed_18S)/colSums(nseqs_removed_18S)
prop_removed_18S

sink(file="04_adjust_files_for_QTAG/summary_of_things_removed.txt")
print("16S removed")
taxa16_nremoved
nseqs_removed_16S
prop_removed_16S
print("18S removed")
taxa18_nremoved
nseqs_removed_18S
prop_removed_18S
sink()



######## REMOVE #######
taxa16_tokeep <- taxa16 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse(length(grep("Chloropl|Mitochon", Taxon))>0 | Domain_NA %in% c("d__Archaea","Unassigned","d__Eukaryota"),"Chloro_Mito_Euk_Archaea_Unassigned","Bacteria")) %>%
  mutate(TREEREMOVE=ifelse(TaxaID %in% tree16, "KEEP","REMOVE")) %>%
  ungroup() %>%
  filter(REMOVE=="Bacteria", TREEREMOVE=="KEEP") %>% pull(TaxaID)

taxa18_tokeep <-  taxa18 %>%
  rowwise() %>%
  mutate(REMOVE=ifelse( Domain_NA %in% c("Unassigned","d__Archaea","d__Bacteria"),"Unassigned","Eukaryote")) %>%
  mutate(TREEREMOVE=ifelse(TaxaID %in% tree18, "KEEP","REMOVE")) %>%
  ungroup()%>%
  filter(REMOVE=="Eukaryote") %>% pull(TaxaID)

otu16_toKeep <- otu16 %>% filter(`X-OTU-ID` %in% c(taxa16_tokeep))
otu16_toKeep <- otu16_toKeep[rowSums(otu16_toKeep[,-1])>0,]
otu18_toKeep <- otu18 %>% filter(`X-OTU-ID` %in% c(taxa18_tokeep))
otu18_toKeep <- otu18_toKeep[rowSums(otu18_toKeep[,-1])>0,]

taxa16 %>% filter(TaxaID %in% otu16_toKeep$`X-OTU-ID`) %>% select(Domain_NA) %>% table()
taxa18 %>% filter(TaxaID %in% otu18_toKeep$`X-OTU-ID`)%>% select(Domain_NA) %>% table()

intersect(otu16_toKeep$`X-OTU-ID`, otu18_toKeep$`X-OTU-ID`)

allOTU <- full_join(otu16_toKeep, otu18_toKeep)

####### Final sample check ###########
### Make NAs zero
allOTU[is.na(allOTU)] <- 0
# Any samples NOT in meta?
any(!colnames(allOTU[,-1]) %in% meta$X.SampleID)
toRemoveOTUpos <- which(!colnames(allOTU[,-1]) %in% meta$X.SampleID)
toRemoveOTU <- c(colnames(allOTU)[toRemoveOTUpos],"CON-MC_chen16")
# Any meta not in OTU?
any(! meta$X.SampleID %in%colnames(allOTU[,-1]))
toRemoveMeta <- which(! meta$X.SampleID %in%colnames(allOTU[,-1]))
meta$X.SampleID[toRemoveMeta]

# Remove targets
meta_final <- meta %>% filter(!X.SampleID %in% meta$X.SampleID[toRemoveMeta]) %>%
  filter(!is.na(Season))

toRemoveSamples <- which(colnames(allOTU)%in% toRemoveOTU)
if (length(toRemoveSamples) > 0 ) {
  allOTU_final <- allOTU[,-toRemoveSamples]
} else {
  allOTU_final <- allOTU
}
# allOTU99_final <- allOTU99[,-which(colnames(allOTU99)%in% toRemoveOTU)]


# Double check OTU names are same
# any(!colnames(allOTU_final[,-1]) %in%  colnames(allOTU99_final[,-1]))
# which(!colnames(allOTU_final[,-1]) %in%  colnames(allOTU99_final[,-1]))

# Any samples NOT in meta?
# any(!colnames(allOTU99[,-1]) %in% meta$X.SampleID)
# any(! meta$X.SampleID %in%colnames(allOTU99[,-1]))

### Filter taxonomy list
# All OTU in tables has taxonomy?
allOTU_final$`X-OTU-ID`[which(!allOTU_final$`X-OTU-ID` %in% allTaxa$TaxaID)]
# allOTU99_final$X.OTU.ID[which(!allOTU99_final$X.OTU.ID%in% allTaxa99$TaxaID)]
# 
# # Check that there are no overlapping OTUs  in 16S and 18S
# otu16_temp <- allOTU_final %>% select(one_of(colnames(otu16)))
# otu16_temp <- otu16_temp[rowSums(otu16_temp[,-1])>0,]
# otu18_temp <- allOTU_final %>% select(one_of(colnames(otu18)))
# otu18_temp <- otu18_temp[rowSums(otu18_temp[,-1])>0,]
# 
# overlapASVs <- intersect(otu16_temp$`X-OTU-ID`,otu18_temp$`X-OTU-ID`)
# # allTaxa %>% filter(TaxaID %in% overlapASVs) %>% View()



######## raw asvs
dir.create(paste0("04_adjust_files_for_QTAG/split_studies"))
# dir.create(paste0("04_adjust_files_for_QTAG/split_studies_clust99"))
for ( s in allStudies) {
  # s="campbell2012"
  meta_filt <- meta_final %>% filter(Study == s) %>% rename(`#SampleID` = X.SampleID) 
  allSamps <- meta_filt$`#SampleID`[which(meta_filt$`#SampleID` %in% colnames(allOTU_final))]
  meta_filt <- meta_filt %>% filter(`#SampleID` %in% allSamps) 
  
  otu_filt <- allOTU_final %>% select("X-OTU-ID", one_of(meta_filt$`#SampleID`)) %>%
    rename(`#OTU ID`=`X-OTU-ID`)
  otu_filt <- otu_filt[rowSums(otu_filt[,-1])>0,]
  # otu_filt99 <- allOTU99_final %>% select("X.OTU.ID", one_of(meta_filt$`#SampleID`)) %>%
    # rename(`#OTU ID`=`X.OTU.ID`)
  # otu_filt99 <- otu_filt99[rowSums(otu_filt99[,-1])>0,]
  
  dir.create(paste0("04_adjust_files_for_QTAG/split_studies/",s))
  # dir.create(paste0("04_adjust_files_for_QTAG/split_studies_clust99/",s))
  # Write out
  write.table(meta_filt, file = paste0("04_adjust_files_for_QTAG/split_studies/",s,"/",s,"_meta.txt"), quote=FALSE, row.names=FALSE, sep="\t")
  write.table(otu_filt, file = paste0("04_adjust_files_for_QTAG/split_studies/",s,"/",s,"_otu.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
  # write.table(otu_filt99, file = paste0("04_adjust_files_for_QTAG/split_studies_clust99/",s,"/",s,"_otu.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
  
  }

## Write whole versions
colnames(allOTU_final)[1] <- "TaxaID"
# colnames(allOTU99_final)[1] <- "TaxaID"
meta_final <- meta_final %>% rename(SampleID = X.SampleID)
write.table(meta_final, file = paste0("04_adjust_files_for_QTAG/all_meta.txt"), quote=FALSE, row.names=FALSE, sep="\t")
write.table(allOTU_final, file = paste0("04_adjust_files_for_QTAG/otu.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
# write.table(allOTU99_final, file = paste0("04_adjust_files_for_QTAG/otu_clust99.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
write.table(allTaxa, file = paste0("04_adjust_files_for_QTAG/taxonomy.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
# write.table(allTaxa99, file = paste0("04_adjust_files_for_QTAG/taxonomy_clust99.txt"), quote=FALSE, row.names=FALSE, sep = "\t")

## Split into 16S and 18S
meta16 <- meta_final %>% filter(target_region=="16S")
allOTU_final_16 <- allOTU_final %>% select(TaxaID, one_of(meta16$SampleID))
allOTU_final_16 <- allOTU_final_16[which(rowSums(allOTU_final_16[,-1])>0),]
colnames(allOTU_final_16)[1] <- "#OTU ID"
write.table(allOTU_final_16, file = paste0("04_adjust_files_for_QTAG/otu16.txt"), quote=FALSE, row.names=FALSE, sep = "\t")

meta18 <- meta_final %>% filter(target_region=="18S")
allOTU_final_18 <- allOTU_final %>% select(TaxaID, one_of(meta18$SampleID))
allOTU_final_18 <- allOTU_final_18[which(rowSums(allOTU_final_18[,-1])>0),]
colnames(allOTU_final_18)[1] <- "#OTU ID"
write.table(allOTU_final_18, file = paste0("04_adjust_files_for_QTAG/otu18.txt"), quote=FALSE, row.names=FALSE, sep = "\t")

########### By study/season #######

dir.create(paste0("04_adjust_files_for_QTAG/split_studies_seasons"))
# dir.create(paste0("04_adjust_files_for_QTAG/split_studies_clust99"))
for ( s in allStudies) {
  # s="chen16"
  allSeasons <- meta_final %>% filter(Study==s) %>% pull(Season) %>% unique()
  for (s2 in allSeasons ) {
    
    meta_filt <- meta_final %>% filter(Study == s, Season==s2) %>% rename(`#SampleID` = SampleID) 
    allSamps <- meta_filt$`#SampleID`[which(meta_filt$`#SampleID` %in% colnames(allOTU_final))]
    meta_filt <- meta_filt %>% filter(`#SampleID` %in% allSamps) 
    
    otu_filt <- allOTU_final %>% select("TaxaID", one_of(meta_filt$`#SampleID`)) %>%
      rename(`#OTU ID`=`TaxaID`)
    otu_filt <- otu_filt[rowSums(otu_filt[,-1])>0,]
    # otu_filt99 <- allOTU99_final %>% select("X.OTU.ID", one_of(meta_filt$`#SampleID`)) %>%
    # rename(`#OTU ID`=`X.OTU.ID`)
    # otu_filt99 <- otu_filt99[rowSums(otu_filt99[,-1])>0,]
    
    dir.create(paste0("04_adjust_files_for_QTAG/split_studies_seasons/",s,s2))
    # dir.create(paste0("04_adjust_files_for_QTAG/split_studies_clust99/",s))
    # Write out
    write.table(meta_filt, file = paste0("04_adjust_files_for_QTAG/split_studies_seasons/",s,s2,"/",s,s2,"_meta.txt"), quote=FALSE, row.names=FALSE, sep="\t")
    write.table(otu_filt, file = paste0("04_adjust_files_for_QTAG/split_studies_seasons/",s,s2,"/",s,s2,"_otu.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
    
  }
}
meta_final %>%select(Study,Season) %>% table()

