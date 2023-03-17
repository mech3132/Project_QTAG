library(tidyverse)
dir.create("02b_adjust_files_for_processing")

load("02a_adjust_files_for_processing/meta.RData")
load("02a_adjust_files_for_processing/otu.RData")
load("02a_adjust_files_for_processing/taxa.RData")
varToKeep <- c("Study","Year","Month","Day","Season","Salinity","Lat","Long","temp_C", "FWprimerName","RVprimerName", "target_region","StudySeason")
# meta_intermediate$StudySeason

## Collapse by CoLRep; need to remove - for R
meta_mat_temp <- meta_intermediate %>% 
  select(`#SampleID`, ColRep)  %>%mutate(presence=1) %>%
  pivot_wider(names_from=ColRep, values_from=presence, values_fill=0) 
meta.mat <- as.matrix(meta_mat_temp[,-1])
rownames(meta.mat) <- pull(meta_mat_temp[,1])

 #any mf sample names not in otu table?
any(!rownames(meta.mat) %in% colnames(otu_temp))
otu_mat <- as.matrix(otu_temp[,rownames(meta.mat)])
# any(is.na(otu_mat))
# otu_mat[is.na(otu_mat)] <- 0
rownames(otu_mat) <- otu_temp$`#OTU ID`
which(duplicated(otu_temp$`#OTU ID`))
# tail(rownames(otu_mat))
# CHECK SAMPLE NAMES MATCH
ncol(otu_mat) == nrow(meta.mat)
otu.collapsed <- as.data.frame(otu_mat %*% meta.mat) %>% rownames_to_column(var="#OTU ID")
# View(otu.collapsed)
meta_final <- meta_intermediate %>% 
  mutate(ColRep = gsub("[-]","_",ColRep)) %>%
  select(ColRep, one_of(varToKeep)) %>% select(-temp_C) %>%
  rename(`#SampleID`=ColRep) %>% 
  distinct() %>%
  filter(`#SampleID` %in% colnames(otu.collapsed)) 

otu_adj <- otu.collapsed %>% select("#OTU ID", one_of(meta_final$`#SampleID`)) 

ncol(otu_adj)
nrow(meta_final)
meta_final %>% select(Study) %>%
  table()

# Make ASV IDs
asvList <- data.frame(TaxaID = paste0("asv",seq(1,length(otu_adj$`#OTU ID`))), asv_sequence = otu_adj$`#OTU ID`)
rename_taxa <- asvList[match(otu_adj$`#OTU ID`, asvList$asv_sequence),"TaxaID"]
otu_adj$`#OTU ID` <- rename_taxa


# asvList %>% View()
# View(otu_adj)
write.table(meta_final, file = "02b_adjust_files_for_processing/metadata_adj.txt", quote=FALSE, row.names=FALSE, sep="\t")
write.table(otu_adj, file = "02b_adjust_files_for_processing/otu_adj.txt", quote=FALSE, row.names=FALSE, sep = "\t")
write.table(asvList, file = "02b_adjust_files_for_processing/asvList.txt", quote = FALSE, row.names=FALSE, sep = "\t")
### Write 16S and 18S tables
all16SSamples <- meta_final %>% filter(target_region=="16S") %>% pull(`#SampleID`)
all18SSamples <- meta_final %>% filter(target_region=="18S") %>% pull(`#SampleID`)

otu16<-otu_adj%>% select(`#OTU ID`,one_of(all16SSamples))
otu16filt <- otu16[rowSums(otu16[,-1])>0,]

otu18<-otu_adj%>% select(`#OTU ID`,one_of(all18SSamples))
otu18filt <- otu18[rowSums(otu18[,-1])>0,]

write.table(otu16filt, file = "02b_adjust_files_for_processing/otu_adj_16S.txt", quote=FALSE, row.names=FALSE, sep = "\t")
write.table(otu18filt, file = "02b_adjust_files_for_processing/otu_adj_18S.txt", quote=FALSE, row.names=FALSE, sep = "\t")

#### Split into studies to do closed reference clustering ######
allStudies <- meta_final %>% pull(Study) %>% unique()
dir.create("02b_adjust_files_for_processing/studies_split")
for ( s in allStudies) {
  dir.create(paste0("02b_adjust_files_for_processing/studies_split/",s))
  samptemp <- meta_final %>% filter(Study==s) %>% select(`#SampleID`) %>% pull()
  metatemp <- meta_final %>% filter(Study==s)
  otu_temp <- otu_adj[,c("#OTU ID",samptemp)]
  otu_temp_filt <- otu_temp[rowSums(otu_temp[,-1])>0,]
  write.table(otu_temp_filt, file = paste0("02b_adjust_files_for_processing/studies_split/",s,"/otu.txt"), quote=FALSE, row.names=FALSE, sep = "\t")
  write.table(metatemp, file = paste0("02b_adjust_files_for_processing/studies_split/",s,"/meta.txt"), quote=FALSE, row.names=FALSE, sep="\t")
  
}
