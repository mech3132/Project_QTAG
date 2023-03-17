# library(phyloseq)
library(vegan)
library(MASS)
library(tidyverse)

## Adjust formatting for QTAG ##
dir.create("02a_adjust_files_for_processing")

meta16 <- read.delim("01_processed_data/chen_metadata_16S.txt", as.is = TRUE) %>%
  mutate(Lat = as.numeric(Lat), Long = as.numeric(Long)) 
meta18 <- read.delim("01_processed_data/chen_metadata_18S.txt", as.is = TRUE)
meta <- full_join(meta16, meta18)
otu16 <- read.delim("01_processed_data/chen_otu_16S_unfiltered_t.txt")
otu18 <- read.delim("01_processed_data/chen_otu_18S_unfiltered_t.txt")

otu_temp <- full_join(otu16, otu18)
remove(otu16)
remove(otu18)
taxa <- read.delim("01_processed_data/chen_128_taxonomy_2022-08-26.txt")
# taxa <- read.delim("01_processed_data/INPUT_DATA_taxonomysilva128_filtered-chen-salinity.txt")
# phylodat <- readRDS("01_processed_data/INPUT_DATA_phyloseq_filtered-chen-salinity.rds")

meta_temp <- meta %>% rename(Salinity=salinity, temp_C=temperature, Study=study)
# otu_temp <- otu
# remove(otu)
# meta_temp$study
# Check number of samples in each
meta_temp %>% select(Study) %>% 
  table()

## Fix dates for some metadata
meta_temp <- meta_temp %>%
  mutate(Season = ifelse(Season=="" | is.na(Season), "",Season)) %>%
  mutate(Month = ifelse(Season=="winter", 02,
                        ifelse(Season=="spring", 05,
                               ifelse(Season=="summer",08,
                                      ifelse(Season=="autumn", 10, Month))))
         , Day = ifelse(Season=="winter", 28,
                          ifelse(Season=="spring", 31,
                                 ifelse(Season=="summer",15,
                                        ifelse(Season=="autumn", 15, Day))))) %>%
  mutate(Season = ifelse(Season!="", Season,
                         ifelse(Month %in% c(12,1,2),"winter",
                                ifelse(Month %in% c(3,4,5), "spring",
                                       ifelse(Month %in% c(6,7,8), "summer"
                                              , ifelse(Month %in% c(9,10,11), "autumn",NA)))))) 

meta_temp %>% select(Study, Season) %>% table()

# meta_temp %>%
#   select(row_names, Year, Month, Day, Season) %>%View()

# samptemp <- meta %>% filter(Study=="chen") %>% pull(row_names)
# samptemp[which(samptemp %in% colnames(otu))]
# colnames(otu)[grep("R.jt.4.28", colnames(otu))]

# Adjust column names so none begin with numeric
newOTUcolnames <- gsub("[.]|-","_",gsub("^X","",colnames(otu_temp)))
colnames(otu_temp) <- paste0("X", newOTUcolnames)
meta_temp$row_names <- gsub("[.]|-","_",paste0("X", meta$row_names))


# Verify that all metadata samples are in OTU and vice versa
shared_samples <- intersect(meta_temp$row_names, colnames(otu_temp))
length(shared_samples)
nrow(meta_temp)
ncol(otu_temp)-1

meta_temp <- meta_temp %>% filter(row_names %in% shared_samples) %>%
  rename('#SampleID'=row_names)
otu_temp <- otu_temp %>% select("Xrow_names", one_of(shared_samples)) %>%
  rename('#OTU ID'=Xrow_names)
taxa_temp <- taxa %>% rename(TaxaID=row_names) %>%
  rename_at(vars(c("domain","phylum","class","order","family","genus","species")), ~c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  mutate(Phylum_NA = ifelse(Phylum==Domain, NA, Phylum)
         , Class_NA = ifelse(Class==Phylum, NA, Class)
         , Order_NA = ifelse(Order==Class, NA, Order)
         , Family_NA = ifelse(Family==Order, NA, Family)
         , Genus_NA = ifelse(Genus==Family, NA, Genus)
         , Species_NA = ifelse(Species==Genus, NA, Species)) 


# View(meta_temp)

## Check that Study is not empty
# meta_temp %>% filter(is.na(Study)) %>% View()

## Check whether ASVs are shared between studies
# otu_mat <- as.matrix(otu_temp %>% select(-"#OTU ID"))
# study_mat <- matrix(nrow=nrow(otu_mat),ncol=length(unique(meta_temp$Study)),  dimnames = list(rownames(otu_mat), unique(meta_temp$Study)))
# sm <- meta_temp %>% select("#SampleID", Study) %>% table() %>% as.matrix()
# study_mat <- sm[match(colnames(otu_mat), rownames(sm)),]
# otu_mat %*% study_mat %>% View()
# ASVs ARE shared between studies. Need to separate out by study

## Look at what to keep in metadata file
# str(meta_temp)
varToKeep <- c("Study","Year","Month","Day","Season","Salinity","Lat","Long","temp_C", "FWprimerName","RVprimerName", "target_region","StudySeason")
# meta_temp %>% select(Study, RVprimerName) %>% table(useNA="ifany")
# The missing 18S are 2014
# meta_temp %>% filter(Study=="chen18", is.na(temp_C)) %>%
#   select(`#SampleID`)
## Adding back some 18S data
meta18 <- read.csv("00_msc/manual_metadata_chen.csv") %>%
  mutate(X.SampleID = paste0("X",X.SampleID)) %>%
  rename_with(~paste0(.,"2")) %>%
  rename(`#SampleID`= X.SampleID2, Study=Study2)
# View(meta_temp)
# table(meta_temp %>% select(Study, target_region))
meta_temp2 <-  left_join(meta_temp, meta18) %>% 
  mutate(Year = ifelse(is.na(Year) & !is.na(Year2), Year2, Year)
         , Month = ifelse(is.na(Month) & !is.na(Month2), Month2, Month)
         , Day = ifelse(is.na(Day) & !is.na(Day2), Day2, Day)
         , Lat = ifelse(is.na(Lat) & !is.na(Lat2), Lat2, Lat)
         , Long = ifelse(is.na(Long) & !is.na(Long2), Long2, Long)
  ) %>% rowwise() %>%
  mutate(BactorEuk = target_region) %>% 
  mutate(Studynew = ifelse(Study == "hu2016", ifelse(BactorEuk == "16S", "hu16", "hu18"),
                           ifelse(Study != "chen", Study,
                                  ifelse(target_region=="16S","chen16", "chen18")))) %>%
  ungroup() %>%
  mutate(Study = Studynew) %>%
  separate(`#SampleID`, into=c("sampleid","suffix"), sep="_chen", fill="right", remove=FALSE) %>% 
  select(-c(sampleid, suffix))
# meta_temp2 %>% filter(Study=="chen18") %>% View()
### Collapse samples from chen ####
# Chen18 has duplicated samples. Need to remove the ones that are not chen_18 (has less metadata)
# Load in CoLRep
colRep16 <- read.delim("00_msc/final_chen_metadata_2022-07-07/chen_16s_mappingfile.txt") %>%
  mutate(sampleid = paste0("X",gsub("-|[.]","_",X.SampleID))) %>%
  mutate(ColRep = paste0(ColRep, "_chen16")) %>%
  select(sampleid, ColRep)
colRep18 <- read.delim("00_msc/final_chen_metadata_2022-07-07/chen_18s_mappingfile.txt") %>%
  mutate(sampleid = paste0("X",gsub("-","_",X.SampleID))) %>%
  mutate(ColRep = paste0(ColRep, "_chen18")) %>%
  select(sampleid, ColRep)

colRep18_merged <- meta_temp2 %>% filter(Study%in% c("chen18")) %>%
  select('#SampleID', Year, Month, Day, Salinity, Lat, Long, replicate) %>% 
  separate('#SampleID', into=c("sampleid","suffix"), sep="_S", remove=FALSE) %>%
  left_join(colRep18) %>% 
  select('#SampleID', ColRep)
colRep16_merged <- meta_temp2 %>% filter(Study%in% c("chen16")) %>%
  select('#SampleID', Year, Month, Day, Salinity, Lat, Long, replicate) %>%
  separate('#SampleID', into=c("sampleid","suffix"), sep="_chen|_[0-9]$|(_[0-9][0-9]$)|(_[0-9][0-9][0-9]$)", remove=FALSE) %>%
  rowwise() %>%
  mutate(sampleid = ifelse(length(grep("_2016$", sampleid))>0, `#SampleID`, sampleid)) %>%
  ungroup() %>%
  left_join(colRep16) %>% 
  select('#SampleID', ColRep)

meta_adj <- meta_temp2 %>%
  select(-ColRep) %>% left_join(rbind(colRep18_merged,colRep16_merged)) %>% 
  mutate(ColRep = ifelse(!is.na(sample_name)&sample_name!="", sample_name, ColRep)) %>% # for herleman
  mutate(ColRep = ifelse(is.na(ColRep), `#SampleID`, ColRep)) %>%
  unite(Study,Season,sep="",remove=FALSE, col="StudySeason") %>%
  # mutate(`#SampleID` = '#SampleID') %>%
  select(`#SampleID`, one_of("ColRep",varToKeep)) %>%
  filter(!is.na(Salinity)) %>%
  filter(`#SampleID`!="XExtrCon_R118")
# View(meta_adj)
# How many samples kept
ncol(otu_temp)
nrow(meta_adj)
colnames(otu_temp)[which(!colnames(otu_temp) %in% meta_adj$`#SampleID`)]
meta_adj$`#SampleID`[which(!meta_adj$`#SampleID` %in%colnames(otu_temp))]
# meta_adj %>% filter(Study%in%c("chen18","chen16")) %>% View()
# Those are just the duplicates I removed
# ### Check for consistency of sample + composition; previously we had issues where there may have been mis-labelling
otu_for_plot <- otu_temp[,which(colnames(otu_temp) %in% meta_adj$`#SampleID`)]
otu_for_plot <- t(t(otu_for_plot)/colSums(otu_for_plot, na.rm = TRUE))
asvIDsForOTU <- otu_temp[,1]
# # View(otu_for_plot)
# # Chen16
otu_for_plotchen <- otu_for_plot[,meta_adj %>% filter(Study=="chen16") %>% pull(`#SampleID`)]
rownames(otu_for_plotchen) <- otu_temp$`#OTU ID`
otu_for_plotchen <- otu_for_plotchen[rowSums(otu_for_plotchen, na.rm = TRUE)>0,]
# Beta diversity map
bc_chen <- vegdist(t(otu_for_plotchen), method="bray")
nmds_bc_chen <- isoMDS(bc_chen, k=2)
gg_chen16 <- nmds_bc_chen$points %>%as.data.frame() %>%
  rownames_to_column(var="#SampleID") %>%
  left_join(meta_adj) %>%
  ggplot() +
  geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity)))+
  geom_text(aes(x=V1, y=V2, label=`#SampleID`), cex=1)+
  scale_color_gradient(low="blue", high="red")
gg_chen16
ggsave("02a_adjust_files_for_processing/chen16_nmds.png"
       ,gg_chen16, height=6, width=8)
# 18S
otu_for_plotchen <- otu_for_plot[,meta_adj %>% filter(Study=="chen18") %>% pull(`#SampleID`)]
rownames(otu_for_plotchen) <- otu_temp$`#OTU ID`
otu_for_plotchen[is.na(otu_for_plotchen)] <- 0
otu_for_plotchen <- otu_for_plotchen[rowSums(otu_for_plotchen, na.rm=TRUE)>0,]
# otu_for_plotchen[1,]
# Beta diversity map
bc_chen <- vegdist(t(otu_for_plotchen), method="bray")
nmds_bc_chen <- isoMDS(bc_chen, k=2)
gg_chen18 <- nmds_bc_chen$points %>%as.data.frame() %>%
  rownames_to_column(var="#SampleID") %>%
  left_join(meta_adj) %>%
  ggplot() +
  geom_point(aes(x=V1, y=V2, col=as.numeric(Salinity)))+
  geom_text(aes(x=V1, y=V2, label=`#SampleID`), cex=1)+
  scale_color_gradient(low="blue", high="red")
gg_chen18
ggsave("02a_adjust_files_for_processing/chen18_nmds.png"
       ,gg_chen18, height=6, width=8)
### TO REMOVE:
chen16_toremove <- c("XE_BA_2016_3_chen16"
,"XE_TB_2016_2_chen16"
,"XE_BL_2016_1_chen16"
,"XR_wp_5_26_r2_111")
# 
# # otu_for_plot
# otu_for_plot_long <- cbind(OTUID = otu_temp$`#OTU ID`, otu_for_plot) %>% as.data.frame() %>%
#   pivot_longer(-OTUID, names_to="#SampleID", values_to="RA") %>%
#   mutate(RA=as.numeric(RA)) %>%
#   filter(RA!=0) %>%
#   left_join(meta_adj %>% select(`#SampleID`, Salinity, Study))
# otu_for_plot_long %>%
#   filter(Study=="chen16") %>%
#   arrange(as.numeric(Salinity)) %>%
#   mutate(`#SampleID` = factor(`#SampleID`, levels=unique(`#SampleID`))) %>%
#   mutate(Salinity2 = factor(Salinity, levels=sort(as.numeric(unique(Salinity))))) %>%
#   # filter(RA>0.1) %>%
#   # filter(`#SampleID`=="XR_mt_6_23_r_2_26") %>%
#   ggplot() + geom_bar(aes(x=Salinity2, group=`#SampleID`,y=RA, fill=OTUID)
#                       , position="dodge"
#                       , stat="identity", show.legend=FALSE) +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
#   facet_grid(.~Study)

#### INTERMEDIATE SAVE TO REDUCE MEMORY LOAD 
meta_intermediate <- meta_adj %>% 
  mutate(ColRep = gsub("[-]","_",ColRep)) %>%
  filter(!`#SampleID` %in% chen16_toremove)

# Remove singletons and zeros to reduce memory load
any(duplicated(colnames(otu_temp)))
# Remove samples not in metadata file
toRemove <- colnames(otu_temp[,-1])[which(!colnames(otu_temp[,-1]) %in% meta_intermediate$`#SampleID`)]
otu_temp <- otu_temp[,-match(toRemove, colnames(otu_temp))]
otu_temp[is.na(otu_temp)] <- 0
otu_temp <- otu_temp[which(rowSums(otu_temp[,-1])>1),]

save(meta_intermediate,file="02a_adjust_files_for_processing/meta.RData")
save(taxa_temp, file="02a_adjust_files_for_processing/taxa.RData")
save(otu_temp, file="02a_adjust_files_for_processing/otu.RData")


