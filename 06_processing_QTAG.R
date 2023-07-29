#!bin/bash Rscript
library(tidyverse)
library(vegan)
library(ape)
library(RColorBrewer)

dir.create("06_processing_QTAG")
############### Processing QTAG output #################
allMeta <- read.delim("04_adjust_files_for_QTAG/all_meta.txt")
taxa <- read.delim("04_adjust_files_for_QTAG/taxonomy.txt")
tre16 <- read.tree("03b_phylogenetic_tree/filtered-tree_16S/exported_16S_tree/tree.nwk")
tre18 <- read.tree("03b_phylogenetic_tree/filtered-tree_18S/exported_18S_tree/tree.nwk")

rare_otus <- read.delim("06_diversity/coverage_based_rarefaction/rarefiedTable_edit.txt")

listStudies16 <- allMeta %>% filter(target_region=="16S") %>% pull(Study) %>% unique()
listStudies18 <- c("chen18")

listStudies <- unique(allMeta$Study)
# Remove Hu18
listStudies <- listStudies[listStudies!="hu18"]
########## both 16S and 18S #############
# Get PA matrix of rarefied table
rare_otu_RA <- as.matrix(rare_otus[,-1])
rare_otu_RA <- apply(rare_otu_RA, MARGIN=2, function(x) x/sum(x))
rownames(rare_otu_RA) <- rare_otus[,1]
rare_otu_RA  <- rare_otu_RA %>% as.data.frame() %>% rownames_to_column(var="TaxaID")
colSums(rare_otu_RA[,-1])
# allDat16 <- data.frame(TaxaID=taxa16$Feature.ID)
allDat <- data.frame(TaxaID="MockID")
for ( s in listStudies) {
  # s = "campbell2013"
  # s = "chen18"
  # Load temp dat
  dat <- read.delim(paste0("05_run_QTAG/",s,"/QTAG_output_quantile.txt")) %>%
    mutate(Study=s)
  # mf <- read.delim(paste0("04_adjust_files_for_QTAG/split_studies/",s,"/",s,"_meta.txt"))
  # get quantile curves to extract HD areas
  quant <- read.delim(paste0("05_run_QTAG/",s,"/QTAG_quantile_curves.txt"))
  SalMax <- apply(quant[,-1],2,function(x) quant[which.max(x),1]) %>%
    as.data.frame() %>% rownames_to_column(var="TaxaID") %>%
    rename(SalMax=".")
  quant_t <- data.frame(t(quant[,-1])) %>%
    rownames_to_column(var="TaxaID")
  colnames(quant_t) <- c("TaxaID",paste0("GradQuant_",quant[,1]))
  # Get OTU table
  otu <- read.delim(paste0("05_run_QTAG/",s,"/OTUTable_RA.txt")) %>% 
    rename(TaxaID = X.OTUID)
  colnames(otu) <- gsub("[.]","-",colnames(otu))
  colnames(otu) <- gsub("Baltic-","Baltic.",colnames(otu))

  # Filter samples out of OTU that were thrown out of rarefaction
  keepNames <- intersect(colnames(otu), colnames(rare_otu_RA))
  otu_rare <- rare_otu_RA[,keepNames]
  # colSums(otu_rare[,-1])
  otu <- otu[,keepNames]
  ## Create collapsed OTU table by salinity
  # Make otu matrix
  otu.mat <- as.matrix(otu[,-1])
  rownames(otu.mat) <- otu[,1]
  # Make otu rare matrix
  otu_rare.mat <- as.matrix(otu_rare[,-1])
  rownames(otu_rare.mat) <- otu_rare[,1]
 # colSums(otu_rare.mat)
  # Make collapsing matrix
  mftemp <- allMeta %>% 
    filter(Study ==s) %>%select(SampleID, Salinity) %>%
    mutate(Salinity=round(Salinity), presence=1) %>%
    pivot_wider(names_from = Salinity, values_from=presence, values_fill = 0 ) 
  mf.mat <- as.matrix(mftemp[,-1])
  rownames(mf.mat) <- pull(mftemp[,1])
  mf.mat <- mf.mat[,order(as.numeric(colnames(mf.mat)))]
  
  # my samples in same order
  # mf.mat <- mf.mat[colnames(otu.mat),]
  mf.mat <- mf.mat[colnames(otu.mat),]
  otu.mat <- otu.mat[, match(rownames(mf.mat),colnames(otu.mat))]
  otu_rare.mat <- otu_rare.mat[, match(rownames(mf.mat),colnames(otu_rare.mat))]
  # colSums(otu_rare.mat)
  # otu.mat %>% colSums()
  otu_col <- t(t(otu.mat %*% mf.mat)/colSums(mf.mat))
  colnames(otu_col) <- paste0("Salinity_",colnames(otu_col))
  otu_col_final <- otu_col %>% as.data.frame() %>% 
    rownames_to_column(var="TaxaID")
  #rarefied version
  otu_rare_col <- t(t(otu_rare.mat %*% mf.mat)/colSums(mf.mat))
  # colSums(otu_rare_col)
  colnames(otu_rare_col) <- paste0("SalinityRare_",colnames(otu_rare_col))
  otu_rare_col_final <- otu_rare_col %>% as.data.frame() %>% 
    rownames_to_column(var="TaxaID") %>%
    mutate(Study=s)
  # change sample names so I can grep them later with "SAMPLE"
  colnames(otu) <- c("TaxaID", paste0("SAMPLE_",colnames(otu)[-1]))
  dat <- dat %>% mutate(Study=s) %>% left_join(SalMax) %>% left_join(quant_t) %>% 
    left_join(otu) %>% left_join(otu_col_final) 
  dat <- dat %>% full_join(otu_rare_col_final)
  # colSums(otu_rare_col_final[,-1])
  # dat %>% select(starts_with("SalinityRare_")) %>% 
  #   colSums(na.rm=TRUE)
  allDat <- full_join(allDat, dat)
  # allMf16 <- rbind(allMf16, mf)
}
# Filter out ones without classification
allDat <- allDat %>% filter(!is.na(class))

## Adjust taxonomy
taxa_adj <- taxa %>% 
  # separate(Taxon, into=c("Kingdom_NA","Phylum_NA","Class_NA","Order_NA","Family_NA","Genus_NA","Species_NA"), sep="; ") %>%
  # rename(TaxaID = Feature.ID) %>%
  # rename(TaxaID = row_names) %>%
  rowwise() %>%
  mutate(Species_NA = ifelse(length(grep("uncultured|bacterium", Species_NA))>0, NA, Species_NA)) %>%
  mutate(Genus_NA = ifelse(length(grep("Clade", Genus_NA))>0 & length(grep("SAR11", Order_NA))>0
                           , paste0("g__"
                                    ,gsub("_clade","",gsub("o__","",Order_NA))
                                    ,gsub("g__","",Genus_NA)) , Genus_NA)) %>%
  ungroup()  %>%
  mutate(Phylum =ifelse(is.na(Phylum_NA), Domain_NA, Phylum_NA)
         , Class =ifelse(is.na(Class_NA), Phylum, Class_NA)
         , Order =ifelse(is.na(Order_NA), Class, Order_NA)
         , Family =ifelse(is.na(Family_NA), Order, Family_NA)
         , Genus =ifelse(is.na(Genus_NA), Family, Genus_NA)
         , Species =ifelse(is.na(Species_NA), Genus, Species_NA)) %>%
  mutate(ScientificName = ifelse(Genus !=Species, gsub("s__","",Species), paste0(Species,"_sp"))
  ) %>%
  filter(TaxaID %in% allDat$TaxaID)

write.table(taxa_adj, "06_processing_QTAG/all_taxonomy.txt", quote=FALSE, sep="\t", row.names = FALSE)
# taxa_adj$ScientificName

## Make data spreadsheet
dupTaxa <- taxa_adj[which(duplicated(taxa_adj$TaxaID)),"TaxaID"] %>% pull()
taxa_adj %>% filter(TaxaID %in% dupTaxa) %>% arrange(TaxaID)

allDat %>%
  left_join(taxa_adj) %>%
  select(Study, Domain_NA) %>% table()
taxa_adj %>% filter(TaxaID %in% allDat$TaxaID)

dat_adj <- allDat %>%
  left_join(taxa_adj) %>%
  mutate(Classification = factor(ifelse(class=="high","Marine"
                                        ,ifelse(class=="intermediate", "Brackish"
                                                , ifelse(class=="low","Fresh", "Unclassified")))
                                 , levels=c("Fresh","Brackish","Marine", "Unclassified"))) %>%
  mutate(aveAbund = (lwrTol+ uprTol)/2) %>%
  # left_join(SalMax) %>%
  arrange(Classification, aveAbund) %>%
  mutate(TaxaID = factor(TaxaID, levels=unique(TaxaID))) %>%
  mutate(Colours = ifelse(Classification == "Marine", "red",
                          ifelse(Classification == "Fresh", "blue",
                                 ifelse(Classification == "Brackish","purple", "grey")))) %>%
  filter(!is.na(Classification))
write.table(dat_adj,file = "06_processing_QTAG/QTAG_output_with_taxa.txt", sep="\t", row.names = FALSE, quote=FALSE)

tre16Filt <- keep.tip(tre16, tip=as.character(dat_adj %>% filter(Domain_NA=="d__Bacteria") %>% pull(TaxaID)))
write.tree(tre16Filt, file="06_processing_QTAG/tree16_filt.tre")
tre18Filt <- keep.tip(tre18, tip=as.character(dat_adj %>% filter(Domain_NA=="d__Eukaryota") %>% pull(TaxaID)))
write.tree(tre18Filt, file="06_processing_QTAG/tree18_filt.tre")
write.table(allMeta, file="06_processing_QTAG/allMeta.txt", quote=FALSE, row.names = FALSE, sep="\t")

## Calculating metrics
dat_adj %>% select(TaxaID, Study) %>%
  group_by(Study) %>% summarise(total_number_ASVs=n())

dat_adj %>% select(Study, one_of(paste0("SAMPLE_",allMeta$SampleID))) %>%
  pivot_longer(-Study, names_to = "SampleID", values_to = "Count") %>%
  drop_na() %>%
  group_by(Study) %>% summarise(maxCount=max(Count), minCount=min(Count))

