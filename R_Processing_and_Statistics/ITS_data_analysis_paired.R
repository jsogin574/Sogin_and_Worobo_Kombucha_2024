#title: "Kombucha-ITS_paired"
#author: "Jonathan Sogin"
#date: "2023"

#Importing libraries
#######################################################
#pre-processing and data handling packages
library("phyloseq"); packageVersion("phyloseq")
library("PERFect"); packageVersion("PERFect")
library("decontam"); packageVersion("decontam")

#visualization packages

library("ggpubr"); packageVersion("ggpubr")
library("ggtext"); packageVersion("ggtext")

#setting seed
addTaskCallback(function(...) {set.seed(02221997);TRUE})
#######################################################

#Custom functions
#######################################################

#Function to create joined table of OTU taxonomy to sample abundances
summary_table <- function(physeq){
  arg_name <- deparse(substitute(physeq))
  table <- as.data.frame(cbind(tax_table(physeq), get_taxa(physeq, sample_names(physeq))))
  var_name <- paste("summary_table", arg_name, sep="_")
  assign(var_name, table, envir = globalenv())
}

#Function to agglomerate taxa, making sure that unreseolved taxa (NA) are not agglomerated.
#This is supposed to be builtin to phyloseq but it was not behaving correctly.
glom_tax <- function(physeq, rank){
  taxa_ranks <- matrix(seq(1,7), ncol=1)
    rownames(taxa_ranks) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rank_number <- taxa_ranks[rank, ]
  tmp_physeq <- physeq
  test <- is.na(tax_table(tmp_physeq)[,rank_number])
  if(!is.na(table(test)["TRUE"])){
    unresolved <- rownames(subset(tax_table(tmp_physeq), test))
    uniqueid_check <- c()
    for(otu in unresolved){
      name_rank=rank_number-1
      if(is.na(tax_table(tmp_physeq)[otu, name_rank])){name_rank=rank_number-2}
      tax_table(tmp_physeq)[otu, rank_number] <- paste0("UR ", tax_table(tmp_physeq)[otu, name_rank], " ", gsub("(.{1})(.{7})", "\\1", otu)) #renaming the taxa with UR higher_tax_rank .{4} from rownames, which are hash ids of the sequences; the chances of the four character code being the same are low
      if(length(setdiff(tax_table(tmp_physeq)[otu, rank_number], uniqueid_check))==0){
        print("some tax hashes were not unique enough to produce a unique unresolved species ID; check hashes")
        stop()
      uniqueid_check <- c(uniqueid_check, tax_table(tmp_physeq)[otu, rank_number])
      }
    }
  }
  tax_glom(tmp_physeq, taxrank=rank, NArm=FALSE)
}
#######################################################


#Importing Data
#######################################################
#file names as environmental variables
biom_file <- "./Qiime_Data/ITS_paired.biom"
sequence_file <- "./Qiime_Data/ITS-dna-sequences_paired.fasta"
metadata_file <- read.csv("./Qiime_Data/metadata.csv", fileEncoding="UTF-8-BOM")

#removing 16S negative controls from metadata_file
metadata <- subset(metadata_file, !is.na(Sample_ID_ITS))
rownames(metadata) <- metadata$Sample_ID_ITS

#converting data to phyloseq objects
fungdata <- import_biom(biom_file, NULL, sequence_file)
fungdata <- merge_phyloseq(fungdata, sample_data(metadata))

#calculating total and relative sugar difference
reported <- sample_data(fungdata)[,"Nutrition.Sugar..g.L."]
sample_data(fungdata)[,"calculated_sugar"] <- sample_data(fungdata)[,"Glucose..g.L."] + sample_data(fungdata)[,"Fructose..g.L."]
sample_data(fungdata)[,"rel_sugar_difference"] <- (sample_data(fungdata)[,"calculated_sugar"] - reported)/reported
#######################################################

#Renaming taxonomy 
#######################################################
#naming taxonomic ranks
colnames(tax_table(fungdata)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Removing extraneous characters from SILVA taxonomy
tax_table(fungdata) = gsub("[[:lower:]]\\_\\_", "", tax_table(fungdata)) 

#Removing underscores formatted as spaces
tax_table(fungdata) = gsub("_{1}", " ", tax_table(fungdata)) 

write.csv(summary_table(fungdata), "Fungdata_raw_paired.csv", quote=F)
#######################################################

#Duplicating a negative control sample to conform with differing extraction batches
#######################################################
duplicate <- merge_phyloseq(otu_table(fungdata)[,"Run2_PcrNeg_1_ITS"], sample_data(fungdata)["Run2_PcrNeg_1_ITS",])
sample_names(duplicate) <- "Run2_PcrNeg_1.2_ITS"
sample_names(fungdata)[50] <- "Run2_PcrNeg_1.1_ITS"
fungdata <- merge_phyloseq(fungdata, duplicate)
sample_data(fungdata)["Run2_PcrNeg_1.1_ITS", "ITS_Sequencing_Run"] <- "2.1"
sample_data(fungdata)["Run2_PcrNeg_1.2_ITS", "ITS_Sequencing_Run"] <- "2.2"
#######################################################

#Filtering for fungi
#######################################################
#filtering features
fungdataFilt <- subset_taxa(fungdata,
    (Kingdom == "Fungi") & #retaining features classified as fungi
    (!is.na(Phylum)) #excluding features unclassified at Phylum level
    )

#table summarizing number of filtered reads
sum_table_filt <- cbind(
    sample_sums(fungdata),
    sample_sums(subset_taxa(fungdata, is.na(Phylum))),
    sample_sums(fungdataFilt)
    )
colnames(sum_table_filt) <- c("total unfiltered","unidentified at phylum", "total passing filters")
#######################################################

#Filtering based on decontamination with fix for taxa believed to be inaccurately classified
#######################################################
#separating extraction and PCR controls for processing
Samples <- subset_samples(fungdataFilt, TRUE.sample.or.CONTROL=="TRUE SAMPLE")
Ext_controls <- subset_samples(fungdataFilt, Sample_ID_ITS=="Run2_ExtNeg_1_ITS" | Sample_ID_ITS=="Run2_ExtNeg_2_ITS")
PCR_controls <- subset_samples(fungdataFilt, Sample_ID_ITS=="Run2_PcrNeg_1_ITS") #this includes the single negative control duplicated

#running decontam with extraction controls
decontam_Ext <- merge_phyloseq(Samples, Ext_controls)
sample_data(decontam_Ext)$is.neg <- sample_data(decontam_Ext)$TRUE.sample.or.CONTROL == "CONTROL SAMPLE"
contam_Ext <- isContaminant(decontam_Ext, batch="ITS_Sequencing_Run", method="prevalence", neg="is.neg", threshold=0.1) #different threshold chosen after looking at data
contam_table_Ext <- as.data.frame(cbind(tax_table(decontam_Ext), contam_Ext, get_taxa(decontam_Ext, sample_names(decontam_Ext))))
#exporting a summary table
contam_table_Ext <- contam_table_Ext[with(contam_table_Ext, order(-contaminant, p)),]
#adding some metadata to the table  
  sample_sums_Ext <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(Ext_controls))))
    colnames(sample_sums_Ext) <- c(colnames(contam_table_Ext)[1:13], names(sample_sums(Samples)), names(sample_sums(Ext_controls)))
    rownames(sample_sums_Ext) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Blinded_Brand", "Blinded_Product", "ITS_Sequencing_Run")])
  Ext_data <- t(sample_data(Ext_controls)[, c("Blinded_Brand", "Blinded_Product", "ITS_Sequencing_Run")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_Ext)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, Ext_data)
contam_table_Ext <- rbind(samples_data, sample_sums_Ext, contam_table_Ext)
write.csv(contam_table_Ext, "ITS_contams_Ext_paired.csv", quote = F)

#running decontam with PCR controls
decontam_PCR <- merge_phyloseq(Samples, PCR_controls)
sample_data(decontam_PCR)$is.neg <- sample_data(decontam_PCR)$TRUE.sample.or.CONTROL == "CONTROL SAMPLE"
contam_PCR <- isContaminant(decontam_PCR, batch="ITS_Sequencing_Run", method="prevalence", neg="is.neg", threshold=0.1) #different threshold chosen after looking at data
contam_table_PCR <- as.data.frame(cbind(tax_table(decontam_PCR), contam_PCR, get_taxa(decontam_PCR, sample_names(decontam_PCR))))
#exporting a summary table
contam_table_PCR <- contam_table_PCR[with(contam_table_PCR, order(-contaminant, p)),]
  #adding some metadata to the table  
  sample_sums_PCR <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(PCR_controls))))
    colnames(sample_sums_PCR) <- c(colnames(contam_table_PCR)[1:13], names(sample_sums(Samples)), names(sample_sums(PCR_controls)))
    rownames(sample_sums_PCR) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Blinded_Brand", "Blinded_Product", "ITS_Sequencing_Run")])
  PCR_data <- t(sample_data(PCR_controls)[, c("Blinded_Brand", "Blinded_Product", "ITS_Sequencing_Run")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_PCR)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, PCR_data)
contam_table_PCR <- rbind(samples_data, sample_sums_PCR, contam_table_PCR)
write.csv(contam_table_PCR, "ITS_contams_PCR_paired.csv", quote = F)

#pooling contaminants from Extraction and PCR and removing OTUs identified as contaminants
contams <- union(row.names(subset(contam_table_Ext, contaminant==T)), row.names(subset(contam_table_PCR, contaminant==T)))
fungdata_decontam <- prune_taxa(setdiff(taxa_names(Samples), contams), fungdataFilt)

#summary table describing reads retained through pipeline
decontam_counts <- as.data.frame(sample_sums(fungdata_decontam))
sum_table <- cbind(sum_table_filt, as.data.frame(sample_sums(fungdata_decontam)))
colnames(sum_table) <- c("total unfiltered", "unidentified at phylum", "passing filters", "contaminants removed/clean")
sum_table

fungdata_decontam <- subset_samples(fungdata_decontam, TRUE.sample.or.CONTROL=="TRUE SAMPLE")
summary(sample_sums(fungdata_decontam))


#removing taxa with zero OTUs because negative controls are not included
fungdata_decontam <- prune_taxa(taxa_sums(fungdata_decontam)!=0, fungdata_decontam)
summary_table(fungdata_decontam)
write.csv(summary_table(fungdata_decontam), "Fungdata_decontamed_table_paired.csv", quote=F)
write.csv(sum_table, "Fungdata_filtering_counts_paired.csv", quote=F)
#######################################################


#Plotting sample relative abundance as is
#######################################################

samples_plotting <- fungdata_decontam

#aglommerating data to Species level for plotting
samples_plottingGlomSpecies <- glom_tax(samples_plotting, "Species")
samples_plottingGlomSpeciesNorm <- transform_sample_counts(samples_plottingGlomSpecies, function(x) x / sum(x))

plotting <- samples_plottingGlomSpeciesNorm

sample_data(plotting)$Blinded_Product <- as.factor(sample_data(plotting)$Blinded_Product)
sample_data(plotting)$Blinded_Brand <- as.factor(sample_data(plotting)$Blinded_Brand)
sample_data(plotting)$Sample_ID <- as.factor(sample_data(plotting)$Sample_ID)
sample_data(plotting) <- sample_data(plotting)[,c("Sample_ID", "Blinded_Product", "Blinded_Brand")]
pruned = filter_taxa(plotting, function(x) mean(x) > 0.01, T)

otu_other = matrix(1-sample_sums(pruned), nrow=1)
  colnames(otu_other) = sample_names(pruned)
  rownames(otu_other) = "other"
  otu_plot = rbind(otu_table(pruned), otu_other)
tax_other = matrix(rep("other",7), nrow=1)
  rownames(tax_other)="other"
  colnames(tax_other) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_plot = rbind(tax_table(pruned), tax_other)
otu_object = otu_table(otu_plot, taxa_are_row=T)
tax_object = tax_table(tax_plot)
plot_object= phyloseq(otu_object, tax_object, sample_data(pruned))

plottable = psmelt(plot_object)
  plottable$Species[plottable$Species!="other"]=gsub("^", "<i>", plottable$Species[plottable$Species!="other"])
  plottable$Species[plottable$Species!="other"]=gsub("$", "</i>", plottable$Species[plottable$Species!="other"])
  plottable$Species[plottable$Species!="other"] <- gsub("^<i>(U[A-Z]) ", "\\1 <i>", plottable$Species[plottable$Species!="other"])
  plottable$Species[plottable$Species!="other"] <- gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", plottable$Species[plottable$Species!="other"])
  plottable$Blinded_Brand=factor(plottable$Blinded_Brand, levels=c("Brand_1", "Brand_2", "Brand_3", "Brand_4", "Brand_5", "Brand_6"))
  plottable$Blinded_Product=factor(plottable$Blinded_Product, levels=c("Product_1", "Product_2", "Product_3", "Product_4", "Product_5", "Product_6", "Product_7", "Product_8", "Product_9", "Product_10", "Product_11", "Product_12", "Product_13", "Product_14", "Product_15", "Product_16", "Product_17", "Product_18", "Product_19", "Product_20", "Product_21", "Product_22", "Product_23", "Product_24"))
  plottable=plottable[order(plottable$Blinded_Product),]
  plottable$Species=factor(plottable$Species, levels=c(setdiff(plottable$Species, "other"), "other"))

rel_abundance_plot <- ggbarplot(plottable, x="Sample", y="Abundance", fill="Species", xlab=F, ylab=F, palette = "simpsons", width=1, legend="bottom")+
  facet_grid(~Blinded_Brand, switch="x", scales="free_x", space="free_x", labeller=as_labeller(function(x){gsub("_", " ", x)}))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  font("legend.title")+
  theme(text=element_text(family="serif"))+
  theme(legend.text=element_markdown())

ggsave(plot=rel_abundance_plot, filename="Fungal_Relative_Abundance_paired.tiff", width=8, height=3, units="in", dpi="print")
#######################################################