#title: "Kombucha-ITS_single"
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
library("ggnewscale"); packageVersion("ggnewscale")

#data analysis packages
library("microbiome"); packageVersion("microbiome")
library("vegan"); packageVersion("vegan")
library("GUniFrac"); packageVersion("GUniFrac")
library("Maaslin2"); packageVersion("Maaslin2")
library("SpiecEasi"); packageVersion("SpiecEasi")

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
biom_file <- "./Qiime_Data/ITS_single.biom"
sequence_file <- "./Qiime_Data/ITS-dna-sequences_single.fasta"
metadata_file <- read.csv("./Qiime_Data/metadata.csv", fileEncoding="UTF-8-BOM")

#removing ITS negative controls from metadata_file
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

write.csv(summary_table(fungdata), "Fungdata_raw_single.csv", quote=F)
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
write.csv(contam_table_Ext, "ITS_contams_Ext_single.csv", quote = F)

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
write.csv(contam_table_PCR, "ITS_contams_PCR_single.csv", quote = F)

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
write.csv(summary_table(fungdata_decontam), "Fungdata_decontamed_table_single.csv", quote=F)
write.csv(sum_table, "Fungdata_filtering_counts_single.csv", quote=F)
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

ggsave(plot=rel_abundance_plot, filename="Fungal_Relative_Abundance_single.tiff", width=8, height=3, units="in", dpi="print")
#######################################################


#Data analysis
#subsetting and denoising data
#######################################################
#subsetting samples
#Samples K07_1, K07_2, and K41 were removed because they contained 'greens', and because they were excluded for 16S analysis
#Samples K37 and K39 were retained
analysis <- prune_samples(setdiff(sample_names(fungdata_decontam), c("Run2_K07_1_ITS", "Run2_K07_2_ITS", "Run2_K41_ITS")), fungdata_decontam)
  analysis <- filter_taxa(analysis, function(x) sum(x)!=0, prune=T)
  
#splitting up analysis object, as the otu table will be modified and the separate parts will need to be merged into a new object after denoising
fungdatadecontam_otu <- otu_table(analysis)
fungdatadecontam_tax <- tax_table(analysis)
fungdatadecontam_sample <- sample_data(analysis)
fungdatadecontam_seqs <- refseq(analysis)

#running PERFseq
#transposing data to use in PERFseq
Counts <- t(fungdatadecontam_otu)
dim(Counts)

res_sim <- PERFect_sim(X = Counts)
  dim(res_sim$filtX)  

simultaneous_filtering_plot <- pvals_Plots(PERFect = res_sim, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
  simultaneous_filtering_plot <- simultaneous_filtering_plot$plot + ggtitle("Simultanenous Filtering")

res_perm <- PERFect_perm(X = Counts, Order = "pvals", pvals_sim = res_sim, algorithm = "full")
  dim(res_perm$filtX)

permutational_filtering_plot <- pvals_Plots(res_perm, Counts)
  permutational_filtering_plot <- permutational_filtering_plot$plot + ggtitle("Full Algorithm")
ggsave(plot=permutational_filtering_plot, filename="Fungal_Permutational_Filtering.tiff", width=8, height=6, units="in", dpi="print")
  
new_otu <- t(res_perm$filtX)

#remerging data to phyloseq object
fungdata_denoised <- phyloseq(new_otu, fungdatadecontam_tax, fungdatadecontam_sample, fungdatadecontam_seqs)

#table summarizing number of filtered reads
sum_table_denoised <- cbind(sum_table[sample_names(fungdata_denoised),], data.frame(sample_sums(fungdata_denoised)))
colnames(sum_table_denoised) <- c("total unfiltered", "unidentified at phylum", "passing filters", "contaminants removed/clean", "denoised")
sum_table_denoised

write.csv(summary_table(fungdata_denoised), "Fungdata_denoised_table.csv", quote=F)
write.csv(sum_table_denoised, "Fungdata_filtering_counts_denoised.csv", quote=F)
#######################################################


#######################################################
analysis <- fungdata_denoised

#Renaming species denoted as Genus sp to NA so that custom naming occuring within glom_tax function will assign a unique code to it
tax_table(analysis)[,"Species"] <- gsub("^[A-Za-z]* sp", NA, tax_table(analysis)[,"Species"])
#Renaming two OTUs that were unidentified at Genus and Species, as the glom_tax function will not properly name them
tax_table(analysis)["010ea95f70031e1f24d47d83842658c8", c("Family", "Genus", "Species")] <- "UR Saccharomycetales 0728"
tax_table(analysis)["278398e4874290b21093b386785a2957", c("Family", "Genus", "Species")] <- "UR Saccharomycetales 2817"
tax_table(analysis)["3621e3dbf9c25eaa08b32a94946ff3be", c("Species")] <- "UR Pichiaceae 3f09"

#exporting analysis object for network analysis between bacteria and fungi
saveRDS(analysis, "phyloseq_analysis_fungi.rds")

#Rarefying data for analysis
#rarefying data for alpha and beta diversity analyses
min_analysis <- summary(sample_sums(analysis))['Min.']
analysis_rare <- rarefy_even_depth(analysis, rngseed=02221997, sample.size=min_analysis, replace=F)
analysis_rare_GlomSpecies <- glom_tax(analysis_rare, "Species")
analysis_rare_GlomSpecies_norm <- transform_sample_counts(analysis_rare_GlomSpecies, function(x) x / sum(x))

#normalizing data and agglomerating to species for differential abundance analysis
analysis_GlomSpecies <- glom_tax(analysis, "Species")
analysis_GlomSpecies_norm <-transform_sample_counts(analysis_GlomSpecies, function(x) x / sum(x))

#######################################################

#Beta diversity ordination
#######################################################
#Plotting beta diversity ordination
phyloseq4ordination <- analysis_GlomSpecies_norm
ordination_NMDS <- ordinate(phyloseq4ordination, "NMDS", "bray")
ordination_plot_data <- plot_ordination(phyloseq4ordination, ordination_NMDS, type="biplot", justDF=T)
ordination_plot_data$Blinded_Brand <- gsub("Brand_", "", ordination_plot_data$Blinded_Brand)
ordination_plot_data$Blinded_Brand=factor(ordination_plot_data$Blinded_Brand, levels=seq(1,6))
ordination_plot_samples <- subset(ordination_plot_data, id.type=="Samples")
ordination_plot_taxa <- subset(ordination_plot_data, id.type=="Taxa")
  taxa_mean_abundance <- rev(sort(taxa_sums(phyloseq4ordination)/length(sample_names(phyloseq4ordination))))
  taxa_to_plot <- tax_table(phyloseq4ordination)[names(taxa_mean_abundance[taxa_mean_abundance>0.01])]
  ordination_plot_taxa <- ordination_plot_taxa[is.element(ordination_plot_taxa[,"Species"], taxa_to_plot[,"Species"]),]
  #plotting taxa that have greater than 0.01 mean abundance onto plot
  ordination_plot_taxa$Species=gsub("^", "<i>", ordination_plot_taxa$Species)
  ordination_plot_taxa$Species=gsub("$", "</i>", ordination_plot_taxa$Species)
  ordination_plot_taxa$Species = gsub("^<i>(U[A-Z]) ", "\\1 <i>", ordination_plot_taxa$Species)
  ordination_plot_taxa$Species = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", ordination_plot_taxa$Species)
ordination_plot <- ggscatter(data=ordination_plot_samples, x="NMDS1", y="NMDS2", color="Blinded_Brand", legend.title="Brand", palette="simpsons", legend="right", point=F, ellipse=T, ellipse.alpha=0.25, ellipse.type="convex")+
  scale_x_continuous(limits=c(-2.25, 1.7), expand=c(0,0))+
  #scale_y_continuous(limits=c(-1.8,1.3), expand=c(0,0))+
  new_scale_color()+
  geom_point(aes(color=Ethanol...., size=4))+
  scale_color_steps2("% Ethanol (v/v)", low="darkgreen", midpoint=0.5, high="darkred")+
  guides(size="none")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
  theme(legend.title=element_text(face="bold"))+
  theme(text=element_text(family="serif"))+
  new_scale_color()+
  geom_point(data=ordination_plot_taxa, shape="square")+
  geom_richtext(data=ordination_plot_taxa, label=ordination_plot_taxa$Species, label.size=NA, alpha=0.8, fill=NA, nudge_x=c(0.205, 0.235, rep(0.225,5)), nudge_y=c(0.09, 0.18, rep(0.09, 5)), family="serif", size=3)
  
#exported as tiff Fungal_Bray_Ordination
ggsave(plot=ordination_plot, filename="Fungal_Bray_Ordination.tiff", width=8, height=4, units="in", dpi="print")
#######################################################

#Beta diversity PERMANOVA
#######################################################
#Conducting PERMANOVA analysis to determine effect of Product and Time_point variables on overall community structure

#Agglomerating data to Species level to reduce the dimentionality of the data
phyloseq4permanova <- analysis_GlomSpecies_norm

otu <- t(abundances(phyloseq4permanova))
dist <- vegdist(otu, method="bray")
meta <- meta(phyloseq4permanova)

#looking at impact of explanatory variables Brand and Product permutated within sequencing run due to unequal sampling
permanova.explanatory <- adonis(otu ~ Blinded_Brand + Blinded_Product, data=meta, method="bray", permutations=10000, strata=meta$ITS_Sequencing_Run)
permanova.explanatory$call; permanova.explanatory$aov.tab

#looking at association of physiochemical variables with community structure after controlling for sequencing run and brand
dbrda.physiochemical <- dbrda(dist ~ ITS_Sequencing_Run + Blinded_Brand + pH + Ethanol.... + Lactic..g.L. + Acetic..g.L. + calculated_sugar + rel_sugar_difference, data=meta, na.action=na.omit)
dbrda.physiochemical.aov <- anova(dbrda.physiochemical, by = 'margin')
dbrda.physiochemical.aov

#looking at impact of microbiological variables
dbrda.microbiological <- dbrda(dist ~ ITS_Sequencing_Run + Blinded_Brand + GYC.N..log10CFU.mL. + MRS.N..log10CFU.mL. + APDA..log10CFU.mL., data=meta, na.action=na.omit)
dbrda.microbiological.aov <- anova(dbrda.microbiological, by = 'margin')
dbrda.microbiological.aov

#checking condition for equal variance around centroid
permdisp2_Sequencing_Run <- betadisper(dist, meta$ITS_Sequencing_Run)
permdisp2_Sequencing_Run.aov <- anova(permdisp2_Sequencing_Run)
permdisp2_Sequencing_Run$call; permdisp2_Sequencing_Run.aov

permdisp2_Brand <- betadisper(dist, meta$Blinded_Brand)
permdisp2_Brand.aov <- anova(permdisp2_Brand)
permdisp2_Brand$call; permdisp2_Brand.aov
#######################################################


#Using MaAsLin2 to look at associations between chemical and microbiological data and the observed community
#######################################################
#setting a filtering cutoff of mean abundance greater than 0.0005 and prevalence greater than 0.25
phyloseq4maaslin_tmp <- analysis_GlomSpecies
  phyloseq4maaslin_tmp_norm <- analysis_GlomSpecies_norm
  otus_to_keep <- rownames(tax_table(filter_taxa(phyloseq4maaslin_tmp_norm, function(x) {mean(x) > 0.0005 & sum(x>0) > length(sample_sums(phyloseq4maaslin_tmp_norm))*0.25}, T)))
phyloseq4maaslin <- prune_taxa(otus_to_keep, phyloseq4maaslin_tmp)
sample_sums(phyloseq4maaslin)/sample_sums(phyloseq4maaslin_tmp)

otu <- t(otu_table(phyloseq4maaslin))
meta <- data.frame(sample_data(phyloseq4maaslin)[,c("ITS_Sequencing_Run", "Blinded_Brand", "pH", "Ethanol....", "Lactic..g.L.", "Acetic..g.L.", "Glucose..g.L.", "Nutrition.Sugar..g.L.", "Fructose..g.L.", "calculated_sugar", "rel_sugar_difference", "GYC.N..log10CFU.mL.", "APDA..log10CFU.mL.", "MRS.N..log10CFU.mL.")])

model_chem <- Maaslin2(input_data=otu, input_metadata=meta, output="Fungi_MaaslinChemResults", fixed_effects=c("pH", "Ethanol....", "Lactic..g.L.", "Acetic..g.L.", "calculated_sugar", "rel_sugar_difference"), random_effects=c("ITS_Sequencing_Run", "Blinded_Brand"))

Maaslin_chem_results <- read.delim(file="./Fungi_MaaslinChemResults/all_results.tsv")
Maaslin_chem_results$metadata <- factor(Maaslin_chem_results$metadata, levels=rev(c("Ethanol....", "pH", "Lactic..g.L.", "Acetic..g.L.", "calculated_sugar", "rel_sugar_difference")))

#'significance' is derived from the Maaslin q value, and is -log10(q)*sign(coef)
Maaslin_chem_results$sig <- (-log10(Maaslin_chem_results$qval)*sign(Maaslin_chem_results$coef))
  Maaslin_chem_results$feature <- gsub("X(.{32})", "\\1", Maaslin_chem_results$feature)
  Maaslin_chem_results$species <- tax_table(phyloseq4maaslin)[Maaslin_chem_results$feature, "Species"]
    Maaslin_chem_results$species <- paste0("<i>", Maaslin_chem_results$species, "</i>")
    Maaslin_chem_results$species <- gsub("^<i>(U[A-Z]) ", "\\1 <i>", Maaslin_chem_results$species)
    Maaslin_chem_results$species <- gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", Maaslin_chem_results$species)

maaslin_chem_cor <- ggscatter(data=Maaslin_chem_results[abs(Maaslin_chem_results$sig)>0.6,], y="metadata", x="species", legend="right", xlab=F, ylab=F)+
  rotate_x_text(angle=20)+
  geom_point(aes(color=sig), size=10)+
  scale_color_steps2(limit = c(-2.5,2.5), low = "blue", high =  "red", mid = "white", midpoint = 0, name="Significance")+
  theme(axis.text.x=element_markdown())+
  theme(text=element_text(family="serif"))+
  scale_y_discrete(labels=rev(c("lactic acid")))
    
model_micro <- Maaslin2(input_data=otu, input_metadata=meta, output="Fungi_MaaslinMicroResults", fixed_effects=c("GYC.N..log10CFU.mL.", "APDA..log10CFU.mL.", "MRS.N..log10CFU.mL."), random_effects=c("ITS_Sequencing_Run", "Blinded_Brand"))

Maaslin_micro_results <- read.delim(file="./Fungi_MaaslinMicroResults/all_results.tsv")
Maaslin_micro_results$metadata <- factor(Maaslin_micro_results$metadata, levels=rev(c("GYC.N..log10CFU.mL.", "MRS.N..log10CFU.mL.", "APDA..log10CFU.mL.")))

Maaslin_micro_results$sig <- (-log10(Maaslin_micro_results$qval)*sign(Maaslin_micro_results$coef))
  Maaslin_micro_results$feature <- gsub("X(.{32})", "\\1", Maaslin_micro_results$feature)
  Maaslin_micro_results$species <- tax_table(phyloseq4maaslin)[Maaslin_micro_results$feature, "Species"]
    Maaslin_micro_results$species <- paste0("<i>", Maaslin_micro_results$species, "</i>")
    Maaslin_micro_results$species <- gsub("^<i>(U[A-Z]) ", "\\1 <i>", Maaslin_micro_results$species)
    Maaslin_micro_results$species <- gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", Maaslin_micro_results$species)

maaslin_micro_cor <- ggscatter(data=Maaslin_micro_results[abs(Maaslin_micro_results$sig)>0.6,], y="metadata", x="species", legend="right", xlab=F, ylab="")+
  rotate_x_text(angle=20)+
  geom_point(aes(color=sig), size=10)+
  scale_color_steps2(limit = c(-2.5,2.5), low = "blue", high =  "red", mid = "white", midpoint = 0, name="Significance")+
  theme(axis.text.x=element_markdown())+
  theme(text=element_text(family="serif"))+
  scale_y_discrete(labels=rev(c("GYC", "MRS", "APDA")))+
  theme(axis.title.y = element_text(margin = margin(t=0, r=0, b=0, l=60)))

combined_maaslin_plot <- ggarrange(maaslin_chem_cor, maaslin_micro_cor, ncol=1, align="v", common.legend=T, legend="right", heights=c(2,3))+
  theme(panel.background=element_rect(fill = "white"), plot.background=element_rect(fill = "white"), legend.background=element_rect(fill = "white"))

ggsave(plot=combined_maaslin_plot, filename="Fungal_Maaslin.tiff", width=6, height=4.5, units="in", dpi="print")
#######################################################


#ZicoSeq analysis for brand
#######################################################

phyloseq4ZicoSeq <- analysis_GlomSpecies
  phyloseq4ZicoSeq_norm <- analysis_GlomSpecies_norm

zico_meta_data <-data.frame(sample_data(phyloseq4ZicoSeq))
colnames(zico_meta_data) <- colnames(sample_data(phyloseq4ZicoSeq))

feature_vals <- as.integer(otu_table(phyloseq4ZicoSeq))
  dims <- dim(otu_table(phyloseq4ZicoSeq))
zico_otu <- matrix(feature_vals, nrow=dims[1])
  rownames(zico_otu) <- as.character(rownames(otu_table(phyloseq4ZicoSeq)))
  colnames(zico_otu) <- as.character(colnames(otu_table(phyloseq4ZicoSeq)))

zico_meta_data <- data.frame(sample_data(phyloseq4ZicoSeq))

#ZicoSeq filtering parameters
  meta.dat=zico_meta_data
  feature.dat=zico_otu
  prev.filter=0.25
  mean.abund.filter=0.0005
  max.abund.filter=0
  min.prop=0
  outlier.pct=0.03
  perm.no=10000

ZicoSeq_Brand <- ZicoSeq(meta.dat=meta.dat, feature.dat=feature.dat, grp.name="Blinded_Brand", adj.name="ITS_Sequencing_Run", feature.dat.type="count", prev.filter=prev.filter, mean.abund.filter=mean.abund.filter, max.abund.filter=max.abund.filter, min.prop=min.prop, winsor.end="top", perm.no=perm.no, return.feature.dat=T)

  otus_tested_Blinded_Brand <- names(ZicoSeq_Brand$p.adj.fdr)
  rel_proportion_tested <- sample_sums(prune_taxa(otus_tested_Blinded_Brand, phyloseq4ZicoSeq_norm))

res_Blinded_Brand <- data.frame(cbind(ZicoSeq_Brand$p.raw, ZicoSeq_Brand$p.adj.fdr))
  colnames(res_Blinded_Brand) <- c("p.raw", "p.adj.fdr")
  res_Blinded_Brand <- res_Blinded_Brand[order(res_Blinded_Brand[,"p.adj.fdr"]),]
  res_Blinded_Brand_tax <- tax_table(phyloseq4ZicoSeq)[rownames(res_Blinded_Brand),]
  res_Blinded_Brand_mean_1 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_1"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_1 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_1"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand_mean_2 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_2"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_2 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_2"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand_mean_3 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_3"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_3 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_3"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand_mean_4 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_4"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_4 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_4"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand_mean_5 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_5"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_5 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_5"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand_mean_6 <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_6"))[rownames(res_Blinded_Brand)])
  res_Blinded_Brand_cv_6 <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Blinded_Brand=="Brand_6"))[rownames(res_Blinded_Brand)], 1, function(x) sd(x)/mean(x))
  res_Blinded_Brand <- cbind(res_Blinded_Brand, res_Blinded_Brand_tax, res_Blinded_Brand_mean_1, res_Blinded_Brand_cv_1, res_Blinded_Brand_mean_2, res_Blinded_Brand_cv_2, res_Blinded_Brand_mean_3, res_Blinded_Brand_cv_3, res_Blinded_Brand_mean_4, res_Blinded_Brand_cv_4, res_Blinded_Brand_mean_5, res_Blinded_Brand_cv_5, res_Blinded_Brand_mean_6, res_Blinded_Brand_cv_6)
  res_Blinded_Brand <- data.frame(res_Blinded_Brand)
write.csv(res_Blinded_Brand, "ITS_ZicoSeq_Blinded_Brand.csv", quote=F)

brand_cor <- res_Blinded_Brand[res_Blinded_Brand$p.adj.fdr<0.05, c("res_Blinded_Brand_mean_1", "res_Blinded_Brand_mean_2", "res_Blinded_Brand_mean_3", "res_Blinded_Brand_mean_4", "res_Blinded_Brand_mean_5", "res_Blinded_Brand_mean_6")]
colSums(brand_cor)
  brand_cor <- brand_cor[rev(order(rowSums(brand_cor))),]
  brand_cor <- -log(brand_cor)
  brand_cor[brand_cor=="Inf"] <- NA
brand_cor_tax <- tax_table(phyloseq4ZicoSeq)[rownames(brand_cor),]
    brand_cor_tax[,"Species"] <- paste0("<i>", brand_cor_tax[,"Species"], "</i>")
    brand_cor_tax[,"Species"] <- gsub("^<i>(U[A-Z]) ", "\\1 <i>", brand_cor_tax[,"Species"])
    brand_cor_tax[,"Species"] <- gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", brand_cor_tax[,"Species"])
brand_cor_plot_data <- cbind(factor(brand_cor_tax[,"Species"], levels=rev(brand_cor_tax[,"Species"])), factor(brand_cor_tax[,"Genus"]), stack(brand_cor))
colnames(brand_cor_plot_data) <- c("Species", "Genus", "abund", "Blinded_Brand")

my_cols <- get_palette(palette="Paired", k=length(unique(brand_cor_tax[,"Genus"])))
names(my_cols) <- unique(brand_cor_tax[,"Genus"])
brand_diff_cor <- ggscatter(data=brand_cor_plot_data, y="Species", x="Blinded_Brand", legend="right", xlab="Brand", ylab=F)+
  geom_point(aes(color=abund, size=10))+
  guides(size="none")+
  scale_color_steps2(limits=c(0,14), low="#b30000", high="#0048b3", midpoint=c(7), mid="#f7ec9b", na.value="#f3f3f3", name="-ln(mean)", breaks=c(seq(0,14,2)))+
  theme(axis.text.y=element_markdown(color=my_cols[rev(brand_cor_tax[,"Genus"])]))+
  theme(text=element_text(family="serif"))+
  scale_x_discrete(labels=seq(1:6))

ggsave(plot=brand_diff_cor, filename="Fungal_ZicoSeq_Brand_Plot.tiff", width=6, height=6, units="in", dpi="print")
#######################################################