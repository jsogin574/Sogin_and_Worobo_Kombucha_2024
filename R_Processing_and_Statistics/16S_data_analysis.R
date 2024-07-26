#title: "Kombucha-16S"
#author: "Jonathan Sogin"
#date: "2024"
#R version 4.2.3
#renv version 0.16.0


#Importing libraries
#######################################################
#pre-processing and data handling packages
library("phyloseq"); packageVersion("phyloseq") # 1.41.1
library("PERFect"); packageVersion("PERFect") #version 1.12.0
library("decontam"); packageVersion("decontam") #version 1.18.0

#visualization packages
library("ggpubr"); packageVersion("ggpubr") #version 0.5.0
library("ggtext"); packageVersion("ggtext") #version 0.1.2
library("ggnewscale"); packageVersion("ggnewscale") #version 0.4.8

#data analysis packages
library("microbiome"); packageVersion("microbiome") #version 1.20.0
library("vegan"); packageVersion("vegan") #version 2.6-6.1
library("GUniFrac"); packageVersion("GUniFrac") #version 1.8

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
biom_file <- "./Qiime_Data/16S.biom"
sequence_file <- "./Qiime_Data/16S-dna-sequences.fasta"
metadata_file <- read.csv("./Qiime_Data/metadata.csv", fileEncoding="UTF-8-BOM")

#removing ITS negative controls from metadata_file
metadata <- subset(metadata_file, !is.na(Sample_ID_16S))
rownames(metadata) <- metadata$Sample_ID_16S

#converting data to phyloseq objects
bacdata <- import_biom(biom_file, NULL, sequence_file)
bacdata <- merge_phyloseq(bacdata, sample_data(metadata))                    

#calculating total and relative sugar difference
reported <- sample_data(bacdata)[,"Nutrition.Sugar..g.L."]
sample_data(bacdata)[,"calculated_sugar"] <- sample_data(bacdata)[,"Glucose..g.L."] + sample_data(bacdata)[,"Fructose..g.L."]
sample_data(bacdata)[,"rel_sugar_difference"] <- (sample_data(bacdata)[,"calculated_sugar"] - reported)/reported
#######################################################

#Renaming taxonomy 
#######################################################
#naming taxonomic ranks
colnames(tax_table(bacdata)) = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Removing extraneous characters from SILVA taxonomy
tax_table(bacdata) = gsub("[[:lower:]]\\_\\_", "", tax_table(bacdata)) 

#Removing underscores formatted as spaces
tax_table(bacdata) = gsub("_{1}", " ", tax_table(bacdata)) 

write.csv(summary_table(bacdata), "Bacdata_raw.csv", quote=F)
#######################################################

#Filtering for bacteria
#######################################################
#filtering features
bacdataFilt <- subset_taxa(bacdata,
    (Kingdom=="Bacteria") & #retaining features classified as bacteria
    (!is.na(Phylum)) & #excluding features unclassified at Phylum level
    (Order!="Chloroplast" | is.na(Order)) & #excluding features classified as chloroplast
    (Family!="Mitochondria" | is.na(Family)) #excluding features classified as mitochondria
    )

#table summarizing number of filtered reads
sum_table_filt <- cbind(
    sample_sums(bacdata),
    sample_sums(subset_taxa(bacdata, Kingdom != "Bacteria")),
    sample_sums(subset_taxa(bacdata, is.na(Phylum))),
    sample_sums(subset_taxa(bacdata, Order == "Chloroplast")),
    sample_sums(subset_taxa(bacdata, Family == "Mitochondria")),
    sample_sums(bacdataFilt)
    )
colnames(sum_table_filt) <- c("total unfiltered", "non-bacteria", "unidentified at phylum", "chloroplast", "mitochondria", "total passing filters")
#######################################################

#Filtering based on decontamination with fix for taxa believed to be inaccurately classified
#######################################################
#separating extraction and PCR controls for processing
Samples <- subset_samples(bacdataFilt, TRUE.sample.or.CONTROL=="TRUE SAMPLE")
Ext_controls <- subset_samples(bacdataFilt, Sample_ID_16S=="Run1_ExtNeg_1_16S" | Sample_ID_16S=="Run2_ExtNeg_2_16S")
PCR_controls <- subset_samples(bacdataFilt, Sample_ID_16S=="Run1_PcrNeg_1_16S" | Sample_ID_16S=="Run2_PcrNeg_1_16S")

#running decontam with extraction controls
decontam_Ext <- merge_phyloseq(Samples, Ext_controls)
sample_data(decontam_Ext)$is.neg <- sample_data(decontam_Ext)$TRUE.sample.or.CONTROL == "CONTROL SAMPLE"
contam_Ext <- isContaminant(decontam_Ext, batch="X16S_Sequencing_Run", method="prevalence", neg="is.neg", threshold=0.45) #different threshold chosen after looking at data
contam_table_Ext <- as.data.frame(cbind(tax_table(decontam_Ext), contam_Ext, get_taxa(decontam_Ext, sample_names(decontam_Ext))))
#exporting a summary table
contam_table_Ext <- contam_table_Ext[with(contam_table_Ext, order(-contaminant, p)),]
#adding some metadata to the table  
  sample_sums_Ext <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(Ext_controls))))
    colnames(sample_sums_Ext) <- c(colnames(contam_table_Ext)[1:13], names(sample_sums(Samples)), names(sample_sums(Ext_controls)))
    rownames(sample_sums_Ext) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Blinded_Brand", "Blinded_Product", "X16S_Sequencing_Run")])
  Ext_data <- t(sample_data(Ext_controls)[, c("Blinded_Brand", "Blinded_Product", "X16S_Sequencing_Run")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_Ext)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, Ext_data)
contam_table_Ext <- rbind(samples_data, sample_sums_Ext, contam_table_Ext)
contam_table_Ext['6002230822167d60320423c9f6025a1e', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
contam_table_Ext['e02a08af27fc7a2ac6847607895ba1c6', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
contam_table_Ext['ca3b957c1734aaf3e8b2af0896250140', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
contam_table_Ext['c8f39f0fa97cb24d81d7b608a058736c', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
write.csv(contam_table_Ext, "16S_contams_Ext.csv", quote = F)

#running decontam with PCR controls
decontam_PCR <- merge_phyloseq(Samples, PCR_controls)
sample_data(decontam_PCR)$is.neg <- sample_data(decontam_PCR)$TRUE.sample.or.CONTROL == "CONTROL SAMPLE"
contam_PCR <- isContaminant(decontam_PCR, batch="X16S_Sequencing_Run", method="prevalence", neg="is.neg", threshold=0.45) #different threshold chosen after looking at data
contam_table_PCR <- as.data.frame(cbind(tax_table(decontam_PCR), contam_PCR, get_taxa(decontam_PCR, sample_names(decontam_PCR))))
#exporting a summary table
contam_table_PCR <- contam_table_PCR[with(contam_table_PCR, order(-contaminant, p)),]
  #adding some metadata to the table  
  sample_sums_PCR <- t(matrix(c(rep(NA, 13), sample_sums(Samples), sample_sums(PCR_controls))))
    colnames(sample_sums_PCR) <- c(colnames(contam_table_PCR)[1:13], names(sample_sums(Samples)), names(sample_sums(PCR_controls)))
    rownames(sample_sums_PCR) <- "Sample_Sums"
  samples_data <- t(sample_data(Samples)[, c("Blinded_Brand", "Blinded_Product", "X16S_Sequencing_Run")])
  PCR_data <- t(sample_data(PCR_controls)[, c("Blinded_Brand", "Blinded_Product", "X16S_Sequencing_Run")])
  blank <- matrix(rep(NA, 13*3), ncol=13)
    colnames(blank) <- colnames(contam_table_PCR)[1:13]
    rownames(blank) <- rownames(samples_data)
  samples_data <- cbind(blank, samples_data, PCR_data)
contam_table_PCR <- rbind(samples_data, sample_sums_PCR, contam_table_PCR)
contam_table_PCR['ca3b957c1734aaf3e8b2af0896250140', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
contam_table_PCR['bfd55440f87b51213f12cb41c9c277c5', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
contam_table_PCR['6002230822167d60320423c9f6025a1e', 'contaminant'] <- F #after looking at data this appeared to be a contaminant
write.csv(contam_table_PCR, "16S_contams_PCR.csv", quote = F)

#pooling contaminants from Extraction and PCR and removing OTUs identified as contaminants
contams <- union(row.names(subset(contam_table_Ext, contaminant==T)), row.names(subset(contam_table_PCR, contaminant==T)))
bacdata_decontam <- prune_taxa(setdiff(taxa_names(Samples), contams), bacdataFilt)

#summary table describing reads retained through pipeline
decontam_counts <- as.data.frame(sample_sums(bacdata_decontam))
sum_table <- cbind(sum_table_filt, as.data.frame(sample_sums(bacdata_decontam)))
colnames(sum_table) <- c("total unfiltered", "non-bacteria", "unidentified at phylum", "chloroplast", "mitochondria", "passing filters", "contaminants removed/clean")
sum_table

bacdata_decontam <- subset_samples(bacdata_decontam, TRUE.sample.or.CONTROL=="TRUE SAMPLE")
summary(sample_sums(bacdata_decontam))

#removing taxa with zero OTUs because negative controls are not included
bacdata_decontam <- prune_taxa(taxa_sums(bacdata_decontam)!=0, bacdata_decontam)
summary_table(bacdata_decontam)
write.csv(summary_table(bacdata_decontam), "Bacdata_decontamed_table.csv", quote=F)
write.csv(sum_table, "Bacdata_filtering_counts.csv", quote=F)
#######################################################


#Plotting sample relative abundance as is
#######################################################
samples_plotting <- bacdata_decontam

#agglomerating data to Genus level for plotting
samples_plottingGlomGenus <- glom_tax(samples_plotting, "Genus")
samples_plottingGlomGenusNorm <- transform_sample_counts(samples_plottingGlomGenus, function(x) x / sum(x))

plotting <- samples_plottingGlomGenusNorm

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

original_plottable = psmelt(plot_object)
  original_plottable$Genus[original_plottable$Genus!="other"]=gsub("^", "<i>", original_plottable$Genus[original_plottable$Genus!="other"])
  original_plottable$Genus[original_plottable$Genus!="other"]=gsub("$", "</i>", original_plottable$Genus[original_plottable$Genus!="other"])
  original_plottable$Blinded_Brand=factor(original_plottable$Blinded_Brand, levels=c("Brand_1", "Brand_2", "Brand_3", "Brand_4", "Brand_5", "Brand_6"))
  original_plottable$Blinded_Product=factor(original_plottable$Blinded_Product, levels=c("Product_1", "Product_2", "Product_3", "Product_4", "Product_5", "Product_6", "Product_7", "Product_8", "Product_9", "Product_10", "Product_11", "Product_12", "Product_13", "Product_14", "Product_15", "Product_16", "Product_17", "Product_18", "Product_19", "Product_20", "Product_21", "Product_22", "Product_23", "Product_24"))
  original_plottable=original_plottable[order(original_plottable$Blinded_Product),]
  original_plottable$Genus=factor(original_plottable$Genus, levels=c(setdiff(original_plottable$Genus, "other"), "other"))

original_rel_abundance_plot <- ggbarplot(original_plottable, x="Sample", y="Abundance", fill="Genus", xlab=F, ylab="<strong>original</strong><br>relative abundance", palette = "simpsons", width=1)+
  facet_grid(~Blinded_Brand, switch="x", scales="free_x", space="free_x", labeller=as_labeller(function(x){gsub("_", " ", x)}))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  font("legend.title")+
  theme(axis.title.y=element_markdown(face="plain", size=15))+
  theme(text=element_text(family="serif"))+
  theme(legend.text=element_markdown())
#######################################################

#Removing added cultures from samples
#######################################################

bacdata_decontam_filt <- subset_taxa(bacdata_decontam,
    (Genus!="Bacillus" | is.na(Genus)) &
    (Genus!="Arthrospira PCC-7345" | is.na(Genus)) &
    (Genus!="Aphanizomenon MDT14a" | is.na(Genus)) &
    (Genus!="Limnohabitans" | is.na(Genus))
    )

sum_table_decontam_filt <- cbind(
    sum_table[sample_names(bacdata_decontam),], 
    sample_sums(subset_taxa(bacdata_decontam, Genus == "Bacillus")),
    sample_sums(subset_taxa(bacdata_decontam, Genus == "Arthrospira PCC-7345")),
    sample_sums(subset_taxa(bacdata_decontam, Genus == "Aphanizomenon MDT14a")),
    sample_sums(subset_taxa(bacdata_decontam, Genus == "Limnohabitans")),
    sample_sums(bacdata_decontam_filt)
    )
colnames(sum_table_decontam_filt) <- c("total unfiltered", "non-bacteria", "unidentified at phylum", "chloroplast", "mitochondria", "total passing filters", "contaminants removed/clean", "bacillus", "arthrospira", "aphanizomenon", "limnohabitans", "passing added culture removal")

write.csv(summary_table(bacdata_decontam_filt), "Bacdata_decontamedfilt_table.csv", quote=F)
write.csv(sum_table_decontam_filt, "Bacdata_filtering_counts_decontamfilt.csv", quote=F)
#######################################################

#Plotting sample relative abundance with added cultures 'removed'
#######################################################

#making colors the same for taxa, even though some taxa are removed
palette="simpsons"
pal_tab <- cbind(levels(original_plottable$Genus), get_palette(palette=palette, k=length(levels(original_plottable$Genus))))
colnames(pal_tab) <- c("Genus", "Hex")
rownames(pal_tab) <- pal_tab[,"Genus"]

samples_plotting <- bacdata_decontam_filt

#agglomerating data to Genus level for plotting
samples_plottingGlomGenus <- glom_tax(samples_plotting, "Genus")
samples_plottingGlomGenusNorm <- transform_sample_counts(samples_plottingGlomGenus, function(x) x / sum(x))

plotting <- samples_plottingGlomGenusNorm

#subsetting same genera as plotted in the original plot
taxa_to_keep <- factor(setdiff(as.vector(tax_table(plot_object)[,"Genus"]), c("other", "Bacillus", "Arthrospira PCC-7345")))
otus_to_keep <- rownames(tax_table(plotting)[,"Genus"][taxa_to_keep])

sample_data(plotting)$Blinded_Product <- as.factor(sample_data(plotting)$Blinded_Product)
sample_data(plotting)$Blinded_Brand <- as.factor(sample_data(plotting)$Blinded_Brand)
sample_data(plotting)$Sample_ID <- as.factor(sample_data(plotting)$Sample_ID)
sample_data(plotting) <- sample_data(plotting)[,c("Sample_ID", "Blinded_Product", "Blinded_Brand")]
pruned = prune_taxa(otus_to_keep, plotting)

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

addedremoved_plottable = psmelt(plot_object)
  addedremoved_plottable$Genus[addedremoved_plottable$Genus!="other"]=gsub("^", "<i>", addedremoved_plottable$Genus[addedremoved_plottable$Genus!="other"])
  addedremoved_plottable$Genus[addedremoved_plottable$Genus!="other"]=gsub("$", "</i>", addedremoved_plottable$Genus[addedremoved_plottable$Genus!="other"])
  addedremoved_plottable$Blinded_Brand=factor(addedremoved_plottable$Blinded_Brand, levels=c("Brand_1", "Brand_2", "Brand_3", "Brand_4", "Brand_5", "Brand_6"))
  addedremoved_plottable$Blinded_Product=factor(addedremoved_plottable$Blinded_Product, levels=c("Product_1", "Product_2", "Product_3", "Product_4", "Product_5", "Product_6", "Product_7", "Product_8", "Product_9", "Product_10", "Product_11", "Product_12", "Product_13", "Product_14", "Product_15", "Product_16", "Product_17", "Product_18", "Product_19", "Product_20", "Product_21", "Product_22", "Product_23", "Product_24"))
  addedremoved_plottable=addedremoved_plottable[order(addedremoved_plottable$Blinded_Product),]
  addedremoved_plottable$Genus=factor(addedremoved_plottable$Genus, levels=c(setdiff(addedremoved_plottable$Genus, "other"), "other"))

addedremoved_rel_abundance_plot <- ggbarplot(addedremoved_plottable, x="Sample", y="Abundance", fill="Genus", xlab=F, ylab="<strong>cultures removed</strong><br>relative abundance", palette=as.vector(pal_tab[levels(addedremoved_plottable$Genus),"Hex"]), width=1)+
  facet_grid(~Blinded_Brand, switch="x", scales="free_x", space="free_x", labeller=as_labeller(function(x){gsub("_", " ", x)}))+
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  font("legend.title")+
  theme(axis.title.y=element_markdown(face="plain", size=15))+
  theme(text=element_text(family="serif"))+
  theme(legend.text=element_markdown())

combined_rel_abundance_plot <- ggarrange(original_rel_abundance_plot, addedremoved_rel_abundance_plot, ncol=1, common.legend = T, legend = "bottom", align="v")
ggsave(plot=combined_rel_abundance_plot, filename="Bacterial_Relative_Abundance.tiff", width=8.5, height=6, units="in", dpi="print", bg="white")
#######################################################

#Data analysis
#subsetting and denoising data
#######################################################
#subsetting data
#Samples K07_1, K07_2, and K41 were removed because they contained a large proportion of reads related to the addition of 'greens' including Arthrospira
#Samples K37 and K39 were removed because they contained very few reads after B. coagulans was removed from the sample
analysis <- prune_samples(setdiff(sample_names(bacdata_decontam_filt), c("Run1_K07_1_16S", "Run2_K07_2_16S", "Run2_K37_16S", "Run2_K39_16S", "Run2_K41_16S")), bacdata_decontam_filt)
  analysis <- filter_taxa(analysis, function(x) sum(x)!=0, prune=T)

#Denoising data
#Using PERFect, which is a statistical means of eliminating taxa that do not contribute to covariance
#this will inherently affect alpha diversity measurements, but will do so with the benefit of greater confidence in preventing artifically high alpha diversity measurements due to sequencing artifacts

#splitting up analysis object, as the otu table will be modified and the separate parts will need to be merged into a new object after denoising
bacdatadecontam_otu <- otu_table(analysis)
bacdatadecontam_tax <- tax_table(analysis)
bacdatadecontam_sample <- sample_data(analysis)
bacdatadecontam_seqs <- refseq(analysis)

#running PERFect
#transposing data to use in PERFect
Counts <- t(bacdatadecontam_otu)
dim(Counts)

res_sim <- PERFect_sim(X = Counts)
  dim(res_sim$filtX)  

simultaneous_filtering_plot <- pvals_Plots(PERFect = res_sim, X = Counts, quantiles = c(0.25, 0.5, 0.8, 0.9), alpha=0.05)
  simultaneous_filtering_plot <- simultaneous_filtering_plot$plot + ggtitle("Simultanenous Filtering")

res_perm <- PERFect_perm(X = Counts, Order = "pvals", pvals_sim = res_sim, algorithm = "full")
  dim(res_perm$filtX)

permutational_filtering_plot <- pvals_Plots(res_perm, Counts)
  permutational_filtering_plot <- permutational_filtering_plot$plot + ggtitle("Full Algorithm")
ggsave(plot=permutational_filtering_plot, filename="Bacterial_Permutational_Filtering.tiff", width=8, height=6, units="in", dpi="print")
  
new_otu <- t(res_perm$filtX)

#remerging data to phyloseq object
bacdata_denoised <- phyloseq(new_otu, bacdatadecontam_tax, bacdatadecontam_sample, bacdatadecontam_seqs)

#table summarizing number of filtered reads
sum_table_denoised <- cbind(sum_table_decontam_filt[sample_names(bacdata_denoised),], data.frame(sample_sums(bacdata_denoised)))
colnames(sum_table_denoised) <- c("total unfiltered", "non-bacteria", "unidentified at phylum", "chloroplast", "mitochondria", "total passing filters", "contaminants removed/clean", "bacillus", "arthrospira", "aphanizomenon", "limnohabitans", "passing added culture removal", "denoised")
sum_table_denoised

write.csv(summary_table(bacdata_denoised), "Bacdata_denoised_table.csv", quote=F)
write.csv(sum_table_denoised, "Bacdata_filtering_counts_denoised.csv", quote=F)
#######################################################


#######################################################
analysis <- bacdata_denoised
#renaming 'uncultured bacterium' Species and 'uncultured' Genus to unique names so they are not collapsed. This was not relevant for previous plotting because the collapsed 'unresolved' category was not present at the given cutoff, meaning none of the individual 'unresolved' taxa would be either.
  test <- tax_table(analysis)[,"Genus"]=="uncultured"
  uncultured <- rownames(subset(tax_table(analysis), test))
  for(otu in uncultured){
    tax_table(analysis)[otu, "Genus"] <- paste0("UC", " ", tax_table(analysis)[otu, "Family"], " ", gsub("(.{1})(.{7})", "\\1", otu)) #renaming the taxa with UR higher_tax_rank .{4} from rownames, which are hash ids of the sequences
  }
  test <- tax_table(analysis)[,"Species"]=="uncultured bacterium"
  uncultured <- rownames(subset(tax_table(analysis), test))
  for(otu in uncultured){
    tax_table(analysis)[otu, "Species"] <- paste0("UC", " ", tax_table(analysis)[otu, "Genus"], " ", gsub("(.{1})(.{7})", "\\1", otu)) #renaming the taxa with UR higher_tax_rank .{4} from rownames, which are hash ids of the sequences
  }
  #previous step makes portions of the name redundant in the species when genus and species were noted as uncultured
  tax_table(analysis)[,"Species"] <- gsub("UC UC", "UC", tax_table(analysis)[,"Species"])
  tax_table(analysis)[,"Species"] <- gsub("( .{4})( .{4})$", "\\1", tax_table(analysis)[,"Species"])

#exporting analysis object for network analysis between bacteria and fungi
saveRDS(analysis, "phyloseq_analysis_bacteria.rds")

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
  scale_x_continuous(limits=c(-1.45, 2.1), expand=c(0,0))+
  scale_y_continuous(limits=c(-1.25,1.35), expand=c(0,0))+
  new_scale_color()+
  geom_point(aes(color=Ethanol...., size=4))+
  scale_color_steps2("% Ethanol (v/v)", low="darkgreen", midpoint=0.5, high="darkred")+
  guides(size="none")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
  theme(legend.title=element_text(face="bold"))+
  theme(text=element_text(family="serif"))+
  #new_scale_color()+
  geom_point(data=ordination_plot_taxa, shape="square")+
  geom_richtext(data=ordination_plot_taxa, label=ordination_plot_taxa$Species, label.size=NA, alpha=0.8, fill=NA, nudge_x=c(rep(0.225,10), 0, 0.225, 0.225), nudge_y=c(rep(0.052,10), 0.1, 0.052, 0.052), family="serif", size=3)
  
#exported as tiff Bacterial_Bray_Ordination
ggsave(plot=ordination_plot, filename="Bacterial_Bray_Ordination.tiff", width=8, height=4, units="in", dpi="print")
#######################################################

#Beta diversity PERMANOVA
#######################################################
#Conducting PERMANOVA analysis to determine effect of Brand and Product variables on overall community structure

#Agglomerating data to Species level to reduce the dimentionality of the data
phyloseq4permanova <- analysis_GlomSpecies_norm

otu <- t(abundances(phyloseq4permanova))
dist <- vegdist(otu, method="bray")
meta <- meta(phyloseq4permanova)

#looking at impact of explanatory variables Brand and Product permutated within sequencing run due to unequal sampling
permanova.explanatory <- adonis(otu ~ Blinded_Brand + Blinded_Product, data=meta, method="bray", permutations=10000, strata=meta$X16S_Sequencing_Run)
permanova.explanatory$call; permanova.explanatory$aov.tab

#looking at association of physiochemical variables with community structure after controlling for sequencing run and brand
dbrda.physiochemical <- dbrda(dist ~ X16S_Sequencing_Run + Blinded_Brand + pH + Ethanol.... + Lactic..g.L. + Acetic..g.L. + calculated_sugar + rel_sugar_difference, data=meta, na.action=na.omit)
dbrda.physiochemical.aov <- anova(dbrda.physiochemical, by = 'margin')
dbrda.physiochemical.aov

permanova.physiochemical <- adonis(otu ~ X16S_Sequencing_Run + Blinded_Brand + Lactic..g.L., data=meta, method="bray", permutations=10000)
permanova.physiochemical$call; permanova.physiochemical$aov.tab

#looking at impact of microbiological variables
dbrda.microbiological <- dbrda(dist ~ X16S_Sequencing_Run + Blinded_Brand + GYC.N..log10CFU.mL. + MRS.N..log10CFU.mL. + APDA..log10CFU.mL., data=meta, na.action=na.omit)
dbrda.microbiological.aov <- anova(dbrda.microbiological, by = 'margin')
dbrda.microbiological.aov


coef_Lactic <- coefficients(permanova.physiochemical)["Lactic..g.L.",]
top_coef_Lactic <- as.data.frame(coef_Lactic[rev(order(abs(coef_Lactic)))[1:10]])
names(top_coef_Lactic) <- "coefficient"
top_coef_Lactic$otu_id <- row.names(top_coef_Lactic)
for(i in 1:nrow(top_coef_Lactic)){
  otuid = top_coef_Lactic$otu_id[i]
  species = tax_table(phyloseq4permanova)[otuid][[7]]
  species_name = paste0("<i>", species, "</i>")
  species_name = gsub("^<i>(U[A-Z]) ", "\\1 <i>", species_name)
  species_name = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", species_name)
  top_coef_Lactic$species_name[i] <- species_name
}
top_coef_Lactic <- top_coef_Lactic[order(top_coef_Lactic$coefficient),]

#plotting the top Lactic species coefficients
plot_coef_Lactic <- ggbarplot(data=top_coef_Lactic, x="species_name", y="coefficient", fill="coefficient", orientation="horiz", legend="none")+
  scale_fill_gradient(low="#FEDD61", high="#9BB8EA", limits=c(-0.15, 0.15))+
  theme(axis.text.y=element_markdown())+
  scale_y_continuous(breaks=seq(-0.15, 0.15, 0.05), limits=c(-0.15, 0.15), labels=c(paste("\U2212", "0.15", sep=""), paste("\U2212", "0.10", sep=""), paste("\U2212", "0.05", sep=""), "0.00", "0.05", "0.10", "0.15"))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(text=element_text(family="serif"))

#exported as tiff Bacterial_Product_PERMANOVA_Coefficients.tiff
ggsave(plot=plot_coef_Lactic, filename="Bacterial_Lactic_PERMANOVA_Coefficients.tiff", width=8, height=3, units="in", dpi="print")
#######################################################

#Amount of lactic acid was associated with observed differences in communities, so plotting the relative abundance of Lactobacillaceae family versus amount of lactic acid
#######################################################
phyloseq4physchem <- analysis
phyloseq4physchem_glomOrder <- glom_tax(phyloseq4physchem, "Order")
phyloseq4physchem_glomOrderNorm <- transform_sample_counts(phyloseq4physchem_glomOrder, function(x) x / sum(x))

phys_chem_plots <- psmelt(phyloseq4physchem_glomOrderNorm)
phys_chem_plots$Blinded_Brand <- gsub("Brand_", "", phys_chem_plots$Blinded_Brand)

#calculating correlation coefficients
cor_data <- subset(phys_chem_plots, Order=="Lactobacillales")
#p-value for kendall correlation test = 5.689E-9
kendall_cor <- cor.test(x=cor_data$Abundance, y=cor_data$Lactic..g.L., method="kendall", exact=F)
  kendall_cor
#r2 value for linear regression = 0.4842
reg_cor <- summary(lm(Lactic..g.L. ~ Abundance, data=cor_data))
  reg_cor

lactic_cor_plot <- ggscatter(data=subset(phys_chem_plots, Order=="Lactobacillales"), x="Abundance", y="Lactic..g.L.", add="reg.line", legend="right", ylab="lactic acid g/L", xlab="<i>Lactobacillales</i> relative abundance")+
  geom_richtext(label=paste("<i>R<sup>2</sup></i> = 0.48, <i>p</i> = 5.7 \U00D7 10<sup>\U2212","9</sup>", sep=""), x=0.25, y=2.25, family="serif", label.color="white")+
  geom_point(aes(color=Blinded_Brand), size=5)+
  scale_color_manual(name="Brand", values=get_palette("simpsons", 6))+
  theme(text=element_text(family="serif"))+
  theme(axis.title.x=element_markdown())

#exported as tiff Bacterial_Product_PERMANOVA_Coefficients.tiff
ggsave(plot=lactic_cor_plot, filename="Bacterial_Lactic_correlation.tiff", width=6, height=4, units="in", dpi="print")
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

ZicoSeq_Brand <- ZicoSeq(meta.dat=meta.dat, feature.dat=feature.dat, grp.name="Blinded_Brand", adj.name="X16S_Sequencing_Run", feature.dat.type="count", prev.filter=prev.filter, mean.abund.filter=mean.abund.filter, max.abund.filter=max.abund.filter, min.prop=min.prop, winsor.end="top", perm.no=perm.no, return.feature.dat=T)

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
write.csv(res_Blinded_Brand, "16S_ZicoSeq_Blinded_Brand.csv", quote=F)

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

ggsave(plot=brand_diff_cor, filename="Bacterial_ZicoSeq_Brand_Plot.tiff", width=6, height=6, units="in", dpi="print")
#######################################################


#A group popped up on the network analysis that was associated with Ginger, so doing zicoseq on that
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

#ZicoSeq filtering parameters, different from Blinded_Brand
  meta.dat=zico_meta_data
  feature.dat=zico_otu
  prev.filter=0.1
  mean.abund.filter=0.0005
  max.abund.filter=0
  min.prop=0
  outlier.pct=0.03
  perm.no=10000

ZicoSeq_Ginger <- ZicoSeq(meta.dat=meta.dat, feature.dat=feature.dat, grp.name="Contains_Ginger", adj.name=c("Blinded_Brand", "ITS_Sequencing_Run"), feature.dat.type="count", prev.filter=prev.filter, mean.abund.filter=mean.abund.filter, max.abund.filter=max.abund.filter, min.prop=min.prop, winsor.end="top", perm.no=perm.no, return.feature.dat=T)

  otus_tested_Ginger <- names(ZicoSeq_Ginger$p.adj.fdr)
  rel_proportion_tested <- sample_sums(prune_taxa(otus_tested_Ginger, phyloseq4ZicoSeq_norm))

res_Ginger <- data.frame(cbind(ZicoSeq_Ginger$p.raw, ZicoSeq_Ginger$p.adj.fdr))
  colnames(res_Ginger) <- c("p.raw", "p.adj.fdr")
  res_Ginger <- res_Ginger[order(res_Ginger[,"p.adj.fdr"]),]
  res_Ginger_tax <- tax_table(phyloseq4ZicoSeq)[rownames(res_Ginger),]
  res_Ginger_mean_T <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Contains_Ginger=="TRUE"))[rownames(res_Ginger)])
  res_Ginger_cv_T <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Contains_Ginger=="TRUE"))[rownames(res_Ginger)], 1, function(x) sd(x)/mean(x))
  res_Ginger_mean_F <- rowMeans(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Contains_Ginger=="FALSE"))[rownames(res_Ginger)])
  res_Ginger_cv_F <- apply(otu_table(subset_samples(phyloseq4ZicoSeq_norm, Contains_Ginger=="FALSE"))[rownames(res_Ginger)], 1, function(x) sd(x)/mean(x))
  res_Ginger <- cbind(res_Ginger, res_Ginger_tax, res_Ginger_mean_T, res_Ginger_cv_T, res_Ginger_mean_F, res_Ginger_cv_F)
  res_Ginger <- data.frame(res_Ginger)
  
write.csv(res_Ginger, "16S_ZicoSeq_Ginger.csv", quote=F)

ginger_cor <- res_Ginger[res_Ginger$p.adj.fdr<0.05, c("res_Ginger_mean_T", "res_Ginger_mean_F")]
colSums(ginger_cor)
  ginger_cor <- ginger_cor[rev(order(rowSums(ginger_cor))),]
  ginger_cor <- -log(ginger_cor)
  ginger_cor[ginger_cor=="Inf"] <- NA
ginger_cor_tax <- tax_table(phyloseq4ZicoSeq)[rownames(ginger_cor),]
ginger_cor_tax["24d821dddfbe031b81a44817e6d6a4c1","Genus"] <- "UR Enterobacteriaceae 2d8e"
ginger_cor_tax["b6e52c74e408bcf8d49afe994fdcf359","Genus"] <- "UR Enterobacteriaceae bed4"
    ginger_cor_tax[,"Species"] <- paste0("<i>", ginger_cor_tax[,"Species"], "</i>")
    ginger_cor_tax[,"Species"] <- gsub("^<i>(U[A-Z]) ", "\\1 <i>", ginger_cor_tax[,"Species"])
    ginger_cor_tax[,"Species"] <- gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", ginger_cor_tax[,"Species"])
ginger_cor_plot_data <- cbind(factor(ginger_cor_tax[,"Species"], levels=rev(ginger_cor_tax[,"Species"])), factor(ginger_cor_tax[,"Genus"]), stack(ginger_cor))
colnames(ginger_cor_plot_data) <- c("Species", "Genus", "abund", "Blinded_Brand")

my_cols <- get_palette(palette="Paired", k=length(unique(ginger_cor_tax[,"Genus"])))
names(my_cols) <- unique(ginger_cor_tax[,"Genus"])
ginger_diff_cor <- ggscatter(data=ginger_cor_plot_data, y="Species", x="Blinded_Brand", legend="right", xlab="Contains Ginger", ylab=F)+
  geom_point(aes(color=abund, size=10))+
  guides(size="none")+
  scale_color_steps2(limits=c(0,14), low="#b30000", high="#0048b3", midpoint=c(7), mid="#f7ec9b", na.value="#f3f3f3", name="-ln(mean)", breaks=c(seq(0,14,2)))+
  theme(axis.text.y=element_markdown(color=my_cols[rev(ginger_cor_tax[,"Genus"])]))+
  theme(text=element_text(family="serif"))+
  scale_x_discrete(labels=c("True", "False"))

ggsave(plot=ginger_diff_cor, filename="Bacterial_ZicoSeq_Ginger_Plot.tiff", width=6, height=6, units="in", dpi="print")
#######################################################