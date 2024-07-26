#title: "Kombucha-Bac_Fung_Correlation"
#author: "Jonathan Sogin"
#date: "2024"
#R version 4.2.3
#renv version 0.16.0


#Importing libraries
#######################################################
#pre-processing and data handling packages
library("phyloseq"); packageVersion("phyloseq") # 1.41.1

#visualization packages
library("ggpubr"); packageVersion("ggpubr") #version 0.5.0
library("ggtext"); packageVersion("ggtext") #version 0.1.2
library("igraph"); packageVersion("igraph") #version 1.6.0
library("ggnetwork"); packageVersion("ggnetwork") #version 0.5.13

#data analysis packages
library("SpiecEasi"); packageVersion("SpiecEasi") #version 1.1.3

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

#function from phylosmith package
set_sample_order <- function(physeq, treatment=NULL){
  metadata <- as(physeq@sam_data, "data.frame")
  metadata <- data.table::data.table(samples = rownames(metadata), metadata)
  if (is.null(treatment)) treatment <- "samples"
  data.table::setkeyv(metadata, treatment)
  phyloseq::otu_table(physeq) <-
    physeq@otu_table[, metadata$samples]
  phyloseq::sample_data(physeq) <- data.frame(metadata, row.names = 1)
  return(physeq)
}

#######################################################


#Importing Data
#######################################################
#these files are output from the 16S and ITS_single scripts; these are required
bacdata <- readRDS("phyloseq_analysis_bacteria.rds")
fungdata <- readRDS("phyloseq_analysis_fungi.rds")

#agglomerating to species level
bacdata_glomSpecies <- glom_tax(bacdata, "Species")
fungdata_glomSpecies <- glom_tax(fungdata, "Species")

#removing taxa with mean abundance less than 0.0005
bac_tmp <- bacdata_glomSpecies
bac_tmp_norm <- transform_sample_counts(bac_tmp, function(x) x / sum(x))
  bac_otus_to_keep <- rownames(tax_table(filter_taxa(bac_tmp_norm, function(x) mean(x) > 0.0005, T)))
bac <- prune_taxa(bac_otus_to_keep, bac_tmp)
sample_sums(bac)/sample_sums(bac_tmp)

fung_tmp <- fungdata_glomSpecies
fung_tmp_norm <- transform_sample_counts(fung_tmp, function(x) x / sum(x))
  fung_otus_to_keep <- rownames(tax_table(filter_taxa(fung_tmp_norm, function(x) mean(x) > 0.0005, T)))
fung <- prune_taxa(fung_otus_to_keep, fung_tmp)
sample_sums(fung)/sample_sums(fung_tmp)

#samples need to be the same for the analysis so removing samples present in ITS that are not in 16S
fung <- prune_samples(setdiff(sample_names(fung), c("Run2_K37_ITS", "Run2_K39_ITS")), fung)

#renaming samples in phyloseq object according to Sample_ID
sample_names(bac) <- as.vector(unlist(sample_data(bac)[,"Sample_ID"]))
sample_names(fung) <- as.vector(unlist(sample_data(fung)[,"Sample_ID"]))

#order matters for this, so reordering according to Sample_ID
bac <- set_sample_order(bac, treatment="Sample_ID")
fung <- set_sample_order(fung, treatment="Sample_ID")

#sanity check on the order
table(sample_names(bac)==sample_names(fung))
#######################################################


#Cooccurrence analysis using Spieceasi
#######################################################
#running spiec.easi
taxadata <- data.frame(rbind(tax_table(bac), tax_table(fung)))
taxadata["24d821dddfbe031b81a44817e6d6a4c1","Genus"] <- "UR Enterobacteriaceae 2d8e"
taxadata["b6e52c74e408bcf8d49afe994fdcf359","Genus"] <- "UR Enterobacteriaceae bed4"
taxadata["32ab1812bd3770f4b7f6aad7273783da","Genus"] <- "UR Enterobacter 3bb2"
taxadata["dab25ac5ca12ad5ad505258e4c8fc7de",c("Genus", "Species")] <- "UR Pichiaceae dcd4"
taxadata["3621e3dbf9c25eaa08b32a94946ff3be",c("Genus", "Species")] <- "UR Pichiaceae 3f09"
taxadata$Genus <- gsub("(UR [A-za-z]*)( .{4})", "\\1", taxadata$Genus)


#spieceasi.net <- readRDS("spieceasi.RDS")

spieceasi.net <- spiec.easi(list(bac, fung), method='mb',lambda.min.ratio=1e-2, nlambda=250, icov.select.params=list(rep.num=1000))
getStability(spieceasi.net)
sum(getRefit(spieceasi.net))/2
betaMat <- as.matrix(symBeta(getOptBeta(spieceasi.net)))
  rownames(betaMat) <- rownames(taxadata)
  colnames(betaMat) <- rownames(taxadata)
stabilityMat <- as.matrix(getOptMerge(spieceasi.net))
  rownames(stabilityMat) <- rownames(taxadata)
  colnames(stabilityMat) <- rownames(taxadata)

ig2.mb <- adj2igraph(getRefit(spieceasi.net),  vertex.attr=list(name=rownames(taxadata)))

Isolated = which(degree(ig2.mb)<=0)
ig2.mb.pruned = delete.vertices(ig2.mb, Isolated)
network_data <- ggnetwork(ig2.mb.pruned)
network_data$degree <- degree(ig2.mb)[network_data$name]
network_data$connecting_node <- NA

#adding useful information to the network table
#name of connecting node
for(i in 1:nrow(network_data)){
  xend_val <- network_data[i,"xend"]
  yend_val <- network_data[i,"yend"]
  connecting_node <- subset(network_data, x==xend_val & y==yend_val)[1,"name"]
  network_data[i,"connecting_node"] <- connecting_node
}
#correspondence value
network_data$correspondence_val <- NA
for(i in 1:nrow(network_data)){
  name <- network_data[i,"name"]
  connecting_node <- network_data[i,"connecting_node"]
  value <- betaMat[name,connecting_node]
  network_data[i,"correspondence_val"] <- value
}
#stability (n/1000)
network_data$stability <- NA
for(i in 1:nrow(network_data)){
  name <- network_data[i,"name"]
  connecting_node <- network_data[i,"connecting_node"]
  value <- stabilityMat[name,connecting_node]
  network_data[i,"stability"] <- value
}

old_colnames <- colnames(network_data)
network_data <- cbind(network_data, taxadata[network_data$name,], taxadata[network_data$connecting_node,])
colnames(network_data) <- c(old_colnames, "Kingdom1", "Phylum1", "Class1", "Order1", "Family1", "Genus1", "Species1", "Kingdom2", "Phylum2", "Class2", "Order2", "Family2", "Genus2", "Species2")

#adding useful information to the output file
network_data$relationship <- NA
try(network_data[network_data$Kingdom1=="Bacteria" & network_data$Kingdom2=="Bacteria",]$relationship <- "Bacteria-Bacteria")
try(network_data[network_data$Kingdom1=="Fungi" & network_data$Kingdom2=="Bacteria",]$relationship <- "Fungi-Bacteria")
try(network_data[network_data$Kingdom1=="Bacteria" & network_data$Kingdom2=="Fungi",]$relationship <- "Fungi-Bacteria")
try(network_data[network_data$Kingdom1=="Fungi" & network_data$Kingdom2=="Fungi",]$relationship <- "Fungi-Fungi")

write.csv(network_data[!is.na(network_data$weight),!is.element(colnames(network_data), c("x", "y", "xend", "yend", "weight"))], "Network_analysis_data.csv", quote=F)

#adding markdown to italicize the names
network_data$Family1 = paste0("<i>", network_data$Family1, "</i>")
network_data$Family1 = gsub("^<i>(U[A-Z]) ", "\\1 <i>", network_data$Family1)
network_data$Family1 = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", network_data$Family1)

network_data$Genus1 = paste0("<i>", network_data$Genus1, "</i>")
network_data$Genus1 = gsub("^<i>(U[A-Z]) ", "\\1 <i>", network_data$Genus1)
network_data$Genus1 = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", network_data$Genus1)

network_data$Species1 = paste0("<i>", network_data$Species1, "</i>")
network_data$Species1 = gsub("^<i>(U[A-Z]) ", "\\1 <i>", network_data$Species1)
network_data$Species1 = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", network_data$Species1)

network_data$Species2 = paste0("<i>", network_data$Species2, "</i>")
network_data$Species2 = gsub("^<i>(U[A-Z]) ", "\\1 <i>", network_data$Species2)
network_data$Species2 = gsub(" ([a-z0-9]{4})</i>$", "</i> \\1", network_data$Species2)

#factoring some of the variables so they show up in the right order on the plot
network_data$Family1 <- factor(network_data$Family1, levels=c(sort(unique(network_data[network_data$Kingdom1=="Bacteria","Family1"])), sort(unique(network_data[network_data$Kingdom1=="Fungi","Family1"]))))
network_data$Genus1 <- factor(network_data$Genus1, levels=c(sort(unique(network_data[network_data$Kingdom1=="Bacteria","Genus1"])), sort(unique(network_data[network_data$Kingdom1=="Fungi","Genus1"]))))
network_data$Species1 <- factor(network_data$Species1, levels=c(sort(unique(network_data[network_data$Kingdom1=="Bacteria","Species1"])), sort(unique(network_data[network_data$Kingdom1=="Fungi","Species1"]))))

network_data <- network_data[order(network_data$Species1),]

positive_weight <- function(x) {x[x$correspondence_val>0,]}
negative_weight <- function(x) {x[x$correspondence_val<0,]}

bacterial_nodes <- subset(network_data, network_data$Kingdom1=="Bacteria")
  bacterial_nodes$node_labels <- gsub("<i>([A-Za-z]{1}).*</i>", "<i>\\1</i>", bacterial_nodes$Genus1)
  bacterial_nodes$node_labels <- gsub("UR ", "", bacterial_nodes$node_labels)
fungal_nodes <- subset(network_data, network_data$Kingdom1=="Fungi")
  fungal_nodes$node_labels <- gsub("<i>([A-Za-z]{1}).*</i>", "<i>\\1</i>", fungal_nodes$Genus1)
  fungal_nodes$node_labels <- gsub("UR ", "", fungal_nodes$node_labels)
  
network_plot <- ggplot(network_data, aes(x=x, y=y, xend=xend, yend=yend))+
  geom_edges(aes(size=abs(correspondence_val), alpha=stability), color="black", data=positive_weight)+
  geom_edges(aes(size=abs(correspondence_val), alpha=stability), color="purple", data=negative_weight, linetype="dashed")+
  geom_nodes(color="black", aes(size=degree*3+4), data=bacterial_nodes, shape="circle")+
  geom_nodes(aes(color=Genus1, size=degree*3), data=bacterial_nodes, shape="circle")+
  geom_richtext(aes(label=node_labels), size=2.5, fill=NA, label.color=NA, data=bacterial_nodes, fontface="bold", family="serif")+
  geom_nodes(color="black", aes(size=degree*3+4), data=fungal_nodes, shape="square")+
  geom_nodes(aes(color=Genus1, size=degree*3), data=fungal_nodes, shape="square")+
  geom_richtext(aes(label=node_labels), size=2.5, fill=NA, label.color=NA, data=fungal_nodes, fontface="bold", family="serif")+
  scale_shape(name="Kingdom")+
  scale_color_manual(name="Genus", values=c(get_palette("startrek", length(unique(network_data[network_data$Kingdom1=="Bacteria","Genus1"]))), get_palette("rickandmorty", length(unique(network_data[network_data$Kingdom1=="Fungi","Genus1"])))))+
  guides(color=guide_legend(override.aes=list(size=3.5)), alpha="none")+
  scale_size(guide="none")+
  theme_blank()+
  theme(legend.position="right", legend.text=element_markdown(family="serif", size=11), legend.title=element_text(family="serif", size=13))+
  theme(panel.background=element_rect(fill = "#f9f9f9"), plot.background=element_rect(fill = "#f9f9f9"), legend.background=element_rect(fill = "#f9f9f9"))

ggsave(plot=network_plot, filename="Network_plot.tiff", width=10, height=6, units="in", dpi="print")
#######################################################
