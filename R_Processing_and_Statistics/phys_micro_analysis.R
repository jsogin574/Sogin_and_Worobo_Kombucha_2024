#title: "Kombucha-chemical-micro"
#author: "Jonathan Sogin"
#date: "2023"


#Importing libraries
#######################################################

#visualization packages

library("ggpubr"); packageVersion("ggpubr")
library("ggtext"); packageVersion("ggtext")
library("ggnewscale"); packageVersion("ggnewscale")

#data analysis packages
library("multcomp"); packageVersion("multcomp")
library("margins"); packageVersion("margins")
library("ecodist"); packageVersion("ecodist")

#setting seed
addTaskCallback(function(...) {set.seed(02221997);TRUE})
#######################################################

#importing data
#######################################################
metadata_file <- read.csv("./Qiime_Data/metadata.csv", fileEncoding="UTF-8-BOM")

reported <- metadata_file$Nutrition.Sugar..g.L.
metadata_file$calculated_sugar <- metadata_file$Glucose..g.L. + metadata_file$Fructose..g.L.
metadata_file$rel_sugar_difference <- (metadata_file$calculated_sugar - reported)/reported

metadata_graphs <- subset(metadata_file, TRUE.sample.or.CONTROL=="TRUE SAMPLE")
metadata_graphs$Blinded_Brand <- gsub("Brand_", "", metadata_graphs$Blinded_Brand)
#######################################################


#######################################################

metadata_graphs$Blinded_Brand=factor(metadata_graphs$Blinded_Brand, levels=seq(1,6))
metadata_graphs$Blinded_Product=factor(metadata_graphs$Blinded_Product, levels=c("Product_1", "Product_2", "Product_3", "Product_4", "Product_5", "Product_6", "Product_7", "Product_8", "Product_9", "Product_10", "Product_11", "Product_12", "Product_13", "Product_14", "Product_15", "Product_16", "Product_17", "Product_18", "Product_19", "Product_20", "Product_21", "Product_22", "Product_23", "Product_24"))

#plotting relevant sample chemistry
#plotting ethanol composition of samples, including those that were excluded for metagenomic analysis
#conducting t-test on each sample to determine whether the mean is greater than 0.5% abv
#no p value correction was done
abv_stats <- c()
for(b in seq(1:6)){
  data <- subset(metadata_graphs, Blinded_Brand==b)
  abv_stats[[b]] <- t.test(data$Ethanol...., alternative="greater", mu=0.5, na.action=na.omit)
}
#plotting ethanol composition
abv_plot <- ggboxplot(metadata_graphs, x="Blinded_Brand", y="Ethanol....", add="jitter", ylab = "% Ethanol (v/v)", xlab="Brand", color="Blinded_Brand", palette="simpsons", legend="none", title="A")+
  geom_hline(yintercept=0.5, col="red", size=0.75, linetype="dashed")+
  scale_y_continuous(limits=c(0,1.4))+
  scale_x_discrete(labels=seq(1,6))+
  theme(plot.title=element_text(face="bold", size=15))+
  theme(text=element_text(family="serif"))+
  geom_text(label="***", y=1.375, x=3, fontface="bold", size=8, family="serif")+
  geom_text(label="**", y=1.375, x=4, fontface="bold", size=8, family="serif")

#conducting t-test on each sample to determine whether the mean is different from zero
#no p value correction was done
sugar_stats <- c()
for(b in seq(1:6)){
  data <- subset(metadata_graphs, Blinded_Brand==b)
  sugar_stats[[b]] <- t.test(data$rel_sugar_difference, alternative="two.sided", mu=0, na.action=na.omit)
}
#plotting sugar difference
sugar_plot <- ggboxplot(metadata_graphs, x="Blinded_Brand", y="rel_sugar_difference", add="jitter", ylab = "Rel Sugar Diff vs Label", xlab="Brand", color="Blinded_Brand", palette="simpsons", legend="none", title="B")+
  geom_hline(yintercept=0, col="black", lwd=0.25, linetype="solid")+
  scale_y_continuous(limits=c(-0.6,0.95), breaks=c(-0.5,0,0.5), expand=c(0,0))+
  scale_x_discrete(labels=seq(1,6))+
  theme(plot.title=element_text(face="bold", size=15))+
  theme(text=element_text(family="serif"))+
  geom_text(label="**", y=0.85, x=2, fontface="bold", size=8, family="serif")+
  geom_text(label="*", y=0.85, x=3, fontface="bold", size=8, family="serif")+
  geom_text(label="*", y=0.85, x=5, fontface="bold", size=8, family="serif")

sugar_eth_cor_plot <- ggscatter(data=metadata_graphs, x="rel_sugar_difference", y="Ethanol....", add="reg.line", legend="right", cor.coef=T, ylab="% Ethanol (v/v)", xlab="Rel Sugar Diff vs Label", cor.coeff.args=list(method="kendall", aes(label=paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.y=1.25, label.x=-0.3, family="serif"))+
  geom_point(aes(color=Blinded_Brand), size=5)+
  scale_color_discrete(name="Brand")+
  theme(text=element_text(family="serif"))

ggsave(plot=sugar_eth_cor_plot, filename="Sugar_Eth_Cor_Plot.tiff", width=6, height=4, units="in", dpi="print")

combined_sample_chemistry_plot <- ggarrange(abv_plot, sugar_plot, ncol=1)
ggsave(plot=combined_sample_chemistry_plot, filename="Sample_Chemistry_Plot.tiff", width=6, height=5, units="in", dpi="print")

#microbiological characteristics of kombucha samples
GYC_data <- cbind(metadata_graphs[,c("Blinded_Brand", "Blinded_Product", "GYC.N..log10CFU.mL.")], rep("GYC", length(metadata_graphs$GYC.N..log10CFU.mL.)))
  colnames(GYC_data)<-c("Blinded_Brand", "Blinded_Product", "log10CFU.mL", "media")
MRS_data <- cbind(metadata_graphs[,c("Blinded_Brand", "Blinded_Product", "MRS.N..log10CFU.mL.")], rep("MRS", length(metadata_graphs$MRS.N..log10CFU.mL.)))
  colnames(MRS_data)<-c("Blinded_Brand", "Blinded_Product", "log10CFU.mL", "media")
APDA_data <- cbind(metadata_graphs[,c("Blinded_Brand", "Blinded_Product", "APDA..log10CFU.mL.")], rep("APDA", length(metadata_graphs$APDA..log10CFU.mL.)))
  colnames(APDA_data)<-c("Blinded_Brand", "Blinded_Product", "log10CFU.mL", "media")
micro_data <- rbind(GYC_data, MRS_data, APDA_data)
micro_data$media=factor(micro_data$media, levels=c("GYC", "MRS", "APDA"))
  micro_data <- subset(micro_data, !is.na(log10CFU.mL))
#microbial load stats
load_stats <- c()
  load_stats[["GYC"]][["model"]] <- aov(lm(data=GYC_data, log10CFU.mL ~ Blinded_Brand))
  load_stats[["MRS"]][["model"]] <- aov(lm(data=MRS_data, log10CFU.mL ~ Blinded_Brand))
  load_stats[["APDA"]][["model"]] <- aov(lm(data=APDA_data, log10CFU.mL ~ Blinded_Brand))
  load_stats[["GYC"]][["summary"]] <- summary(load_stats[["GYC"]][["model"]])
  load_stats[["MRS"]][["summary"]] <- summary(load_stats[["MRS"]][["model"]])
  load_stats[["APDA"]][["summary"]] <- summary(load_stats[["APDA"]][["model"]])
  load_stats[["GYC"]][["TukeyHSD"]] <- glht(load_stats[["GYC"]][["model"]], linfct=mcp(Blinded_Brand="Tukey"))
  load_stats[["MRS"]][["TukeyHSD"]] <- glht(load_stats[["MRS"]][["model"]], linfct=mcp(Blinded_Brand="Tukey"))
  load_stats[["APDA"]][["TukeyHSD"]] <- glht(load_stats[["APDA"]][["model"]], linfct=mcp(Blinded_Brand="Tukey"))
  load_stats[["GYC"]][["cld"]] <- cld(load_stats[["GYC"]][["TukeyHSD"]])
  load_stats[["MRS"]][["cld"]] <- cld(load_stats[["MRS"]][["TukeyHSD"]])
  load_stats[["APDA"]][["cld"]] <- cld(load_stats[["APDA"]][["TukeyHSD"]])

#exporting compact letter display to use with the plot
GYC_lettering <- cbind(load_stats[["GYC"]][["cld"]][["mcletters"]][["Letters"]], rep("GYC", 6), names(load_stats[["GYC"]][["cld"]][["mcletters"]][["Letters"]]))
MRS_lettering <- cbind(load_stats[["MRS"]][["cld"]][["mcletters"]][["Letters"]], rep("MRS", 6), names(load_stats[["MRS"]][["cld"]][["mcletters"]][["Letters"]]))
APDA_lettering <- cbind(load_stats[["APDA"]][["cld"]][["mcletters"]][["Letters"]], rep("APDA", 6), names(load_stats[["APDA"]][["cld"]][["mcletters"]][["Letters"]]))
stats_lettering <- rbind(GYC_lettering, MRS_lettering, APDA_lettering)
  colnames(stats_lettering) <- c("letter", "media", "Blinded_Brand")
  stats_lettering <- data.frame(stats_lettering)
  #substituting letters for the different media types to not confuse comparisons between media
  stats_lettering$letter[stats_lettering$media=="MRS"] <- gsub("a", "c", stats_lettering$letter[stats_lettering$media=="MRS"])
  stats_lettering$letter[stats_lettering$media=="MRS"] <- gsub("b", "d", stats_lettering$letter[stats_lettering$media=="MRS"])
  stats_lettering$letter[stats_lettering$media=="APDA"] <- gsub("a", "e", stats_lettering$letter[stats_lettering$media=="APDA"])
  stats_lettering$letter[stats_lettering$media=="APDA"] <- gsub("b", "f", stats_lettering$letter[stats_lettering$media=="APDA"])
  stats_lettering$media=factor(stats_lettering$media, levels=c("GYC", "MRS", "APDA"))
  
micro_plot <- ggboxplot(micro_data, x="Blinded_Brand", y="log10CFU.mL", color="media", palette="simpsons", add="jitter", xlab="Brand", ylab="log<sub>10</sub>(CFU/mL)", legend="none")+
  facet_wrap(vars(media), scales="free_x", strip.position="top")+
  scale_y_continuous(limits=c(0,8), breaks=c(2,4,6), expand=c(0,0), labels=c("2.0", "4.0", "6.0"))+
  scale_x_discrete(labels=as_labeller(function(x){gsub("Brand_", " ", x)}))+
  theme(axis.title.y=element_markdown())+
  theme(text=element_text(family="serif"))+
  geom_text(data=stats_lettering, aes(x=Blinded_Brand, y=7.7, label=letter), family="serif", fontface="italic")
ggsave(plot=micro_plot, filename="Microbial_Load_Plot.tiff", width=8, height=3, units="in", dpi="print")

#summary statistics for chemical variables
summary(metadata_graphs$pH)
summary(metadata_graphs$Ethanol....)
summary(metadata_graphs$Acetic..g.L.)
summary(metadata_graphs$Lactic..g.L.)
summary(metadata_graphs$calculated_sugar)
summary(metadata_graphs$rel_sugar_difference)

#summary statistics for microbial load variables
summary(metadata_graphs$GYC.N..log10CFU.mL.)
summary(metadata_graphs$MRS.N..log10CFU.mL.)
summary(metadata_graphs$APDA..log10CFU.mL.)

#multiple linear regressesions for explanatory variables, chemical variables, and microbial load variables
explanatory_fit <- lm(Ethanol.... ~ Blinded_Brand + Blinded_Product, data=data.frame(metadata_graphs), na.action=na.omit)
summary(explanatory_fit)

physiochemical_fit <- lm(Ethanol.... ~ Blinded_Brand + pH + Lactic..g.L. + Acetic..g.L. + calculated_sugar + rel_sugar_difference, data=data.frame(metadata_graphs), na.action=na.omit)
summary(physiochemical_fit)

microbiological_fit <- lm(Ethanol.... ~ Blinded_Brand + GYC.N..log10CFU.mL. + MRS.N..log10CFU.mL. + APDA..log10CFU.mL., data=data.frame(metadata_graphs), na.action=na.omit)
summary(microbiological_fit)
#######################################################

#plotting pH vs ethanol content
#######################################################
pH_cor_plot <- ggscatter(data=metadata_graphs, x="pH", y="Ethanol....", add="reg.line", legend="right", cor.coef=T, ylab="% Ethanol (v/v)", xlab="pH", cor.coeff.args=list(method="kendall", aes(label=paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), label.x=2.9, family="serif"))+
  geom_point(aes(color=Blinded_Brand), size=5)+
  scale_color_manual(name="Brand", values=get_palette("simpsons", 6))+
  theme(text=element_text(family="serif"))

ggsave(plot=pH_cor_plot, filename="pH_Eth_Cor_Plot.tiff", width=6, height=4, units="in", dpi="print")
#######################################################

#physiochemical ordination
#######################################################

phys_ord <- metadata_graphs

sample_data <- data.frame(metadata_graphs[,c("pH", "Ethanol....", "Lactic..g.L.", "Acetic..g.L.", "calculated_sugar", "rel_sugar_difference")])
mahal_sq <- distance(sample_data, method="mahal")
mahal <- sqrt(mahal_sq)
nmds.out <- nmds(mahal)
scores <- nmds.min(nmds.out)
chem_ordination <- cbind(metadata_graphs, scores)

chem_ordination_plot <- ggscatter(data=chem_ordination, x="X1", y="X2", xlab="NMDS1", ylab="NMDS2", color="Blinded_Brand", legend.title="Brand", palette="simpsons", legend="right", point=F, ellipse=T, ellipse.alpha=0.25, ellipse.type="convex")+
  new_scale_color()+
  geom_point(aes(color=Ethanol...., size=4))+
  scale_color_steps2("% Ethanol (v/v)", low="darkgreen", midpoint=0.5, high="darkred")+
  guides(size="none")+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank())+
  theme(legend.title=element_text(face="bold"))+
  theme(text=element_text(family="serif"))

ggsave(plot=chem_ordination_plot, filename="Chemical_Ordination.tiff", width=8, height=4, units="in", dpi="print")
#######################################################


