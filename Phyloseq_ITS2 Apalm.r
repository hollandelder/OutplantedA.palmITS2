#####################################

# installing packages: execute just once when first using a script:
#Second bioconductor install method is commented out because it is for older than 3.5.0 versions of R. The first version of the install is for recent versions of R
#you only need to install once
source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")
#OR
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
#
# install.packages("~/Downloads/dada2-1.4", #to install from source, just indicate pkg download location
                 # repos = NULL,
                 # type = "source",
                 # dependencies = c("Depends", "Suggests","Imports"))

#####################################

library(ggplot2); packageVersion("ggplot2") #3.3.2
library(phyloseq); packageVersion("phyloseq") #1.32.0
library("Rmisc")
library("ggpubr")
library(cowplot)
library("vegan")
library(writexl)
library(readr)
library(gridExtra)

# All trimming, filtering, merging, pcr duplication removal, and pcr contaminate removal was handled by phyloseq
# See NotesOnSymportalCodeHistory.txt for record of Symportal code

############ Begin Phyloseq anaylsis here ###############

#Set path to absolute values of collapsed symb types and the sample meta data
path <- "C:/Users/elder/OneDrive/Documents/UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis" # CHANGE ME to the main working directory with the absolute values of symb types and sample meta data

###########################################################
# import files to make the collapsed type phyloseq object #
###########################################################

# import the collapsed absolute abundance and meta data file from SymPortal that has been adjusted for phyloseq. 
# Change the SymPortal number ID to the sample ID these must match the meta data file.
# Remove the top seven rows of the absolute abundance and meta data file from SymPortal.
counts <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/Rready_ProfilesAbsoluteAbund.csv", header = TRUE, row.names=1, check.names = FALSE)
View(counts)
# plot counts checking for samples with zero.
plot(rowSums(counts))
# there are none :)

#import metadata file
meta <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/samplemetadata.csv", header = TRUE)
rownames(meta) <- meta$Sample
View(meta)

#import the taxa file made from the 
taxa <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/symportaltaxa.csv", header = TRUE)
rownames(taxa) <- as.factor(taxa$ITS2typeprofile)
mtaxa <- as.matrix(taxa)
View(mtaxa)

# import all the types absolute abundance and meta data file from SymPortal for phyloseq. 
# Change the SymPortal number ID to the sample ID these must match the meta data file.
all_counts <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/All_RreadyProfilesAbsoluteAbund.csv", header = TRUE, row.names=1, check.names = FALSE)
View(all_counts)
# plot counts checking for samples with zero.
plot(rowSums(all_counts))
# there are none :)

#import metadata file
all_meta <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/samplemetadata.csv", header = TRUE)
rownames(all_meta) <- all_meta$Sample
View(all_meta)

# import the taxa file
# Made this taxa file from the top row of the absolute abundance file added some clade and ID info, but it is not as complete as the post analysis taxa file
all_taxa <- read.csv("UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/All_symportaltaxa.csv", header = TRUE)
rownames(all_taxa) <- as.factor(all_taxa$ITS2typeprofile)
all_mtaxa <- as.matrix(all_taxa)
View(all_mtaxa)


#########################################################
####### handoff 2 phyloseq from SymPortal outputs #######
#########################################################

ps=phyloseq(otu_table(counts, taxa_are_rows=FALSE), sample_data(meta), tax_table(mtaxa))
ps

tax_table(ps)

psall=phyloseq(otu_table(all_counts, taxa_are_rows=FALSE), sample_data(all_meta), tax_table(all_mtaxa))
psall

tax_table(psall)


#Visualize alpha-diversity - ***Should be done on raw, untrimmed dataset***
#total species diversity in a landscape (gamma diversity) is determined by two different things, the mean species diversity in sites or habitats at a more local scale (alpha diversity) and the differentiation among those habitats (beta diversity)
#Shannon:Shannon entropy quantifies the uncertainty (entropy or degree of surprise) associated with correctly predicting which letter will be the next in a diverse string. Based on the weighted geometric mean of the proportional abundances of the types, and equals the logarithm of true diversity. When all types in the dataset of interest are equally common, the Shannon index hence takes the value ln(actual # of types). The more unequal the abundances of the types, the smaller the corresponding Shannon entropy. If practically all abundance is concentrated to one type, and the other types are very rare (even if there are many of them), Shannon entropy approaches zero. When there is only one type in the dataset, Shannon entropy exactly equals zero (there is no uncertainty in predicting the type of the next randomly chosen entity).
#Simpson:equals the probability that two entities taken at random from the dataset of interest represent the same type. equal to the weighted arithmetic mean of the proportional abundances pi of the types of interest, with the proportional abundances themselves being used as the weights. Since mean proportional abundance of the types increases with decreasing number of types and increasing abundance of the most abundant type, λ obtains small values in datasets of high diversity and large values in datasets of low diversity. This is counterintuitive behavior for a diversity index, so often such transformations of λ that increase with increasing diversity have been used instead. The most popular of such indices have been the inverse Simpson index (1/λ) and the Gini–Simpson index (1 − λ).

# make both into an RDS and only run this once
#saveRDS(ps,"UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/ps.its2.RDS")
#saveRDS(psall,"UniversitySouthernCa/CoralProjects/Apalm_TimeCourse/Processing/AmpliconSeq/ITS2/ITS23rdAnalysis/ps.all.its2.RDS")

# read in RDS if you need to
ps <- readRDS("ps.its2.RDS")
ps.all <- readRDS("ps.all.its2.RDS")


#######################################
########## Look at Diversity ##########
#######################################

# plot Shannon and Simpson diversity
plot_richness(ps, x="Genotype", measures=c("Shannon", "Simpson"), color="Site") + theme_bw()
#plot all Shannon and Simpson diversity
plot_richness(psall, x="Genotype", measures=c("Shannon", "Simpson"), color="Site") + theme_bw()


# first looks at data
#all pre-collapsed
plot_bar(psall, "Sample", fill="ITS2typeprofile")
#post-collapsed
plot_bar(ps, "Sample", fill="ITS2typeprofile")


#plots separated by genotype add site panels
plot_bar(ps,"Genotype", fill="ITS2typeprofile",facet_grid=~Site)


# subset for later exploration of different time points
ps_T24 = subset_samples(ps, sample_data(ps)$TimePoint != "T0")
ps_T24

ps_T0 = subset_samples(ps, sample_data(ps)$TimePoint != "T24")
ps_T0


# Absolute abundance plot
T24.sitebygenotype_absolbar <- plot_bar(ps_T24,"Genotype", fill="ITS2typeprofile",facet_grid=~Site) +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T24.sitebygenotype_absolbar




# Get Relative Abundance of symbs by timepoint
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
plot_bar(ps.rel,"ITS2typeprofile", fill="ITS2typeprofile",facet_grid=~Site)

psT0.rel = subset_samples(ps.rel, sample_data(ps.rel)$TimePoint != "T24")
plot_bar(psT0.rel,"ITS2typeprofile", fill="ITS2typeprofile",facet_grid=~Site)

psT24.rel = subset_samples(ps.rel, sample_data(ps.rel)$TimePoint != "T0")
plot_bar(psT24.rel,"ITS2typeprofile", fill="ITS2typeprofile",facet_grid=~Site)


# subset ps_T24 into Sites
ps_T24BP.rel = subset_samples(psT24.rel, SiteInitials == "BP")
ps_T24BP.rel

ps_T24DL.rel = subset_samples(psT24.rel, SiteInitials == "DL")
ps_T24DL.rel

ps_T24EDR.rel = subset_samples(psT24.rel, SiteInitials == "EDR")
ps_T24EDR.rel

ps_T24ES.rel = subset_samples(psT24.rel, SiteInitials == "ES")
ps_T24ES.rel

ps_T24LK.rel = subset_samples(psT24.rel, SiteInitials == "LK")
ps_T24LK.rel

ps_T24M32.rel = subset_samples(psT24.rel, SiteInitials == "M32")
ps_T24M32.rel

ps_T24MS.rel = subset_samples(psT24.rel, SiteInitials == "MS")
ps_T24MS.rel

ps_T24WS.rel = subset_samples(psT24.rel, SiteInitials == "WS")
ps_T24WS.rel


#Relative abundance bar plots 
all.bar <- plot_bar(ps.rel,"Sample",fill="ITS2typeprofile") +
  geom_bar(stat="identity")+
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

all.bar

T0.bar <- plot_bar(psT0.rel,"Sample",fill="ITS2typeprofile") +
  geom_bar(stat="identity")+
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T0.bar

T0.genotype_bar <- plot_bar(psT0.rel, "Genotype", fill = "ITS2typeprofile") +
  geom_bar(stat = "identity")+
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T0.genotype_bar

T24.bar <- plot_bar(psT24.rel,"Sample",fill="ITS2typeprofile")+
  geom_bar(stat="identity")+  
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T24.bar

T24.genotype_bar <- plot_bar(psT24.rel,"Genotype",fill="ITS2typeprofile")+
  geom_bar(stat="identity")+  
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T24.genotype_bar

T24.site_bar <- plot_bar(psT24.rel,"Site",fill="ITS2typeprofile")+
  geom_bar(stat="identity")+  
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T24.site_bar

T24.sitebygenotype_bar <- plot_bar(psT24.rel, x="Genotype", fill="ITS2typeprofile") + 
  facet_wrap(~Site, scales="free_y") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))

T24.sitebygenotype_bar

####### Finally subset by sample
T24.sitebysample_bar <- plot_bar(psT24.rel, x="Sample", fill="ITS2typeprofile", facet_grid=~Site) +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))
T24.sitebysample_bar
# This is a mess. Need to make the panels individually and then put them together.

# Subset by sample plot all samples and then stitch back together.
BP_RelAbundT24 <- plot_bar(ps_T24BP.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
BP_RelAbundT24

DL_RelAbundT24 <- plot_bar(ps_T24DL.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
DL_RelAbundT24

EDR_RelAbundT24 <- plot_bar(ps_T24EDR.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
EDR_RelAbundT24

ES_RelAbundT24 <- plot_bar(ps_T24ES.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
ES_RelAbundT24

LK_RelAbundT24 <- plot_bar(ps_T24LK.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
LK_RelAbundT24

M32_RelAbundT24 <- plot_bar(ps_T24M32.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
M32_RelAbundT24

MS_RelAbundT24 <- plot_bar(ps_T24MS.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral")) +
  theme(legend.position="none")
MS_RelAbundT24

WS_RelAbundT24 <- plot_bar(ps_T24WS.rel, x="Sample", fill = "ITS2typeprofile") +
  scale_fill_manual(name="ITS2 type profiles", values=c("paleturquoise2","cyan3","aquamarine2","darkgoldenrod1","darkgoldenrod3","mediumseagreen","lightcoral"))
WS_RelAbundT24

grid.arrange(BP_RelAbundT24, DL_RelAbundT24, EDR_RelAbundT24, ES_RelAbundT24, LK_RelAbundT24, M32_RelAbundT24, MS_RelAbundT24, WS_RelAbundT24,  nrow = 2)

# Ended up using the grid extra graph and rearranging, because the other one is ugly.

#There is D present in 13-XK nursery samples
#Suggests that the frags were capable of reshuffling to D in the field. They are capable of picking up new symbs too and all genotypes are able to do this. However, they do not seem to do this in large numbers. 
#In other words, few individuals switch. Only 5 individuals actually completely switch to predominately and A type. While 9 individuals take up new types, but the new symbiont does not become the predominate symbiont type.


#############################################
##################  PCoA  ###################
#############################################

#Bray Curtis PCoA All collapsed types with PCA plot
ps.ord <- ordinate(ps,"PCoA",distance="bray")
ps.ord

plot_ordination(ps, ps.ord, color ="Genotype", shape="Site")+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"))+
  scale_shape_manual(values=c(15,16,17,18,20,8,11,0,6))+
  geom_point(alpha=0.5, size=5)+
  #stat_ellipse(aes(linetype=Genotype))+
  theme_cowplot()

ps.bp <- subset_samples(ps,SiteInitials=="BP")
gg.pc.bp.types <- plot_ordination(ps.bp, ordinate(ps.bp,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Big Pine")
gg.pc.bp.types

ps.dl <- subset_samples(ps,SiteInitials=="DL")
gg.pc.dl.types <- plot_ordination(ps.dl, ordinate(ps.dl,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Dave's Ledge")
gg.pc.dl.types

ps.edr <- subset_samples(ps,SiteInitials=="EDR")
gg.pc.edr.types <- plot_ordination(ps.edr, ordinate(ps.edr,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Eastern Dry Rocks")
gg.pc.edr.types

ps.es <- subset_samples(ps,SiteInitials=="ES")
gg.pc.es.types <- plot_ordination(ps.es, ordinate(ps.es,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Eastern Sambo")
gg.pc.es.types

ps.lk <- subset_samples(ps,SiteInitials=="LK")
gg.pc.lk.types <- plot_ordination(ps.lk, ordinate(ps.lk,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Looe Key")
gg.pc.lk.types

ps.m32 <- subset_samples(ps,SiteInitials=="M32")
gg.pc.m32.types <- plot_ordination(ps.m32, ordinate(ps.m32,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Marker 32")
gg.pc.m32.types

ps.ms <- subset_samples(ps,SiteInitials=="MS")
gg.pc.ms.types <- plot_ordination(ps.ms, ordinate(ps.ms,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Maryland Shoals")
gg.pc.ms.types

ps.ws <- subset_samples(ps,SiteInitials=="WS")
gg.pc.ws.types <- plot_ordination(ps.ws, ordinate(ps.ws,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Western Sambo")
gg.pc.ws.types

ps.rw <- subset_samples(ps,SiteInitials=="RW")
gg.pc.rw.types <- plot_ordination(ps.rw, ordinate(ps.rw,"PCoA",distance="bray"), color ="Genotype")+
  geom_point(size = 10)+
  stat_ellipse()+
  scale_color_manual(values=c("rosybrown", "paleturquoise4", "seagreen", "goldenrod", "salmon"),labels=c("13-X5","13-X7", "13-X9", "13-XK", "14-5"))+
  guides(color=guide_legend(title="Genotype"))+
  theme_cowplot()+
  #xlab("Axis 1 (61.4%)")+
  #ylab("Axis 2 (12.2%)")+
  #theme(legend.position="none")+
  ggtitle("Race Way")
gg.pc.rw.types


#####################   Stats   ########################
library(vegan)
install.packages("remotes")
remotes::install_github("Jtrachsel/funfuns")
library(funfuns)
library(dplyr)
library(edgeR)


##################### All #########################
seq.ps <- data.frame(psall@otu_table)
samdf.ps <- data.frame(psall@sam_data)
dist.ps <- vegdist(seq.ps)
bet.ps <- betadisper(dist.ps,samdf.ps$Site)
anova(bet.ps) #ns
permutest(bet.ps,pairwise=TRUE,permutations=999)
plot(bet.ps)
adonis(seq.ps ~ Site, data=samdf.ps, permutations=999) #p<001***
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)    
# Site        8    1.9145 0.239312  2.8279 0.14631  0.001 ***
# Residuals 132   11.1706 0.084626         0.85369           
# Total     140   13.0851                  1.00000           
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.adonis(seq.ps, factors=samdf.ps$Site, permutations=999) # The Raceway versus all the sites with A, and teh sites with A versus other sites with no A.
#                                pairs    F.Model          R2 p.value p.adjusted
# 1                 Raceway vs BigPine 0.96460200 0.017235931   0.433      0.433
# 2              Raceway vs DavesLedge 9.28395692 0.128437309   0.002      0.002
# 3         Raceway vs EasternDryRocks 2.82028550 0.044894503   0.042      0.042
# 4            Raceway vs EasternSambo 3.69917979 0.054641427   0.036      0.036
# 5                 Raceway vs LooeKey 4.32540950 0.064246316   0.018      0.018
# 6                Raceway vs Marker32 1.76568870 0.029543518   0.189      0.189
# 7          Raceway vs MarylandShoals 2.42394878 0.037049870   0.114      0.114
# 8            Raceway vs WesternSambo 2.43657382 0.037235657   0.119      0.119
# 9              BigPine vs DavesLedge 1.90732724 0.106510995   0.172      0.172
# 10        BigPine vs EasternDryRocks 0.56089835 0.041361445   0.751      0.751
# 11           BigPine vs EasternSambo 0.27953602 0.016177287   0.777      0.777
# 12                BigPine vs LooeKey 0.68986443 0.041334334   0.578      0.578
# 13               BigPine vs Marker32 0.15755551 0.014120970   0.943      0.943
# 14         BigPine vs MarylandShoals 0.01629302 0.001017278   0.962      0.962
# 15           BigPine vs WesternSambo 0.75431663 0.045022226   0.425      0.425
# 16     DavesLedge vs EasternDryRocks 1.00034935 0.045469703   0.365      0.365
# 17        DavesLedge vs EasternSambo 4.60597334 0.155575812   0.021      0.021
# 18             DavesLedge vs LooeKey 1.25414414 0.049660924   0.226      0.226
# 19            DavesLedge vs Marker32 3.14584848 0.142051386   0.059      0.059
# 20      DavesLedge vs MarylandShoals 4.64029706 0.162019865   0.024      0.024
# 21        DavesLedge vs WesternSambo 5.88777097 0.196995988   0.002      0.002
# 22   EasternDryRocks vs EasternSambo 1.17138902 0.050553250   0.319      0.319
# 23        EasternDryRocks vs LooeKey 0.09843059 0.004665304   0.962      0.962
# 24       EasternDryRocks vs Marker32 0.90906622 0.053762059   0.480      0.480
# 25 EasternDryRocks vs MarylandShoals 1.32988422 0.059556252   0.208      0.208
# 26   EasternDryRocks vs WesternSambo 2.23726597 0.096279226   0.031      0.031
# 27           EasternSambo vs LooeKey 1.13502590 0.043429301   0.442      0.442
# 28          EasternSambo vs Marker32 0.43999342 0.021526104   0.592      0.592
# 29    EasternSambo vs MarylandShoals 1.13821590 0.043546044   0.288      0.288
# 30      EasternSambo vs WesternSambo 2.74710052 0.099004958   0.076      0.076
# 31               LooeKey vs Marker32 1.13983672 0.056596125   0.282      0.282
# 32         LooeKey vs MarylandShoals 1.74039464 0.067613363   0.122      0.122
# 33           LooeKey vs WesternSambo 2.86656191 0.106696269   0.011      0.011
# 34        Marker32 vs MarylandShoals 0.25258610 0.013119593   0.941      0.941
# 35          Marker32 vs WesternSambo 0.85520636 0.043072147   0.388      0.388
# 36    MarylandShoals vs WesternSambo 0.91169055 0.036596896   0.376      0.376


#############################  T0  #############################
# subset all T0
ps_T0_stats = subset_samples(psall, SiteInitials == "RW")
ps_T0_stats
# do stats within T0
seq.psT0 <- data.frame(ps_T0_stats@otu_table)
samdf.psT0 <- data.frame(ps_T0_stats@sam_data)
dist.psT0 <- vegdist(seq.psT0)
bet.psT0 <- betadisper(dist.psT0,samdf.psT0$Genotype)
anova(bet.psT0) #ns
permutest(bet.psT0,pairwise=TRUE,permutations=999)
plot(bet.psT0)
adonis(seq.psT0 ~ Genotype, data=samdf.psT0, permutations=999) #p<001*** 13-XK seems to be driving this pattern
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   4    3.8277 0.95692  94.902 0.88983  0.001 ***
# Residuals 47    0.4739 0.01008         0.11017           
# Total     51    4.3016                 1.00000           
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.adonis(seq.psT0, factors=samdf.psT0$Genotype, permutations=999) # not sure what is going on here
#             pairs     F.Model         R2 p.value p.adjusted
# 1  13-X5 vs 13-X7   2.4583140 0.10946120   0.130      0.130
# 2  13-X5 vs 13-X9   8.1375699 0.28920656   0.006      0.006
# 3  13-X5 vs 13-XK 253.0905235 0.92676421   0.001      0.001
# 4   13-X5 vs 14-5   0.7652043 0.03685031   0.417      0.417
# 5  13-X7 vs 13-X9  10.9082862 0.37734116   0.003      0.003
# 6  13-X7 vs 13-XK 177.8636797 0.90809935   0.001      0.001
# 7   13-X7 vs 14-5   0.5763405 0.03102551   0.484      0.484
# 8  13-X9 vs 13-XK 204.2800586 0.91902108   0.001      0.001
# 9   13-X9 vs 14-5   5.8283787 0.24459821   0.018      0.018
# 10  13-XK vs 14-5 152.4087595 0.89437163   0.001      0.001


#############################  T24  #############################
# subset all T24
ps_T24_stats = subset_samples(psall, TimePoint == "T24")
ps_T24_stats
# do stats within T24 Genotype
seq.psT24 <- data.frame(ps_T24_stats@otu_table)
samdf.psT24 <- data.frame(ps_T24_stats@sam_data)
dist.psT24 <- vegdist(seq.psT24)
bet.psT24 <- betadisper(dist.psT24,samdf.psT24$Genotype)
anova(bet.psT24) #ns
permutest(bet.psT24,pairwise=TRUE,permutations=999)
plot(bet.psT24)
adonis(seq.psT24 ~ Genotype, data=samdf.psT24, permutations=999) #p<001***
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Genotype   4    0.2441 0.061015 0.65513 0.03025  0.764
# Residuals 84    7.8233 0.093135         0.96975       
# Total     88    8.0674                  1.00000        
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.adonis(seq.psT24, factors=samdf.psT24$Genotype, permutations=999) #No significant differenecs between genotypes
#             pairs    F.Model          R2 p.value p.adjusted
# 1   14-5 vs 13-X7 1.17022737 0.030658119   0.310      0.310
# 2   14-5 vs 13-X9 1.54529701 0.042284429   0.209      0.209
# 3   14-5 vs 13-XK 0.71063859 0.019899913   0.449      0.449
# 4   14-5 vs 13-X5 0.60053446 0.015971434   0.495      0.495
# 5  13-X7 vs 13-X9 0.50444530 0.015519271   0.591      0.591
# 6  13-X7 vs 13-XK 0.05512782 0.001719782   0.962      0.962
# 7  13-X7 vs 13-X5 0.28072048 0.008188873   0.758      0.758
# 8  13-X9 vs 13-XK 0.36930162 0.012160359   0.722      0.722
# 9  13-X9 vs 13-X5 0.78230468 0.023863627   0.443      0.443
# 10 13-XK vs 13-X5 0.21350716 0.006627877   0.852      0.852

# do stats within T24 Site
bet.psT24site <- betadisper(dist.psT24,samdf.psT24$Site)
anova(bet.psT24site) #ns
permutest(bet.psT24site,pairwise=TRUE,permutations=999)
plot(bet.psT24site)
adonis(seq.psT24 ~ Site, data=samdf.psT24, permutations=999) #p<031*
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# Site      7    1.1984 0.171199  2.0188 0.14855  0.031 *
# Residuals 81    6.8690 0.084803         0.85145         
# Total     88    8.0674                  1.00000         
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

pairwise.adonis(seq.psT24, factors=samdf.psT24$Site, permutations=999) # only DL vs. ES, MS, and WS as well as EDR vs. WS as well as LK vs. WS are different
These are all the A sites versus the D dominated sites.
#             pairs    F.Model          R2 p.value p.adjusted
# 1              BigPine vs DavesLedge 1.90732724 0.106510995   0.154      0.154
# 2         BigPine vs EasternDryRocks 0.56089835 0.041361445   0.754      0.754
# 3            BigPine vs EasternSambo 0.27953602 0.016177287   0.760      0.760
# 4                 BigPine vs LooeKey 0.68986443 0.041334334   0.615      0.615
# 5                BigPine vs Marker32 0.15755551 0.014120970   0.936      0.936
# 6          BigPine vs MarylandShoals 0.01629302 0.001017278   0.955      0.955
# 7            BigPine vs WesternSambo 0.75431663 0.045022226   0.426      0.426
# 8      DavesLedge vs EasternDryRocks 1.00034935 0.045469703   0.329      0.329
# 9         DavesLedge vs EasternSambo 4.60597334 0.155575812   0.025      0.025
# 10             DavesLedge vs LooeKey 1.25414414 0.049660924   0.216      0.216
# 11            DavesLedge vs Marker32 3.14584848 0.142051386   0.073      0.073
# 12      DavesLedge vs MarylandShoals 4.64029706 0.162019865   0.016      0.016
# 13        DavesLedge vs WesternSambo 5.88777097 0.196995988   0.001      0.001
# 14   EasternDryRocks vs EasternSambo 1.17138902 0.050553250   0.303      0.303
# 15        EasternDryRocks vs LooeKey 0.09843059 0.004665304   0.960      0.960
# 16       EasternDryRocks vs Marker32 0.90906622 0.053762059   0.495      0.495
# 17 EasternDryRocks vs MarylandShoals 1.32988422 0.059556252   0.219      0.219
# 18   EasternDryRocks vs WesternSambo 2.23726597 0.096279226   0.044      0.044
# 19           EasternSambo vs LooeKey 1.13502590 0.043429301   0.447      0.447
# 20          EasternSambo vs Marker32 0.43999342 0.021526104   0.624      0.624
# 21    EasternSambo vs MarylandShoals 1.13821590 0.043546044   0.264      0.264
# 22      EasternSambo vs WesternSambo 2.74710052 0.099004958   0.073      0.073
# 23               LooeKey vs Marker32 1.13983672 0.056596125   0.262      0.262
# 24         LooeKey vs MarylandShoals 1.74039464 0.067613363   0.106      0.106
# 25           LooeKey vs WesternSambo 2.86656191 0.106696269   0.013      0.013
# 26        Marker32 vs MarylandShoals 0.25258610 0.013119593   0.933      0.933
# 27          Marker32 vs WesternSambo 0.85520636 0.043072147   0.392      0.392
# 28    MarylandShoals vs WesternSambo 0.91169055 0.036596896   0.367      0.367




