####SMB2024_roughVisualization_2025-05-22###

#Load packages
library(phyloseq)
library(ggplot2)
library(pairwiseAdonis)
library(vegan)

#set working directory 
setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/SMB/SMB_Data/SMB2024_PhyloseqObjects/")

#read in data
SMB2024rare<- readRDS("rarefied_SMB2024_default_phyloseq.RDS")

#remove positive controls and one batilaria sample that clustered with the positive controls because they squew the ordination
noPC1 = subset_samples(SMB2024rare, sampleID != c("SMBPC1")
noPC2= subset_samples(noPC1, sampleID != "SMBPC2")
noPC3= subset_samples(noPC2, sampleID != "SMBPC3")
noPC4 = subset_samples(noPC3, sampleID != "SMBPC4")


#Need to figure out what bat sample it is
# Extract NMDS site scores
site_scores <- as.data.frame(scores(SMB.ord, display = "sites"))
site_scores$SampleID <- rownames(site_scores) # Add sample IDs
outlier_sample <- site_scores[site_scores$NMDS1 > 40, ]# Find outlier sample with NMDS1 > 40
print(outlier_sample) # Print the problematic sample(s)
# NMDS1     NMDS2 SampleID
# SMB041 46.16863 0.4954586   SMB041

#remove this sample
noPC_rare = subset_samples(noPC4, sampleID != "SMB041")

##Make NMDS
SMB.ord = ordinate(noPC_rare, "NMDS","bray") #ordinate data for NMDS
stressplot(SMB.ord) ## see what the fit looks like
SMB.ord #see what the stress is --> 0.163718
View(noPC_rare@sam_data) # metadata

# make the NMDS plot colour by species and site by shape
plot_ordination(noPC_rare, SMB.ord, color="species", shape = "location")+
  geom_point(size=3) + theme_classic()

## make the NMDS plot colour by species and facetwrap by site
plot_ordination(noPC_rare, SMB.ord, color="species")+
  geom_point(size=3)+facet_wrap("location")



#######Subset just the Littorina scutulata###################################
scut = subset_samples(SMB2024rare, species == "L. scutulata")

##Make NMDS
scut.ord = ordinate(scut, "NMDS","bray") #ordinate data for NMDS
stressplot(scut.ord)## see what the fit looks like
scut.ord#see what the stress is --> 0.1836462 

# metadata
View(scut@sam_data)

## make the NMDS plot colour by species and site by shape
plot_ordination(scut, SMB.ord, color="location", shape = "timepoint") + theme_classic() 




#######Subset just the Quadra###################################
qb = subset_samples(SMB2024rare, location == "QB")

##Make NMDS
qb.ord = ordinate(qb, "NMDS","bray") #ordinate data for NMDS
stressplot(qb.ord)## see what the fit looks like
scut.ord#see what the stress is --> 0.1836462 

# metadata
View(qb@sam_data)

## make the NMDS plot colour by species and site by shape
plot_ordination(qb, SMB.ord, color="species", shape = "timepoint") + theme_classic() +geom_point(size = 3)

