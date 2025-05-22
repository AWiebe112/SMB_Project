####SMB_DataVis_2025-05-22###

#Load packages
library(phyloseq)
library(vegan)
library(tidyverse)
#library(tidyr)

#switch working directory to where you want the file saved
setwd("C:/Users/allyw/Desktop/UBC_Master_Zoology/SMB/SMB_Data/SMB2024_PhyloseqObjects/")

#Save the final phyloseq object as an RDS file
nonrare<-readRDS("SMB2024_phyloseq.RDS")

###INVESTIGATE PHYOSEQ OBJECT
## View the 3 different parts of the phyloseq object (use @ to view specific tables)
# metadata
View(nonrare@sam_data)
# otu table
View(nonrare@otu_table)
# view taxonomy
View(nonrare@tax_table)
#view phy_tree
View(nonrare@phy_tree)


###FILTERING 
# remove off target taxa --> already done in the sockeye deblur script
#nonrare = subset_taxa(nonrare,
                      # domain!="Unassigned"&
                      #   species!="Chloroplast" &
                      #   species!="Mitochondria" &
                      #   domain!="Eukaryota")

## look at your taxonomy table to check everyhting got removed with View()
#after first round of filtering there were still mitochondrial categories in rank 5 so added additional restriction to subset function

## print the sample sums (number of reads) of each sample in the console
sample_sums(nonrare)

## make a plot of the sample sums
plot(sort(sample_sums(nonrare)))

## add total number of reads in each sample to the metadata
nonrare@sam_data$sample_sums_unfiltered = as.numeric(sample_sums(nonrare))

##Invertigate min and max reads per sample
min(nonrare@sam_data$sample_sums_unfiltered)
#Output:20
max(nonrare@sam_data$sample_sums_unfiltered)
#Output:79699

## remove the samples by your threshold (equal or greater to 13000)
nonrare.high <- prune_samples(sample_sums(nonrare) >= 13000, nonrare)

## make another object with only the samples you lost (less than 13000 reads)
nonrare.below6000 <- prune_samples(sample_sums(nonrare) < 13000, nonrare)

## get the metadata out of phyloseq from your low reads object --> 3 samples below cutoff (blank, NC1, SMB124)
nonrare.below6000 = as.matrix(nonrare.below6000@sam_data)

## write file to know which samples were lost here. This is important for the methods section.
write.csv(nonrare.below6000, "SMB2024_samples_less_than_13000.csv")


## extract otu dataframe (asv table) from phyloseq object
# notice the t before as matrix here! This is pivoting the otu table 90 degrees.
otutab <- as.data.frame((as.matrix(otu_table(nonrare.high@otu_table))))

## check your asvtab dataframe visually to see if your ASVs are the ROWS an the sample names are the COLUMNS.
# This is critical otherwise you will not filter to correct thing.
## calculate the sum of each row in your otu table
otutab$asv_abundance = rowSums(otutab)

## get the minimum asv_abundance value
min(otutab$asv_abundance)
#Output: 0, low abundance ASVs indicates polymerase mistake --> remove low ab.

## remove low frequency asvs
otu.pruned = subset(otutab, otutab$asv_abundance>=100)

## now get the minimum asv abundance from your cleaned up otu table.
# this should be at or above your filtering treshold.
min(otu.pruned$asv_abundance)
#Output: 100

## remove asv_abundance column from your OTU table since we don't want to analyse it
# our asv_abundance column gets tacked onto the end when you make it, so you just need to delete the last ## how many columns in dataframe?
widthotu = ncol(otu.pruned)
#Output: 125 columns in otu.prume
#Output after removing last column: 124 columns in otu.prume

## see how widthotu appears as a Value in the environement?
## keep everything in the otu.pruned dataset except the last columns
otu.pruned = otu.pruned[,-c(widthotu)] 
#run line 89 again to confirm that a column was removed


## create a function to count the number of occurence along rows where the number in the cell (ASV occurence)
ASVoccur = function(x){return(sum(x>0))}

## calculate the occurence of each ASV in your dataframe
otu.pruned$asv_occur_count = apply(otu.pruned, 1, ASVoccur)

## see what the ASV occurence counts look like
summary(otu.pruned$asv_occur_count)
# #Output:
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    7.00   12.00   17.07   22.00  116.00 

## lets remove ASVs found two or less times in the dataset
otu.highfreq = subset(otu.pruned, otu.pruned$asv_occur_count>2)

## see if the filtering worked.The minimum should now be above your threshold (2 in this case)
summary(otu.highfreq$asv_occur_count)
# #Output: 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3.00    8.00   13.00   18.05   23.00  116.00 

## remove the asv_occur_count column
otu.highfreq = otu.highfreq[,-c(widthotu)]


##De-noising the data
otu.clean <- mutate_all(otu.highfreq, funs(ifelse(. < 5, 0, .)))


## make the cleaned phyloseq object
cleanSMB2024 <- phyloseq(sample_data(nonrare.high@sam_data),
                         tax_table(nonrare.high@tax_table),
                         otu_table(as.matrix(otu.clean), taxa_are_rows = TRUE),
                         nonrare.high@phy_tree)  # Directly pass the phylogenetic tree


## basic check to see that the object was created
cleanSMB2024

## add the number of reads after filtering
cleanSMB2024@sam_data$sample_sums_filtered = sample_sums(cleanSMB2024)

## save clean filtered data as RDS 
write_rds(cleanSMB2024, "filtered_SMB2024_default_phyloseq.RDS")


###RAREFACTION
## use the rarecurve function in the package vegan to plot the rarefaction
## there are a lot of things going on here because you need your SampleID to be the row name and the ASV ## If you open the otu.clean dataframe, you will see that it's not in that format, so we need to fix this.
## R works like math, you read the inner most () first and work outwards
rarecurve(
  as.data.frame( ## 3. Turn the matrix back into a dataframe
    ( ## 2. turn the otu.clean matrix 90 degrees to have the correct orientation
      as.matrix( ## 1. turn the otu.clean dataframe into a matrix
        otu.clean))), ## 0. the dataframe we are doing stuff to
  step=50, cex=0.5, label=FALSE) ## now that the data are formatted, sample the reads (these are just parameters


## open the cleanseagrass metadata
View(cleanSMB2024@sam_data)

## set seed to tell R which set of random numbers to use.
# This is important because it will allow you to sample randomly the same way every time
set.seed(5)
## rarefy every sample to a set number of reads here
raresg <- rarefy_even_depth(cleanSMB2024, sample.size = 11988)
  #this is the lowest number of read a sample has 11988 --> I will keep everything to start

## calculate the rarefied sample sums
raresg@sam_data$rare_sample_sums = sample_sums(raresg)
## see distribution of rare_sample_sums
summary(raresg@sam_data$rare_sample_sums)
#Output: 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 11988   11988   11988   11988   11988   11988 

## save rarefied dataframe
write_rds(raresg, "rarefied_SMB2024_default_phyloseq.RDS")



