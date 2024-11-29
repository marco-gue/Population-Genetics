# Full script to analyse selections signatures based on variation of linkage disequilibrium
# Software: varLD
# Presented by: Gabor Meszaros, see Genomics Boot Camp YouTube channel

# Workflow:

# 1) Install Java
# 2) Aquire the data for two populations, for example: http://widde.toulouse.inra.fr/widde/widde/main.do?module=cattle
# 3) Aquire the program: https://blog.nus.edu.sg/sshsphphg/varld/
#    Link in the publication does not work
# 4) Check the manual, i.e. the README for instructions
# 5) Follow the script below to figure out how it is done


# Clear workspace and load packages
rm(list = ls())
library(tidyverse)

# Set the location of the working directory
PATH <- "d:/analysis/2020_GenomicsBootCamp_Demo/varLD/"
setwd(PATH)
# Set the number of autosomes for the loops
autosomeNr <- 29


####
# Quality control 
####
# breed 1
system("plink --file Charolais --cow --autosome --geno 0.1 --mind 0.1 --hwe 0.000001 --nonfounders --allow-no-sex --make-bed --out breed1_QC")

# breed 2 
system("plink --file Holstein --cow --autosome --geno 0.1 --mind 0.1 --hwe 0.000001 --nonfounders --allow-no-sex --make-bed --out breed2_QC")


####
# Find the SNP appearing in both populations 
####
map1 <- read_tsv("breed1_QC.bim",col_names = F)
map2 <- read_tsv("breed2_QC.bim",col_names = F)

mapOverlap <- inner_join(map1, map2, by="X2") %>%
  select(X2) %>%
  write_tsv("commonSnpPop1Pop2.txt", col_names = F)

# START MAIN LOOP
for (chrNr in 1:autosomeNr) {
  
  #chrNr <- 1
  
  ####
  # Population 1 
  ####
  
  # change PLINK file to format close to varLD requirement: A-transpose
  system(str_c("plink --bfile breed1_QC --cow --allow-no-sex --nonfounders --chr ", chrNr ,
               " --extract commonSnpPop1Pop2.txt --recode A-transpose --out pop1_forVarLD"))
  
  
  # read in and modify the traw file to varLD format
  traw <- read_tsv("pop1_forVarLD.traw")
  
  # add 1 to all genotypes, to meet the requirements of varLD for AA=1, AB=2, BB=3 - missing values remain NA
  traw[,7:ncol(traw)] <- traw[,7:ncol(traw)]+1
  
  # replace missing NA values with 4, according to varLD requirements and keep only selected columns
  # very unfortunate column names, therefore the index references with numbers
  pop1_VarLD <- traw %>%
    mutate(across(starts_with("CHL_"), ~replace_na(.x, 4))) %>%
    select(-c(1,3,5,6)) %>%
    write_tsv(paste0("pop1_chr",chrNr,"_varLD.txt"), col_names = F)
  
  
  ####
  # Population 2
  ####
  
  # change PLINK file to format close to varLD requirement: A-transpose
  system(str_c("plink --bfile breed2_QC --cow --allow-no-sex --nonfounders --chr ", chrNr ,
               " --extract commonSnpPop1Pop2.txt --recode A-transpose --out pop2_forVarLD"))
  
  
  # read in and modify the traw file to varLD format
  traw <- read_tsv("pop2_forVarLD.traw")
  
  # add 1 to all genotypes, to meet the requirements of varLD for AA=1, AB=2, BB=3 - missing values remain NA
  traw[,7:ncol(traw)] <- traw[,7:ncol(traw)]+1
  
  # replace missing NA values with 4, according to varLD requirements and keep only selected columns
  # very unfortunate column names, therefore the index references with numbers
  pop2_VarLD <- traw %>%
    mutate(across(starts_with("HOL_"), ~replace_na(.x, 4))) %>%
    select(-c(1,3,5,6)) %>%
    write_tsv(paste0("pop2_chr",chrNr,"_varLD.txt"), col_names = F)
  
  
  ####
  # run varLD
  ####
  system(paste0("java -jar rgenetics-1.0.jar -p VarLD pop1_chr",chrNr,"_varLD.txt pop2_chr",chrNr,"_varLD.txt -o varLD_chr",chrNr,".txt"))
  
  ####
  # Change decimal comma to decimal point - regional setting Windows
  # Use if necessary...
  ####
  #outVarLD <- read.table(str_c("varLD_chr",chrNr,".txt"), sep="\t", header = T)
  #outVarLD$position <- gsub(',', '.', outVarLD$position)
  #outVarLD$raw_score <- gsub(',', '.', outVarLD$raw_score)
  #write.table(outVarLD,str_c("pointvarLD_chr",chrNr,".txt"), sep="\t", quote = F)
  
  # END MAIN LOOP
}


####
# Standardize varLD values
####


# This software is supplied without any warranty or guaranteed support whatsoever.
# NUS CME can not be responsible for its use, misuse, or functionality.
#
# This R script can be used for standardizing the genome-wide scores to have 
# a mean of 0 and a variance of 1, and to output files with the standardized scores. 
# It assumes that varld scores are available for all 22 autosomal chromosomes. 
# This code also identifies the varLD scores that correspond to the stated percentiles.
#
# Author: rick and yy teo
# Version: 1.0
# Last Updated : 07 March 2010
######################################################################################
## global variables to modify for own use 
# specifies the folder containing all varld output for 22 autosomal chromosomes
#PATH = "XXXXXXXX"
# filename prefix before the chromosome number, i.e. CHS_INS_chr1, CHS_INS_chr2 
FILENAME = "varLD_chr" 
# percentile of the genomewide distribution to highlight, default to 95%, 99%, 99.9% and 99.99%.
percentile.out = c(0.95, 0.99, 0.999, 0.9999)


varLD.out <- {}
chr.store <- {}
for (chr in 1:autosomeNr){
  varLD.temp <- read.table(paste(PATH, FILENAME, chr, ".txt", sep=""), sep="\t", header = T)
  varLD.out <- rbind(varLD.out, varLD.temp)
  chr.store <- c(chr.store, rep(chr, dim(varLD.temp)[1]))
  print(paste("completed reading in unstandardized varLD output file for chromosome ", chr, sep=""))
}
varLD.mean <- mean(varLD.out[,"raw_score"])
varLD.sd <- sd(varLD.out[,"raw_score"])
standardized_score <- (varLD.out[,"raw_score"] - varLD.mean)/varLD.sd
varLD.out <- cbind(varLD.out, standardized_score)
varLD.threshold <- quantile(standardized_score, probs = percentile.out)

for (chr in 1:autosomeNr){
  chr.flag <- which(chr.store == chr)
  write.table(varLD.out[chr.flag,], paste(PATH, FILENAME, chr, "_standardized.out", sep=""), sep="\t", quote=F, row.names=F)
  print(paste("completed writing out standardized varLD output file for chromosome ", chr, sep=""))
}
n.length.percentile <- length(percentile.out)
for (i in 1:n.length.percentile){
  print(paste("varLD threshold for ", percentile.out[i], " = ", varLD.threshold[i], sep=""))
} 





####
# plot VarLD values
####
# This software is supplied without any warranty or guaranteed support whatsoever.
# NUS CME can not be responsible for its use, misuse, or functionality.
#
# This is example code for producing the varLD region plot of VKORC1 that is similar to 
# Figure 1 in the Bioinformatics Application Note. 
# Anybody who is moderately proficient in R should be able to modify the codes here
# to customise the plot required. 
#
# Author: rick and yy teo
# Version: 1.0
# Last Updated : 07 March 2010
######################################################################################
# Loop for all chromosomes
for (chrPlot in 1:autosomeNr) {  
  ## global variables to modify for own use 
  # specifies the folder where the output files from varLD is stored
  #PATH = "XXXXXXXX"
  # filename prefix before the chromosome number, i.e. CHS_INS_chr1, CHS_INS_chr2
  FILENAME = "varLD_chr" 
  # chromosome of region to plot
  chr = chrPlot
  # percentile of the genomewide distribution to highlight 
  percentile.out 	= c(0.95, 0.99, 0.999, 0.9999)
  # varLD score corresponding to the stated percentiles in percentile.out
  # !!! skip this, as already computed by the previous script
  #varLD.threshold 	= c(1.8, 3.4, 6.0, 7.0)
  # start and end base-pair coordinates of the region to plot
  region.start = 0
  region.end = 320000000
  
  
  temp.in <- read.table(paste(PATH, FILENAME, chr, "_standardized.out", sep=""), sep="\t", header = T)
  region.flag <- which(temp.in[,"position"] >= region.start & temp.in[,"position"] <= region.end)
  
  plot(0, 1, xlab = paste0("Physical position (Mb) - ",chrPlot), ylab = "Standardized score", las = 1, type = "n", xlim = c(region.start, region.end)/10^6, ylim = c(min(varLD.threshold), max(temp.in[region.flag,"standardized_score"])), bty = "n", cex.lab = 1.4, cex.main = 1.4)
  points(temp.in[region.flag, "position"]/10^6, temp.in[region.flag, "standardized_score"], pch = 16, col = "red")
  n.length.threshold <- length(varLD.threshold)
  if (n.length.threshold >= 2){
    for (i in 2:n.length.threshold){
      abline(h = varLD.threshold[i], lty = 2)
      text(region.start/10^6 + 0.2, varLD.threshold[i] + 0.15, labels = paste("Top ", round((1 - percentile.out[i]) * 100, 2), "%", sep=""))
    }
  }
  
  
}



####
# Full manhattan plot
####
# read all files and join to one big file
fullData <- NULL
for (tmpCounter in 1:autosomeNr){
  tmp <- read_tsv(str_c("varLD_chr",tmpCounter,"_standardized.out"), col_names = T)
  tmp$chr <- tmpCounter
  fullData <- rbind(fullData,tmp)
  
}

# select only the required columns
forManhattan <- fullData %>%
  select(chr, position, pop1, standardized_score)

# draw the manhattan plot with qqman
library(qqman)
manhattan(forManhattan, chr = "chr", bp = "position", p = "standardized_score", snp = "pop1", logp = F,
          xlab = "Chromosomes", ylab="Standardized score", suggestiveline=varLD.threshold[3],
          genomewideline = varLD.threshold[4])
