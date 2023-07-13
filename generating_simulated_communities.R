


##################################################

# Project: Uncovering the hidden: Leveraging sub-pixel spectral diversity to estimate plant diversity from space

# Script Purpose: Simulates plant communities using 17 abundant species, 2 soil condition and their spectral signal

# Date: 21.02.2023

# Author: Christian Rossi

##################################################



library(tidyverse)
library(stringr)
library(vegan)
library(phytools)
library(V.PhyloMaker)



#load spectral signature of 17 dominant grasses and two soil spectra---------
endmember_grass <- read.table(
  "gis_work/spectral_libraries/TGPP_insitu_average_spectral_signatures.csv",
  sep = ",",
  header = TRUE
)


#load species data to generate phylogentetic tree ---------------------

species.df <- read.table(
  "code/input/TGPP_2021_AbunMatrix.csv",
  #includes 253 species
  sep = ";",
  quote = "\"*",
  header = TRUE,
  fileEncoding = "UTF-8-BOM"
) %>% mutate(species_names = str_trim(paste(genus, species, sep = " "))) %>%
  dplyr::select(common.name:species) %>%
  dplyr::filter(!row_number() %in% c(254, 255))  %>% #remove soil and rock
  mutate(species = str_trim(paste(genus, species))) %>%
  dplyr::select(-common.name) %>% relocate(species) %>% relocate(family, .after = genus) %>%
  dplyr::filter(species  %in% endmember_grass$Endmember) #select only the 17 dominant species used for the simulation


#phylogenetic tree generation--------------------------------

phylo.tree <-
  phylo.maker(
    sp.list = species.df,
    tree = GBOTB.extended,
    nodes = nodes.info.1,
    scenarios = "S1"
  )$scenario.1

phylo.tree$tip.label <- gsub("_", " ", phylo.tree$tip.label)


#function to compute simpson index without normalization----------------

simpson_nonorm <- function(x, index = "invsimpson")
{
  x <- x * x
  if (length(dim(x)) > 1)
    H <- apply(x, 1, sum, na.rm = TRUE)
  else
    H <- sum(x, na.rm = TRUE)
  if (index == "simpson")
    H <- 1 - H
  else if (index == "invsimpson")
    H <- 1 / H
  
  H
}




#generate simulated communities with dominant species and two soil spectra (at least three per community) -------------------------


info_simu <- NULL

for (s in 3:19)
{
  for (i in 1:900)
  {
    x <- runif(s, 0.03, 1)
    abundance <- data.frame(t(x / sum(x)))
    
    #sample s random species and soil
    random_species <- endmember_grass %>% sample_n(s)
    
    #mix spectra based on abundance
    spectrum_mixed <-
      colSums(random_species[, 2:length(random_species)] * t(abundance))
    
    colnames(abundance) <- random_species$Endmember
    
    #diversity
    
    empty_species <- which(colSums(abundance) == 0)
    
    if (length(empty_species > 0))
    {
      abundance <- abundance[,-empty_species]
    }
    
    # remove soil for diversity calculation
    
    if (any(colnames(abundance) %in% "Soil_1"))
    {
      abundance <- abundance %>% dplyr::select(-Soil_1)
    }
    
    if (any(colnames(abundance) %in% "Soil_2"))
    {
      abundance <- abundance %>% dplyr::select(-Soil_2)
    }
    
    soil_abundance <- 1 - sum(abundance)
    
    
    if (length(abundance) > 0)
    {
      #taxonomic diversity
      Richness <- specnumber(abundance)
      Simpson_no <- simpson_nonorm(abundance, index = "invsimpson")
      
      
    } else{
      Richness <- specnumber(abundance)
      Simpson_no <- NA
      
    }
    
    
    if (length(abundance) > 1)
    {
      #phylogenetic diversity
      phylo.species.rich <- psr(abundance, phylo.tree)$PSR
      phylo.species.eve <- pse(abundance, phylo.tree)$PSEs
    } else{
      phylo.species.rich <- NA
      phylo.species.eve <- NA
      
      
    }
    
    
    info_simu <-
      info_simu %>% rbind(
        data.frame(
          soil_abundance,
          Richness,
          Simpson_no,
          phylo.species.rich,
          phylo.species.eve,
          t(spectrum_mixed)
        )
      )
    
  }
  
  
}



#add noise
dim_data_spec <- dim(info_simu[, 6:2006])
B <- dim_data_spec[2]
N <- dim_data_spec[1]

SNR <- 60 # dB
variance <- sum(info_simu[, 6:2006] ^ 2) / (10 ^ (SNR / 10)) / N / B
n <- matrix(rnorm(B * N, mean = 0, sd = sqrt(variance)), ncol = N)
comm_data_spec_n <- t(info_simu[, 6:2006]) + n

info_simu_noise <- cbind(info_simu[, 1:5], t(comm_data_spec_n))

write.table(
  info_simu_noise,
  "simulated_communities_field_red100_soil2_3speciesup03_SNR60.csv",
  row.names = FALSE,
  sep = ";"
)
