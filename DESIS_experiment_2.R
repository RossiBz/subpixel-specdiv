

##################################################

# Project: Uncovering the hidden: Leveraging sub-pixel spectral diversity to estimate plant diversity from space

# Script Purpose: Calculate endmember diversity for spaceborne hyperspectral data, as done in Rossi and Gholizadeh (2023) Experiment 2

# Date: 01.08.2022

# Author: Christian Rossi

##################################################

#load packages

#library(devtools)
#install_github("RossiBz/unmix")
library(raster)
library(unmix)
library(sf)

#load data------------------------------------------------------------------

#load preprocessed DESIS data
data_desis <-
  stack(
    "DESIS-HSI-L2A-DT0618506904_016-20210806T174855-V0213-SPECTRAL_IMAGE_reduced_187bands_mask_smooth.dat"
  )
data_desis <- stack(calc(
  data_desis,
  fun = function(x) {
    x / 10000
  }
))

#load polygon with the nine equally-sized non-overlapping subimages
gridPolygon <-
  st_read("grid_tiles/gridpolygon_desis.shp")

#load mask including homogeneous regions based on the Sobel-Felman operator
mask_sd <- raster("NDVI_edges_desis_mask.tif")
mask_sd[is.na(mask_sd)] <- 0

#add info on homogeneous area to the desis data
data_desis <- stack(data_desis, mask_sd)


#findining the number of endmembers in the DESIS image using the noise-whitened Harsanyi–Farrand–Chang method  ----------------

num_end <- unmix::hfcvd(data_desis[[-dim(data_desis)[3]]])


#extracting the endmembers from each subimages with the pixel purity index -----------------------------------

endmembers_ppi <- NULL

for (l in 1:nrow(gridPolygon))
{
  tile <- crop(data_desis, gridPolygon[l, ])
  tile <- raster::mask(tile, gridPolygon[l, ])
  tile[tile < 0] <- NA
  
  matrix_tile <- t(na.omit(as.matrix(tile)))
  
  #endmember extraction
  endmembers_ppi_prov <-
    unmix::ppi(matrix_tile[-nrow(matrix_tile), ], q = num_end, reductionmethod =
                 "MNF")
  
  
  #check if endmembers area in homogenous area
  index_pos <- sapply(1:ncol(endmembers_ppi_prov), function(z) {
    which(apply(matrix_tile[-dim(data_desis)[3], ], 2, function(x)
      return(all(
        x == endmembers_ppi_prov[, z]
      ))))
  })
  
  # if not they are removed
  endmembers_ppi <- cbind(endmembers_ppi,
                          endmembers_ppi_prov[, as.logical(matrix_tile[nrow(matrix_tile), index_pos])])
  
  
}

#averaging spectrally similar endmembers ------------------------------------------------
matrix = endmembers_ppi
range = max(matrix) - min(matrix)
threshold = 0.02

# Get the number of rows and columns in the matrix
m = nrow(matrix)
n = ncol(matrix)

# Iterate through the columns
for (i in 1:n) {
  # Initialize a variable to store the indices of similar measurements
  similar = c()
  # Iterate through the remaining columns
  if (i < n)
  {
    for (j in (i + 1):n) {
      # Compare the measurements
      diff = sqrt(mean((matrix[, i] - matrix[, j]) ^ 2)) / range
      if (diff < threshold) {
        # Add the column indices to the similar measurements
        similar = c(similar, j)
      }
    }
  }
  # Average the similar measurements
  if (length(similar) > 0) {
    matrix[, i] = rowMeans(matrix[, c(i, similar)], na.rm = TRUE)
    # Remove the similar measurements
    matrix = matrix[,-similar]
    n = n - length(similar)
  }
}

endmember_ppi_red <- matrix

#plot endmembers
matplot(matrix, type = "l")

#spectral unmixing ------------------------------------------------------

end_abund <-
  estimateabundanceLS(data_desis[[-dim(data_desis)[3]]], matrix, method =
                        "fcls")

#disregarded endmembers with abundance values of less than 1%
end_abund[end_abund < 0.01 | is.na(end_abund)] <- 0


#endmember diversity calculation for each pixel-------------------------

div_endmembers <- unmix::endDiv(end_abund)
