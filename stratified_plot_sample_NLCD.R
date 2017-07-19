# this script generates a stratified random sample of points for each NLCD class
# output is a raster image clipped to the OSBS extent (square) and a list of random x, y locations and the class
# NLCD class legend: https://www.mrlc.gov/nlcd11_leg.php

# IMPROVEMENTS TO MAKE
# contrain points to be within property boundary (need to get one), or HS image extent rather than a crude OSBS extent - NOW WITHIN HS FLIGHT PATHS, BUT BEYOND OSBS BOUNDARY
# ensure points fall near the center of flight paths
# adjust the number of random points to be depended on the size of the class (proportional stratification) - DONE

library(raster)
library(rgdal)

if(file.exists("../../data/NEON_NLCD_data/NLCD_raster_subsets/NLCD_2011_OSBS_2014.tif")==T){
  
  r <- raster("../../data/NEON_NLCD_data/NLCD_raster_subsets/NLCD_2011_OSBS_2014.tif")
  
} else{
  
  # read in image file
  r <- raster("C:/Users/Sarah/Downloads/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")
  
  extent(r)
  crs(r)
  
  p <- readOGR(dsn="../../data/NEON_AOP_supplementary_data",layer="OSBS_2014_flight_boundary_NAD83")
  
  extent(p)
  crs(r)
  plot(r)
  plot(p,add=T)
  
  # first crop (more simple funcitn to restrict area by the extent)
  rc <- crop(r,p)
  
  # now mask by polygon
  rcm <- mask(rc,p)
  
  summary(rcm)
  
  
  plot(rcm)
  
  r <- rcm
  
  rm(rc,rcm)
  
  # where are the data stored?
  r@data@values <- r[]
  
  # save clipped raster
  writeRaster(r,filename = "../../data/NEON_NLCD_data/NLCD_raster_subsets/NLCD_2011_OSBS_2014.tif",format="GTiff")
  
}

plot(r)

# want to loop through all classes and create random point

# class list
r_classes <- sort(unique(r[]))
my_classes <- c("41","42","43","52","90")
class_match <- r_classes[r_classes %in% my_classes]

# calculate proportion for each class in class_match
class_pix <- table(r[])
class_pix_match <- class_pix[names(class_pix) %in% class_match]
class_prop <- class_pix_match/sum(class_pix_match)

# select random pixels
nsample <- 30

# calculate number of samples
# how many points to sample?
class_sample <- round(class_prop*nsample,0)



# for each class, subset raster, create raster with random value for these pixels, select random pixels and extract coordinates
# store coordinates
rand_coords <- NA

set.seed(19)
for(i in 1:length(class_match)){
  
  # mask raster - all nonclass values turn to NA
  r_class_i <- mask(r, mask = r, maskvalue = class_match[i], inverse=T)
  plot(r_class_i)

  # generate random values in non NA pixels for this subset
  pix_i <- length(r[r==class_match[i]])
  
  # for the subset raster, only with values equal to the class
  # assign a random value - non-repeating, all values are unique from 1 to the number of pixels in this class
  r_class_i_rand <- r_class_i
  r_class_i_rand[r_class_i_rand==class_match[i]] <- sample.int(n=pix_i)

  # get coordinates of pixels with these values
  coords_i <- rasterToPoints(r_class_i_rand, function(x) {x <= class_sample[i]})
  
  rand_coords <- rbind(rand_coords,coords_i)
  
}

# remove first row
rand_coords <- rand_coords[-1,]

# change final column to be the class
tmp <- data.frame(class_sample,stringsAsFactors = F)
rand_coords[,3] <- as.character(rep(tmp$Var1,tmp$Freq))
colnames(rand_coords) <- c("x","y","class")

# save
write.csv(rand_coords,"../../data/NEON_NLCD_data/stratified_random_points_NAD83.csv",row.names = F)

# NOW IN A GIS DO THE FOLLOWING
# select which points you want
# add attribute for which flight path the points fall
