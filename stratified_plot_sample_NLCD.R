# this script generates a stratified random sample of points for each NLCD class
# output is a raster image clipped to the OSBS extent (square) and a list of random x, y locations and the class

# IMPROVEMENTS TO MAKE
# contrain points to be within property boundary (need to get one), or HS image extent rather than a crude OSBS extent
# ensure points fall near the center of flight paths
# adjust the number of random points to be depended on the size of the class (proportional stratification)

library(raster)
library(rgdal)

# read in image file
r <- raster("C:/Users/Sarah/Downloads/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img")

extent(r)
crs(r)

p <- readOGR(dsn="C:/Users/Sarah/Downloads/nlcd_2011_landcover_2011_edition_2014_10_10",
             layer="OSBS_extent")

extent(p)
crs(r)
plot(p)

rc <- crop(r,p)

plot(rc)

# where are the data stored?
rc@data@values <- rc[]
hist(rc)

# save clipped raster
#writeRaster(rc,filename = "../Downloads/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10_OSBSclip_reclass.tif",format="GTiff")

# clean up environment
r <- rc
rm(rc,p)

# want to loop through all classes and create random point

# class list
r_classes <- sort(unique(r[]))
my_classes <- c("41","42","43","52","90","95")
class_match <- r_classes[r_classes %in% my_classes]

# select random pixels
nsample <- 10

# for each class, subset raster, create raster with random value for these pixels, select random pixels and extract coordinates
# store coordinates
rand_coords <- NA

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
  coords_i <- rasterToPoints(r_class_i_rand, function(x) {x <= nsample})
  
  rand_coords <- rbind(rand_coords,coords_i)
  
}

# remove first row
rand_coords <- rand_coords[-1,]

# change final column to be the class
rand_coords[,3] <- sort(rep(class_match,nsample))
colnames(rand_coords) <- c("x","y","class")

# save
write.csv(rand_coords,"../../data/NEON_NLCD_data/stratified_random_points.csv",row.names = F)
