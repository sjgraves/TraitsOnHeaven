# script to read .asd files
# Sarah Graves, November 15, 2016

# files stored in separate folder that is not part of the repo
# in total there are 20,005 files at 700 MB

library(dplyr)

# set reflectance file names
asd_files <- list.files("../../data/dimensions_field_data/2015/OSBS/spectra/sample_refl/")

# NEXT STEP
# want to average across all files for 1 leaf
# process is; open file, store in matrix, calculate average, write average to file
# need to somehow know which files to group together

