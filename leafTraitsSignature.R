# script for calibrating the leaf lv spectr PLSR models and ncertainty estimates.
#--------------------------------------------------------------------------------------------------#
# Close all devices and delete all variables.
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- Load required libraries ---------------------------------------------------------#
# Info: Loads required R libraries and warns if package is not availible.
ok = require(pls) ; if (! ok) 
  stop("*** Package pls is not available.  This is needed for model optimization ***")

ok = require(plotrix) ; if (! ok) 
  stop("*** Package plotrix is not available.  This is needed visualization of results ***")

ok = require(prospectr) ; if (! ok) 
  stop("*** Package prospectr is not available.  This is needed to perform standardNormalVariate ***")

# Script options
pls.options(plsralg = "oscorespls")
pls.options("plsralg")
rm(ok)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Set ouput directory
out.dir = paste(getwd(), "/Outputs/LMA/", sep="") 
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Import dry spectra dataset
spec.dir <- paste(getwd(), "/Inputs/", sep="") 
dry.spectra <- read.table(paste(spec.dir,'D17_ASD_Chem.csv',sep=""), header=TRUE,sep=",")
spec <- dry.spectra[,11:dims[2]]
head(dry.spectra)
rm(spec.dir)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Generate diagnostic spectra figure?

Create.fig <- FALSE

if (Create.fig) {
  ### Create dry spec diagnostic figure
  dims <- dim(dry.spectra)
  spec <- dry.spectra[,11:dims[2]]
  mean.spec <- colMeans(spec)
  spec.quant <- apply(spec,2,quantile,probs=c(0.05,0.95))
  
  # View/Plot spectra
  pdf(paste(out.dir,'FFT_Raw_Dry_Spectra_QC.pdf',sep=""),height=7,width=8)
  par(mar=c(4,4,1,1.2)) #B,L,T,R
  matplot(t(spec), type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n")
  ind <- pretty(seq(from = 350, to = 2500, by = 1)) # Using pretty to standardize the axis
  ind <- ind[ind >= 350 & ind <= 2500]
  ind <- (ind - 349) / 1
  axis(1, ind, colnames(spec)[ind]) # add column names to wavelengths
  # Mean spectra
  lines(mean.spec,lwd=9)
  # lines(mean.spec+(sd.spec*1.96),lty=3,lwd=5,col="black")
  # lines(mean.spec-(sd.spec*1.96),lty=3,lwd=5,col="black")
  # CIs
  lines(spec.quant[1,],lty=1,lwd=7,col="dark grey")
  lines(spec.quant[2,],lty=1,lwd=7,col="dark grey")
  legend("bottomright",legend=c("Mean","95% CI"),lty=c(1,1),
         col=c("black","dark grey"),lwd=3)
  box(lwd=2.2)
  dev.off()
  rm(mean.spec,spec.quant,ind)
}
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#


# Raw spectra
spec <- as.matrix(spec)
# Raw
spec_corr = cor(as.matrix(spec), dry.spectra$mass_per_area)
#Log 1/R transform

pdf(paste(out.dir,'FFT_Leaf_LMA_Spectra_Correlations.pdf',sep=""),height=12,width=8)
par(mfrow=c(2,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
matplot(t(spec), type = "l", lty = 1, ylab = "Reflectance (%)", xaxt = "n",ylim=c(0,0.9))
ind <- pretty(seq(from = 350, to = 2500, by = 1)) # Using pretty to standardize the axis
ind <- ind[ind >= 350 & ind <= 2500]
ind <- (ind - 349) / 1
axis(1, ind, colnames(spec)[ind]) # add column names to wavelengths

plot(waves[,2],spec_corr,xlab="WAVELENGTH (nm)",ylab="CORRELATION",
     main="Leaf LMA (gDW/m2)", cex=0.01)
lines(waves[,2],spec_corr,lwd=4)
abline(h=0,lty=2,lwd=1.5,col="grey80")
box(lwd=2)
dev.off()

# Output correlation data
spec_corr <- data.frame(spec_corr)
names(spec_corr) <- c("Correlation")
write.csv(spec_corr,paste(out.dir,'FFT_Leaf_LMA_Spectra_Correlations.csv',sep=""),
          row.names=TRUE)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Subset spectra for analysis 

dry.spectra2 <- as.matrix(spec[,151:2051]) #151:2101 (500-2400nm), 551:2101 (900 - 2400nm),
# 851:2051 (1200 - 2400nm)
# dry.spectra2[1:5,1:10]
# dry.spectra2[1:5,1899:1901]
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Build full PLSR dataset 
full.plsr.data <- data.frame(dry.spectra$Site,dry.spectra[,6:10],dry.spectra2)
full.plsr.data[1:5,1:16]
dim(full.plsr.data)
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
# Subset data into cal/val by site
sites <- unique(full.plsr.data$dry.spectra.Site)
create.seed <- TRUE  #TRUE/FALSE

# Sample proportion for cal data
prop <- 0.8

# Random seed
if (create.seed){
  set.seed(as.vector(round(runif(5,min=5,max=9))))
  ### Write out seed
  .Random.seed[1:6]
  seed.save <- .Random.seed
  write.table(seed.save,paste(out.dir,"random.seed",sep=""));
}

### Read in previous random seed
seed <- read.table(paste(out.dir,"random.seed",sep=""))[,1];
.Random.seed <- seed

cal.plsr.data <- 0
val.plsr.data <- 0
j <- 1
for (i in sites){
  print(paste("Site: ",i,sep=""))
  temp.data <- full.plsr.data[which(full.plsr.data$dry.spectra.Site==i),]
  rows <- sample(1:nrow(temp.data),floor(prop*nrow(temp.data)))
  cal_data = droplevels(temp.data[rows,])
  val_data = droplevels(temp.data[-rows,])
  
  if(j==1){
    cal.plsr.data <- cal_data
    val.plsr.data <- val_data
  } else {
    cal.plsr.data <- rbind(cal.plsr.data,cal_data)
    val.plsr.data <- rbind(val.plsr.data,val_data)
  }
  
  j <- j+1
}
rm(temp.data)

# Datasets:
# cal.plsr.data -- For building PLSR model
print(paste("Cal observations: ",dim(cal.plsr.data)[1],sep=""))
# val.plsr.data -- Independent (external) PLSR model validation data ~20 of data
print(paste("Val observations: ",dim(val.plsr.data)[1],sep=""))

pdf(paste(out.dir,'FFT_Leaf_LMA_Cal_Val_Histograms.pdf',sep=""),height=12,width=8)
par(mfrow=c(2,1),mar=c(4,4.6,1,1.4)) #B, L, T, R
hist(cal.plsr.data$mass_per_area)
hist(val.plsr.data$mass_per_area)
dev.off()
write.csv(full.plsr.data,file=paste(out.dir,"FFT_Leaf_LMA_Full_PLSR_Dataset.csv",sep=""),row.names=FALSE)
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Run calibration PLSR analysis to select optimal number of components
dims = dim(cal.plsr.data)
k <- round(dims[1]/10)
segs = cvsegments(dims[1],k = k, type="random")
leaf.lma <- as.matrix(cal.plsr.data$mass_per_area)
spectra <- as.matrix(cal.plsr.data[,7:length(cal.plsr.data[1,])])
#standard normal variate transform [log(first derivative)]
spectra_log_dif_snv <- standardNormalVariate(X = t(diff(t(log(spectra)),differences=1, lag=3)))
train.PLS = data.frame(Y = I(leaf.lma), spec=I(spectra_log_dif_snv))
LeafLMA.pls = plsr(Y ~ spec ,scale=FALSE, ncomp=15,validation="LOO", trace=TRUE, method = "oscorespls", data = train.PLS) #segments if CV

# Examine raw PLSR output
summary(LeafLMA.pls)
plot(RMSEP(LeafLMA.pls), legendpos = "topright")
names(LeafLMA.pls)

# Examine raw the R2 
R2(LeafLMA.pls)
plot(R2(LeafLMA.pls), legendpos = "topright")

#assuming RMSEP is our method to define the components number
LeafLMA.ncomp = ceiling(which(RMSEP(LeafLMA.pls)$val==min(RMSEP(LeafLMA.pls)$val))/2)-1
# plot of prediction over measurement
plot(LeafLMA.pls,ncomp=LeafLMA.ncomp, asp=1, line=TRUE)

### Output cal data for VIP models
write.csv(cal.plsr.data,file=paste(out.dir,'FFT_Leaf_LMA_Calibration_Dataset.csv', sep=""), row.names=FALSE)
#--------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------#
# External validation - Use to test initial and final PLSR models
leaf.lma.test <- as.matrix(val.plsr.data$mass_per_area)
spectra.test <- as.matrix(val.plsr.data[,7:length(val.plsr.data[1,])])
#standard normal variate transform [log(first derivative)]
spectra_log_dif_snv.test <- standardNormalVariate(X = t(diff(t(log(spectra.test)),differences=1, lag=3)))
test.PLS = data.frame(Y = I(leaf.lma.test), spec=I(spectra_log_dif_snv.test))

plot(RMSEP(LeafLMA.pls, newdata = test.PLS))
plot(RMSEP(LeafLMA.pls,estimate=c("test"),newdata = test.PLS), main="MODEL RMSEP",
     xlab="NUM OF COMPONENTS",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)

R2(LeafLMA.pls, newdata = test.PLS)
plot(R2(LeafLMA.pls,estimate=c("test"),newdata = test.PLS), main="MODEL R2",
     xlab="NUM OF COMPONENTS",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)

# Quick validation diagnostic plots
predplot(LeafLMA.pls, ncomp = 9:11, newdata = val.plsr.data, asp = 1, line = TRUE,
         which = c("train","validation", "test"),
         xlim=c(5,300),ylim=c(5,300))

### Output validation data for VIP model
write.csv(val.plsr.data,file=paste(out.dir,'FFT_Leaf_LMA_Validation_Dataset.csv',
                                   sep=""), row.names=FALSE)
#--------------------------------------------------------------------------------------------------#

