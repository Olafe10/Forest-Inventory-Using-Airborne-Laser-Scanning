
# ALS-BASED FOREST INVENTORY WITH LIDR


#==================================
# Initialize R
#==================================

# Clear memory from old stuff

rm(list = ls())

# Only run the following four installations when starting an a new computer
# Comment when not needed

install.packages("terra")
install.packages("rgdal")
install.packages("raster")
install.packages("lidR")
install.packages("microbenchmark")

# Load installed libraries into active use

require(raster)
require(lidR)
require(sf)
require(microbenchmark)
library('terra')
# Check if the help page of the following command shows up. 
# If yes, the installation is ok.

?polygon_metrics

# Set working directory path - all your files should be here

setwd("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab1") # Edit your own folder path here

# Load a set of custom functions

source("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab1/lidr_functions.r")


#==================================
# Read data sets into R
#==================================

# Read background and plot aerial image

ortho_image <- brick("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab1/kiihtelys_ortho.tif")

plotRGB(ortho_image)

# Read field plot shapefile from the working directory

ars <- st_read("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab1/ars_plots.shp")

# Show field plot locations over the aerial image

plot(ars, add=TRUE, col="red")

# View the attributes of the field plots

View(ars)

# Read laser data file into R's LAS object

Klas <- lidR::readLAS("C:/Users/Olamidayo Fasalejo/OneDrive - University of Eastern Finland/Documents/Remote sensing/Lab1/kiihtelys.laz")

# View a summary of the LAS object

summary(Klas)

# View the coordinate system of the field plots

crs(ars)

# Set the LAS object's coordinate reference system as TM35FIN 

crs(Klas) <- 3067

# Review the summary of the LAS object

summary(Klas)
plot(Klas)
# Show the records of the LAS object

tail(Klas@data)

# View the a table of classifications of the las object

table(Klas$Classification)

# Add echo type to User Data, overwrite the previous las object

Klas <- add_echo_type(Klas)

# View the updated User Data field in tabular form

table(Klas$UserData)

# Check laser data integrity

las_check(Klas)


#========================================
# Visualization 
#========================================

# Clip a smaller portion of the total area for a better view

Klas_clip<- clip_rectangle(Klas, 665000, 6934700, 665300, 6935000)

# Plot the smaller area in different ways

plot(Klas_clip, color="Z") # Elevation
plot(Klas_clip, color="UserData") # User Data
plot(Klas_clip, color="Classification") # Classification
plot(Klas_clip, color="Intensity") # Intensity

# Re-plot Intensity with rescaled colors
plot(Klas_clip, color="Intensity", breaks = "quantile", nbreaks = 50)

# Clip point clouds of plots


Klas_clip_poicld <- clip_roi(Klas, ars)

# View point clouds of plots

plot(Klas_clip_poicld[[1]], colour='UserData')

plot(Klas_clip_poicld[[1]], colour='Z')

#=================================================
# Height normalization and removal of outliers
#=================================================

# Construct a Digital Terrain Model
# Set output file, input file, spatial resolution, and interpolation method

nKlas <- grid_terrain(Klas, res = 1, tin())

plot(nKlas)

# View the DTM in 2D and 3D hillshade

plot(nKlas)
plot_dtm3d(nKlas)

# Height normalize LAS object by subtracting the DTM

htnormlas <- Klas - nKlas

# Remove the non-normalized LAS object and the DTM to free memory

rm(...)

# Reset graphical parameters

dev.off()

# View normalized Z and intensity histograms

hist(htnormlas$Z)
hist(htnormlas$Intensity, nclass=2000, xlim=range(0,500))

# Overwrite the LAS object by keeping only the echoes with 
# realistic Z & Intensity

htnormlas <- filter_poi(htnormlas, Z < 50 & Intensity < 250)

# Review echo histograms and repeat previous step if needed

hist(htnormlas$Z)
hist(htnormlas$Intensity)

# Set negative elevations as zero

htnormlas$Z[htnormlas$Z < 0] <- 0


max(htnormlas$Z) # View Z maximum
min(htnormlas$Z) # View Z minimum

max(htnormlas$Intensity) # View Intensity maximum
min(htnormlas$Intensity) # View Intensity minimum

#============================================
# Computation of predictor variables
#============================================

# Compute preditor variables for each plot into data frame d
# Inputs: LAS file, function that computes the variables, plot polygons

d <- polygon_metrics(htnormlas, ~uef_metrics(Z, Intensity, UserData), ars)

tail(d)      # View structure of the predictor file
class(d)     # Note that it can be treated as a georeferenced data frame

View(d)
# Convert into a normal data frame and remove the last column with geometry
d <- as.data.frame(d[, 1:ncol(d)-1])

# Merge with field plots using column bind

dbind <- cbind(ars, d)

tail(dbind)      # View structure of the predictor file
class(dbind)     # Note that it can be treated as a georeferenced data frame

# Make a scatter plot of max echo height (f_hmax) and dominant height (hdom)
# If the predictors are OK, there should be a high correlation

plot(dbind$f_hmax, dbind$hdom, xlim=c(0,30), ylim=c(0,30)); abline(0,1)
cor(dbind$f_hmax, dbind$hdom)
View(dbind)
#=================================================
# Construct a regression model for plot volume
#=================================================

# Make your own plot level volume model 

modvol<-(lm(v~f_hmean+f_hstd, data = dbind))
summary(modvol)

modvol1<-(lm(v~f_hmean+f_hmax+f_hstd, data = dbind))
summary(modvol1)
plot(fitted(modvol1),resid(modvol1), xlab='Fitted Values', ylab='Residuals')
abline(h=0)

rmse(dbind$v, fitted(modvol1))

exh_var_search(dbind, 3, "v", 12, 68)

newmodel<-lm(v~f_hmean+l_p80+l_p40, data = dbind)
summary(newmodel)
plot(fitted(newmodel),resid(newmodel), xlab='Fitted Values', ylab='Residuals')
rmse(dbind$v, fitted(newmodel))
bias(dbind$v, fitted(newmodel))

dbind$lgf_hmean<- log(dbind$f_hmean)
modvol2<-(lm(v~lgf_hmean+f_hmax+f_hstd, data = dbind))
summary(modvol1)
plot(fitted(modvol2),resid(modvol2), xlab='Fitted Values', ylab='Residuals')
abline(h=0)

rmse(dbind$v, fitted(modvol2))


#==================================================
# Wall-to-wall prediction
#==================================================

# Compute a raster stack of ALS metrics for the whole area
# Parameters: LAS object, function to compute, resolution

abagrid <- grid_metrics(htnormlas, ~uef_metrics(Z, Intensity, UserData), res = 15)

# Plot the f_hmean from the raster stack

plot(abagrid$f_hmean)

# Predict volume using the raster brick and your model object
# Save as v_pred

v_pred <- 10557.025 + 27.758*abagrid$f_hmean - 15344.057*abagrid$l_p80 + 4773.840*abagrid$l_p40

v_pred1 <- -32.290 + 25.659*abagrid$f_hmean - 4.110*abagrid$f_hstd

# Force areas where top echo is height is < 2 m as 0
v_pred[v_pred < 0] <-0
v_pred1[v_pred1 < 0] <- 0
v_pred[abagrid$f_hmax < 2] <- 0

rmse(abagrid$, predict(v_pred))

# Reset graphics and plot volume grid

dev.off()
plot(v_pred, main="Volume m3/ha")
plot(v_pred1, main="Volume m3/ha")

# Save volume grid as a .tif file and view it in a GIS

writeRaster(v_pred, "volume", overwrite=T)
writeRaster(v_pred1, "volume1", overwrite=T)

# Free memory

rm(abagrid)

#============================================
# Construction of a canopy height model
#============================================

# Interpolate a CHM
# Parameters: LAS object, spatial resolution, interpolation algorithm

CHM <- rasterize_canopy(htnormlas, res = 0.5, algorithm = pitfree(thresholds = c(0,5,10,15,20,25)))

# Plot CHM

plot(CHM)

# Add field plot polygons 

plot(ars, add=T, col="black")

# Crop a smaller area from the CHM for closer inspection and plot it

CHMcp <-crop(CHM, extent(664030, 664090, 6935380, 6935440))

plot(CHMcp)


#============================================
# Individual tree detection
#============================================

# Define a function for computing a height-dependent moving window size

window_size <- function(height){ 1.2 + 0.003*height^2 } 

# Compute a georeferenced data frame with tree coordinates and heights

trees_chm <- find_trees(CHM, lmf(window_size))

tail(trees_chm) # View output


# Visualize tree locations

plot(CHMcp) # Re-plot the small-area CHM

plot(ars, add=T, lwd=2, col=rgb(0,0,0,0)) # Re-add plot polygons as transparent

plot(trees_chm, add=T) # Add detected tree locations on top


# Apply height correction model with the tree top data frame

trees_chm$corrected_height <- 1.01 + trees_chm$Z*1.2

tail(trees_chm) # View the result

#============================================
# Construction of a tree level volume model
#============================================

dev.off() # Reset graphics

# Read in tree file

ars_trees <-read.csv("ars_trees.csv", header=T, sep=",", dec=".")

# Make your own tree level volume model 
tree_model<- lm(ars_trees$v~ars_trees$h)
summary(tree_model)
rmse(ars_trees$v, fitted(tree_model))

plot(fitted(tree_model),resid(tree_model), xlab='Fitted Values', ylab='Residuals')

ars_trees$sqh <- sqrt(ars_trees$h)
ars_trees$sqv <- sqrt(ars_trees$v)
tree_model_sq<- lm(ars_trees$sqv~ars_trees$sqh)
summary(tree_model_sq)
rmse(ars_trees$v, fitted(tree_model_sq))

plot(fitted(tree_model_sq),resid(tree_model_sq), xlab='Fitted Values', ylab='Residuals')

ars_trees$lgh <- log(ars_trees$h)
ars_trees$lgv <- log(ars_trees$v)
tree_model_lg<- lm(ars_trees$lgv~ars_trees$lgh)
summary(tree_model_lg)
rmse(ars_trees$v, fitted(tree_model_lg))

plot(fitted(tree_model_lg),resid(tree_model_lg), xlab='Fitted Values', ylab='Residuals')
abline(h=0)
head(ars_trees)
plot(ars_trees$h, ars_trees$v)

x <- seq(0,30,1) # Creates numeric sequence 0, 1, 2, ., 30
y <- exp(-11.56894+3.52736*log(x)) # Apply your model coefficients to the sequence x
lines(x, y, col="red")
rm(x)
rm(y)
# Apply your model to predict a new column in the tree top file

ars_trees$newvol <- exp(-11.56894 + 3.52736*ars_trees$lgh)
trees_chm$v_itd <- exp(-11.56894 + 3.52736*log(trees_chm$corrected_height)+ 0.447^2/2)
tail(ars_trees)
head(trees_chm)

rmse(ars_trees$v, fitted(tree_model_lg))
bias(ars_trees$v, fitted(tree_model_lg))
#============================================
# Plot level accuracy assessment for ITD
#============================================

# Convert the detected trees to a different spatial format

topsf <- st_as_sf(trees_chm)         
  
# Set the coordinate system of this object to be the same as the field plots

st_crs(topsf) <- 3067 

# Make a spatial overlay of field plots and detected trees 
View(ars)
topsf_ars <- st_join(ars,topsf) 

# View the last rows of the new object with plot-labeled itd trees

head(topsf_ars)

remove(topsf_ars)
# Compute a new column called v_itd in your FIELD PLOT TABLE using function 
# tapply with your plot-labeled itd trees. Tapply applies a function (here sum())
# to the values of a given column (v_itd) using a given grouping variable (plot).

ars$v_itd <- tapply(topsf_ars$v_itd, topsf_ars$plot_id, sum)

# Repeat, but now compute the plot level stem number by using the length() function
# instead of sum() and save it as n_itd

ars$n_itd <- tapply(topsf_ars$v_itd, topsf_ars$plot_id, length)

# Scale v_itd and n_itd to per hectare using plot radius
head(ars)
ars$v_itd_hectare <- ars$v_itd * (10000 / (pi*10^2))
ars$n_itd_hectare <- ars$n_itd * (10000/ (pi*10^2))

# View the last rows of the field plot table

tail(ars)
rm(topsf, topsf_ars, topsf_ars1)
# Make scatter plots of plot level predicted vs observed volume and stem numbe

plot(ars$v_itd_hectare, ars$v)
abline(a=0, b=1)

# Compute plot level RMSE-% and bias-% of both volume and stem number

relrmse(ars$v, ars$v_itd_hectare)
relbias(ars$v, ars$n_itd_hectare)

