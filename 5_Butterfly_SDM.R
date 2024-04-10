# DESCRIPTION:

# building SDMs using variables from thinned 2020 LiDAR data
# Analyse colinearity patterns among LiDAR metrics.
# Run SDM models (glm, rf, maxent) per species.
# Derive model performance statistics.
# Produce species probability distributions for 2020.
# Apply model to 2010 LiDAR metrics to produce species probability distribution for 2010.
# Creation: 09-2023
# Author: Martha Ledger

#### Loading data & Packages ####

# Install packages:
install.packages("devtools")
install.packages("ecospat")
install.packages("usdm")
install.packages("sdm")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("corrplot")
install.packages("rgeos")
install.packages("tidyverse")
install.packages("GGally")
install.packages("FactoMineR")
install.packages("factoextra")

# Read packages:
library(raster)
library(dismo)
library(devtools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(usdm)
library(sdm)
library(tidyverse)
library(GGally)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(gridExtra)
library(grid)
#installAll() # run just once to ensure that all methods are available for modelling in 'sdm' (sourced from other packages)


##### 1. Lethe chandica #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Lethe_chandica/2020_thinned/"
setwd(workingdirectory)

# Import data frame of occurences data + LiDAR metrics

lidarmetrics_LC <- read.csv("Lethe_chandica_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics
lidarmetrics_LC <- lidarmetrics_LC %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_LC <- na.exclude(lidarmetrics_LC)

# count presence and absence points in data frame
lidarmetrics_LC %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_LC[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_LC[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop: Open area, 5-20m density, 1-5m density (colinearity problem),
#.    Based on VIF, drop any above 3: <0.2m density, Low.veg.roughness..m., Slope
#     Used as input for species distribution modelling
lidarmetrics_LC <- lidarmetrics_LC %>%
  subset(select = c(Lon, Lat, presence, X0.2.1.m.density...., X.20.m.density...., Total.veg.roughness..m., Elevation..m., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

LC_sdm_metrics <- sdmData(formula=presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_LC, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelLC <- sdm(presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Open.patches..count.+Land.cover, data=LC_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelLC, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Lethe_chandica/2020_thinned/modelLC.rds")


#Derive model performance

# AUC, COR, TSS, Deviance over 100 runs
modelLC

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelLC)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Install packages:
# install.packages("devtools")
# install.packages("ecospat")
# install.packages("usdm")
# install.packages("sdm")
# install.packages("rJava")
# install.packages("tidyverse")
# install.packages("RColorBrewer")
# # Updates of usdm & sdm:
# devtools::install_github('babaknaimi/usdm', force = TRUE)
# devtools::install_github('babaknaimi/sdm', force = TRUE)

# Read packages:
library(raster)
library(dismo)
library(devtools)
library(rgeos)
library(rgdal)
library(ggplot2)
library(usdm)
library(sdm)
library(tidyverse)
library(GGally)
library(grid)
library(ggpubr)
library(gridExtra)
library(RColorBrewer)


# Read model object:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Lethe_chandica/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelLC <- readRDS("modelLC.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelLC, method = 'rf', main = expression(italic('Lethe chandica')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))



#Variable importances

fea_imp=getVarImp(modelLC, method = 'rf')
varimp_LC <- plot(fea_imp, 'cor', main = 'Lethe chandica', col=c('#E69F00', '#009E73', '#009E73', '#0072B2', '#0072B2', '#0072B2')) + 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('land cover', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_LC


#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Lethe chandica - greyed out variables are not included in the model (high colinearity)
dens.02_1=rcurve(modelLC,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.a20=rcurve(modelLC,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
# height=rcurve(modelLC,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# height.df=height$data
total.veg.rough=rcurve(modelLC,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
# low.veg.rough=rcurve(modelLC,n=c("Low.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# low.veg.rough.df=low.veg.rough$data
elevation=rcurve(modelLC,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
# slope=rcurve(modelLC,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# slope.df=slope$data
# openarea=rcurve(modelLC,n=c("Open.area..ha"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelLC,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
# edge.ext=rcurve(modelLC,n=c("Edge.extent..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# edge.ext.df=edge.ext$data
landcover=rcurve(modelLC,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots - greyed out variables are not included in the model (high colinearity)
# p1=ggplot(data=dens.b02.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("0.2-1 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 15, 30, 45)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
# p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 20), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 50), breaks = c(0, 10, 20, 30, 40, 50)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 15500), breaks = c(0, 5000, 10000, 15000, 20000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p13=ggplot(data=edge.ext.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("edge extent (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 750), breaks = c(0, 250, 500, 750)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

# Lethe chandica
LCresponse <- grid.arrange(
  p2,
  p5,
  p7,
  p9,
  p12,
  p14,
  nrow = 1, ncol = 6,
  top = textGrob("Lethe chandica", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("LCplots.tiff", width = 6, height = 4, units = 'in', res = 600)
LC_plots <- grid.arrange(
  varimp_LC,
  LCresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS 2020

library(terra)

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)

DTM <- raster(list[2])
open_patches <- raster(list[7])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[16])


# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- raster::resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Open.patches..count.", "Land.cover")

library(tictoc)
tic()
LC_2020_presence <- predict(modelLC, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Lethe_chandica_distribution", overwrite = T)
toc()

plot(LC_2020_presence)


# PREDICTIONS 2010

library(terra)

# raster stack of 2010 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[14])


# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- raster::resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Open.patches..count.", "Land.cover")

workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

library(tictoc)
tic()
LC_2020_presence <- predict(modelLC, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Lethe_chandica_distribution_2010", overwrite = T)
toc()

plot(LC_2020_presence)


##### 2. Prosotas dubiosa #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_dubiosa/2020_thinned/"
setwd(workingdirectory)

# Import Lidar metrics

lidarmetrics_PD <- read.csv("Prosotas_dubiosa_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics
lidarmetrics_PD <- lidarmetrics_PD %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_PD <- na.exclude(lidarmetrics_PD)

# count presence and absence points in data frame
lidarmetrics_PD %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_PD[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_PD[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop Open area, 5-20m density (colinearity problem)
#.    Based on VIF, drop any above 3: <0.2m density, 1-5m density, height, low veg roughness.
#     Used as input for species distribution modelling
lidarmetrics_PD <- lidarmetrics_PD %>%
  subset(select = c(Lon, Lat, presence, X0.2.1.m.density...., X.20.m.density...., Total.veg.roughness..m., Elevation..m., Slope..degrees., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

PD_sdm_metrics <- sdmData(formula=presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_PD, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelPD <- sdm(presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover, data=PD_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelPD, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_dubiosa/2020_thinned/modelPD.rds")


#Derive model performance

# AUC, COR TSS, Deviance over 100 runs
modelPD

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelPD)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_dubiosa/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelPD <- readRDS("modelPD.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelPD, method = 'rf', main = expression(italic('Prosotas dubiosa')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))



#Variable importances

fea_imp=getVarImp(modelPD, method = 'rf')
varimp_PD <- plot(fea_imp, 'cor', main = 'Prosotas dubiosa', col=c('#E69F00', '#009E73', '#009E73', '#0072B2', '#0072B2', '#0072B2', '#0072B2')) + 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('urban proximity', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_PD


#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Prosotas dubiosa

dens.02_1=rcurve(modelPD,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.a20=rcurve(modelPD,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
# height=rcurve(modelPD,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# height.df=height$data
total.veg.rough=rcurve(modelPD,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
elevation=rcurve(modelPD,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
slope=rcurve(modelPD,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
slope.df=slope$data
# openarea=rcurve(modelPD,n=c("Open.area..ha."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelPD,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
landcover=rcurve(modelPD,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots
# p1=ggplot(data=dens.b02.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) + 
  geom_line(linewidth=1, colour='#E69F00', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("0.2-1 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 15, 30, 45)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
# p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total veg. roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 950), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

PDresponse <- grid.arrange(
  p2,
  p5,
  p7,
  p9,
  p10,
  p12,
  p14,
  nrow = 1, ncol = 7,
  top = textGrob("Prosotas dubiosa", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("PDplots.tiff", width = 6, height = 4, units = 'in', res = 600)
PD_plots <- grid.arrange(
  varimp_PD,
  PDresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[16])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
PD_2020_presence <- predict(modelPD, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Prosotas_dubiosa_distribution", overwrite = T) # 10 cores for PC, 6 for iMac
toc()

plot(PD_2020_presence)


# raster stack of 2010 metrics

workingdirectory <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[14])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
PD_2010_presence <- predict(modelPD, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Prosotas_dubiosa_distribution_2010", overwrite = T) # 10 cores for PC, 6 for iMac
toc()

plot(PD_2010_presence)


##### 3. Zizula hylax #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Zizula_hylax/2020_thinned/"
setwd(workingdirectory)

# Import Lidar metrics

lidarmetrics_ZH <- read.csv("Zizula_hylax_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics
lidarmetrics_ZH <- lidarmetrics_ZH %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_ZH <- na.exclude(lidarmetrics_ZH)

# count presence and absence points in data frame
lidarmetrics_ZH %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_ZH[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_ZH[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop: Open.area..ha. X5.20.m.density....(colinearity problem)
#.    Based on VIF, drop any above 3: <0.2m density, X1.5.m.m.density...., low veg roughness.
#     Used as input for species distribution modelling
lidarmetrics_ZH <- lidarmetrics_ZH %>%
  subset(select = c(Lon, Lat, presence, X0.2.1.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Elevation..m., Slope..degrees., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

ZH_sdm_metrics <- sdmData(formula=presence~X0.2.1.m.density....+X.20.m.density....+Height..m.+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_ZH, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelZH <- sdm(presence~X0.2.1.m.density....+X.20.m.density....+Height..m.+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover, data=ZH_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelZH, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Zizula_hylax/2020_thinned/modelZH.rds")


#Derive model performance

# AUC, COR TSS, Deviance over 100 runs
modelZH

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelZH)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Zizula_hylax/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelZH <- readRDS("modelZH.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelZH, method = 'rf', main = expression(italic('Zizula hylax')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))



#Variable importances

# Species combined in 1 figure

fea_imp=getVarImp(modelZH, method = 'rf')
varimp_ZH <- plot(fea_imp, 'cor', main = 'Zizula hylax', col=c('#E69F00', '#009E73', '#009E73', '#009E73', '#0072B2', '#0072B2', '#0072B2','#0072B2'))+ 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('land cover', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_ZH


#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Zizula hylax

dens.02_1=rcurve(modelZH,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.a20=rcurve(modelZH,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
height=rcurve(modelZH,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
height.df=height$data
# low.veg.rough=rcurve(modelZH,n=c("Low.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# low.veg.rough.df=low.veg.rough$data
total.veg.rough=rcurve(modelZH,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
elevation=rcurve(modelZH,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
slope=rcurve(modelZH,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
slope.df=slope$data
# openarea=rcurve(modelZH,n=c("Open.area..ha."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelZH,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
# edge.ext=rcurve(modelZH,n=c("Edge.extent..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# edge.ext.df=edge.ext$data
landcover=rcurve(modelZH,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots
# p1=ggplot(data=dens.b02.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) + 
  geom_line(linewidth=1, colour='#E69F00', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("0.2-1 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 15, 30, 45)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
# p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total veg. roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p13=ggplot(data=edge.ext.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("edge extent (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 750), breaks = c(0, 250, 500, 750)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

ZHresponse <- grid.arrange(
  p2,
  p5,
  p6,
  p7,
  p9,
  p10,
  p12,
  p14,
  nrow = 1, ncol = 8,
  top = textGrob("Zizula hylax", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("ZHplots.tiff", width = 6, height = 4, units = 'in', res = 600)
ZH_plots <- grid.arrange(
  varimp_ZH,
  ZHresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

CHM <- raster(list[1])
DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[16])


# check that resolutions and extents now match (output should be TRUE)
compareRaster(CHM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
CHM_10 <- aggregate(CHM, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, CHM_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Height..m.", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
ZH_2020_presence <- predict(modelZH, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Zizula_hylax_distribution", overwrite = T)
toc()

plot(ZH_2020_presence)


# raster stack of 2010 metrics

workingdirectory <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

CHM <- raster(list[1])
DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[11])
above20 <- raster(list[12])
betw021 <- raster(list[14])


# check that resolutions and extents now match (output should be TRUE)
compareRaster(CHM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing
slope_10 <- aggregate(slope, fact = 2, fun = "mean")
CHM_10 <- aggregate(CHM, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, CHM_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Height..m.", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
ZH_2010_presence <- predict(modelZH, newdata = rasters_10, method = 'rf', mean = T, nc = 6, filename = "Zizula_hylax_distribution_2010", overwrite = T)
toc()

plot(ZH_2010_presence)


##### 4. Notocrypta paralysos #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Notocrypta_paralysos/2020_thinned/"
setwd(workingdirectory)

# Import Lidar metrics

lidarmetrics_NP <- read.csv("Notocrypta_paralysos_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics
lidarmetrics_NP <- lidarmetrics_NP %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_NP <- na.exclude(lidarmetrics_NP)

# count presence and absence points in data frame
lidarmetrics_NP %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_NP[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_NP[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop: Open.area..ha., X5.20.m.density...., Height..m. (colinearity problem)
#.    Based on VIF, drop any above 3: none
#     Used as input for species distribution modelling
lidarmetrics_NP <- lidarmetrics_NP %>%
  subset(select = c(Lon, Lat, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X.20.m.density...., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

NP_sdm_metrics <- sdmData(formula=presence~X.0.2.m.density....+X0.2.1.m.density....+X1.5.m.m.density....+X.20.m.density....+Total.veg.roughness..m.+Low.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_NP, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelNP <- sdm(presence~X.0.2.m.density....+X0.2.1.m.density....+X1.5.m.m.density....+X.20.m.density....+Total.veg.roughness..m.+Low.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover, data=NP_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelNP, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Notocrypta_paralysos/2020_thinned/modelNP.rds")


#Derive model performance

# AUC, COR TSS, Deviance over 100 runs
modelNP

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelNP)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Notocrypta_paralysos/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelNP <- readRDS("modelNP.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelNP, method = 'rf', main = expression(italic('Notocrypta paralysos')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))



#Variable importances

# Species combined in 1 figure

fea_imp=getVarImp(modelNP, method = 'rf')
varimp_NP <- plot(fea_imp, 'cor', main = 'Notocrypta paralysos', col=c('#E69F00', '#E69F00', '#009E73', '#009E73', '#009E73', '#E69F00', '#0072B2', '#0072B2', '#0072B2','#0072B2'))+ 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('land cover', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_NP

#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Notocrypta paralysos

dens.02=rcurve(modelNP,n=c("X.0.2.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02.df=dens.02_1$data
dens.02_1=rcurve(modelNP,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.1_5=rcurve(modelNP,n=c("X1.5.m.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.1_5.df=dens.1_5$data
dens.a20=rcurve(modelNP,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
# height=rcurve(modelNP,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# height.df=height$data
low.veg.rough=rcurve(modelNP,n=c("Low.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
low.veg.rough.df=low.veg.rough$data
total.veg.rough=rcurve(modelNP,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
elevation=rcurve(modelNP,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
slope=rcurve(modelNP,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
slope.df=slope$data
# openarea=rcurve(modelNP,n=c("Open.area..ha."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelNP,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
landcover=rcurve(modelNP,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots
p1=ggplot(data=dens.02.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total veg. roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p13=ggplot(data=edge.ext.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("edge extent (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 750), breaks = c(0, 250, 500, 750)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

NPresponse <- grid.arrange(
  p1,
  p2,
  p3,
  p5,
  p7,
  p8,
  p9,
  p10,
  p12,
  p14,
  nrow = 1, ncol = 10,
  top = textGrob("Notocrypta paralysos", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("NPplots.tiff", width = 6, height = 4, units = 'in', res = 600)
NP_plots <- grid.arrange(
  varimp_NP,
  NPresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

below02 <- raster(list[13])
betw15 <- raster(list[14])
DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
low_rough <- raster(list[3])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
betw021 <- raster(list[16])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(below02, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw15, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(low_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
low_rough_10 <- aggregate(low_rough, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw15_10 <- aggregate(betw15, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
below02_10 <- aggregate(below02, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(below02_10, betw021_10, betw15_10, above20_10, total_rough_10, low_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X.0.2.m.density....", "X0.2.1.m.density....", "X1.5.m.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Low.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
NP_2020_presence <- predict(modelNP, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Notocrypta_paralysos_distribution", overwrite = T)
toc()

plot(NP_2020_presence)


# raster stack of 2010 metrics

workingdirectory <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

below02 <- raster(list[13])
betw15 <- raster(list[15])
DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
low_rough <- raster(list[3])
urbanicity <- raster(list[11])
above20 <- raster(list[12])
betw021 <- raster(list[14])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(below02, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw15, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(low_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)


# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
low_rough_10 <- aggregate(low_rough, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw15_10 <- aggregate(betw15, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
below02_10 <- aggregate(below02, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(below02_10, betw021_10, betw15_10, above20_10, total_rough_10, low_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X.0.2.m.density....", "X0.2.1.m.density....", "X1.5.m.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Low.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
NP_2010_presence <- predict(modelNP, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Notocrypta_paralysos_distribution_2010", overwrite = T)
toc()

plot(NP_2010_presence)


##### 5. Prosotas nora #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_nora/2020_thinned/"
setwd(workingdirectory)

# Import Lidar metrics

lidarmetrics_PN <- read.csv("Prosotas_nora_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics
lidarmetrics_PN <- lidarmetrics_PN %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_PN <- na.exclude(lidarmetrics_PN)

# count presence and absence points in data frame
lidarmetrics_PN %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_PN[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_PN[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop: Open.area..ha., X1.5.m.m.density...., X5.20.m.density...., Height..m. (colinearity problem)
#.    Based on VIF, drop any above 3: Low veg roughness.
#     Used as input for species distribution modelling
lidarmetrics_PN <- lidarmetrics_PN %>%
  subset(select = c(Lon, Lat, presence, X.0.2.m.density...., X0.2.1.m.density...., X.20.m.density...., Total.veg.roughness..m., Elevation..m., Slope..degrees., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

PN_sdm_metrics <- sdmData(formula=presence~X.0.2.m.density....+X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_PN, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelPN <- sdm(presence~X.0.2.m.density....+X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Slope..degrees.+Open.patches..count.+Land.cover, data=PN_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelPN, file = "/Users/marthaledger/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_nora/2020_thinned/modelPN.rds")
saveRDS(modelPN, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_nora/2020_thinned/modelPN.rds")


#Derive model performance

# AUC, COR TSS, Deviance over 100 runs
modelPN

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelPN)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Prosotas_nora/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelPN <- readRDS("modelPN.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelPN, method = 'rf', main = expression(italic('Prosotas nora')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))



#Variable importances

# Species combined in 1 figure

fea_imp=getVarImp(modelPN, method = 'rf')
varimp_PN <- plot(fea_imp, 'cor', main = 'Prosotas nora', col=c('#E69F00', '#E69F00', '#009E73', '#009E73', '#0072B2', '#0072B2', '#0072B2', '#0072B2'))+ 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('land cover', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_PN


#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Prosotas nora

dens.02=rcurve(modelPN,n=c("X.0.2.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02.df=dens.02$data
dens.02_1=rcurve(modelPN,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.a20=rcurve(modelPN,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
# height=rcurve(modelPN,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# height.df=height$data
# low.veg.rough=rcurve(modelPN,n=c("Low.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# low.veg.rough.df=low.veg.rough$data
total.veg.rough=rcurve(modelPN,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
elevation=rcurve(modelPN,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
slope=rcurve(modelPN,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
slope.df=slope$data
# openarea=rcurve(modelPN,n=c("Open.area..ha."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelPN,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
landcover=rcurve(modelPN,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots
p1=ggplot(data=dens.02.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) +
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("0.2-1 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 15, 30, 45)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total veg. roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p13=ggplot(data=edge.ext.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("edge extent (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 750), breaks = c(0, 250, 500, 750)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

PNresponse <- grid.arrange(
  p1,
  p2,
  p5,
  p7,
  p9,
  p10,
  p12,
  p14,
  nrow = 1, ncol = 8,
  top = textGrob("Prosotas nora", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("PNplots.tiff", width = 6, height = 4, units = 'in', res = 600)
PN_plots <- grid.arrange(
  varimp_PN,
  PNresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[10])
above20 <- raster(list[12])
below02 <- raster(list[13])
betw021 <- raster(list[16])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(below02, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)

# first, resample rasters to 5m so that spatial resolutions match the open_patches raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
below02_10 <- aggregate(below02, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(below02_10, betw021_10, above20_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X.0.2.m.density....", "X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
PN_2020_presence <- predict(modelPN, newdata = rasters_10, method = 'rf', mean = T, nc = 6, filename = "Prosotas_nora_distribution", overwrite = T)
toc()

plot(PN_2020_presence)


# raster stack of 2010 metrics

workingdirectory <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
slope <- raster(list[8])
total_rough <- raster(list[9])
urbanicity <- raster(list[11])
above20 <- raster(list[12])
below02 <- raster(list[13])
betw021 <- raster(list[14])

# check that resolutions and extents now match (output should be TRUE)
compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(slope, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(below02, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)

# first, resample rasters to 5m so that spatial resolutions match the open_patches raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

slope_10 <- aggregate(slope, fact = 2, fun = "mean")
DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
below02_10 <- aggregate(below02, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(below02_10, betw021_10, above20_10, total_rough_10, DTM_10, slope_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X.0.2.m.density....", "X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Slope..degrees.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
PN_2010_presence <- predict(modelPN, newdata = rasters_10, method = 'rf', mean = T, nc = 6, filename = "Prosotas_nora_distribution_2010", overwrite = T)
toc()

plot(PN_2010_presence)


##### 6. Neptis nata #####

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Neptis_nata/2020_thinned/"
setwd(workingdirectory)

# Import Lidar metrics

lidarmetrics_NN <- read.csv("Neptis_nata_point_metrics_2020_thinned.csv", sep = ",", stringsAsFactors = F)

#Data preparation

# Order metrics

lidarmetrics_NN <- lidarmetrics_NN %>%
  subset(select = c(ID, Lat, Lon, presence, X.0.2.m.density...., X0.2.1.m.density...., X1.5.m.m.density...., X5.20.m.density...., X.20.m.density...., Height..m., Total.veg.roughness..m., Low.veg.roughness..m., Elevation..m., Slope..degrees., Open.area..ha., Open.patches..count., Land.cover))

lidarmetrics_NN <- na.exclude(lidarmetrics_NN)

# count presence and absence points in data frame
lidarmetrics_NN %>% count(presence, 0)

#Colinearity assessment

# Spearman rank correlation plots
#tiff("C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/presence/corBplot.tiff", width = 6, height = 5, units = 'in', res = 600)
p <- ggcorr(lidarmetrics_NN[,c(5:17)], c("pairwise", "spearman"), name = expression(italic("Spearman's r")), label=TRUE, label_alpha=TRUE, label_size=2.5, hjust=0.85, size=3, layout.exp=3)
p
#dev.off()

# Pairwise VIF variable selection with usdm (spearman-based)
usdm::vifcor(lidarmetrics_NN[,c(5:17)], th=0.7)

# VIF variable selections

#     Filter metrics selected by VIF procedure
#     Drop: Open.area..ha., Height..m., X1.5.m.m.density...., X5.20.m.density.... Slope..degrees. (colinearity problem)
#.    Based on VIF, drop any above 3: X.0.2.m.density...., Low veg roughness.
#     Used as input for species distribution modelling
lidarmetrics_NN <- lidarmetrics_NN %>%
  subset(select = c(Lon, Lat, presence, X0.2.1.m.density...., X.20.m.density...., Total.veg.roughness..m., Elevation..m., Open.patches..count., Land.cover)) 


#SDM Models

# Inputs for sdm
crs_HK1980 <- CRS("+proj=tmerc +lat_0=22.3121333333333 +lon_0=114.178555555556 +k=1 +x_0=836694.05 +y_0=819069.8 +ellps=intl +towgs84=-162.619,-276.959,-161.764,-0.067753,2.243648,1.158828,-1.094246 +units=m +no_defs +type=crs")

NN_sdm_metrics <- sdmData(formula=presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Open.patches..count.+Land.cover + coords(Lon+Lat), train=lidarmetrics_NN, crs=crs_HK1980)

# Run models, default settings (500 trees, 20 nodes)
set.seed(111)
modelNN <- sdm(presence~X0.2.1.m.density....+X.20.m.density....+Total.veg.roughness..m.+Elevation..m.+Open.patches..count.+Land.cover, data=NN_sdm_metrics, methods=c('glm', 'rf', 'maxent'), replication=c('boot'), n=100, test.percent=30)

# Save models for visualisation script
saveRDS(modelNN, file = "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Neptis_nata/2020_thinned/modelNN.rds")


#Derive model performance

# AUC, COR TSS, Deviance over 100 runs
modelNN

# Derive standard deviation of statistics

# Select model
ModelStats <- getEvaluation(modelNN)

ModelStats_glm <- ModelStats %>%
  slice(1:100)
ModelStats_rf <- ModelStats %>%
  slice(101:200)
ModelStats_max <- ModelStats %>%
  slice(201:300)

apply(ModelStats_glm, 2, sd)
apply(ModelStats_rf, 2, sd)
apply(ModelStats_max, 2, sd)

# SDM visualisation

# Read data:
workingdirectory="C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/Neptis_nata/2020_thinned/"
setwd(workingdirectory)

# Import data:

# Read final models
#     3 methods, ran 100x
modelNN <- readRDS("modelNN.rds")

#Model Performance

# ROC (Receiver Operating Characteristic)
#par(mar = c(4, 4, 3, 2), mgp = c(2.5,1,0))

#tiff("rocHplot.tiff", width = 6, height = 5, units = 'in', res = 600)
roc(modelNN, method = 'rf', main = expression(italic('Neptis nata')), 
    xlab = 'False positive rate', ylab = 'True positive rate', 
    cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.9) 
#dev.off()

# return to default
#par(mar = c(5.1,4.1,4.1,2.1), mgp = c(3,1,0))


#Variable importances

# Species combined in 1 figure

fea_imp=getVarImp(modelNN, method = 'rf')
varimp_NN <- plot(fea_imp, 'cor', main = 'Neptis nata', col=c('#E69F00', '#009E73', '#009E73', '#0072B2', '#0072B2', '#0072B2')) + 
  scale_y_continuous(name ="Relative variable importance", limits = c(0,0.5), expand=c(0,0)) +
  scale_x_discrete(limits=c('Land.cover', 'Open.patches..count.', 'Open.area..ha.', 'Slope..degrees.', 'Elevation..m.', 'Low.veg.roughness..m.', 'Total.veg.roughness..m.',
                            'Height..m.', 'X.20.m.density....', 'X5.20.m.density....', 'X1.5.m.m.density....', 'X0.2.1.m.density....', 'X.0.2.m.density....'),
                   labels=c('land cover', 'open patches', 'open area', 'slope', 'elevation', 'low veg. roughness', 'total veg. roughness',
                            'height', '>20m density', '5-20m density', '1-5m density', '0.2-1m density', '<0.2m density')) +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 8), axis.text = element_text(size = 7),
        plot.title = element_text(size = 9, face = 'italic'))
varimp_NN


# Combined plot
#varimp_all <- ggarrange(varimp_B, varimp_LC, ncol = 2, nrow = 1, widths = c(1, 1))


#Response curves

# 1: Derive all feature responses per species. Range 101-200 = RF model.
# 2: Run all response plots (p1-p10); deny if absent for focal species. 
# 3: Derive response plots per species (mean + confidence interval)


# Neptis nata

dens.02_1=rcurve(modelNN,n=c("X0.2.1.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.02_1.df=dens.02_1$data
dens.a20=rcurve(modelNN,n=c("X.20.m.density...."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
dens.a20.df=dens.a20$data
# height=rcurve(modelNN,n=c("Height..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# height.df=height$data
# low.veg.rough=rcurve(modelNN,n=c("Low.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# low.veg.rough.df=low.veg.rough$data
total.veg.rough=rcurve(modelNN,n=c("Total.veg.roughness..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
total.veg.rough.df=total.veg.rough$data
elevation=rcurve(modelNN,n=c("Elevation..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
elevation.df=elevation$data
# slope=rcurve(modelNN,n=c("Slope..degrees."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# slope.df=slope$data
# openarea=rcurve(modelNN,n=c("Open.area..ha."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# openarea.df=openarea$data
patches_low=rcurve(modelNN,n=c("Open.patches..count."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
patches_low.df=patches_low$data
# edge.ext=rcurve(modelNN,n=c("Edge.extent..m."),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
# edge.ext.df=edge.ext$data
landcover=rcurve(modelNN,n=c("Land.cover"),id= c(seq(from=101,to=200,by=1)),mean=T,confidence=T)
landcover.df=landcover$data

# Response plots
# p1=ggplot(data=dens.b02.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("<0.2 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
p2=ggplot(data=dens.02_1.df, aes(x=Value, y=Response)) + 
  geom_line(linewidth=1, colour='#E69F00', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value, ymin=lower, ymax=upper), linetype=2, alpha=0.3, show.legend = FALSE) +
  xlab("0.2-1 m density (%)") + ylab("Probability of occurrence") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 15, 30, 45)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 8))
# p3=ggplot(data=dens.1_5.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#E69F00', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("1-5 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p4=ggplot(data=dens.5_20.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("5-20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 60), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p5=ggplot(data=dens.a20.df, aes(x=Value, y=Response)) + 
  geom_line(size=1, colour='#009E73', show.legend = FALSE) + 
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab(">20 m density (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p6=ggplot(data=height.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#009E73', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("height (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p7=ggplot(data=total.veg.rough.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#009E73', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("total veg. roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 25), breaks = c(0, 10, 20)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p8=ggplot(data=low.veg.rough.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("low vegetation roughness (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 0.75), breaks = c(0, 0.25, 0.5, 0.75)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p9=ggplot(data=elevation.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("elevation (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 800), breaks = c(0, 200, 400, 600, 800)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p10=ggplot(data=slope.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("slope (degrees)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 45), breaks = c(0, 10, 20, 30, 40)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p11=ggplot(data=openarea.df, aes(x=Value, y=Response)) +
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("open area (ha)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 2.5), breaks = c(0, 1, 2)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p12=ggplot(data=patches_low.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("open patches (#)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 16000), breaks = c(0, 5000, 10000, 15000)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
# p13=ggplot(data=edge.ext.df, aes(x=Value, y=Response)) + 
#   geom_line(size=1, colour='#0072B2', show.legend = FALSE) + 
#   geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
#   xlab("edge extent (m)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 750), breaks = c(0, 250, 500, 750)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
#   theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())
p14=ggplot(data=landcover.df, aes(x=Value, y=Response)) +
  geom_line(size=1, colour='#0072B2', show.legend = FALSE) +
  geom_ribbon(aes(x=Value,ymin=lower, ymax=upper), linetype=2, alpha=0.3,show.legend = FALSE) +
  xlab("urban land cover (%)") + theme_grey(base_size = 8) + scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100)) + scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 7), axis.title.y = element_blank(), axis.text.y = element_blank())

NNresponse <- grid.arrange(
  p2,
  p5,
  p7,
  p9,
  p12,
  p14,
  nrow = 1, ncol = 6,
  top = textGrob("Neptis nata", gp=gpar(fontsize=9,font=3)),
  widths = c(2.5,2,2,2,2,2)
)


#Final plots

# Built from variable importance & response curves
tiff("NNplots.tiff", width = 6, height = 4, units = 'in', res = 600)
NN_plots <- grid.arrange(
  varimp_NN,
  NNresponse,
  nrow = 2, ncol = 1
)
dev.off()


# PREDICTIONS

# raster stack of 2020 metrics

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2020/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
total_rough <- raster(list[9])
urbanicity <- raster(list[11])
above20 <- raster(list[12])
betw021 <- raster(list[16])

# check that resolutions and extents now match (output should be TRUE)

compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)

# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
NN_2020_presence <- predict(modelNN, newdata = rasters_10, method = 'rf', mean = T, nc = 6, filename = "Neptis_nata_distribution", overwrite = T) # 8 cores for PC, 6 for Mac
toc()

plot(NN_2020_presence)


# raster stack of 2010 metrics

workingdirectory <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
setwd(workingdirectory)

path <- "C:/Users/TCB-Martha/OneDrive - The University Of Hong Kong/HKU/Butterflies/Data/HK_RS_data/2010/L2Processed_Martha/H4/"
list <- list.files(path, "*.tif", full=TRUE)
list

DTM <- raster(list[2])
open_patches <- raster(list[7])
total_rough <- raster(list[9])
urbanicity <- raster(list[11])
above20 <- raster(list[12])
betw021 <- raster(list[14])

# check that resolutions and extents now match (output should be TRUE)

compareRaster(DTM, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(total_rough, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(urbanicity, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(above20, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)
compareRaster(betw021, open_patches, extent=TRUE, rowcol=TRUE, crs=TRUE, res=FALSE, orig=FALSE,
              rotation=TRUE, values=FALSE, stopiffalse=TRUE, showwarning=FALSE)

# first, resample rasters to 5m so that spatial resolutions match the open_area raster
urbanicity <- resample(urbanicity, open_patches, method="bilinear")

# PREDICTIONS
# use RF method only to simplify processing

DTM_10 <- aggregate(DTM, fact = 2, fun = "mean")
total_rough_10 <- aggregate(total_rough, fact = 2, fun = "mean")
urbanicity_10 <- aggregate(urbanicity, fact = 2, fun = "mean")
above20_10 <- aggregate(above20, fact = 2, fun = "mean")
betw021_10 <- aggregate(betw021, fact = 2, fun = "mean")
open_patches_10 <- aggregate(open_patches, fact = 2, fun = "sum")

rasters_10 <- raster::stack(betw021_10, above20_10, total_rough_10, DTM_10, open_patches_10, urbanicity_10)
names(rasters_10) <- c("X0.2.1.m.density....", "X.20.m.density....", "Total.veg.roughness..m.", "Elevation..m.", "Open.patches..count.", "Land.cover")


library(tictoc)
tic()
NN_2010_presence <- predict(modelNN, newdata = rasters_10, method = 'rf', mean = T, nc = 10, filename = "Neptis_nata_distribution_2010", overwrite = T) # 10 cores for PC, 6 for Mac
toc()

plot(NN_2010_presence)

