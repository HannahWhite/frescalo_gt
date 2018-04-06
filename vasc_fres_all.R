##################################################################################
### Estimated species richness for all vascular plants (native and non-native) ###
##################################################################################

## Hannah White


### FRESCALO Analysis

library(sparta)

vasc_full <- read.csv('G:/Postdoc Grassland Resilience/Biological recordings datasets/vasc_full.csv', header = TRUE)


## get lat and long for vasc full

## Run os2eastnorth function
source('G:\\Postdoc Grassland Resilience\\R functions\\os2eastnorth.R')

# Convert Irish Grid to Lat Long
LatLong <- os2eastnorth(vasc_full$GR)

# Add easting and northing to frescalo estimates
vasc_full$east <- LatLong$en[LatLong$index, 1] # extracts index of unique location from list
vasc_full$north <- LatLong$en[LatLong$index, 2]

## calculate weights for FRESCALO

sites <- unique(data.frame(vasc_full$GR, vasc_full$east, vasc_full$north))
names(sites) <- c('GR', 'east', 'north')
site_dist <- dist(sites[,2:3])
site_dist <-as.matrix(site_dist, labels=TRUE)
rownames(site_dist) <- colnames(site_dist) <- sites$GR

# reshape matrix to df for frescalo
library(reshape2)
my_dist <- subset(melt(site_dist))

## Resample CORINE data

library(raster)
library(rgdal)
library(sp)


corine <- raster('G:\\Postdoc Grassland Resilience\\HannahData\\CORINE_IE.grd')

test = rasterToPoints(corine, fun=function(x){x!=0}, spatial=T)
test_grid = spTransform(test, CRS("+init=epsg:29903"))  # Plotting these coords looks like all of Ireland
corine_sp = as.data.frame(test_grid)  # Plotting these coords looks like all of Ireland

# An array to hold the number of land use classes in each hectad
corine_sites = array(0, dim=c(nrow(sites),44))
colnames(corine_sites) = paste('clc',c(1:44),sep='_')

for (s in 1:nrow(sites)) {
  # Find centres of corine raster cells that lie within each hectad (some cells may cross hectad boundary but this is ignored)
  corine_sub = subset(corine_sp, x>=sites$east[s] & x<sites$east[s]+10^4 & y>=sites$north[s] & y<sites$north[s]+10^4)
  
  # Count number of corine classes within this hectad
  corine_sites[s,] = sapply(c(1:44),FUN=function(x) {sum(x==corine_sub$g100_clc12_V18_5)}) 
}

corine_sites_GR = cbind(sites, corine_sites)

#### Use CORINE to create weights

corine10 <- as.data.frame(corine_sites_GR)

# rename CORINE land classes
names(corine10) <- c('GR', 'east', 'north', 'Urban1', 'Urban2',
                     'Industrial', 'Road',
                     'Port', 'Airport',
                     'Quarry', 'Dump',
                     'Construction', 'Green.Urban',
                     'Sport', 'Arable1',
                     'Arable2', 'Rice',
                     'Vineyard', 'Orchard',
                     'Olive', 'Pasture',
                     'Annual.Crops', 'Cultivation',
                     'Natural', 'AgroForestry',
                     'BroadLeaf', 'Conferous',
                     'Mixed', 'Grasslands',
                     'Moors', 'Sclero',
                     'Trans.Woodland', 'Dunes',
                     'Rocks', 'Sparse',
                     'Burnt', 'Glacier',
                     'Marsh', 'Peat',
                     'Saltmarsh', 'Salines',
                     'Intertidal', 'Water.Course',
                     'Water.Body', 'Coastal.Lagoon',
                     'Estuary', 'Sea')

corine10[,4:47] <- corine10[,4:47]/10000


### Aggregate corine10 to broader land classes

Urban <- rowSums(data.frame(corine10$Urban1, corine10$Urban2))

Industrial <- rowSums(data.frame(corine10$Industrial, corine10$Road, 
                                 corine10$Port, corine10$Airport))

Quarry <- rowSums(data.frame(corine10$Quarry, corine10$Dump,
                             corine10$Construction))

ArtificialVeg <- rowSums(data.frame(corine10$Green.Urban, corine10$Sport))

Arable <- rowSums(data.frame(corine10$Arable1, corine10$Arable2,
                             corine10$Rice))

PermanentCrops <- rowSums(data.frame(corine10$Vineyard, corine10$Orchard,
                                     corine10$Olive))

Pasture <- corine10$Pasture

HeteroAg <- rowSums(data.frame(corine10$Annual.Crops, corine10$Cultivation,
                               corine10$Natural, corine10$AgroForestry))

Forest <- rowSums(data.frame(corine10$BroadLeaf, corine10$Conferous,
                             corine10$Mixed))

Scrub <- rowSums(data.frame(corine10$Grasslands, corine10$Moors,
                            corine10$Sclero, corine10$Trans.Woodland))

OpenGround <- rowSums(data.frame(corine10$Dunes, corine10$Rocks,
                                 corine10$Sparse, corine10$Burnt, corine10$Glacier))

InlandWet <- rowSums(data.frame(corine10$Marsh, corine10$Peat))

MaritimeWet <- rowSums(data.frame(corine10$Saltmarsh, corine10$Salines,
                                  corine10$Intertidal))

InlandWater <- rowSums(data.frame(corine10$Water.Course, corine10$Water.Body))

Marine <- rowSums(data.frame(corine10$Coastal.Lagoon, corine10$Estuary, corine10$Sea))


my_corine <- data.frame(corine10$GR, Urban, Industrial, Quarry, ArtificialVeg,
                        Arable, PermanentCrops, Pasture, HeteroAg, Forest,
                        Scrub, OpenGround, InlandWet, MaritimeWet, 
                        InlandWater, Marine)

rm(corine_sites, cedar_use, corine_sites_GR,nbdc_use, nbdc_vasc, corine_sub, corine10, site_dist,sites,  Arable, ArtificialVeg, Forest, gr.length, HeteroAg, Industrial,
   InlandWater, InlandWet, Marine, MaritimeWet, OpenGround, Pasture, PermanentCrops, Quarry)
rm(test, test_grid, Urban, Scrub)

# create site weights for frescalo
my_weights <- createWeights(distances = my_dist,
                            attributes = my_corine)

my_fresc_path <- 'C:\\Users\\Hannah\\Documents\\Frescalo_3a_windows.exe'

### Run FRESCALO

my_time <- data.frame(start = 1970, end = 2017)


# specify where to save results
my_folder <- 'G:/Postdoc Grassland Resilience/Frescalo_out/vasc_ire'

frescalo_results <- frescalo(Data = vasc_full,
                             frespath = my_fresc_path,
                             time_periods = my_time,
                             site_col = 'GR',
                             sp_col = 'SpecName',
                             year = 'SampleYear',
                             Fres_weights = my_weights,
                             sinkdir = my_folder,
                             phi = 0.78)




