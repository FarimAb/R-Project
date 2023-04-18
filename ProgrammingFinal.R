#############
### START ###
#############

install.packages("sp")
install.packages("raster")
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("tidyr")
install.packages("dplyr")
install.packages("tibble")
install.packages("sf")
install.packages("terra")
install.packages("mapview")
install.packages("mapedit")
install.packages("caret")
install.packages("Rstoolbox")
install.packages("devtools")
install.packages("vctrs")
install.packages("scales")


library(sp)
library(raster)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(dplyr)
library(tibble)
library(sf)
library(terra)
library(mapview)
library(mapedit)
library(caret)
library(RStoolbox)
library(devtools)
library(vctrs)
library(scales)


#set seed to get comparable results
set.seed(1404)

#set working directory 
setwd("C:/Users/farim/Documents/SemeseterOne/ProgrammingAssignment")

#raster data
SARVH <- raster("Input/2020/SARVH_2020_s.tif")
SARVV <- raster("Input/2020/SARVV_2020_s.tif")
SDVV <- raster("Input/2020/SARVV_2020_s_stddev.tif")
SDVH <- raster("Input/2020/SARVH_2020_s_stddev.tif")
true_rgb <- brick("Input/2020/true_color_2020.tif")
ndvi <- raster("Input/2020/summer_2020_ndvi.tif")
ndwi <- raster("Input/2020/summer_2020_ndwi.tif")

#rename bands depending on rasters and bands used
sen_stack <- brick(SARVH, SARVV, SDVV, SDVH, true_rgb, ndvi, ndwi)
bands <- c("VV_mean", "VH_mean", "VH_sd","VH_sd", "B4", "B3", "B2", "ndvi", "ndwi")
names(sen_stack) <- bands

#create samples 
file_samples <- "Input/samples_winter_TUK.gpkg"
if(!file.exists(file_samples)){
  # sampling:
  #land
  land <- drawFeatures(
    map = viewRGB(sen_stack, r=1, maxpixels = 1000000)
  )
  
  # water
  shoreline <- drawFeatures(
    map = viewRGB(sen_stack, r=1, maxpixels = 1000000)
  )
  
  # ice
  water <- drawFeatures(
    map = viewRGB(sen_stack, r=2, maxpixels = 1000000)
  )
  
  
  # add landuse attribute
  land$landcover <- "land"
  shoreline$landcover <- "shoreline"
  water$landcover <- "water"
  
  
  # codify landcover as integers using mapvalues
  labeled_cover <- rbind(land, shoreline, water)
  labeled_cover$classid <- as.numeric(
    plyr::mapvalues(labeled_cover$landcover, 
                    from = unique(labeled_cover$landcover), 
                    to = 1:length(unique(labeled_cover$landcover)))
  )
  
  # save
  st_write(labeled_cover, "Input/samples_winter_TUK.gpkg")
}else{
  labeled_cover <- st_read("Input/samples_summer_TUK.gpkg")
}

#load polygone
# load polygons
labeled_cover <- st_transform(labeled_cover, st_crs(sen_stack))
labeled_cover$resp_var <- labeled_cover$classid


# let's randomly select some points in our polygons and save there labels to them
labeled_points <- list()
for(i in unique(labeled_cover$resp_var)){
  message(paste0("Sampling points from polygons with resp_var=", i))
  
  # sample points for polygons of resp_var = i
  labeled_points[[i]] <- st_sample(
    x = labeled_cover[labeled_cover$resp_var == i,], 
    size = 100
  )
  labeled_points[[i]] <- st_as_sf(labeled_points[[i]])
  labeled_points[[i]]$resp_var <- i
}
labeled_points <- do.call(rbind, labeled_points)

# extract features from points
# we can use our RS data for this since we do not have any ground-truth from site
sen_stack <- normImage(sen_stack)
sen_stack <- rescaleImage(sen_stack, ymin = 0, ymax = 1)

# extract features and label them with our response variable!
unlabeled_features <- raster::extract(sen_stack, labeled_points, df = T)
unlabeled_features <- unlabeled_features[,-1] # no ID column needed
labeled_features <- cbind(
  resp_var = labeled_points$resp_var,
  unlabeled_features
)

# remove duplicates (in case multiple points fall into the same pixel)
dupl <- duplicated(labeled_features)
which(dupl)
length(which(dupl)) # number of duplicates in labeled features that we need to remove!
labeled_features <- labeled_features[!dupl,]

# x = features
x <- labeled_features[,2:ncol(labeled_features)] # remove ID column
y <- as.factor(labeled_features$resp_var) #we want caret to treat this as categories, thus factor
levels(y) <- paste0("class_", levels(y))


# fit the model, here a ranodm Forest
model <- train(
  x = x,
  y = y, 
  trControl = trainControl(
    p = 0.75, # percentage of samples used for training, rest for validation
    method  = "cv", #cross validation
    number  = 5, # 5-fold
    verboseIter = TRUE, # progress update per iteration
    classProbs = TRUE # probabilities for each example
  ),
  #preProcess = c("center", "scale"), #center/scale if not done by you on the raster (see previous code rescl)
  method = "rf" # used algorithm
) 

# performance
model
confusionMatrix(model)

# predict classes for raster stack
summer_2020 <- predict(sen_stack, model, type='raw')


# write 
writeRaster(summer_2020, filename = "output/summer_2020.tif",
            datatype = "INT1U")

cols <- c("green", "black", "blue")

# map
mapview(summer_2020, col.regions = cols)

# plot for visualization
ggplot() + ggR(summer_2018, geom_raster = T, ggLayer = T) +
  scale_fill_manual(values = cols, 
                    name = "Classses") +
  coord_sf(crs = st_crs(sen_stack), datum = st_crs(4326))


#list for class names, might need to be changed
classes <- list("land", "shoreline", "water", "NA")

#add Year, Month and class names 
df_winter22 <- as.data.frame(winter_2022) %>%
  group_by(value) %>%
  tally() %>%
  mutate(area = n * res(winter_2022)[1] * res(winter_2022)[2])

#divide by 1mio
df_winter22$Areakm2 <- with(df_winter22, area/1000000)

#df_winter17 <- data.frame(freq(winter_2017))
df_winter22$Class <- classes
df_winter22['Year'] <- '2022'
df_winter22['Season'] <- 'Winter'

#save df
df_winter22 <- apply(df_winter22,2,as.character)
write.csv(df_total, "output/df_total.csv")

#bind the dataframes
df_total <- rbind(df_summer17, df_summer18, df_summer19, df_summer20, df_summer21,
                  df_summer22, df_winter16, df_winter17, df_winter18, df_winter19,
                  df_winter20, df_winter21, df_winter22)

#remove NA's from df
df_total <- na.omit(df_total)










#set woring directory
setwd("C:/Users/ellys/Desktop/EAGLE/Permafrost/output")

#load df
df_agg <- read.csv("df_total.csv")

#split df by season
df_test <- split(df_agg, df_agg$Season)

#df for winter and summer
df_summer <- as.data.frame(df_test$Summer)
df_winter <- as.data.frame(df_test$Winter)

#get area difference for every class
summer_area_land <- df_summer$Areakm2[16] - df_summer$Areakm2[1]
summer_area_shoreline <- df_summer$Areakm2[17] -  df_summer$Areakm2[2]  
summer_area_water <- df_summer$Areakm2[18] - df_summer$Areakm2[3]         

#same shit for winter
winter_area_land <- df_winter$Areakm2[16] - df_winter$Areakm2[1]
winter_area_shoreline <- df_winter$Areakm2[17] -  df_winter$Areakm2[2]  
winter_area_water <- df_winter$Areakm2[18] - df_winter$Areakm2[3] 

#bind new dataframes
df_summer_change <- rbind(summer_area_land, summer_area_shoreline, summer_area_water)
df_winter_change <- rbind(winter_area_land, winter_area_shoreline, winter_area_water)



##Thermokarst lakes in permafrost area##





#set the directory of dataframe and call it #

setwd("C:/Users/farim/Documents/SemeseterOne/ProgrammingAssignment/Output")
PR <- read.csv("df_total.csv")
PR

#Adding the Area Column
#Perma$Area_m2 <- with(Perma,(count * 100))
#Perma$Area_k2 <- with(Perma, Area_m2/1000)
#Perma = na.omit(Perma)

#Set the directory again, save the png and plot the final results

setwd("C:/Users/farim/Documents/SemeseterOne/ProgrammingAssignment")
png("Thermokarst Lake.png", width = 1800 , height = 1800)
par(mar = c(5,5,5,5))
ggplot(data = Perma, aes(x = Year , y = Areakm2 , linetype = Season, color = Class), group=1) +
  geom_line(size= 1 )+
  labs(title="Thermokarst Lake Dynamics 2016-2023", 
       subtitle="Tuktoyaktuk, Northwest Territories, CA") +
  theme_classic() +
  theme(plot.margin = margin(5, 5, 5, 5)) +
  scale_color_manual(values = c("#ae7346", "#7fc097", "#4680ae"))
dev.off() 

plot

library(raster)
library(tmap)

