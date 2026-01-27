library(terra)
library(readr)
library(sf)
library(dplyr)
library(sswids)
library(landscapemetrics)
library(purrr)
library(tidyverse)
library(daymetr)
library(readxl)
library(MODISTools)
library(leaflet)
library(lubridate)
library(mgcv)
library(rlandfire)
library(data.table)
library(CropScapeR)


buffers <- c(100,500,1000,2500,5000) # buffer sizes in meters, buffer size 5000 not enough memory
df <- read_rds(paste(getwd(), "BearCPUESpring2019-2024.rds", sep="/")) #4/15-6/30
BEARCPUE <- df%>%filter(ppn_classified >= 0.95  & year > 2018)#& days_active > 4
BEARCPUE2 <- BEARCPUE%>%group_by(cam_site_id, year)%>%add_tally(name = "n.occ")%>%
  group_by(cam_site_id)%>%add_tally(name="n.occ.total")
BEARCPUE3 <- BEARCPUE2%>%filter(n.occ > 5)%>%select(1:3,4,5,8,13,14,15:18)#1:3,5,8,14,15,16,17,18
camsites <- BEARCPUE3%>%distinct(cam_site_id, geometry)
BEARCPUE3%>%group_by(year)%>%summarise(n=n_distinct(cam_site_id))
camsiteyears <- BEARCPUE3%>%ungroup()%>%distinct(cam_site_id, year, geometry)%>%group_split(year)%>%setNames(2019:2024)

# camsite_buffers <- st_buffer(camsites, dist = 500,endCapStyle = "SQUARE")
# st_write(camsite_buffers, "./Bearcamsite_buffers.shp")
# st_write(camsites, "./BearCamsites.shp")
leaflet() %>% 
  addProviderTiles('Esri.WorldImagery') %>%
  addCircleMarkers(data=st_transform(camsites, 4326), fillColor = "white", fillOpacity = 1, stroke=F, radius=3, label=camsites$cam_site_id)


# set min/max years for obtaining data
min_year <- 2019
max_year <- 2024

# season start/end dates (in -MM-DD format)
# start with this and consider early spring and fall for getting newborns and yearlings 
min_date <- '-04-15'
max_date <- '-06-30'


# create data frame of seasons for data filtering
seasons_df <-
  create_season_dates(
    min_date = min_date,
    max_date = max_date,
    years = c(2019,2024)
  )
check_season_dates(seasons_df)

num_occasions <- 11
day_occasion_df <-
  seasons_df %>%mutate(year=season+2018)%>%
  group_by(year) %>%
  nest() %>%
  # create date sequence for each year
  mutate(date = map(data, date_sequence)) %>%
  unnest(date) %>%
  select(-data) %>%
  # using row_number give day of season starting with day 1
  mutate(day_of_season = row_number()) %>%
  ungroup() %>%
  # split season into equal intervals (1-day, 3-day, ...1 week)
  # ntile() assigns each day into a sampling occasion
  mutate(
    occ = ntile(day_of_season, num_occasions)
  )

occdates <- day_occasion_df%>%group_by(year, occ)%>%summarise("start_date"=min(date), "end_date"=max(date))

BEARCPUE3 <- left_join(BEARCPUE3, occdates, by=c("year", "occ"))
################################################################################################
####                            Grid                                                         ###
################################################################################################
#will need to pull environmental covariates for all grid cells in Wisconsin
Wisconsin <- st_read("C:/Users/wildeefb/Documents/GeoSpatial/VectorLayers/Wisconsin_State_Boundary_24K.shp")
counties <- get_spatial_data("counties")
plot(st_geometry(Wisconsin))
cell_area <- units::as_units(8.5, "km^2") # target grid size
cell_size <- units::as_units(8.5, "km")
grid_spacing <- sqrt(2*5000/sqrt(3)) # size of hexagon calculated from area
WiscGrid <- st_make_grid(Wisconsin, cellsize = cell_area, square=FALSE)%>%st_as_sf()%>%mutate("cell"=row_number())
WiscGridsquare <- st_make_grid(Wisconsin, cellsize = cell_size)%>%st_as_sf()%>%mutate("cell"=row_number())
plot(st_geometry(Wisconsin))
plot(st_geometry(WiscGrid), add=TRUE)
WiscGrid2 <- st_intersection(Wisconsin, WiscGrid)
WiscGridSquare2 <- st_intersection(Wisconsin, WiscGridsquare)
plot(st_geometry(Wisconsin))
plot(st_geometry(WiscGrid2), add=TRUE)
plot(st_geometry(WiscGridSquare2), add=TRUE)
countygrids <- st_join(counties, WiscGridSquare2)
camsites <- st_join(camsites, WiscGrid2[,7:8])

#######################################################################################################################
####                                       Land Cover                                                              ####
#######################################################################################################################
Wiscland1 <- get_spatial_data(layer_name = 'wiscland2', level = 1)
Wiscland2 <- get_spatial_data(layer_name = 'wiscland2', level = 2)
Wiscland3 <- get_spatial_data(layer_name = 'wiscland2', level = 3)
crs(Wiscland2) #EPSG\",9001
Wiscland1.3071 <- terra::project(Wiscland1, "EPSG:3071")
Wiscland2.3071 <- project(Wiscland2, "EPSG:3071")
Wiscland3.3071 <- project(Wiscland3, "EPSG:3071")
WisclandGuide <- read_xlsx("C:/Users/wildeefb/Documents/GeoSpatial/wiscland2/user_guide/Wiscland2 Color Scheme.xlsx",
                           range = "B2:C70", .name_repair = make.names)



# lc <- extract(project(Wiscland2, "EPSG:3071"), st_coordinates(dfAll$geometry))
# st_crs(st_as_sf(dfAll, wkt="geometry", crs=3071))
# st_as_sf(dfAll, wkt="geometry", crs=3071)
# st_coordinates(dfAll$geometry)
# table(dfAll$lc)
# camsites$buffer <- st_buffer(camsites$geometry, 500)
# extract(project(Wiscland2, "EPSG:3071"), camsites$buffer, fun=table)
# st_crs(Wiscland2)



# Wiscland2 has 30m resolution


# Wiscland 3 prop land cover
lm_output <- 
  buffers %>% 
  set_names() %>% 
  # produce a dataframe after this is all done
  map_dfr( 
    ~sample_lsm(
      # raster layer
      landscape = Wiscland3,
      # camera locations
      y = camsites,
      # get landcover class level metrics
      level = "class",
      # return NA values for classes not in buffer
      # all_classes = TRUE, 
      # camera site IDs here
      plot_id = camsites$cam_site_id,
      # can do multiple metrics at once
      what = 'lsm_c_pland',
      # buffer sizes to use
      size = ., 
      # default is square buffer
      shape = "circle", 
      # turn warnings on or off
      verbose = FALSE 
    ), 
    # get buffer size column in the output
    .id = "buffer_size"
  )

lm_output <- left_join(lm_output, WisclandGuide, by= join_by(class == dn.label))
lm_output$label <- gsub(pattern = "\\W", replacement = "", x = lm_output$label)

# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
lm_output <- 
  lm_output %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label, buffer_size),
    values_from = c(value),
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  ) %>%
  # clean up names
  rename(cam_site_id = plot_id)


# join pland back to camera sf object
forestLCs <- c("AspenPaperBirch", "RedMaple", "Oak", "CentralHardwoods", "NorthernHardwoods","AspenForestedWetland", "BottomlandHardwoods", "SwampHardwoods",
               "MixedDeciduousConiferousForest", "MixedDeciduousConiferousForestedWetland")
developedLCs <- c("DevelopedHighIntensity","DevelopedLowIntensity")

lm_output2 <- cbind(lm_output, lapply(buffers, function(x) rowSums(lm_output[,paste0(forestLCs, "_", x)])))
forestcols <- grep(pattern = "c\\(", x = colnames(lm_output2))
colnames(lm_output2)[forestcols] <- paste0("Forest", "_",buffers)

ForestProp <- lm_output2[,c(1,forestcols)]
camsites <- left_join(camsites, ForestProp, by="cam_site_id")

lm_output2 <- cbind(lm_output, lapply(buffers, function(x) rowSums(lm_output[,paste0(developedLCs, "_", x)])))
devcols <- grep(pattern = "c\\(", x = colnames(lm_output2))
colnames(lm_output2)[devcols] <- paste0("Developed", "_",buffers)

DevelopedProp <- lm_output2[,c(1,devcols)]
camsites <- left_join(camsites, DevelopedProp, by="cam_site_id")



#shared Wiscland level 2-3 landcover class: c(DevelopedHighIntensity, DevelopedLowIntensity, Cranberries, 
#                                             MixedDeciduousConiferousForest, OPENWATER, FloatingAquaticHerbaceousVegetation,
#                                             BARREN, SHRUBLAND)

# lcvars <- grep(pattern = "^[A-Z].*_\\d+",x = colnames(RUGRCPUESaved),value=TRUE)
# lcvars2 <- as.data.frame(str_split_fixed(lcvars, pattern = "_", n=2))
# lcvars3 <- left_join(lcvars2, WisclandGuide, by=join_by(V1==label))

#Wiscland 2 prop land cover
lm_output <- 
  buffers[1:3] %>% 
  set_names() %>% 
  # produce a dataframe after this is all done
  map_dfr( 
    ~sample_lsm(
      # raster layer
      landscape = Wiscland2,
      # camera locations
      y = camsites,
      # get landcover class level metrics
      level = "class",
      # return NA values for classes not in buffer
      # all_classes = TRUE, 
      # camera site IDs here
      plot_id = camsites$cam_site_id,
      # can do multiple metrics at once
      what = 'lsm_c_pland',
      # buffer sizes to use
      size = ., 
      # default is square buffer
      shape = "circle", 
      # turn warnings on or off
      verbose = FALSE 
    ), 
    # get buffer size column in the output
    .id = "buffer_size"
  )

lm_output <- left_join(lm_output, WisclandGuide, by= join_by(class == dn.label))
lm_output$label <- gsub(pattern = "\\W", replacement = "", x = lm_output$label)

# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
lm_output <- 
  lm_output %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label, buffer_size),
    values_from = c(value),
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  ) %>%
  # clean up names
  rename(cam_site_id = plot_id)





#mean patch area in hectares
lsm_l_condent
cats(Wiscland3)
# Define the lookup table (original values to new values)
lookup <- data.frame(id=c(1, 2, 3), new_cover=c("forest", "water", "urban"))

# Reclassify using subst
r_reclassified <- subst(r, lookup$id, lookup$new_cover)

patcharea <- 
  buffers %>% 
  set_names() %>% 
  # produce a dataframe after this is all done
  map_dfr( 
    ~sample_lsm(
      # raster layer
      landscape = Wiscland1,
      # camera locations
      y = camsites,
      # get landcover class level metrics
      level = "landscape",
      # return NA values for classes not in buffer
      # all_classes = TRUE, 
      # camera site IDs here
      plot_id = camsites$cam_site_id,
      # can do multiple metrics at once
      what = 'lsm_c_area_mn',
      # buffer sizes to use
      size = ., 
      # default is square buffer
      shape = "circle", 
      # turn warnings on or off
      verbose = FALSE 
    ), 
    # get buffer size column in the output
    .id = "buffer_size"
  )

patcharea <- left_join(patcharea, WisclandGuide, by= join_by(class == dn.label))
patcharea$label <- gsub(pattern = "\\W", replacement = ".", x = patcharea$label)

# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
patcharea <- 
  patcharea %>%
  # MAY NEED TO ADD distinct() HERE???
  distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label, buffer_size),
    values_from = c(value),
    names_prefix = "meanpatch",
    values_fill = 0
  ) %>%
  # clean up names
  rename(cam_site_id = plot_id)

camsites <- left_join(camsites, patcharea, by="cam_site_id")


########################################################################################################
####                                        DayMet                                                  ####
########################################################################################################

campts <- data.frame("site_id"=camsites$cam_site_id,st_coordinates(st_transform(camsites, crs=4326)))
get_daymet <- function(i, locationsdf, starttime, endtime){
  
  temp_lat <- locationsdf[i, ] %>% pull(Y)
  temp_lon <- locationsdf[i, ] %>% pull(X)
  temp_site <- locationsdf[i, ] %>% pull(site_id)
  
  temp_daymet <- download_daymet(
    lat = temp_lat,
    lon = temp_lon,
    start = starttime,
    end = endtime
  ) %>% 
    #--- just get the data part ---#
    .$data %>% 
    #--- convert to tibble (not strictly necessary) ---#
    as_tibble() %>% 
    #--- assign site_id so you know which record is for which site_id ---#
    mutate(site_id = temp_site) %>% 
    #--- get date from day of the year ---#
    mutate(date = as.Date(paste(year, yday, sep = "-"), "%Y-%j"))
  
  return(temp_daymet)
}  

daymet_all_points <- lapply(1:nrow(campts), function (i) {get_daymet(i, locationsdf = campts, starttime = 2019, endtime = 2024) %>% 
  #--- need to combine the list of data.frames into a single data.frame ---#
  bind_rows()
  print(i/nrow(campts)*100)})

daymet_all_points$GDD <- sapply(1:nrow(daymet_all_points), function (i) 
       {if (tmin < 10) {tmin = 10} 
       if(tmax > 30) {tmax = 30} 
       (tmin + tmax /2) - 10})


###########################################################################################################
####                                            GEE NDVI                                               ####
###########################################################################################################
EVIbuffers <- st_buffer(camsites, dist=1000, endCapStyle = 'SQUARE')
st_write(EVIbuffers, "./EVIBuffers.shp")
EVIbuffers <- st_read("./EVIBuffers.shp")
#Modis 16 day EVI, not available for 2023 or is it 463m pixel size 
#https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
EVI_MODIS_16DAY <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/MODIS16DayGlobal250_2019-2024Ext.csv") #named for start time
EVI_MODIS_16DAY$year <- year(as.Date(EVI_MODIS_16DAY$date))
EVI_MODIS_16DAY$DOY <- yday(EVI_MODIS_16DAY$date)
EVI_MODIS_16DAY$camsiteyear <- paste0(EVI_MODIS_16DAY$camsiteid, EVI_MODIS_16DAY$year)
EVI_MODIS_16DAYpercamyear <- EVI_MODIS_16DAY%>%group_by(camsiteid, year)%>%summarise(n=n())
  
  #try a gam instead of logistic equation for more flexibility in fits to data?
siteyear <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == unique(EVI_MODIS_16DAY$camsiteyear)[1],]  
gam <- gam(meanEVI ~ s(DOY, k = 5), data = siteyear, method ="REML")
year <- siteyear$year[1]  
predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
  predsiteyear$DOY <- yday(predsiteyear$date)
  predsiteyear$preds <- predict(gam, newdata=predsiteyear)
  siteyear$preds <- predict(gam)
  ggplot(siteyear, aes(x=DOY, y=meanEVI)) + geom_point(aes(color="observed")) + 
    geom_point(aes(y=preds, color="predicted")) + geom_smooth(aes(y=preds, color="predicted")) +
    scale_color_manual(values = c("orange", "deepskyblue3"))
  
#remove old EVIModissiteyearocc before running each time if necessary
  rm(EVIModissiteyearocc)
  for(i in 1:length(unique(EVI_MODIS_16DAY$camsiteyear))){
    siteyear <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == unique(EVI_MODIS_16DAY$camsiteyear)[i],]
    curvefit <- try(gam(meanEVI ~ s(DOY, k = 5), data = siteyear, method ="REML"))
    if(inherits(curvefit, "try-error")){
        predsiteocc <- data.frame(occ=1:11, meanEVI=rep(NA,11), cam_site_id=siteyear$camsiteid[1], year=siteyear$year[1])
      if(!exists("EVIModissiteyearocc")){
        EVIModissiteyearocc <- predsiteocc
      }
      else{
        EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
      }
      next
    } else {
      year <- siteyear$year[1]
      predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
      predsiteyear$DOY <- yday(predsiteyear$date)
      predsiteyear$preds <- predict(curvefit, newdata=predsiteyear)
      predsiteyear <- left_join(predsiteyear, day_occasion_df, by="date")
      predsiteocc <- predsiteyear%>%group_by(occ)%>%summarise(meanEVI=mean(preds))%>%mutate("cam_site_id"=siteyear$camsiteid[1], "year"=year)
      
      if(!exists("EVIModissiteyearocc")){
        EVIModissiteyearocc <- predsiteocc
      }
      else{
        EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
      }
    }
  }
  table(is.na(EVIModissiteyearocc$meanEVI))
  what <- EVIModissiteyearocc[is.na(EVIModissiteyearocc$meanEVI),]
  saveRDS(EVIModissiteyearocc, "./EVIsiteyearocc.rds")
  EVIModisBearCPUE <- left_join(BEARCPUE3, EVIModissiteyearocc, by=c("cam_site_id", "year", "occ"))
  table(is.na(EVIModisBearCPUE$meanEVI)) #2250 false
  BEARCPUE4 <- left_join(BEARCPUE3, EVIModissiteyearocc, by=c("cam_site_id", "year", "occ"))
  
##########################################################################################
####                               Human Footprint Index                              ####
##########################################################################################
HFI <- rast("C:/Users/wildeefb/Documents/GeoSpatial/HFI2020.tif") #300m resolution
crs(HFI)
Wisconsin <- st_read("C:/Users/wildeefb/Documents/GeoSpatial/VectorLayers/Wisconsin_State_Boundary_24K.shp")
WisconsinHFI <- crop(HFI, st_transform(Wisconsin, 4326), mask=TRUE)
plot(WisconsinHFI)
WisconsinHFI <- focal(WisconsinHFI, w=9, fun=mean, na.policy="only", na.rm=T)
  
WisconsinHFI3071 <- project(WisconsinHFI, "epsg:3071")

camsite_buffers <- lapply(buffers, function(x) vect(st_buffer(camsites,  dist=x)))
names(camsite_buffers) <- buffers
HFIcamsites <- lapply(camsite_buffers, function (x) terra::extract(x=WisconsinHFI3071, y=x, fun=mean, na.rm=TRUE))
HFIcamsites <- lapply(HFIcamsites, function (x) x%>%mutate(cam_site_id=camsites$cam_site_id))
HFIcamsites <- rbindlist(HFIcamsites, idcol = TRUE)
HFIcamsites <- pivot_wider(data = HFIcamsites, names_from = .id, values_from=focal_mean)
colnames(HFIcamsites)[3:7] <- paste("HFI",colnames(HFIcamsites)[3:7], sep="_")

leaflet() %>% 
  addProviderTiles('Esri.WorldImagery') %>%
  addCircleMarkers(data=st_transform(camsites, 4326), fillColor = "black", fillOpacity = 1,   stroke=F, radius=5, popup=camsites$cam_site_id)%>%
  addRasterImage(WisconsinHFI, opacity = 0.8, group="HFI", maxBytes = Inf)

camsites <- left_join(camsites, HFIcamsites, by="cam_site_id")

BEARCPUESaved3 <- left_join(BEARCPUESaved3, HFIcamsites, by="cam_site_id")
NA.HFI <- HFIcamsites[which(is.na(HFIcamsites$HFI2020)),]
saveRDS(BEARCPUESaved3, "ModelingDF.rds")
  
#################################################################################################
####                             Roads                                                       ####
#################################################################################################
majorroads <- get_spatial_data('major_roads')
majorroads3071 <- st_transform(majorroads, crs=3071)
countyroads <- get_spatial_data('county_local_roads')
countyroads3071 <- st_transform(countyroads, crs=3071)


  
##########################################################################################
#####                               Latitude                                          ####
##########################################################################################
camsites$Lat <- st_coordinates(st_transform(camsites, crs=4326))[,2]

BEARCPUE4 <- left_join(BEARCPUE3, st_drop_geometry(camsites), by="cam_site_id")

saveRDS(BEARCPUE4, "ModelingDF.rds")

#####################################################################################################
###                                   R Landfire                                                 ####
#####################################################################################################
aoi <- getAOI(Wisconsin)
products <- "HDIST2023"
email <- "eli.wildey@wisconsin.gov"
projection <- 3071
resolution <- 90
path <- tempfile(fileext = ".zip")
lf_dir <- "C:/Users/wildeefb/Documents/GeopSpatial/LANDFIRE/HDist2023.zip"
hdist2023 <-landfireAPIv2(products = products,
                          aoi = aoi, 
                          email = email,
                          projection = projection, 
                          resolution = resolution,
                          path = path,
                          verbose = TRUE)
lf_dir <- file.path(tempdir(), "lf")
utils::unzip(path, exdir = lf_dir)
hdist <- terra::rast(list.files(lf_dir, pattern = ".tif$", 
                                full.names = TRUE, 
                                recursive = TRUE))
dbf <- list.files(lf_dir, pattern = ".dbf$",
                  full.names = TRUE,
                  recursive = TRUE)
dbf_tbl  <- foreign::read.dbf(dbf)
HDistLU <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/LANDFIRE/LF2024_HDist24.csv")
dbf2 <- left_join(dbf_tbl, HDistLU, by=join_by(Value ==VALUE))%>%mutate(EarlySuccess=ifelse(Value > 0, "0-10", "10+"))
levels(hdist) <- dbf2[,c(1,20)]
hdist <- addCats(hdist, value=dbf2[,c(5,7:10)])
cats(hdist)
activeCat(hdist)
levels(hdist)
plot(hdist)
activeCat(hdist) <- 1


hdist2 <- ifel(hdist > 1, 1, hdist)
plot(hdist2)
# Wiscland2 has 30m resolution

# Wiscland 3 prop land cover
lm_output <- 
  buffers %>% 
  set_names() %>% 
  # produce a dataframe after this is all done
  map_dfr( 
    ~sample_lsm(
      # raster layer
      landscape = hdist2,
      # camera locations
      y = camsites,
      # get landcover class level metrics
      level = "class",
      # return NA values for classes not in buffer
      # all_classes = TRUE, 
      # camera site IDs here
      plot_id = camsites$cam_site_id,
      # can do multiple metrics at once
      what = 'lsm_c_pland',
      # buffer sizes to use
      size = ., 
      # default is square buffer
      shape = "circle", 
      # turn warnings on or off
      verbose = FALSE 
    ), 
    # get buffer size column in the output
    .id = "buffer_size"
  )


lm_output$label <- ifelse(lm_output$class == 0, "10+", "0-10")

# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
lm_output <- 
  lm_output %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label, buffer_size),
    values_from = c(value),
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  ) %>%
  # clean up names
  rename(cam_site_id = plot_id)

lm_output <- lm_output%>%select(-matches("\\+"))
colnames(lm_output)[2:6] <- gsub(x = colnames(lm_output)[2:6], pattern = "0-10", "Dist")

#########################################################################################
####                             CropScape                                          #####
#########################################################################################
lapply(unique(BEARCPUE$year), function (x) GetCDLData(aoi = 55, year = x, type = 'f', format = 'raster', crs = '+init=epsg:4326',
                                                      save_path=paste0("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/Crops", x, ".tif")))#aoi=FIPS code for WI
GetCDLData(aoi = 55, year = 2019, type = 'f', format = 'raster', crs = '+init=epsg:4326',
           save_path=paste0("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/Crops", "2019", ".tif"))

CropRasts <- list.files("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/", full.names = TRUE)
CropRasts2 <- lapply(CropRasts, function (x) rast(x))
data("linkdata")
camsiteyearsCDL <- lapply(camsiteyears, function (x) st_transform(x, crs(CropRasts2[[1]])))

# Wiscland 3 prop land cover
buffers2=buffers %>% 
  set_names()
lm_output <- 
  mapply(x= CropRasts2, y=camsiteyearsCDL,
         # produce a dataframe after this is all done
         \(x,y) lapply(buffers2[2], 
                       \(z) sample_lsm(
                         # raster layer
                         landscape = x,
                         # camera locations
                         y = y,
                         # get landcover class level metrics
                         level = "class",
                         # return NA values for classes not in buffer
                         # all_classes = TRUE, 
                         # camera site IDs here
                         plot_id = y$cam_site_id,
                         # can do multiple metrics at once
                         what = 'lsm_c_pland',
                         # buffer sizes to use
                         size = z, 
                         # default is square buffer
                         shape = "circle", 
                         # turn warnings on or off
                         verbose = FALSE 
                       )
         ), SIMPLIFY=FALSE
  )
names(lm_output) <- 2019:2024

lm_output2 <- rbindlist(lapply(lm_output, function(x) rbindlist(x, idcol = "buffer_size")), idcol="year")
lm_output2 <- left_join(lm_output2, linkdata, by=join_by(class==MasterCat))
Corn <- lm_output2[grep(pattern = "Corn", x = lm_output2$Crop, ignore.case = TRUE),]
Corn2 <- Corn%>%group_by(year, buffer_size, plot_id)%>%summarise(CornProp=sum(value))
Corn2$buffer_size <- as.numeric(Corn2$buffer_size)
# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
Corn3 <- 
  Corn2 %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    names_from = buffer_size,
    names_prefix="Corn",
    values_from = CornProp,
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  ) %>%
  # clean up names
  rename(cam_site_id = plot_id)
Corn3$season <- as.numeric(Corn3$year)-2018
Corn3$year <- as.numeric(Corn3$year)
ModelingDFNew <- left_join(ModelingDF, st_drop_geometry(camsites[,c(1,8:12)]))
ModelingDFNew2 <- left_join(ModelingDFNew, Corn3)
colSums(is.na(ModelingDFNew2))
CornNA <- ModelingDFNew2[is.na(ModelingDFNew2$Corn1000),]
CornNA2 <- lm_output2[paste0(lm_output2$plot_id, lm_output2$year) %in% paste0(CornNA$cam_site_id, CornNA$year),]
ModelingDFNew2 <- ModelingDFNew2%>%mutate(across(matches("Corn"), ~replace_na(.x, 0)))
saveRDS(ModelingDFNew2, "./ModelingDFSpring.rds")
ModelingDFNew3 <- readRDS("./ModelingDFSpring.rds")
ModelingDFNew3 <- left_join(ModelingDFNew3, Corn3)%>%relocate(Corn500, .before = Corn1000)%>%
  mutate(across(matches("Corn500"), ~replace_na(.x, 0)))
saveRDS(ModelingDFNew3, "./ModelingDFSpring.rds")

BEARCPUE5 <- left_join(BEARCPUE4, Corn3)
BEARCPUE5 <- BEARCPUE5%>%mutate(across(matches("Corn"), ~replace_na(.x, 0)))

##########################################################################################
####                   save sitecovs and detcovs  for use in modeling                #####
##########################################################################################
camsites <- left_join(camsites, lm_output, by="cam_site_id")
camsites <- camsites[,-9]
saveRDS(camsites, "./sitevariables.rds")
BEARCPUE5 <- left_join(BEARCPUE4,  st_drop_geometry(camsites) ,by="cam_site_id")
saveRDS(BEARCPUE5, "./ModelingDF.rds")
############################################################################################
#####                             scrap                                                 ####
############################################################################################
activity_plot = df %>%
  filter(days_active>0) %>%
  group_by(occ) %>%
  summarise(bear_adult = sum(BEAR_ADULT_AMT), nsites=n_distinct(cam_site_id, season), amt.cam=bear_adult/nsites) %>%
  pivot_longer(cols=c(bear_adult))

#vlines <- c(17, 29)
ggplot(activity_plot, aes(x=occ,y=value)) +
  facet_wrap(~name, scales="free_y") +
  #geom_vline(xintercept = vlines, col = "forestgreen", size=2) +
  geom_point() +
  geom_line() +
  # geom_text(
  #       aes(x = 17, y = 200, label = "4/15-6/30"),
  #       size = 3, vjust = 0, hjust = 0, color = "black"
  #     ) +
  labs(title="Total triggers by occasion")





# join pland back to camera sf object, don't know what this is for
dfAll2 <- left_join(dfAll, lm_output, by="cam_site_id")
dfAll3 <- dfAll2[!is.na(dfAll2$BEAR_ADULT_AMT),]
dfAll3 <- filter(dfAll3, ppn_classified >= 0.95)
dfAll3$mean_date2 <- scale(dfAll3$mean_date^2)
dfAll3$mean_date <- scale(dfAll3$mean_date)
dfAll3$days_active <- scale(dfAll3$days_active)
ggplot(dfAll3, aes(x=mean_date, y=BEAR_ADULT_AMT)) + geom_point()
gsub(pattern="_100$", replacement="",grep(pattern="_100$",x=colnames(lm_output),value=TRUE))
Wiscland2key <- read_csv("C:/Users/wildeefb/Documents/GIT/sswids/vignettes/data_demo/Wiscland2_Key_Level2.csv")
left_join(Wiscland2)


#################################NDVI testing############################

#Landsat-8 8-day NDVI
NDVI <- rast("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/20190415NDVI.tif")
plot(NDVI)
values(NDVI)
NDVI <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/NDVI_reduceRegions.csv") #named for start time
NDVI$year <- year(as.Date(NDVI$date))
NDVI$DOY <- yday(NDVI$date)
NDVI$camsiteyear <- paste0(NDVI$camsiteid, NDVI$year)
cam0316_2019 <- NDVI[NDVI$camsiteid == "cam_0316" & NDVI$year == 2019,]
ggplot(cam0316_2019, aes(x=DOY, y=meanNDVI)) + geom_point() + geom_smooth(method = "lm", se = FALSE)
NDVIpercamyear <- NDVI%>%group_by(camsiteid, year)%>%summarise(n=n())
hist(NDVIpercamyear$n)
NDVI10 <- NDVIpercamyear[NDVIpercamyear$n >= 10,]
NDVI10year <- left_join(NDVI10, NDVI, by=c("camsiteid", "year"))
ggplot(NDVI10year, aes(x=DOY, y=meanNDVI)) + facet_wrap(.~camsiteid) + 
  geom_point() + geom_smooth(method = "lm", formula = y~x)

for(i in 1:length(unique(NDVI$camsiteyear))){
  siteyear <- NDVI[NDVI$camsiteyear == unique(NDVI$camsiteyear)[i],]
  curvefit <- try(nls(meanNDVI ~ SSlogis(DOY, Asym, xmid, scal), siteyear))
  if(inherits(curvefit, "try-error")){
    next
  } else {
    year <- siteyear$year[1]
    predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
    predsiteyear$DOY <- yday(predsiteyear$date)
    predsiteyear$preds <- predict(curvefit, newdata=predsiteyear)
    predsiteyear <- left_join(predsiteyear, day_occasion_df, by="date")
    predsiteocc <- predsiteyear%>%group_by(occ)%>%summarise(meanNDVI=mean(preds))%>%mutate("cam_site_id"=siteyear$camsiteid[1], "year"=siteyear$year[1])
    if(!exists("NDVIsiteyearocc")){
      NDVIsiteyearocc <- predsiteocc
    }
    else{
      NDVIsiteyearocc <- rbind(NDVIsiteyearocc, predsiteocc)
    }
  }
}

NDVIBearCPUE <- left_join(BEARCPUE3, NDVIsiteyearocc, by=c("cam_site_id", "year", "occ"))
table(is.na(NDVIBearCPUE$meanNDVI)) #27415 NAs

#Landsat 8 8-day EVI
#https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_COMPOSITES_C02_T1_L2_8DAY_EVI
EVI_LS8_8DAY <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/EVI_LS8DAY2019-2024Ext.csv") #named for start time
EVI_LS8_8DAY$year <- year(as.Date(EVI_LS8_8DAY$date))
EVI_LS8_8DAY$DOY <- yday(EVI_LS8_8DAY$date)
EVI_LS8_8DAY$camsiteyear <- paste0(EVI_LS8_8DAY$camsiteid, EVI_LS8_8DAY$year)
EVI_LS8_8DAYpercamyear <- EVI_LS8_8DAY%>%group_by(camsiteid, year)%>%summarise(n=n())
for(i in 1:length(unique(EVI_LS8_8DAY$camsiteyear))){
  siteyear <- EVI_LS8_8DAY[EVI_LS8_8DAY$camsiteyear == unique(EVI_LS8_8DAY$camsiteyear)[i],]
  curvefit <- try(nls(meanNDVI ~ SSlogis(DOY, Asym, xmid, scal), siteyear))
  if(inherits(curvefit, "try-error")){
    next
  } else {
    year <- siteyear$year[1]
    predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
    predsiteyear$DOY <- yday(predsiteyear$date)
    predsiteyear$preds <- predict(curvefit, newdata=predsiteyear)
    predsiteyear <- left_join(predsiteyear, day_occasion_df, by="date")
    predsiteocc <- predsiteyear%>%group_by(occ)%>%summarise(meanEVI=mean(preds))%>%mutate("cam_site_id"=siteyear$camsiteid[1], "year"=siteyear$year[1])
    if(!exists("EVILS8siteyearocc")){
      EVILS8siteyearocc <- predsiteocc
    }
    else{
      EVILS8siteyearocc <- rbind(EVILS8siteyearocc, predsiteocc)
    }
  }
}

EVILS8BearCPUE <- left_join(BEARCPUE3, EVILS8siteyearocc, by=c("cam_site_id", "year", "occ"))
table(is.na(EVILS8BearCPUE$meanEVI)) #26847
NAEVI8Day <- EVILS8BearCPUE[is.na(EVILS8BearCPUE$meanEVI),]
hist(table(NAEVI8Day$cam_site_id))




#Modis 16 day EVI, not available for 2023 or is it 463m pixel size 
#https://developers.google.com/earth-engine/datasets/catalog/MODIS_MCD43A4_006_EVI#bands
EVI_MODIS_16DAY <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/MODIS16DayGlobal250_2019-2024.csv") #named for start time
EVI_MODIS_16DAY$year <- year(as.Date(EVI_MODIS_16DAY$date))
EVI_MODIS_16DAY$DOY <- yday(EVI_MODIS_16DAY$date)
EVI_MODIS_16DAY$camsiteyear <- paste0(EVI_MODIS_16DAY$camsiteid, EVI_MODIS_16DAY$year)
testsite <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == "cam_00012019",]
ggplot(testsite, aes(x=DOY, y=meanNDVI)) + geom_point() + geom_smooth(method = "lm", se = FALSE)
EVI_MODIS_16DAYpercamyear <- EVI_MODIS_16DAY%>%group_by(camsiteid, year)%>%summarise(n=n())
for(i in 1:length(unique(EVI_MODIS_16DAY$camsiteyear))){
  siteyear <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == unique(EVI_MODIS_16DAY$camsiteyear)[i],]
  curvefit <- try(nls(meanNDVI ~ SSlogis(DOY, Asym, xmid, scal), siteyear))
  if(inherits(curvefit, "try-error")){
    predsiteocc <- data.frame(occ=1:11, meanEVI=rep(NA,11), cam_site_id=siteyear$camsiteid[1], year=siteyear$year[1])
    if(!exists("EVIModissiteyearocc")){
      EVIModissiteyearocc <- predsiteocc
    }
    else{
      EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
    }
    next
  } else {
    year <- siteyear$year[1]
    predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
    predsiteyear$DOY <- yday(predsiteyear$date)
    predsiteyear$preds <- predict(curvefit, newdata=predsiteyear)
    predsiteyear <- left_join(predsiteyear, day_occasion_df, by="date")
    predsiteocc <- predsiteyear%>%group_by(occ)%>%summarise(meanEVI=mean(preds))%>%mutate("cam_site_id"=siteyear$camsiteid[1], "year"=siteyear$year[1])
    
    if(!exists("EVIModissiteyearocc")){
      EVIModissiteyearocc <- predsiteocc
    }
    else{
      EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
    }
  }
}

EVIModisBearCPUE <- left_join(BEARCPUE3, EVIModissiteyearocc, by=c("cam_site_id", "year", "occ"))
table(is.na(EVIModisBearCPUE$meanEVI)) #2250 false
NAEVI16Day <- EVIModisBearCPUE[is.na(EVIModisBearCPUE$meanEVI),]
hist(table(NAEVI16Day$cam_site_id))
NAEVI16DaybyCam <- as.data.frame(table(NAEVI16Day$cam_site_id))
NAEVI16DaybyCam <- right_join(camsites, NAEVI16DaybyCam, by=join_by("cam_site_id" == "Var1"))%>%st_transform(., 4326)
NAEVI16DaybyCamBig <- NAEVI16DaybyCam[NAEVI16DaybyCam$Freq > 12,]
NAEVI <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear %in% unique(paste0(NAEVI16Day$cam_site_id, NAEVI16Day$year)),]
NAEVI%>%group_by(camsiteyear)





leaflet() %>% 
  addProviderTiles('Esri.WorldImagery') %>%
  addCircleMarkers(data=NAEVI16DaybyCamBig, fillColor = "white", fillOpacity = 1, stroke=F, radius=3, label=NAEVI16DaybyCamBig$Freq)

EVIModisBearCPUE%>%group_by(cam_site_id, year)%>%summarise()

which(!(paste0(BEARCPUE3$cam_site_id, BEARCPUE3$year) %in% EVI_MODIS_16DAY$camsiteyear))



#Modis Daily EVI, not available for 2023
#https://developers.google.com/earth-engine/datasets/catalog/MODIS_MOD09GA_006_EVI
EVI_MODIS_Daily <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/EVI_MODISDaily.csv") #named for start time
EVI_MODIS_Daily$year <- year(as.Date(EVI_MODIS_Daily$date))
EVI_MODIS_Daily$DOY <- yday(EVI_MODIS_Daily$date)
EVI_MODIS_Daily$camsiteyear <- paste0(EVI_MODIS_Daily$camsiteid, EVI_MODIS_Daily$year)
EVI_MODIS_Daily$date <- as.Date(EVI_MODIS_Daily$date)
EVI_MODIS_Dailypercamyear <- EVI_MODIS_Daily%>%group_by(camsiteid, year)%>%summarise(n=n())
EVI_MODIS_Daily <- left_join(EVI_MODIS_Daily, day_occasion_df, by=c("date", "year"))
colnames(EVI_MODIS_Daily)[2] <- "cam_site_id"
EVI_MODIS_Daily2 <- EVI_MODIS_Daily%>%group_by(cam_site_id, year, occ)%>%summarise(meanEVI=mean(meanNDVI))
EVIDailyCPUE3 <- left_join(BEARCPUE3, EVI_MODIS_Daily2, by=c("cam_site_id", "year", "occ"))
table(is.na(EVIDailyCPUE3$meanEVI)) #10619 NAs
NAEVIDaily <- EVIDailyCPUE3[is.na(EVIDailyCPUE3$meanEVI),]
hist(table(NAEVIDaily$cam_site_id))
NAEVItable(NAEVIDaily$cam_site_id)

## fit without uncertainty estimation
data(bartlett2009.filtered)
bartlett2009.filtered <- bartlett2009.filtered[-1]
index(bartlett2009.filtered) <- lubridate::yday(index(bartlett2009.filtered))
fitted.beck <- BeckFit(bartlett2009.filtered)
days <- index(bartlett2009.filtered)
plot(days, bartlett2009.filtered)
lines(fitted.beck$fit$predicted, col='red')



agepal <- colorNumeric(palette = "Spectral",domain=values(NDVI) , na.color="transparent")
lcpal(1:11)
previewColors(colorFactor("Spectral", domain = NULL), values)
brewer.pal(17,"Spectral")
camsites3857 <- st_transform(camsites, 3857)

leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addCircleMarkers(data=camsites3857, fillColor = "black", fillOpacity = 1,   stroke=F, radius=5, popup=camsites3857$cam_site_id) %>%
  addRasterImage(NDVI,colors = agepal, opacity = 0.8, )# %>%
addLegend(pal=agepal, position = "bottomright", values=values(NDVI))



#EVI  MODIS 16day 250m 04-01 to 07-15
#https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1
EVI_MODIS_16DAY <- read.csv("C:/Users/wildeefb/Documents/GeoSpatial/BearNDVI/EVI_MODIS16Day2019-2024.csv") #named for start time
EVI_MODIS_16DAY$year <- year(as.Date(EVI_MODIS_16DAY$date))
EVI_MODIS_16DAY$DOY <- yday(EVI_MODIS_16DAY$date)
EVI_MODIS_16DAY$camsiteyear <- paste0(EVI_MODIS_16DAY$camsiteid, EVI_MODIS_16DAY$year)
EVI_MODIS_16DAYpercamyear <- EVI_MODIS_16DAY%>%group_by(camsiteid, year)%>%summarise(n=n())
testsite <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == "cam_00012019",]
ggplot(testsite, aes(x=DOY, y=meanNDVI)) + geom_point() + geom_smooth(method = "lm", se = FALSE)
EVI_MODIS_16DAYpercamyear <- EVI_MODIS_16DAY%>%group_by(camsiteid, year)%>%summarise(n=n())
for(i in 1:length(unique(EVI_MODIS_16DAY$camsiteyear))){
  siteyear <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteyear == unique(EVI_MODIS_16DAY$camsiteyear)[i],]
  curvefit <- try(nls(meanEVI ~ SSlogis(DOY, Asym, xmid, scal), siteyear))
  if(inherits(curvefit, "try-error")){
    meanEVI <- 
      predsiteocc <- data.frame(occ=1:11, meanEVI=rep(NA,11), cam_site_id=siteyear$camsiteid[1], year=siteyear$year[1])
    if(!exists("EVIModissiteyearocc")){
      EVIModissiteyearocc <- predsiteocc
    }
    else{
      EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
    }
    next
  } else {
    year <- siteyear$year[1]
    predsiteyear= data.frame("date"=(seq.Date(from=as.Date(paste(year, "04-15", sep="-")),to = as.Date(paste(year, "06-30", sep="-")), by = 1)))
    predsiteyear$DOY <- yday(predsiteyear$date)
    predsiteyear$preds <- predict(curvefit, newdata=predsiteyear)
    predsiteyear <- left_join(predsiteyear, day_occasion_df, by="date")
    predsiteocc <- predsiteyear%>%group_by(occ)%>%summarise(meanEVI=mean(preds))%>%mutate("cam_site_id"=siteyear$camsiteid[1], "year"=year)
    
    if(!exists("EVIModissiteyearocc")){
      EVIModissiteyearocc <- predsiteocc
    }
    else{
      EVIModissiteyearocc <- rbind(EVIModissiteyearocc, predsiteocc)
    }
  }
}

EVIModissiteyearocc$camsiteyear <- paste(EVIModissiteyearocc$cam_site_id, EVIModissiteyearocc$year)
EVIModisBearCPUE <- left_join(BEARCPUE3, EVIModissiteyearocc, by=c("cam_site_id", "year", "occ"))
table(is.na(EVIModisBearCPUE$meanEVI)) #2250 false
NAEVI <- EVIModisBearCPUE[is.na(EVIModisBearCPUE$meanEVI),]
NAEVIbycam <- NAEVI%>%group_by(cam_site_id)%>%summarise(n=n())
NAEVIbycam3857 <- st_transform(NAEVIbycam, crs=4326)
cam_2231 <- EVI_MODIS_16DAY[EVI_MODIS_16DAY$camsiteid == "cam_2228",]
ggplot(cam_2231, aes(x=DOY, y=meanEVI)) + facet_wrap(~year) + geom_point() + geom_smooth()


#################################  NDVI-Modistools  ###############################
products <- mt_products()
bands <- mt_bands(product = "MOD13Q1")
dates <- mt_dates(product = "MOD13Q1", lat = 44.7, lon = -92.5)
subset <- mt_subset(product = "MOD13Q1",
                    lat = 40,
                    lon = -110,
                    band = "250m_16_days_NDVI",
                    start = "2019-01-01",
                    end = "2020-01-01",
                    km_lr = 1,
                    km_ab = 1,
                    site_name = "testsite",
                    internal = TRUE,
                    progress = FALSE)

NDVI <- apply(camlocdf[1303:3022,],MARGIN=1, function (x) mt_subset(product = "MOD13Q1",
                                                                    lat = x[2],
                                                                    lon = x[3],
                                                                    band = "250m_16_days_NDVI",
                                                                    start = "2019-03-01",
                                                                    end = "2023-06-01",
                                                                    km_lr = 1,
                                                                    km_ab = 1,
                                                                    site_name = x[1],
                                                                    out_dir="C:/Users/wildeefb/Documents/GeoSpatial/NDVI2",
                                                                    internal = FALSE,
                                                                    progress = TRUE))

# test batch download
camlocdf <- data.frame("site_name"=camsites$cam_site_id, "lat"=st_coordinates(st_transform(camsites, crs=4326))[,2],
                       "lon"=st_coordinates(st_transform(camsites, crs=4326))[,1])
NDVI <- mt_batch_subset(df = camlocdf,
                        product = "MOD13Q1",
                        band = "250m_16_days_NDVI",
                        internal = FALSE,
                        out_dir= "C:/Users/wildeefb/Documents/GeoSpatial/NDVI",
                        start = "2019-01-01",
                        end = "2019-01-03",
                        km_lr = 1,
                        km_ab = 1)

apply(camlocdf[1:10,],MARGIN=1, function (x) x[1])

# convert the coordinates
lat_lon <- sin_to_ll(arcachon_lc$xllcorner, arcachon_lc$yllcorner)

# bind with the original dataframe
subset <- cbind(arcachon_lc, lat_lon)

# convert to bounding box
bb <- apply(arcachon_lc, 1, function(x){
  mt_bbox(xllcorner = x['xllcorner'],
          yllcorner = x['yllcorner'],
          cellsize = x['cellsize'],
          nrows = x['nrows'],
          ncols = x['ncols'])
})

LC_r <- mt_to_terra(df = arcachon_lc, reproject = TRUE)


# Make a new data.frame that will contain binarized VIQ values.
VIQbin = VIQ

# Solve for VI Quality
# Source: https://gis.stackexchange.com/questions/144441/how-can-i-parse-modis-mod13q1-quality-layers-in-r
first_k_bits <- function(int, k=16, reverse=T) {
  integer_vector <- as.integer(intToBits(int))[1:k]
  if(reverse) integer_vector <- rev(integer_vector)
  return(paste(as.character(integer_vector), collapse=""))
}
# Binarize the VIQ values:
VIQbin_list = lapply(VIQ$value,
                     FUN = first_k_bits)
VIQbin_vector = unlist(VIQbin_list)
VIQbin$value = as.character(VIQbin_vector)

