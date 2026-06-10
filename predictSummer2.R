library(CropScapeR)
library(sf)
library(terra)
library(readxl)
library(rlandfire)
library(ggplot2)
library(dplyr)
library(sswids)
library(MCMCvis)
library(tidyterra)

samplesSummer <- readRDS("./RNsamplesFullModelBearRangeSummerGOF.rds")
bearzones <- get_spatial_data("bear_zones")%>%rename(zone=bear_mgmt_zone_id)
b_Fixed <- MCMCsummary(samplesSummer[1:3], 
                       params = c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Wetland"),
                       ISB = TRUE,
                       round=2)

#get scale attributes from model site covariates
selectedscales <- c("Corn2500", "Forest_500", "Developed_2500", "Dist_1000", "Wetland_500")
scaleattr <- lapply(selectedscales, function (x) {scaleattr <- c(attr(sitecovs%>%pull(x), "scaled:center"), attr(sitecovs%>%pull(x), "scaled:scale"))})
names(scaleattr) <- selectedscales

#Corn
CropRasts <- list.files("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/", full.names = TRUE)
CropRasts2 <- lapply(CropRasts, function (x) rast(x))
data("linkdata")


CropRast2019 <- CropRasts2[[1]]

cormat <- cbind(c(1, 12, 13, 225, 226, 228, 237, 241), rep(1, 8))
CornLC <- lapply(CropRasts2, function (x) classify(x, rcl=cormat, others=0))
plot(CornLC[[1]])
CornLC <- lapply(CornLC, function (x) crop(x, st_transform(bearrange2, crs=crs(x))))
CornLC <- lapply(CornLC, function (x) mask(x, st_transform(bearrange2, crs=crs(x))))
plot(CornLC[[1]])
temp <- CropRast2019
res(temp) <- 5000
CornLC5000 <- lapply(CornLC, function (x) resample(x, temp, method=mean))
CornLC5000 <- lapply(CornLC5000, function (x) {values(x) <- (values(x)*100- scaleattr[["Corn2500"]][1])/scaleattr[["Corn2500"]][2];x})#scaled sitecov variables on percentage level
plot(CornLC5000[[1]])
names(CornLC5000) <- paste0("Corn", 2019:2024)



Wiscland3 <- get_spatial_data(layer_name = 'wiscland2', level = 3)
Wiscland3 <- project(Wiscland3, crs(CropRast2019))
WisclandGuide <- read_xlsx("C:/Users/wildeefb/Documents/GeoSpatial/wiscland2/user_guide/Wiscland2 Color Scheme.xlsx",
                           range = "B2:C70", .name_repair = make.names)

forestLCs <- c("AspenPaperBirch", "RedMaple", "Oak", "CentralHardwoods", "NorthernHardwoods",
               "MixedDeciduousConiferousForest")
wetlandLCs <- c("AspenForestedWetland", "BottomlandHardwoods", "SwampHardwoods", "OtherEmergentWetMeadow",
                "MixedDeciduousConiferousForestedWetland", "ConiferousForestedWetland", "BroadleavedDeciduousScrubShrub",
                "BroadleavedEvergreenScrubShrub", "NeedleleavedScrubShrub")
developedLCs <- c("DevelopedHighIntensity","DevelopedLowIntensity")

forestLCs <- WisclandGuide$dn.label[which(gsub(pattern = "\\W", replacement = "", x = WisclandGuide$label) %in% forestLCs)]
wetlandLCs <-WisclandGuide$dn.label[which(gsub(pattern = "\\W", replacement = "", x = WisclandGuide$label) %in% wetlandLCs)]
developedLCs <-WisclandGuide$dn.label[which(gsub(pattern = "\\W", replacement = "", x = WisclandGuide$label) %in% developedLCs)]

forestmat <- cbind(forestLCs, rep(1, length(forestLCs)))
wetlandmat <- cbind(wetlandLCs, rep(1, length(wetlandLCs)))
developedmat <- cbind(developedLCs, rep(1, length(developedLCs)))

WisclandRC <- lapply(list(forestmat, wetlandmat, developedmat), function (x) classify(Wiscland3, x, others=0))
WisclandRC <- lapply(WisclandRC, function (x) crop(x, st_transform(bearrange2, crs=crs(x))))
WisclandRC <- lapply(WisclandRC, function (x) mask(x, st_transform(bearrange2, crs=crs(x))))
names(WisclandRC) <- c("Forest", "Wetland", "Developed")

#Forest-0.5, Developed-2.5, Wetland-0.5
temp <- Wiscland3
res(temp) <- 1000
ForestLC1000 <- resample(WisclandRC[[1]], temp, method=mean)
values(ForestLC1000) <- (values(ForestLC1000)*100-scaleattr[["Forest_500"]][1])/scaleattr[["Forest_500"]][2]
WetlandLC1000 <- resample(WisclandRC[[2]], temp, method=mean)
values(WetlandLC1000) <- (values(WetlandLC1000)*100-scaleattr[["Wetland_500"]][1])/scaleattr[["Wetland_500"]][2]
temp <- Wiscland3
res(temp) <- 5000
DevelopedLC5000 <- resample(WisclandRC[[3]], temp, method=mean)
values(DevelopedLC5000) <- (values(DevelopedLC5000)*100-scaleattr[["Developed_2500"]][1])/scaleattr[["Developed_2500"]][2]
plot(DevelopedLC5000)


path <- "C:/Users/wildeefb/Documents/GeoSpatial"
lf_dir <- file.path(path, "lf")
hdist <- terra::rast(list.files(lf_dir, pattern = ".tif$", 
                                full.names = TRUE, 
                                recursive = TRUE))
dbf <- list.files(lf_dir, pattern = ".dbf$",
                  full.names = TRUE,
                  recursive = TRUE)
dbf_tbl  <- foreign::read.dbf(dbf)
levels(hdist) <- dbf_tbl[,c(1,5)]
hdist <- addCats(hdist, value=dbf_tbl[,c(3:4)])
cats(hdist)
activeCat(hdist)
levels(hdist)
plot(hdist)
activeCat(hdist) <- 1


hdist2 <- ifel(hdist>0, 1, hdist)
hdist2 <- project(hdist2, crs(CropRast2019))
hdist2 <- crop(hdist2, st_transform(bearrange2, crs=crs(hdist2)))
hdist2 <- mask(hdist2, st_transform(bearrange2, crs=crs(hdist2)))
temp <- hdist2
res(temp) <- 2000
hdistLC2000 <- resample(hdist2, temp, method=mean)
values(hdistLC2000) <- (values(hdistLC2000)*100-scaleattr[["Dist_1000"]][1])/scaleattr[["Dist_1000"]][2]
plot(hdistLC2000)

#Get all the covariates on the bear home range scale with terra::resample, method mean for finer than HR resolution and method bilinear for coarser than HR resolution
cellside <- sqrt(5.24)*(1.60934)*1000 #5.24 sq mi home range estimate from Erin
temp <- CornLC5000[[1]]
res(temp) <- cellside
CornLCHR <- lapply(CornLC5000, function (x) resample(x, temp, method="bilinear"))
ForestLCHR <- resample(ForestLC1000, temp, method="mean")
WetlandsLCHR <- resample(WetlandLC1000, temp, method="mean")
DevelopLCHR <- resample(DevelopedLC5000, temp, method="bilinear")
hdistLCHR <- resample(hdistLC2000, temp, method="mean")

#Get spatial spline at HR scale
#do this at HR scale
hdistLCHR.3071 <- project(hdistLCHR, "EPSG:3071")
cellcoords2 <- as.data.frame(xyFromCell(hdistLCHR.3071, cell=which(terra::values(!is.na(hdistLCHR.3071)))))
#get Zs for spatial spline
cellcoords2 <- cellcoords2 %>%
  mutate(X.scale = (x-mean_x)/sd_x,    #need mean_xy and sd_xy from original data frame of camera points
         Y.scale = (y-mean_y)/sd_y)
bearrangeknots2 <- read.csv("./bearrangeknots.csv")
bearrangeknots2 <- bearrangeknots2 %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y) %>%
  dplyr::select(X.scale,Y.scale)
sp.cov.dist_predall = fields::rdist(x1=cbind(cellcoords2$X.scale,cellcoords2$Y.scale),x2=bearrangeknots2)
sp.Z_K.predall = sp.cov.dist_predall^2*log(sp.cov.dist_predall) # basis
sp.Z.predall <- t(solve(sp.sqrt.omega_all,t(sp.Z_K.predall)))
sp.Z.predall <- (sp.Z.predall - sp.meanZ)/sp.sdZ
spatspline.bs <- MCMCsummary(samplesSummer[1:3],  
                             params = c('spat.spline.b'),
                             ISB = TRUE,
                             round=2)
#this gets the mean spatial spline abundance but can't use for CIs because of use of posterior beta means
inprodspline2 <- sp.Z.predall%*%spatspline.bs$mean
inprodspline3 <- terra::values(hdistLCHR.3071)
inprodspline3[!is.na(inprodspline3)] <- inprodspline2
Zrast <- rast(hdistLCHR.3071, nlyr=1)
terra::values(Zrast) <- inprodspline3
Zrast <- crop(Zrast, st_transform(bearrange2, crs=crs(Zrast)))
Zrast <- mask(Zrast, st_transform(bearrange2, crs=crs(Zrast)))
Zrast <- project(Zrast, crs(ForestLCHR))
Zrast2 <- resample(Zrast, ForestLCHR)
plot(Zrast2[[1]])

#for raw Z values
sp.Zrast <- rast(hdistLCHR.3071, nlyr=50)
#make template raster for 50 column Z matrix
sp.Zrast_values <- matrix(data=rep(terra::values(hdistLCHR.3071), 50), nrow = length(terra::values(hdistLCHR.3071)),ncol=50)
sp.Zrast_values[!is.na(sp.Zrast_values)] <- sp.Z.predall
terra::values(sp.Zrast) <- sp.Zrast_values
sp.Zrast <- crop(sp.Zrast, st_transform(bearrange2, crs=crs(sp.Zrast)))
sp.Zrast <- mask(sp.Zrast, st_transform(bearrange2, crs=crs(sp.Zrast)))
sp.Zrast <- project(sp.Zrast, crs(ForestLCHR))
sp.Zrast <- resample(sp.Zrast, ForestLCHR)
plot(sp.Zrast[[10]])
#



bearzones <- st_transform(bearzones, crs(ForestLCHR))
bearzonesrast <- rasterize(vect(bearzones), temp, field="zone")
bearzonesrast <- crop(bearzonesrast, st_transform(bearrange2, crs=crs(hdist2)))
bearzonesrast <- mask(bearzonesrast, st_transform(bearrange2, crs=crs(hdist2)))
bearzonesrast <- resample(bearzonesrast, ForestLCHR)
levels(bearzonesrast)
plot(bearzones[[1]])
# Combine into a list
raster_list <- c(CornLCHR, ForestLCHR, WetlandsLCHR, DevelopLCHR, hdistLCHR, Zrast2, sp.Zrast, bearzonesrast)
names(raster_list) <- c(names(CornLC5000), names(WisclandRC), "Disturbance", "SpatialSplineMean", "Zmatrix", "Zone")
# Create the SpatRaster stack
stacked_raster <- rast(raster_list)
covariate_matrix <- as.matrix(stacked_raster)
covariate_matrix <- covariate_matrix[complete.cases(covariate_matrix),]
covariate_matrix <- left_join(as.data.frame(covariate_matrix), as.data.frame(levels(bearzonesrast)), by=join_by("Zone"=="ID"))




b_ZoneEff <- MCMCsummary(samplesSummer[1:3], 
                         params = c("b_ZoneYr", "b_Zone"),
                         ISB = TRUE,
                         round=9)
b_ZoneEff2 <- data.frame("beta.Zone"=c(b_ZoneEff$mean[7:12]), 
                         "beta.ZoneYr"= c(b_ZoneEff$mean[1:6]) , 
                         "zone"=LETTERS[1:6])
covariate_matrix <- left_join(covariate_matrix, b_ZoneEff2)

yearX <- unique(sitecovs$year)
lambda.predicted.grid <- lapply(1:length(yearX), function (i) 
  lambda <- exp(covariate_matrix$beta.ZoneYr*yearX[i] + covariate_matrix$beta.Zone  + 
                  b_Fixed$mean[1]*covariate_matrix$Developed + 
                  b_Fixed$mean[2]*covariate_matrix$Forest + b_Fixed$mean[3]*covariate_matrix[,i] +
                  b_Fixed$mean[4]*covariate_matrix$Disturbance + b_Fixed$mean[5]*covariate_matrix$Wetland +
                  covariate_matrix$SpatialSplineMean
                  ))
names(lambda.predicted.grid) <- 2019:2024
completecells <- which(complete.cases(as.matrix(stacked_raster)))
lambda.rast <- lapply(lambda.predicted.grid, function (x) {
  rasttemp <- rast(ForestLCHR, nlyr=1)
  terra::values(rasttemp)[completecells] <- x
  return(rasttemp)
  })
log.lambda.rast <- lapply(lambda.rast, function (x) log(x))

lambda.predicted.grid2 <- do.call(cbind, lambda.predicted.grid)
predict.lambda4 <- cbind(covariate_matrix, lambda.predicted.grid2)
colnames(predict.lambda4)[66:71] <- paste0("TotalLambda", 2019:2024)
saveRDS(predict.lambda4, "./predict.lambda4Zones.rds")

#####################################################################################
#                             abundance by zone                                     #
#####################################################################################
popbyzone <- predict.lambda4%>%group_by(zone)%>%summarise(across(matches("TotalLambda"), ~sum(.x)))
total <- data.frame("zone"="Total", t(colSums(popbyzone[,2:7])))
popbyzone <- rbind(popbyzone, total)
popbyzone2 <- pivot_longer(popbyzone, cols = -c(zone), names_to = "Year", values_to = "lambda")
popbyzone2$Year <- as.numeric(gsub(x = popbyzone2$Year, pattern = "TotalLambda", replacement = ""))
ggplot(filter(popbyzone2, zone != "Total"), aes(x=Year, y=lambda, colour = zone)) + geom_point(size=3) + geom_line(lwd=1)

#variation
#combine MCMC chains into one
allchains <- MCMCchains(samplesSummer[1:3], params =c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Wetland", "b_Zone", "b_ZoneYr", "spat.spline.b"), ISB = TRUE)
#loop through estimated parameter at each iteration

allchainssample <- allchains[sample(nrow(allchains), 10000), ]

predict.df <- predict.lambda4[complete.cases(predict.lambda4),]
CIsforZones <- function(predict.df, nZones, year, yearX) {
  for (z in 1:nZones){
    ZonePredict <- predict.df[!is.na(predict.df$zone) & predict.df$zone == LETTERS[z],]
    print(cat(LETTERS[z], nrow(ZonePredict)))
    #zfac <- ifelse(z==6, 5, z) #need to change this with E,F merge/un-merge
    lambda.CIs <- array(dim = c(nrow(ZonePredict), 10000))
    yr <- year-2018
    for(j in 1:10000){
      lambda.CIs[,j] <- exp(allchainssample[j, paste0("b_Zone[", z,"]")] +
                              allchainssample[j, paste0("b_ZoneYr[", z,"]")]*yearX[yr] +
                              allchainssample[j,"b_Dev"]*ZonePredict$Developed + 
                              allchainssample[j,"b_Forest"]*ZonePredict$Forest + 
                              allchainssample[j,"b_Corn"]*ZonePredict[,paste0("Corn", year)] +
                              allchainssample[j,"b_Dist"]*ZonePredict$Disturbance +
                              allchainssample[j,"b_Wetland"]*ZonePredict$Wetland +
                              (as.matrix(ZonePredict[,grep(x = colnames(ZonePredict), pattern = "Zmatrix")])%*%
                             allchainssample[j, grep(x = colnames(allchainssample), pattern = "spat.spline.b")])
                            ) #this latter column specification needs to change depending on # of zones
    }
    CLcells <- apply(X = lambda.CIs, MARGIN = 1, FUN = quantile(0.975)-quantile(0.025))
    zone.lambda.distr <- colSums(lambda.CIs)
    CL <- c(mean(zone.lambda.distr),quantile(zone.lambda.distr, c(0.025, 0.5, 0.975)))
    if(z==1){
      CL2 <- data.frame("Zone"=LETTERS[z], "Year"=yr, "mean"=CL[1],"CL2.5"=CL[2], "CL0.5"=CL[3],"CL97.5"=CL[4])
    }else{
      CL2 <- rbind(CL2, data.frame("Zone"=LETTERS[z], "Year"=yr, "mean"=CL[1], "CL2.5"=CL[2], "CL0.5"=CL[3],"CL97.5"=CL[4]))
    }
  }
  return(CL2)
}
CIs2019 <- CIsforZones(predict.df = predict.lambda4, nZones = 6, year = 2019, yearX)
CIs2020 <- CIsforZones(predict.df = predict.lambda4, nZones = 6, year = 2020, yearX)
CIsZone <- lapply(2019:2024, function (i) CIsforZones(predict.df = predict.lambda4, nZones = 6, year = i, yearX))
CIsZone2 <- rbindlist(CIsZone)
CIsZone2$Year <- CIsZone2$Year+2018



popbyzone3 <- left_join(popbyzone2, CIsZone2, by=join_by("zone"=="Zone", "Year"))
cbPalette <- c("#F0E442","#009E73", "#D55E00","#E69F00", "#0072B2","#56B4E9")
ggplot(filter(popbyzone3, zone != "Total"), aes(x=Year, y=lambda, colour = zone)) + 
  geom_point(size=2) +
  #geom_pointrange(aes(ymin = CL2.5, ymax = CL97.5,color= Zone), size=1, lineend = "square", lwd=1) + 
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5,color= zone, fill=zone), alpha=0.3) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(y=expression("Predicted Relative Abundance ( "*lambda* ")"), color="Zone", fill="Zone") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
  )

#calculate interval
lambda.CIs2019.2 <- bind_cols(lambda.CIs2019, predict.lambda.sf2$bear_mgmt_zone_id)
colnames(lambda.CIs2019.2[10001]) <- "Zone"
zonepop.posterior <- lambda.CIs2019.2%>%group_by(Zone)
totalpop.posterior <- lambda.CIs2019
CL2019A <- apply(lambda.CIs, 1, function(x){quantile(x, prob = c(0.025, 0.975))})




#####################################################################################
#                   spatial prediction across years                                 #
#####################################################################################
lambda.rast.stack <- rast(lambda.rast)
Wisconsin.RastCRS <- st_transform(Wisconsin2, crs=crs(lambda.rast.stack))
ggplot() +
  geom_spatraster(data = lambda.rast.stack) +
  geom_sf(data=Wisconsin.RastCRS, fill=NA) +
  scale_fill_viridis_c(na.value = "transparent", option="C") +
  facet_wrap(~lyr) + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),        # Removes x and y axis tick numbers/text
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(fill=expression("Predicted\n Relative\n Abundance ("*lambda* ")"),
  )



##################################### scrap ###################################
library(leaflet)
cellcoordsSF <- st_as_sf(cellcoords2, coords=c("x", "y"), crs=3071)
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  #addCircleMarkers(data=st_transform(notbadinits, 4326), popup=notbadinits$cam_site_id, fillColor = "yellow", fillOpacity = 1,   stroke=F, radius=5)%>%
  addCircleMarkers(data=st_transform(cellcoordsSF, 4326), fillColor = "black", fillOpacity = 1,   stroke=F, radius=5)%>%
  addCircleMarkers(data=st_transform(predict.lambda.sf, 4326), fillColor = "red", fillOpacity = 1,   stroke=F, radius=5)%>%
  addPolygons(data=st_transform(Wisconsin, crs=4326))

z <- cbind(inprodspline2, cellcoords2)
ggplot(z, aes(x=x, y=y, colour = inprodspline2)) + geom_point() + scale_colour_viridis(option="C")


whichagg <- data.frame(summarise5=as.numeric(terra::values(hdistLC5000)), aggfrom2=as.numeric(terra::values(hdistLCagg)))%>%drop_na()
ggplot(whichagg, aes(x=summarise5, y=aggfrom2)) + geom_point()


for (z in 1:6){cat(paste0("b_Zone[", z,"]"), LETTERS[z])}



lambda.predicted.grid <- lapply(1:length(yearX), function (i) 
  lambda <- exp(b_Fixed$mean[1]*covariate_matrix$Developed + 
                  b_Fixed$mean[2]*covariate_matrix$Forest + b_Fixed$mean[3]*covariate_matrix[,i] +
                  b_Fixed$mean[4]*covariate_matrix$Disturbance + b_Fixed$mean[5]*covariate_matrix$Wetland))
names(lambda.predicted.grid) <- 2019:2024
lambda.predicted.grid2 <- do.call(cbind, lambda.predicted.grid)
predict.lambda4 <- cbind(covariate_matrix, lambda.predicted.grid2)
colnames(predict.lambda4)[66:71] <- paste0("TotalLambda", 2019:2024)
popbyzone <- predict.lambda4%>%group_by(zone)%>%summarise(across(matches("TotalLambda"), ~sum(.x)))
total <- data.frame("zone"="Total", t(colSums(popbyzone[,2:7])))
popbyzone <- rbind(popbyzone, total)
popbyzone2 <- pivot_longer(popbyzone, cols = -c(zone), names_to = "Year", values_to = "lambda")
popbyzone2$Year <- as.numeric(gsub(x = popbyzone2$Year, pattern = "TotalLambda", replacement = ""))


CIsZone <- function(predict.df, nZones, year, yearX) {
  for (z in 1:nZones){
    ZonePredict <- predict.df[!is.na(predict.df$zone) & predict.df$zone == LETTERS[z],]
    print(cat(LETTERS[z], nrow(ZonePredict)))
    #zfac <- ifelse(z==6, 5, z) #need to change this with E,F merge/un-merge
    lambda.CIs <- array(dim = c(nrow(ZonePredict), 10000))
    yr <- year-2018
    for(j in 1:10000){
      lambda.CIs[,j] <- exp(allchainssample[j,"b_Dev"]*ZonePredict$Developed + 
                              allchainssample[j,"b_Forest"]*ZonePredict$Forest + 
                              allchainssample[j,"b_Corn"]*ZonePredict[,paste0("Corn", year)] +
                              allchainssample[j,"b_Dist"]*ZonePredict$Disturbance) #this latter column specification needs to change depending on # of zones
    }
    zone.lambda.distr <- colSums(lambda.CIs)
    CL <- c(mean(zone.lambda.distr),quantile(zone.lambda.distr, c(0.025, 0.5, 0.975)))
    if(z==1){
      CL2 <- data.frame("Zone"=LETTERS[z], "Year"=yr, "mean"=CL[1],"CL2.5"=CL[2], "CL0.5"=CL[3],"CL97.5"=CL[4])
    }else{
      CL2 <- rbind(CL2, data.frame("Zone"=LETTERS[z], "Year"=yr, "mean"=CL[1], "CL2.5"=CL[2], "CL0.5"=CL[3],"CL97.5"=CL[4]))
    }
  }
  return(CL2)
}
CIsZone <- lapply(2019:2024, function (i) CIsZone(predict.df = predict.lambda4, nZones = 6, year = i, yearX))
CIsZone2 <- rbindlist(CIsZone)
CIsZone2$Year <- CIsZone2$Year+2018

lapply(2019:2024, function(i) print(i-2018))

popbyzone3 <- left_join(popbyzone2, CIsZone2, by=join_by("zone"=="Zone", "Year"))
cbPalette <- c("#F0E442","#009E73", "#D55E00","#E69F00", "#0072B2","#56B4E9")
ggplot(filter(popbyzone3, zone != "Total"), aes(x=Year, y=lambda, colour = zone)) + 
  geom_point(size=2) +
  #geom_pointrange(aes(ymin = CL2.5, ymax = CL97.5,color= Zone), size=1, lineend = "square", lwd=1) + 
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5,color= zone, fill=zone), alpha=0.3) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(y="Mean Relative Abundance")


lambda.predicted.fixed <- lapply(1:length(yearX), function (i) 
  lambdas <- data.frame(ZoneYr=covariate_matrix$beta.ZoneYr*yearX[i] , Zone=covariate_matrix$beta.Zone  , 
                  Developed=b_Fixed$mean[1]*covariate_matrix$Developed ,
                  Forest=b_Fixed$mean[2]*covariate_matrix$Forest , Corn=b_Fixed$mean[3]*covariate_matrix[,i] ,
                  Disturbance=b_Fixed$mean[4]*covariate_matrix$Disturbance , Wetland=b_Fixed$mean[5]*covariate_matrix$Wetland
  ))
lambda.predicted.fixed <-lapply(lambda.predicted.fixed, function (x) cbind (x, covariate_matrix$SpatialSplineMean))
