library(landscapemetrics)
library(stars)
library(sf)
library(dplyr)
library(terra)
library(readxl)
library(rlandfire)
library(CropScapeR)
library(purrr)
library(tidyr)
library(leaflet)
library(ggplot2)




####################################predict EVI spline############################################
#predicting EVI spline
EVI.pred <- seq(from=range(EVI2$meanEVI)[1], to=range(EVI2$meanEVI)[2], length.out=100)
pred.data <- expand.grid(camversion=1, daysactive=0.2211643, EVI=EVI.pred)
# prediction dataset for spatial smoothing
pred.EVI.Z_K <- (abs(outer(as.numeric(pred.data$EVI),EVI.knots,"-")))^3

pred.EVI.Z <- t(solve(EVI.sqrt.OMEGA_all,t(pred.EVI.Z_K)))

# standardize for better performance
pred.EVI.Z <- (pred.EVI.Z-meanZ)/sdZ
pred.data2 <- cbind(pred.data, pred.EVI.Z)

#get mean response of spline
MCMCsummary(samples[[1]],params = "b", round = 3, ISB = T)
bmeans <- MCMCsummary(samples[[1]],params = "b", round = 3, ISB = T)
camversionmeans <- MCMCsummary(samples[[1]],params = "a_version", round = 3, ISB = T)
daysactivemean <- MCMCsummary(samples[[1]],params = "a_daysactive", round = 3, ISB = T)
EVImean <- MCMCsummary(samples[[1]],params = "a_EVI", round = 3, ISB = T)
rho.mean <- plogis(camversionmeans$mean[1] + daysactivemean$mean*0.2211643 + EVImean$mean*pred.data$EVI +bmeans$mean[1]*pred.EVI.Z[,1] + bmeans$mean[2]*pred.EVI.Z[,2] + bmeans$mean[3]*pred.EVI.Z[,3] +
                     bmeans$mean[4]*pred.EVI.Z[,4] + bmeans$mean[5]*pred.EVI.Z[,5])
#get CRIs for spline
#combine MCMC chains into one
allchains <- MCMCchains(samples[[1]], params =c("a_version", "a_daysactive", "a_EVI", "b"))
#loop through estimated parameter at each iteration
rho.preds <- array(dim = c(length(pred.data$EVI), length(allchains[,"a_version[1]"])))
for(j in 1:length(allchains[,"a_version[1]"])){
  rho.preds[,j] <- plogis(allchains[,"a_version[1]"][j] + allchains[,"a_daysactive"][j]*0.2211643 + allchains[,"a_EVI"][j]*pred.data$EVI + 
                            allchains[,"b[1]"][j]*pred.EVI.Z[,1] + allchains[,"b[2]"][j]*pred.EVI.Z[,2] + allchains[,"b[3]"][j]*pred.EVI.Z[,3] +
                            allchains[,"b[4]"][j]*pred.EVI.Z[,4] + allchains[,"b[5]"][j]*pred.EVI.Z[,5])
}
#calculate interval
CL <- apply(rho.preds, 1, function(x){quantile(x, prob = c(0.025, 0.975))})

rhoplotdf <- cbind(pred.data2, rho.mean, t(CL))
rhoplotdf$EVI <- rhoplotdf$EVI * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
colnames(rhoplotdf)[10:11] <- c("CL2.5", "CL97.5")
ggplot(rhoplotdf, aes(x=EVI, y=rho.mean)) + geom_line(color="blue") + 
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5), alpha=0.3,fill="blue", color="blue")

######################################  predict spatial spline #####################################################
#make prediction grid for spatial spline
cellsize <- rep(sqrt(5.18e7), 2)#20mi^2
predict.grid <- st_make_grid(st_union(bearrange2), cellsize, what="centers")
cellsize <- rep(sqrt(2.59e7), 2)#10mi^2
knots.grid <- st_make_grid(st_union(bearrange2), cellsize, what="centers")
knots.grid2 <- st_coordinates(knots.grid)
bearrangeknots <- cover.design(knots.grid2, 50)
predict.grid.polys <- st_make_grid(st_union(bearrange2), cellsize, what="polygons")
predict.grid2 <- as.data.frame(do.call(rbind, st_intersection(predict.grid, bearrange2)))%>%rename("X"="V1", "Y"="V2")
predict.grid2.polys <- st_intersection(predict.grid.polys, bearrange2)
predict.grid2.polys <- as.data.frame(do.call(rbind, st_intersection(predict.grid.polys, bearrange2)))
predict.grid3 <- predict.grid2 %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y)
sp.cov.dist_pred = fields::rdist(x1=cbind(predict.grid3$X.scale,predict.grid3$Y.scale),x2=spknots)
sp.Z_K.pred = sp.cov.dist_pred^2*log(sp.cov.dist_pred) # basis
sp.Z.pred <- t(solve(sp.sqrt.omega_all,t(sp.Z_K.pred)))
sp.Z.pred <- (sp.Z.pred - sp.meanZ)/sp.sdZ #for prediction?
bXY <- MCMCsummary(samples[[1]], 
                   params = c('b_X', 'b_Y'),
                   ISB = FALSE,
                   round=2)
spatspline.bs <- MCMCsummary(samples[[1]], 
                             params = c('spat.spline.b'),
                             ISB = TRUE,
                             round=2)
lambda.pred.spline <- bXY$mean[1]*predict.grid3$X.scale + bXY$mean[2]*predict.grid3$Y.scale

lambda.pred.spline2 <- exp(lambda.pred.spline + inprodsplineman)
inprodspline <- sp.Z.pred%*%spatspline.bs$mean
predict.grid3$lambda.pred <- lambda.pred.spline2
predict.grid4 <- st_buffer(st_as_sf(predict.grid3, coords=c("X", "Y"), crs=3071), dist = cellsize/2)
plot(predict.grid4["lambda.pred"])
plot(st_geometry(bearrange2), add=TRUE)

numpal <- colorNumeric(
  palette = "viridis",                          # Use a ColorBrewer palette name
  domain = predict.grid4$lambda.pred                     # The numeric range of values
)
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addPolygons(data=st_transform(bearrange2, 4326), color="black", fillColor=NA)%>%
  addPolygons(data=st_transform(predict.grid4, crs=4326), fillColor=~numpal(lambda.pred), fillOpacity = 1)


#---------------------------------------------------------------------
#                       predict overall lambda                        
#---------------------------------------------------------------------
predict.centers <- st_make_grid(st_union(bearrange2), cellsize=c(2000,2000), what="centers")
plot(st_geometry(predict.centers))
predict.lambda <- as.data.frame(do.call(rbind, st_intersection(predict.centers, bearrange2)))%>%rename("X"="V1", "Y"="V2")
predict.lambda.sf <- st_as_sf(predict.lambda, coords=c("X", "Y"), crs=3071)
predict.lambda.sf2 <- st_join(predict.lambda.sf, st_transform(st_make_valid(get_spatial_data("bear_zones")), 3071))
#get Zs for spatial spline
predict.lambda2 <- predict.lambda %>%
  mutate(X.scale = (X-mean_x)/sd_x,    #need mean_xy and sd_xy from original data frame of camera points
         Y.scale = (Y-mean_y)/sd_y)
sp.cov.dist_predall = fields::rdist(x1=cbind(predict.lambda2$X.scale,predict.lambda2$Y.scale),x2=bearrangeknots2)
sp.Z_K.predall = sp.cov.dist_predall^2*log(sp.cov.dist_predall) # basis
sp.Z.predall <- t(solve(sp.sqrt.omega_all,t(sp.Z_K.predall)))
sp.Z.predall <- (sp.Z.predall - sp.meanZ)/sp.sdZ  #standardize on same scale as camera data
spatspline.bs <- MCMCsummary(samples,  
                             params = c('spat.spline.b'),
                             ISB = TRUE,
                             round=2)
inprodspline <- sp.Z.predall%*%spatspline.bs$mean

#-------------------------fixed effects
# Forest and Developed land cover
yearX <- unique(sitecovs$year)
b_Fixed <- MCMCsummary(samples, 
                     params = c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Yr"),
                     ISB = TRUE,
                     round=2)
MCMCtrace(samples, 
          params = c('abundance_scale'),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
Wiscland3 <- get_spatial_data(layer_name = 'wiscland2', level = 3)
Wiscland3.3071 <- project(Wiscland3, "EPSG:3071")
WisclandGuide <- read_xlsx("C:/Users/wildeefb/Documents/GeoSpatial/wiscland2/user_guide/Wiscland2 Color Scheme.xlsx",
                           range = "B2:C70", .name_repair = make.names)
lm_pred <- sample_lsm(
  # raster layer
  landscape = Wiscland3.3071,
  # camera locations
  y = predict.lambda.sf,
  # get landcover class level metrics
  level = "class",
  # return NA values for classes not in buffer
  # all_classes = TRUE, 
  # can do multiple metrics at once
  what = 'lsm_c_pland',
  # buffer sizes to use
  size = 500, 
  # default is square buffer
  shape = "square", 
  # turn warnings on or off
  verbose = FALSE 
)
lm_output <- left_join(lm_pred, WisclandGuide, by= join_by(class == dn.label))
lm_output$label <- gsub(pattern = "\\W", replacement = "", x = lm_output$label)

lm_output <- 
  lm_output %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label),
    values_from = c(value),
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  )

forestLCs <- c("AspenPaperBirch", "RedMaple", "Oak", "CentralHardwoods", "NorthernHardwoods","AspenForestedWetland", "BottomlandHardwoods", "SwampHardwoods",
               "MixedDeciduousConiferousForest", "MixedDeciduousConiferousForestedWetland")
developedLCs <- c("DevelopedHighIntensity","DevelopedLowIntensity")

forest.pred <- rowSums(lm_output[,forestLCs])
forest.pred2 <- (forest.pred - attr(sitecovs$Forest_500, "scaled:center"))/attr(sitecovs$Forest_500, "scaled:scale")
predict.lambda3 <- cbind(predict.lambda2, forest.pred2)

dev.pred <- rowSums(lm_output[,developedLCs])
dev.pred2 <- (dev.pred - attr(sitecovs$Developed_500, "scaled:center"))/attr(sitecovs$Developed_500, "scaled:scale")
predict.lambda3 <- cbind(predict.lambda3, dev.pred2)

#Disturbance
aoi <- getAOI(Wisconsin)
products <- "HDIST2023"
email <- "eli.wildey@wisconsin.gov"
projection <- 3071
resolution <- 90
path <- tempfile(fileext = ".zip")#"C:/Users/wildeefb/Documents/GeopSpatial/LANDFIRE/HDist2023.zip"
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
levels(hdist) <- dbf_tbl[,c(1,10)]
hdist <- addCats(hdist, value=dbf_tbl[,c(4:9)])
cats(hdist)
activeCat(hdist)
levels(hdist)
plot(hdist)
activeCat(hdist) <- 1


hdist2 <- ifel(hdist > 1, 1, hdist)
# Wiscland2 has 30m resolution

# Wiscland 3 prop land cover
lm_dist <- sample_lsm(
      # raster layer
      landscape = hdist2,
      # camera locations
      y = predict.lambda.sf,
      # get landcover class level metrics
      level = "class",
      # return NA values for classes not in buffer
      # all_classes = TRUE, 
      # can do multiple metrics at once
      what = 'lsm_c_pland',
      # buffer sizes to use
      size = 1000, 
      # default is square buffer
      shape = "square", 
      # turn warnings on or off
      verbose = FALSE 
    )


lm_dist$label <- ifelse(lm_dist$class == 0, "10+", "0-10")

# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
lm_dist <- 
  lm_dist %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    id_cols = plot_id,
    names_from = c(label),
    values_from = c(value),
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  )

lm_dist <- lm_dist%>%select(-matches("\\+"))
colnames(lm_dist)[2] <- gsub(x = colnames(lm_dist)[2], pattern = "0-10", "Dist1000")
lm_dist2 <- (lm_dist$Dist1000 - attr(sitecovs$Dist_1000, "scaled:center"))/attr(sitecovs$Dist_1000, "scaled:scale")

predict.lambda3 <- cbind(predict.lambda3, lm_dist2)

#Corn
CropRasts <- list.files("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/", full.names = TRUE)
CropRasts2 <- lapply(CropRasts, function (x) rast(x))
data("linkdata")
predict.lambdaCDL <- st_coordinates(st_transform(predict.lambda.sf, crs=crs(CropRasts2[[1]])))


lm_corn <- 
  map_dfr(.x= CropRasts2, 
         # produce a dataframe after this is all done
         ~sample_lsm(
                         # raster layer
                         landscape = .x,
                         # camera locations
                         y = predict.lambdaCDL,
                         # get landcover class level metrics
                         level = "class",
                         # return NA values for classes not in buffer
                         # all_classes = TRUE, 
                         # can do multiple metrics at once
                         what = 'lsm_c_pland',
                         # buffer sizes to use
                         size = 500, 
                         # default is square buffer
                         shape = "square", 
                         # turn warnings on or off
                         verbose = FALSE 
                       ), .id = "year"
         )



lm_corn2 <- left_join(lm_corn, linkdata, by=join_by(class==MasterCat))
Corn.pred <- lm_corn2[grep(pattern = "Corn", x = lm_corn2$Crop, ignore.case = TRUE),]
Corn.pred2 <- Corn.pred%>%group_by(year,  plot_id)%>%summarise(CornProp=sum(value))
Corn.pred2$plot_id <- factor(Corn.pred2$plot_id, levels=1:nrow(predict.lambda))
Corn.pred3<- Corn.pred2%>%group_by(year)%>%complete(., plot_id, fill = list(CornProp=0))%>%ungroup()
# in this data frame plot_id = camera ID
# class = landcover type
# value = % of that landcover type in the buffer
# make each landcover type x buffer into a column
Corn.pred4 <- 
  Corn.pred3 %>%
  # MAY NEED TO ADD distinct() HERE???
  #distinct() %>% # this removes duplicate rows before pivot. Not sure why there are duplicate rows in the first place
  pivot_wider(
    names_from = year,
    names_prefix="Corn",
    values_from = CornProp,
    # give class 0 if it doesn't exist in buffer
    values_fill = 0
  ) 

Corn.pred5 <- as.data.frame(Corn.pred4%>%mutate(across(matches("Corn"), ~(. - attr(sitecovs$Corn500, "scaled:center"))/attr(sitecovs$Corn500, "scaled:scale"))))



lambda.predicted.grid <- lapply(1:length(yearX), function (i) 
  lambda <- exp(b_Fixed$mean[5]*yearX[i] + b_Fixed$mean[1]*predict.lambda3$dev.pred2 + 
                b_Fixed$mean[2]*predict.lambda3$forest.pred2 + b_Fixed$mean[3]*as.numeric(Corn.pred5[,i+1]) +
                b_Fixed$mean[4]*predict.lambda3$lm_dist2 + sp.Z.predall%*%spatspline.bs$mean))
names(lambda.predicted.grid) <- 2019:2024
lambda.predicted.grid2 <- do.call(cbind, lambda.predicted.grid)
predict.lambda4 <- cbind(predict.lambda3, lambda.predicted.grid2)
colnames(predict.lambda4)[10:15] <- paste0("TotalLambda", 2019:2024)

#####################################################################################
#                             abundance by zone                                     #
#####################################################################################
#add zone information to predictions
predict.lambda5 <- cbind(predict.lambda4, predict.lambda.sf2$bear_mgmt_zone_id)
colnames(predict.lambda5)[16] <- "Zone"
popbyzone <- predict.lambda5%>%group_by(Zone)%>%summarise(across(matches("TotalLambda"), ~sum(.x)))%>%drop_na()
total <- data.frame("Zone"="Total", t(colSums(popbyzone[,2:7])))
popbyzone <- rbind(popbyzone, total)
popbyzone2 <- pivot_longer(popbyzone, cols = -c(Zone), names_to = "Year", values_to = "lambda")
popbyzone2$Year <- as.numeric(gsub(x = popbyzone2$Year, pattern = "TotalLambda", replacement = ""))
ggplot(filter(popbyzone2, Zone != "Total"), aes(x=Year, y=lambda, colour = Zone)) + geom_point() + geom_line()

#####################################################################################
#                   spatial prediction across years                                 #
#####################################################################################
predict.lambda4SF <- st_as_sf(predict.lambda4, coords=c("X", "Y"), crs=3071)
predict.lambda4.polys <- st_buffer(st_as_sf(predict.lambda4, coords=c("X", "Y"), crs=3071), dist = 1000)
predict.lambda4stars <- st_as_stars(predict.lambda4.polys)
plot(predict.lambda4.polys["TotalLambda2019"], key.pos = 1)
plot(st_geometry(bearrange2))
plot(predict.lambda4stars["TotalLambda2024"])
saveRDS(lambda.predicted.grid, "./predictedlambdasSummer.rds")



############################ scrap ######################################
#spatial spline
coordsdf <- ModelingDF1%>%select("cam_site_id", "X", "Y")%>%distinct()
coordsmatrix <- coordsdf%>%st_drop_geometry()%>%select("X", "Y")%>%as.matrix()
#make grid of potential knots based on bear range
cellsize <- rep(sqrt(2.59e7), 2)#10mi^2
knots.grid <- st_make_grid(st_union(bearrange2), cellsize, what="centers")
knots.grid2 <- as.data.frame(do.call(rbind, st_intersection(knots.grid, bearrange2)))%>%rename("X"="V1", "Y"="V2")
#calculate knots from potential knots
bearrangeknots <- cover.design(knots.grid2, 50)
bearrangeknots2 <- as.data.frame(bearrangeknots$design)
# scale coordinates 
mean_x <- mean(coordsdf$X)
sd_x <- sd(coordsdf$X)
mean_y <- mean(coordsdf$Y)
sd_y <- sd(coordsdf$Y)

bearrangeknots2 <- bearrangeknots2 %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y) %>%
  dplyr::select(X.scale,Y.scale)

dat <- coordsdf %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y)

# get matrix ready for spatial smoothing
spknots.dist <- dist(bearrangeknots2,"euclidean",diag=T,upper=T)
sp.omega_all = spknots.dist^2*log(spknots.dist) # basis
sp.svd.omega_all <- svd(sp.omega_all)
sp.sqrt.omega_all <- t(sp.svd.omega_all$v %*%
                         (t(sp.svd.omega_all$u)*sqrt(sp.svd.omega_all$d)))

# now for spline for data
sp.cov.dist_all = fields::rdist(x1=cbind(dat$X.scale,dat$Y.scale),x2=bearrangeknots2)
sp.Z_K = sp.cov.dist_all^2*log(sp.cov.dist_all) # basis
sp.Z <- t(solve(sp.sqrt.omega_all,t(sp.Z_K)))
sp.meanZ <- mean(sp.Z)
sp.sdZ <- sd(sp.Z)
sp.Z <- (sp.Z - sp.meanZ)/sp.sdZ



library(leaflet)
Wisc.BearRange4326 <- st_transform(bearrange2, crs=4326)
numpal <- colorNumeric(
  palette = "YlOrRd",                          # Use a ColorBrewer palette name
  domain = Wisconsin$INSIDE_WI_                     # The numeric range of values
)
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addCircleMarkers(data=st_transform(notbadinits, 4326), popup=notbadinits$cam_site_id, fillColor = "yellow", fillOpacity = 1,   stroke=F, radius=5)%>%
  addCircleMarkers(data=st_transform(badinits3, 4326), popup=badinits3$cam_site_id, fillColor = "black", fillOpacity = 1,   stroke=F, radius=5)%>%
  addPolygons(data=st_transform(counties, crs=4326))

leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addCircleMarkers(data=st_transform(predictgrid2, 4326), fillColor = "yellow", fillOpacity = 1,   stroke=F, radius=3)%>%
  addPolygons(data=st_transform(bearrange2, crs=4326))

predict.lambdaNAzone <- predict.lambda.sf2[is.na(predict.lambda.sf2$bear_mgmt_zone_id),]
bearzones <- get_spatial_data("bear_zones")
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addCircleMarkers(data=st_transform(predict.lambdaNAzone, 4326), fillColor = "yellow", fillOpacity = 1,   stroke=F, radius=3)%>%
  addPolygons(data=st_transform(bearzones, crs=4326), label = bearzones$bear_mgmt_zone_id)



inprodsplineman <- spatspline.bs$mean[1]*sp.Z.pred[,1] + spatspline.bs$mean[2]*sp.Z.pred[,2] + 
  spatspline.bs$mean[3]*sp.Z.pred[,3] + spatspline.bs$mean[4]*sp.Z.pred[,4] + 
  spatspline.bs$mean[5]*sp.Z.pred[,5] + spatspline.bs$mean[6]*sp.Z.pred[,6] + 
  spatspline.bs$mean[7]*sp.Z.pred[,7] + spatspline.bs$mean[8]*sp.Z.pred[,8] + 
  spatspline.bs$mean[9]*sp.Z.pred[,9] + spatspline.bs$mean[10]*sp.Z.pred[,10] + 
  spatspline.bs$mean[11]*sp.Z.pred[,11] + spatspline.bs$mean[12]*sp.Z.pred[,12] + 
  spatspline.bs$mean[13]*sp.Z.pred[,13] + spatspline.bs$mean[14]*sp.Z.pred[,14] + 
  spatspline.bs$mean[15]*sp.Z.pred[,15] + spatspline.bs$mean[16]*sp.Z.pred[,16] + 
  spatspline.bs$mean[17]*sp.Z.pred[,17] + spatspline.bs$mean[18]*sp.Z.pred[,18] + 
  spatspline.bs$mean[19]*sp.Z.pred[,19] + spatspline.bs$mean[20]*sp.Z.pred[,20] + 
  spatspline.bs$mean[21]*sp.Z.pred[,21] + spatspline.bs$mean[22]*sp.Z.pred[,22] + 
  spatspline.bs$mean[23]*sp.Z.pred[,23] + spatspline.bs$mean[24]*sp.Z.pred[,24] + 
  spatspline.bs$mean[25]*sp.Z.pred[,25] + spatspline.bs$mean[26]*sp.Z.pred[,26] + 
  spatspline.bs$mean[27]*sp.Z.pred[,27] + spatspline.bs$mean[28]*sp.Z.pred[,28] + 
  spatspline.bs$mean[29]*sp.Z.pred[,29] + spatspline.bs$mean[30]*sp.Z.pred[,30] + 
  spatspline.bs$mean[31]*sp.Z.pred[,31] + spatspline.bs$mean[32]*sp.Z.pred[,32] + 
  spatspline.bs$mean[33]*sp.Z.pred[,33] + spatspline.bs$mean[34]*sp.Z.pred[,34] + 
  spatspline.bs$mean[35]*sp.Z.pred[,35] + spatspline.bs$mean[36]*sp.Z.pred[,36] + 
  spatspline.bs$mean[37]*sp.Z.pred[,37] + spatspline.bs$mean[38]*sp.Z.pred[,38] + 
  spatspline.bs$mean[39]*sp.Z.pred[,39] + spatspline.bs$mean[40]*sp.Z.pred[,40] + 
  spatspline.bs$mean[41]*sp.Z.pred[,41] + spatspline.bs$mean[42]*sp.Z.pred[,42] + 
  spatspline.bs$mean[43]*sp.Z.pred[,43] + spatspline.bs$mean[44]*sp.Z.pred[,44] + 
  spatspline.bs$mean[45]*sp.Z.pred[,45] + spatspline.bs$mean[46]*sp.Z.pred[,46] + 
  spatspline.bs$mean[47]*sp.Z.pred[,47] + spatspline.bs$mean[48]*sp.Z.pred[,48] + 
  spatspline.bs$mean[49]*sp.Z.pred[,49] + spatspline.bs$mean[50]*sp.Z.pred[,50]
