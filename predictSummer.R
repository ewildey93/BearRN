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
library(MCMCvis)
library(data.table)
library(viridis)


samplesSummer <- readRDS("./RNsamplesFullModelBearRange5.rds")
bearrangeknots2 <- read.csv("./bearrangeknots.csv")
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

####################################predict EVI spline############################################
#predicting EVI spline
EVI.pred <- seq(from=quantile(EVI2$meanEVI, probs=0.01), to=quantile(EVI2$meanEVI, probs=0.99), length.out=100)
pred.data <- expand.grid(camversion=1, daysactive=0.2211643, EVI=EVI.pred)
# prediction dataset for spatial smoothing
pred.EVI.Z_K <- (abs(outer(as.numeric(pred.data$EVI),EVI.knots,"-")))^3

pred.EVI.Z <- t(solve(EVI.sqrt.OMEGA_all,t(pred.EVI.Z_K)))

# standardize for better performance
pred.EVI.Z <- (pred.EVI.Z-meanZ)/sdZ
pred.data2 <- cbind(pred.data, pred.EVI.Z)

#get mean response of spline
bmeans <- MCMCsummary(samplesSummer[1:3],params = "b", round = 3, ISB = T)
camversionmeans <- MCMCsummary(samplesSummer[1:3],params = "a_version", round = 3, ISB = T)
daysactivemean <- MCMCsummary(samplesSummer[1:3],params = "a_daysactive", round = 3, ISB = T)
EVImean <- MCMCsummary(samplesSummer[1:3],params = "a_EVI", round = 3, ISB = T)
rho.mean <- plogis(camversionmeans$mean[1] + daysactivemean$mean*0.2211643 + EVImean$mean*pred.data$EVI +bmeans$mean[1]*pred.EVI.Z[,1] + bmeans$mean[2]*pred.EVI.Z[,2] + bmeans$mean[3]*pred.EVI.Z[,3] +
                     bmeans$mean[4]*pred.EVI.Z[,4] + bmeans$mean[5]*pred.EVI.Z[,5])
#get CRIs for spline
#combine MCMC chains into one
allchains <- MCMCchains(samplesSummer[1:3], params =c("a_version", "a_daysactive", "a_EVI", "b"))
allchainssample <- allchains[sample(nrow(allchains), 10000), ]
#loop through estimated parameter at each iteration
rho.preds <- array(dim = c(length(pred.data$EVI), 10000))
for(j in 1:10000){
  rho.preds[,j] <- plogis(allchainssample[,"a_version[1]"][j] + allchainssample[,"a_daysactive"][j]*0.2211643 + allchainssample[,"a_EVI"][j]*pred.data$EVI + 
                            allchainssample[,"b[1]"][j]*pred.EVI.Z[,1] + allchainssample[,"b[2]"][j]*pred.EVI.Z[,2] + allchainssample[,"b[3]"][j]*pred.EVI.Z[,3] +
                            allchainssample[,"b[4]"][j]*pred.EVI.Z[,4] + allchainssample[,"b[5]"][j]*pred.EVI.Z[,5])
}
#calculate interval
CL <- apply(rho.preds, 1, function(x){quantile(x, prob = c(0.025, 0.975))})

rhoplotdf <- cbind(pred.data2, rho.mean, t(CL))
rhoplotdf$EVI <- rhoplotdf$EVI * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
colnames(rhoplotdf)[10:11] <- c("CL2.5", "CL97.5")
ggplot(rhoplotdf, aes(x=EVI, y=rho.mean)) + geom_line(color="#0072B2") + 
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5), alpha=0.3,fill="#0072B2", color="#0072B2") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none"
  ) +
  labs(
    title="Effect of Green-Up on\n Detection",
    y=expression("Individual Detection Probability"),
    x="Greenness of Vegetation"
  )

#--------------------------------------------------------------------
#                      detection betas
#--------------------------------------------------------------------
a_Fixed <- MCMCsummary(samples, 
                        params = c("a_version", "a_daysactive", "a_EVI"),
                        ISB = TRUE,
                        round=2)
a_Fixed <- data.frame("Var"=c("Version 3", "Version 4", "Effort", "NDVI"),a_Fixed)
a_Fixed$Var <- factor(a_Fixed$Var, levels=c("Version 3", "Version 4", "Effort", "NDVI"))
a_Fixed <- a_Fixed%>%arrange(Var)

ggplot(a_Fixed, aes(x=Var, y=mean)) + geom_pointrange(aes(ymin=X2.5., ymax=X97.5.), size=1, lwd=1) +
  coord_flip() + geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(
          face = ifelse(a_Fixed$Var %in% c("Version 3", "Version 4", "Effort"), "bold", "plain"),
          hjust = 0.5
        )
  ) +
  labs(
    y=expression("Posterior " * beta *""),
  ) #+
  #scale_color_manual(values=rev(cbPalette)) +
  #scale_x_discrete(labels = c("Forest Disturbance" = "Forest\nDisturbance"))



#---------------------------------------------------------------------
#                       predict overall lambda                        
#---------------------------------------------------------------------
predict.centers <- st_make_grid(st_union(bearrange2), cellsize=c(2000,2000), what="centers")
plot(st_geometry(predict.centers))
predict.lambda <- as.data.frame(do.call(rbind, st_intersection(predict.centers, bearrange2)))%>%rename("X"="V1", "Y"="V2")
predict.lambda.sf <- st_as_sf(predict.lambda, coords=c("X", "Y"), crs=3071)
plot(st_geometry(predict.lambda.sf))
predict.lambda.sf2 <- st_join(predict.lambda.sf, st_transform(st_make_valid(get_spatial_data("bear_zones")), 3071))
#get Zs for spatial spline
predict.lambda2 <- predict.lambda %>%
  mutate(X.scale = (X-mean_x)/sd_x,    #need mean_xy and sd_xy from original data frame of camera points
         Y.scale = (Y-mean_y)/sd_y)
bearrangeknots2 <- read.csv("./bearrangeknots.csv")
bearrangeknots2 <- bearrangeknots2 %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y) %>%
  dplyr::select(X.scale,Y.scale)
sp.cov.dist_predall = fields::rdist(x1=cbind(predict.lambda2$X.scale,predict.lambda2$Y.scale),x2=bearrangeknots2)
sp.Z_K.predall = sp.cov.dist_predall^2*log(sp.cov.dist_predall) # basis
sp.Z.predall <- t(solve(sp.sqrt.omega_all,t(sp.Z_K.predall)))
sp.Z.predall <- (sp.Z.predall - sp.meanZ)/sp.sdZ  #standardize on same scale as camera data
spatspline.bs <- MCMCsummary(samplesSummer[1:3],  
                             params = c('spat.spline.b'),
                             ISB = TRUE,
                             round=2)
inprodsplineSummer <- exp(sp.Z.predall%*%spatspline.bs$mean)
inprodspline2Summer <- st_buffer(predict.lambda.sf, dist = 1000, endCapStyle = "SQUARE")
inprodspline2Summer <- cbind(inprodspline2Summer, inprodsplineSummer)

# 1. Determine global min and max
global_min <- min(c(inprodspline2Spring$inprodsplineSpring, inprodspline2Summer$inprodsplineSummer), na.rm = TRUE)
global_max <- max(c(inprodspline2Spring$inprodsplineSpring, inprodspline2Summer$inprodsplineSummer), na.rm = TRUE)
common_limits <- c(global_min, global_max)
# 0.007108554 6.018643190



ggplot() + geom_sf(data=inprodspline2Summer, aes(fill=inprodsplineSummer), color=NA) +
   geom_sf(data=Wisconsin2, fill=NA) +
   geom_sf(data=bearzones, color="black", fill=NA, linewidth=1.25) +
  scale_fill_viridis(option="C") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),        # Removes x and y axis tick numbers/text
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(fill=expression("Mean site\nabundance"),
       title="Effect of Spatial Spline on Abundance"
       )
colorpal <- colorNumeric(palette = "magma", 
                        domain=inprodspline2Summer$inprodsplineSummer)
bearzones <- st_make_valid(get_spatial_data("bear_zones"))
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addPolygons(data=st_transform(inprodspline2Summer, crs=4326), color = ~colorpal(inprodspline2Summer$inprodsplineSummer), group="lambda")%>%
  addPolygons(data=st_transform(bearzones, 4326), fill=NA, color="black", group="zones")%>%
  addLayersControl(overlayGroups = c("lambda", "zones"))
#-------------------------fixed effectinprodspline2#-------------------------fixed effects
# Forest and Developed land cover
yearX <- unique(sitecovs$year)
b_Fixed <- MCMCsummary(samplesSummer[1:3], 
                     params = c('b_Dev', 'b_Forest', "b_Corn", "b_Dist"),
                     ISB = TRUE,
                     round=2)
b_ZoneEff <- MCMCsummary(samplesSummer[1:3], 
                        params = c("b_ZoneYr", "b_Zone"),
                        ISB = TRUE,
                        round=2)
b_ZoneEff2 <- data.frame("beta.Zone"=c(b_ZoneEff$mean[7:12]), 
                         "beta.ZoneYr"= c(b_ZoneEff$mean[1:6]) , 
                         "Zone"=LETTERS[1:6])
colnames(predict.lambda.sf2)[1] <- "Zone"
betas.ZoneYr <- left_join(predict.lambda.sf2[,1], b_ZoneEff2, by="Zone")
predict.lambda3[, 64:66] <- st_drop_geometry(betas.ZoneYr)
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
predict.lambda3 <- cbind(predict.lambda3, forest.pred2)

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
#dbf2 <- left_join(dbf_tbl, HDistLU, by=join_by(Value ==VALUE))%>%mutate(EarlySuccess=ifelse(Value > 0, "0-10", "10+"))
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

predict.lambda3 <- cbind(predict.lambda3, Corn.pred5[,2:7])
predict.lambda3 <- cbind(predict.lambda3, sp.Z.predall)

saveRDS(predict.lambda3, "./predict.grid.summer.rds")
predict.lambda3 <- readRDS("./predict.grid.summer.rds")
sp.Z.predall2 <- as.matrix(predict.lambda3[,14:63])
lambda.predicted.grid <- lapply(1:length(yearX), function (i) 
  lambda <- exp(predict.lambda3$beta.Zone  + predict.lambda3$beta.ZoneYr*yearX[i] + 
                b_Fixed$mean[1]*predict.lambda3$dev.pred2 + 
                b_Fixed$mean[2]*predict.lambda3$forest.pred2 + b_Fixed$mean[3]*predict.lambda3[,i+7] +
                b_Fixed$mean[4]*predict.lambda3$lm_dist2 + sp.Z.predall%*%spatspline.bs$mean))
names(lambda.predicted.grid) <- 2019:2024
lambda.predicted.grid2 <- do.call(cbind, lambda.predicted.grid)
predict.lambda4 <- cbind(predict.lambda3, lambda.predicted.grid2)
colnames(predict.lambda4)[67:72] <- paste0("TotalLambda", 2019:2024)
saveRDS(predict.lambda4, "./predict.lambda4Zones.rds")
#####################################################################################
#                             abundance by zone                                     #
#####################################################################################
#add zone information to predictions
predict.lambda5 <- cbind(predict.lambda4, predict.lambda.sf2$bear_mgmt_zone_id)
colnames(predict.lambda5)[70] <- "Zone"
popbyzone <- predict.lambda4%>%group_by(Zone)%>%summarise(across(matches("TotalLambda"), ~sum(.x)))%>%drop_na()
total <- data.frame("Zone"="Total", t(colSums(popbyzone[,2:7])))
popbyzone <- rbind(popbyzone, total)
popbyzone2 <- pivot_longer(popbyzone, cols = -c(Zone), names_to = "Year", values_to = "lambda")
popbyzone2$Year <- as.numeric(gsub(x = popbyzone2$Year, pattern = "TotalLambda", replacement = ""))
ggplot(filter(popbyzone2, Zone != "Total"), aes(x=Year, y=lambda, colour = Zone)) + geom_point(size=3) + geom_line(lwd=1) +
  scale_color_manual(values=cbPalette)

#variation
#combine MCMC chains into one
allchains <- MCMCchains(samplesSummer[1:3], params =c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Zone", "b_ZoneYr", "spat.spline.b"), ISB = TRUE)
#loop through estimated parameter at each iteration

allchainssample <- allchains[sample(nrow(allchains), 10000), ]

predict.df <- predict.lambda4
CIsforZones <- function(predict.df, nZones, year, yearX) {
  for (z in 1:nZones){
ZonePredict <- predict.df[!is.na(predict.df$Zone) & predict.df$Zone == LETTERS[z],]
#zfac <- ifelse(z==6, 5, z) #need to change this with E,F merge/un-merge
lambda.CIs <- array(dim = c(nrow(ZonePredict), 10000))
yr <- year-2018
  for(j in 1:10000){
      lambda.CIs[,j] <- exp(allchainssample[j, paste0("b_Zone[", z,"]")] +
                               allchainssample[j, paste0("b_ZoneYr[", z,"]")]*yearX[yr] +
                               allchainssample[j,"b_Dev"]*ZonePredict$dev.pred2 + 
                               allchainssample[j,"b_Forest"]*ZonePredict$forest.pred2 + 
                               allchainssample[j,"b_Corn"]*ZonePredict[,paste0("Corn", yr)] +
                               allchainssample[j,"b_Dist"]*ZonePredict$lm_dist2 +
                               as.matrix(ZonePredict[,14:63])%*%allchainssample[j, 17:66]) #this latter column specification needs to change depending on # of zones
                }
    zone.lambda.distr <- colSums(lambda.CIs)
    CL <- quantile(zone.lambda.distr, c(0.025, 0.975))
    if(z==1){
      CL2 <- data.frame("Zone"=LETTERS[z], "Year"=yr, "CL2.5"=CL[1], "CL97.5"=CL[2])
    }else{
      CL2 <- rbind(CL2, data.frame("Zone"=LETTERS[z], "Year"=yr, "CL2.5"=CL[1], "CL97.5"=CL[2]))
    }
  }
  return(CL2)
}
CIs2019 <- CIsforZones(predict.df = predict.lambda4, nZones = 6, year = 2019, yearX)
CIs2020 <- CIsforZones(predict.df = predict.lambda4, nZones = 6, year = 2020, yearX)
CIsZone <- lapply(2019:2024, function (i) CIsforZones(predict.df = predict.lambda3, nZones = 6, year = i, yearX))
CIsZone2 <- rbindlist(CIsZone)
CIsZone2$Year <- CIsZone2$Year+2018


popbyzone3 <- left_join(popbyzone2, CIsZone2, by=c("Zone", "Year"))
ggplot(filter(popbyzone3, Zone != "Total"), aes(x=Year, y=lambda, colour = Zone)) + 
  geom_point(size=2) +
  #geom_pointrange(aes(ymin = CL2.5, ymax = CL97.5,color= Zone), size=1, lineend = "square", lwd=1) + 
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5,color= Zone, fill=Zone), alpha=0.3) +
  scale_color_manual(values=cbPalette) +
  scale_fill_manual(values=cbPalette) +
  labs(y="Mean Relative Abundance")

#calculate interval
lambda.CIs2019.2 <- bind_cols(lambda.CIs2019, predict.lambda.sf2$bear_mgmt_zone_id)
colnames(lambda.CIs2019.2[10001]) <- "Zone"
zonepop.posterior <- lambda.CIs2019.2%>%group_by(Zone)
totalpop.posterior <- lambda.CIs2019
CL2019A <- apply(lambda.CIs, 1, function(x){quantile(x, prob = c(0.025, 0.975))})
#####################################################################################
#                   spatial prediction across years                                 #
#####################################################################################
predict.lambda4SF <- st_as_sf(predict.lambda4, coords=c("X", "Y"), crs=3071)
predict.lambda4.polys <- st_buffer(st_as_sf(predict.lambda4, coords=c("X", "Y"), crs=3071), dist = 1000, , endCapStyle = "SQUARE")
predict.lambda4stars <- st_as_stars(predict.lambda4.polys)
plot(predict.lambda4.polys["TotalLambda2019"], key.pos = 1)
plot(st_geometry(bearrange2))
plot(predict.lambda4stars["TotalLambda2024"])
saveRDS(lambda.predicted.grid, "./predictedlambdasSummer.rds")

ggplot() + geom_sf(data=predict.lambda4.polys, aes(fill=TotalLambda2024), color=NA) +
  geom_sf(data=Wisconsin2, fill=NA) +
  scale_fill_viridis(option="C") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),        # Removes x and y axis tick numbers/text
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(fill="Mean site\nabundance",
       title=expression("Mean abundance 2024 (" *lambda* ")"),
       subtitle = "Full model prediction"
  )

predict.lambda4.polyslong <- pivot_longer(predict.lambda4.polys, 
                                          cols = matches("TotalLambda"), 
                                          names_to = "Year",
                                          names_pattern = "(\\d{4})",
                                          values_to = "Lambda")
ggplot() + geom_sf(data=predict.lambda4.polyslong, aes(fill=Lambda), color=NA) +
  facet_wrap(~Year) +
  geom_sf(data=Wisconsin2, fill=NA) +
  scale_fill_viridis(option="D") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),        # Removes x and y axis tick numbers/text
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(fill=expression("Mean \nabundance")
  )


#########################################################################
#                        Fixed effects betas                            #
#########################################################################
b_Fixed2 <- MCMCsummary(samplesSummer[1:3], 
                        params = c('b_Dev', 'b_Forest', "b_Corn", "b_Dist"),
                        ISB = TRUE,
                        round=2)
b_Fixed2 <- data.frame("Var"=c("Developed", "Deciduous", "Corn", "Forest Disturbance"),b_Fixed2)
b_Fixed2$Var <- factor(b_Fixed2$Var, levels=rev(c("Corn",  "Deciduous", "Developed", "Forest Disturbance")))
b_Fixed2 <- b_Fixed2%>%arrange(Var)

FixedPalette <- c("#E69F00",  "#56B4E9","#009E73","#F0E442", "#0072B2")
ggplot(b_Fixed2, aes(x=Var, y=mean, colour = Var)) + geom_pointrange(aes(ymin=X2.5., ymax=X97.5.), size=1, lwd=1) +
  coord_flip() + geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(
          face = ifelse(b_Fixed2$Var %in% c("Developed", "Forest Disturbance"), "bold", "plain"),
          hjust = 0.5
        )
  ) +
    labs(y=expression("Posterior " * beta *"")) +
    scale_color_manual(values=FixedPalette) +
  scale_x_discrete(labels = c("Forest Disturbance" = "Forest\nDisturbance"))

b_ZoneEff3 <- MCMCsummary(samplesSummer[1:3], 
                          params = c('b_Zone', 'b_ZoneYr'),
                          ISB = TRUE,
                          round=2)
b_ZoneEff3 <- data.frame("Var"=c(paste0("Zone ", LETTERS[1:6]), paste0("ZoneYr ", LETTERS[1:6])), b_ZoneEff3, 
                         group=as.character(rep(1:6,2)))
b_ZoneEff3 <- b_ZoneEff3%>%arrange(desc(Var))
b_ZoneEff3$Var <- factor(b_ZoneEff3$Var, levels=b_ZoneEff3$Var)

ggplot(b_ZoneEff3, aes(x=Var, y=mean, group=group, color=group)) + geom_pointrange(aes(ymin=X2.5., ymax=X97.5.), size=1, lwd=1) +
  coord_flip() + geom_hline(yintercept = 0) +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(
          face = ifelse(b_ZoneEff3$Var %in% c("ZoneYr D"), "bold", "plain"),
          hjust = 0.5
        )
  ) +
  labs(y=expression("Posterior " * beta *"")) +
  scale_color_manual(values=cbPalette) +
  scale_x_discrete(labels = c("Forest Disturbance" = "Forest\nDisturbance"))


########################################################################
#                      Zone and Year Effects                           #
########################################################################
lambdaZoneYr <- lapply(1:length(yearX), function (i)
lambda <- exp(b_Fixed$mean[5]*yearX[i]  + predict.lambda3$beta.Zone + 
      predict.lambda3$beta.ZoneYr*yearX[i]*predict.lambda3$beta.Zone))
names(lambdaZoneYr) <- 2019:2024
lambdaZoneYr2 <- do.call(cbind, lambdaZoneYr)
lambdaZoneYr3 <- cbind(predict.lambda3[,c(1,2,64:66)], lambdaZoneYr2)
lambdaZoneYr3 <- st_as_sf(lambdaZoneYr3, coords=c("X", "Y"), crs=3071)
lambdaZoneYr3long <- lambdaZoneYr3%>%pivot_longer(., cols=c(matches("\\d{4}")), names_to = "year", values_to = "lambda")
lambdaZoneYr3long.polys <- st_buffer(lambdaZoneYr3long, dist = 1000, , endCapStyle = "SQUARE")
ggplot() + geom_sf(data=lambdaZoneYr3long.polys, aes(fill=lambda), color=NA) +
  facet_wrap(~year) +
  geom_sf(data=Wisconsin2, fill=NA) +
  scale_fill_viridis(option="D") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),        # Removes x and y axis tick numbers/text
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = "bottom") +
  labs(fill=expression("Mean \nabundance")
  )
########################################################################
#                        Scale Selection                               #
########################################################################
#combine MCMC chains into one
allchains.scalevars <- MCMCchains(samplesSummer[1:3], params =c("abundance_scale"), ISB = TRUE)
abundancescale1 <- data.frame(table(allchains.scalevars[,1]))
abundancescale1 <- rbind(abundancescale1, data.frame("Var1"="4", "Freq"=0))
abundancescale24 <- sapply(2:4, function (i) table(allchains.scalevars[,i]))
abundancescales <- cbind(abundancescale1, abundancescale24)
colnames(abundancescales) <- c("Buffer", "Developed", "Forest Disturbance", "Deciduous", "Corn")
buffers <- c(100,500,1000,2500,5000)/1000
abundancescales <- abundancescales%>%mutate(Buffer=case_when(Buffer==1 ~ buffers[2],
                                                             Buffer==2 ~ buffers[3],
                                                             Buffer==3 ~ buffers[4],
                                                             Buffer==4 ~ buffers[5]
                                                            ))
scalevars <- colnames(abundancescales[2:5])
abundancescales <- abundancescales%>%mutate(across(scalevars, ~./sum(.)))
abundancescaleslong <- pivot_longer(abundancescales, cols = scalevars, names_to = "Var", values_to = "Probability")
abundancescaleslong$Buffer <- as.factor(abundancescaleslong$Buffer)

cbPalette <- c("#F0E442", "#56B4E9", "#009E73", "#E69F00" , "#0072B2", "#D55E00")
ggplot(abundancescaleslong, aes(x=Buffer, y=Probability, fill=Var)) + 
  facet_wrap(~Var) +
  geom_col() +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position="none",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)
        ) +
  labs(
       title="Scale Selection of Environmental Variables",
       y="Posterior Probability",
       x="Spatial scale (km)"
  ) +
  scale_fill_manual(values=cbPalette)

#########################################################################
#             Marginal effect disturbance                               #
#########################################################################
range(sitecovs$Dist_1000)
Dist1000.pred <- seq(-0.25, 4, length.out=100)
Dist1000.pred.ogscale <- round(Dist1000.pred * attr(sitecovs$Dist_1000, "scaled:scale") + attr(sitecovs$Dist_1000, "scaled:center"), 2)
lambda.dist <- exp(b_ZoneEff$mean[7]  + 
                  b_ZoneEff$mean[1]*yearX[1] + 
                  b_Fixed$mean[1]*mean(sitecovs$Developed_500) + 
                  b_Fixed$mean[2]*mean(sitecovs$Forest_500) + b_Fixed$mean[3]*mean(as.matrix(sitecovs[,7:11])) +
                  b_Fixed$mean[4]*Dist1000.pred)
#get CRIs for spline
#combine MCMC chains into one
fixedvars <- c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Zone[1]", "b_ZoneYr[1]")
allchains <- MCMCchains(samplesSummer[1:3], params =fixedvars, ISB=FALSE)
#loop through estimated parameter at each iteration
allchainssample <- allchains[sample(nrow(allchains), 10000), ]
lambda.preds.dist <- array(dim = c(length(Dist1000.pred), 10000))
for(j in 1:10000){
  lambda.preds.dist[,j] <- exp(allchainssample[j,"b_Zone[1]"] + allchainssample[j,"b_ZoneYr[1]"]*yearX[1] + 
                                 allchainssample[j,"b_Dev"]*mean(sitecovs$Developed_500) + 
                                 allchainssample[j,"b_Forest"]*mean(sitecovs$Forest_500) + 
                                 allchainssample[j,"b_Corn"]*mean(as.matrix(sitecovs[,7:11])) +
                                 allchainssample[j,"b_Dist"]*Dist1000.pred) #+  sp.Z.predall%*%allchainssample[j, 6:55]
}
#calculate interval
CL.lambda.dist <- apply(lambda.preds.dist, 1, function(x){quantile(x, prob = c(0.025, 0.975))})

dist.plot.df <- cbind.data.frame(Dist1000.pred, Dist1000.pred.ogscale, lambda.dist, t(CL.lambda.dist))
colnames(dist.plot.df) <- c("DistScale", "Dist", "lambda", "CL2.5", "CL97.5")
ggplot(dist.plot.df, aes(x=Dist, y=lambda)) + geom_line(color="#E69F00") + 
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5), alpha=0.5,fill="#E69F00", color="#E69F00") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none"
  ) +
  labs(
    title="Marginal Effect of Forest Disturbance on\n Site Abundance",
    y=expression("Mean Site Abundance ( " * lambda * ")"),
    x="Proportion of Forest Disturbance (1km)"
  )

########################################################################
#                Marginal effect developed                             #
########################################################################
range(sitecovs$Developed_500)
Dev500.pred <- seq(-0.25, 6, length.out=100)
Dev500.pred.ogscale <- round(Dev500.pred * attr(sitecovs$Developed_500, "scaled:scale") + attr(sitecovs$Developed_500, "scaled:center"), 2)
lambda.dev <- exp(b_ZoneEff$mean[7]  + 
                    b_ZoneEff$mean[1]*yearX[1] +  + b_Fixed$mean[1]*Dev500.pred + 
                     b_Fixed$mean[2]*mean(sitecovs$Forest_500) + b_Fixed$mean[3]*mean(as.matrix(sitecovs[,7:11])) +
                     b_Fixed$mean[4]*mean(sitecovs$Dist_1000))
#get CRIs for spline
#loop through estimated parameter at each iteration
lambda.preds.dev <- array(dim = c(length(Dev500.pred), 10000))
for(j in 1:10000){
  lambda.preds.dev[,j] <- exp(allchainssample[j,"b_Zone[1]"] + allchainssample[j,"b_ZoneYr[1]"]*yearX[1]  + allchainssample[j,"b_Dev"]*Dev500.pred + 
                                 allchainssample[j,"b_Forest"]*mean(sitecovs$Forest_500) + 
                                 allchainssample[j,"b_Corn"]*mean(as.matrix(sitecovs[,7:11])) +
                                 allchainssample[j,"b_Dist"]*mean(sitecovs$Dist_1000)) #+  sp.Z.predall%*%allchainssample[j, 6:55]
}
#calculate interval
CL.lambda.dev<- apply(lambda.preds.dev, 1, function(x){quantile(x, prob = c(0.025, 0.975))})

dev.plot.df <- cbind.data.frame(Dev500.pred, Dev500.pred.ogscale, lambda.dev, t(CL.lambda.dev))
colnames(dev.plot.df) <- c("DevScale", "Dev", "lambda", "CL2.5", "CL97.5")
ggplot(dev.plot.df, aes(x=Dev, y=lambda)) + geom_line(color="#56B4E9") + 
  geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5), alpha=0.5,fill="#56B4E9", color="#56B4E9") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.5),
        legend.position="none"
  ) +
  labs(
    title="Marginal Effect of Human Development on\n Site Abundance",
    y=expression("Mean Site Abundance( " * lambda *")"),
    x="Proportion of Developed Land Cover (0.5km)"
  )

#######################################################################################################
############################ scrap ####################################################################
#######################################################################################################
predict.lambda4[,14:63] <- sp.Z.predall


plot(inprodspline2["inprodspline"])
plot(st_geometry(inprodspline2))
splinestars <- st_as_stars(inprodspline2)
plot(splinestars)
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
