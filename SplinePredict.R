library(fields)
library(leaflet)
library(stars)

samples <- readRDS("./RNsamplesFullModelBearRangeNoXY.rds")
samples2 <- samples[[1]]

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
sp.cov.dist_pred = fields::rdist(x1=cbind(predict.grid3$X.scale,predict.grid3$Y.scale),x2=bearrangeknots2)
sp.Z_K.pred = sp.cov.dist_pred^2*log(sp.cov.dist_pred) # basis
sp.Z.pred <- t(solve(sp.sqrt.omega_all,t(sp.Z_K.pred)))
sp.Z.pred <- (sp.Z.pred - sp.meanZ)/sp.sdZ #for prediction?
# bXY <- MCMCsummary(samples2, 
#                    params = c('b_X', 'b_Y'),
#                    ISB = FALSE,
#                    round=2)
spatspline.bs <- MCMCsummary(samples, 
                             params = c('spat.spline.b'),
                             ISB = TRUE,
                             round=2)
lambda.pred.spline <- bXY$mean[1]*predict.grid3$X.scale + bXY$mean[2]*predict.grid3$Y.scale

lambda.pred.spline2 <- exp(lambda.pred.spline + inprodsplineman)
inprodspline <- sp.Z.pred%*%spatspline.bs$mean
lambda.pred.spline <- exp(inprodspline)
predict.grid3$lambda.pred <- lambda.pred.spline
predict.grid4 <- st_buffer(st_as_sf(predict.grid3, coords=c("X", "Y"), crs=3071), dist = cellsize/2)
predict.grid5 <- st_as_stars(predict.grid4)
plot(predict.grid4["lambda.pred"])
plot(st_geometry(bearrange2))

PR <- rnorm(15000, 0, 2)
MCMCtrace(samples2, 
          params = c( 'b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Yr"),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples2, 
          params = c('spat.spline.b[1]', 'spat.spline.b[2]', 'spat.spline.b[3]', 'spat.spline.b[4]', 'spat.spline.b[5]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples2, 
          params = c('abundance_scale'),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples2, 
          params = c("a_version", "a_daysactive", "a_EVI"),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)



############################ scrap ######################################
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

numpal <- colorNumeric(
  palette = "viridis",                          # Use a ColorBrewer palette name
  domain = predict.grid4$lambda.pred                     # The numeric range of values
)
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addPolygons(data=st_transform(bearrange2, 4326), color="black", fillColor=NA)%>%
  addPolygons(data=st_transform(predict.grid4, crs=4326), fillColor=~numpal(lambda.pred), fillOpacity = 1)


#make bear range from counties layer
Wisconsin <- st_read("C:/Users/wildeefb/Documents/GeoSpatial/VectorLayers/Wisconsin_State_Boundary_24K.shp")
Wisconsin2 <- Wisconsin%>%filter(INSIDE_WI_ == 1)
counties <- get_spatial_data("counties")
counties$bearrange <- ifelse(counties$COUNTY_NAM %in% 
                               c("Vernon", "Crawford", "Richland", "Sauk", "Iowa", "Grant",
                                 "Lafayette", "Dane", "Green", "Rock", "Walworth", "Racine",
                                 "Kenosha", "Milwaukee", "Waukesha", "Jefferson", "Dodge",
                                 "Washington", "Ozaukee", "Sheboygan", "Fond du Lac",
                                 "Green Lake", "Winnebago", "Calumet", "Manitowoc", 
                                 "Kewaunee", "Door"), 0, 1)
table(counties$bearrange)
counties2 <- counties%>%select(COUNTY_NAM, bearrange)
bearrange <- counties2%>%filter(bearrange == 1)
bearrange <- st_cast(bearrange, "POLYGON")
#merge with Wisconsin sf layer to get rid of all the islands
bearrange2 <- st_intersection(bearrange, Wisconsin2)
#join camsites to bear range to get rid of camera sites outside of bear range
camsitessf <- st_as_sf(ModelingDF, coords=c("X","Y"),crs=3071)%>%distinct(cam_site_id, geometry)
camsites.bearrange <- st_join(camsitessf, bearrange2)%>%drop_na(bearrange)

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


#see where camsites with bad initial values are
badinits <- c('lambda[138, 1]', 'lambda[1255, 1]', 'lambda[1309, 1]', 'lambda[1107, 2]', 'lambda[1108, 2]',
              'lambda[1209, 3]', 'lambda[1210, 3]', 'lambda[1303, 4]', 'lambda[1289, 5]', 'lambda[1290, 5]',
              'lambda[1492, 6]')
badinits <- data.frame(siteID=str_extract(string = badinits, pattern = "(lambda\\[)(\\d+)", group = 2),
                       yearID=str_extract(string = badinits, pattern = "(lambda\\[\\d+, )(\\d+)", group = 2))
badinits2 <- camsiteskey%>%filter(paste0(siteID, yearID) %in% paste0(badinits$siteID, badinits$yearID))
#3/4 bad initis down in kenosha county, other one in superior, WI and has lambda of 0.58
badinits3 <- camsitessf%>%filter(cam_site_id %in% badinits2$cam_site_id)
notbadinits <- camsitessf%>%filter(!(cam_site_id %in% badinits2$cam_site_id))
badinitssitecovs <- sitecovs%>%filter(paste0(siteID, yearID) %in% paste0(badinits$siteID, badinits$yearID))