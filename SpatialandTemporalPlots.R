library(ggplot2)
library(tidyterra)

load("path to file/SpatialandTemporalPlots.RData")

#####################################################################################
#                             abundance by zone                                     #
#####################################################################################

yearX <- unique(as.numeric(data$Year[!is.na(data$Year)]))

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


lambda.predicted.grid2 <- do.call(cbind, lambda.predicted.grid)
predict.lambda4 <- cbind(covariate_matrix, lambda.predicted.grid2)
colnames(predict.lambda4)[66:71] <- paste0("TotalLambda", 2019:2024)





popbyzone <- predict.lambda4%>%group_by(zone)%>%summarise(across(matches("TotalLambda"), ~sum(.x)))
total <- data.frame("zone"="Total", t(colSums(popbyzone[,2:7])))
popbyzone <- rbind(popbyzone, total)
popbyzone2 <- pivot_longer(popbyzone, cols = -c(zone), names_to = "Year", values_to = "lambda")
popbyzone2$Year <- as.numeric(gsub(x = popbyzone2$Year, pattern = "TotalLambda", replacement = ""))
ggplot(filter(popbyzone2, zone != "Total"), aes(x=Year, y=lambda, colour = zone)) + geom_point(size=3) + geom_line(lwd=1)

#variation
#combine MCMC chains into one
allchains <- MCMCchains(samples, params =c('b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Wetland", "b_Zone", "b_ZoneYr", "spat.spline.b"), ISB = TRUE)
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
      ) 
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



#####################################################################################
#                   spatial prediction across years                                 #
#####################################################################################
Wisconsin <- st_read("C:/Users/wildeefb/Documents/GeoSpatial/VectorLayers/Wisconsin_State_Boundary_24K.shp")
Wisconsin2 <- Wisconsin%>%filter(INSIDE_WI_ == 1)

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
