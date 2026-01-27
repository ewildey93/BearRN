library(dplyr)
library(sf)
library(data.table)
library(tidyr)
library(nimble)
library(MCMCvis)
library(sswids)
library(stringr)
library(fields)
library(sswids)



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
                                 "Kewaunee", "Door", "Columbia"), 0, 1)
table(counties$bearrange)
counties2 <- counties%>%select(COUNTY_NAM, bearrange)
bearrange <- counties2%>%filter(bearrange == 1)
bearrange <- st_cast(bearrange, "POLYGON")
#merge with Wisconsin sf layer to get rid of all the islands
bearrange2 <- st_intersection(bearrange, Wisconsin2)


#5/21-8/26 this data frame doesn't have NAs for cam site-year-occs that dont have effort its just missing those rows
#results in n.occs column being wrong but don't think that matters
ModelingDF <- readRDS("./ModelingDFSummer.rds")%>%st_transform(., 3071)
ModelingDF1 <- st_join(ModelingDF, bearrange2)%>%drop_na(bearrange)%>%select(-c(37:44))
ModelingDF1 <- ModelingDF1%>%cbind(., st_coordinates(.))%>%
  st_drop_geometry()%>%ungroup()%>%filter(occ > 3 & occ < 18)%>%filter(BEAR_ADULT_AMT < 200) #5/21-8/26
#remove occasions which have multiple camera versions
ModelingDF1 <- ModelingDF1[-which(ModelingDF1$camera_version %in% c("V2,V4","V2,V3")),] 
only1occ <- ModelingDF1%>%group_by(cam_site_id, year)%>%summarise(N=n())%>%filter(N == 1)
ModelingDF1 <- ModelingDF1%>%filter(
  !(paste0(cam_site_id, year) %in% paste0(only1occ$cam_site_id,only1occ$year)))%>%relocate(Corn500, .before=Corn1000)
#split up data frame by year and reshape to wide format for detections histories for each year
dethist <- split(ModelingDF1, ModelingDF1$year)
dethist <- lapply(dethist, function (x) tidyr::pivot_wider(x, id_cols = c(cam_site_id, year), names_from = occ, values_from = BEAR_ADULT_AMT, names_sort = TRUE)%>%ungroup())
dethistall <- rbindlist(dethist)%>%group_by(year)%>%mutate(siteID=row_number())%>%ungroup()
Dethistlong <- dethistall%>%pivot_longer(cols=3:16, names_to="occ", values_to="det")%>%
  mutate(yearID=year-2018, det=ifelse(det > 0,1,0), occ=as.numeric(occ))%>%
  arrange(year, cam_site_id, occ)%>%ungroup()
Dethistlong2 <- Dethistlong%>%group_by(yearID,siteID)%>%arrange(yearID, siteID, is.na(det), occ)%>%mutate(occre=row_number())%>%ungroup()
nsite <- sapply(dethist,nrow)

maxsites <- max(nsite)
no.occs <- length(unique(ModelingDF1$occ))
no.years <- length(unique(ModelingDF1$year))
# create a 3D matrix for the counts (dim1 = site; dim2 = repeated visit; dim3 = year).
# Important note: The rows in the different array slices can be different sites among years.
# Array stores data in multiple dimensions and "nsite" is the number of sites in a
# given year (e.g. year 2 has 1078 sites).
Bear_All <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years



# In the "all" array (above), loop through each year of data (yearID) "t", each weekly
# count (rep) "j", each site (siteID) "i" (number of sites depends on year; addressed
# above and incorporated below), and assign counts.
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      Bear_All[ i, j, t] <- Dethistlong2[ c( Dethistlong2$siteID == i & Dethistlong2$occre == j & Dethistlong2$yearID == t), "det"]$det
    }
  }
}

#scale continuous variables and rename season to yearID
ModelingDF2 <- ModelingDF1%>%mutate(across(c(10, 12, 14:37), scale))%>%rename(yearID=season)
#create vector of selection scale variables
scalevars <- grep(x = colnames(ModelingDF2), pattern = "\\d+$", value = TRUE)
#create vector of site covariate variables
sitecovcols <- c("cam_site_id", "yearID", "year", "X", "Y", scalevars)
#make dataframe of site covariates and add site ID
sitecovs <- ModelingDF2%>%select(all_of(sitecovcols))%>%distinct()%>%arrange(year, cam_site_id)%>%
  group_by(year)%>%mutate(siteID = row_number())%>%ungroup()%>%relocate(siteID)
#make a siteID x year matrix of camera site IDs
camsiteskey <- sitecovs%>%select(yearID, siteID, cam_site_id)%>%mutate(cam_site_id2 = as.numeric(as.factor(cam_site_id)))
camsites <- camsiteskey%>%select(-cam_site_id)%>%pivot_wider(values_from = cam_site_id2, names_from = yearID)%>%select(-siteID)%>%as.matrix()
#make a siteID x year matrix of number of occasions
nsurveys <- ModelingDF2%>%group_by(year)%>%mutate(siteID=as.numeric(factor(cam_site_id)))%>%
  ungroup()%>%group_by(cam_site_id,siteID,year,yearID)%>%summarise(nsurveys=n())%>%arrange(yearID, siteID)%>%ungroup()
nsurveys2 <- nsurveys%>%pivot_wider(id_cols=siteID,names_from=yearID, values_from=nsurveys)%>%select(-siteID)%>%
  as.matrix()%>%unname() 

# make a siteID x year matrix of scaled year
yr <- sitecovs%>%select(yearID, year, siteID)%>%
  mutate(year=year-2019)%>% # Reformatted (0, 1, 2, 3, 4).
  mutate(year = as.numeric(scale(year)))%>%
  pivot_wider(names_from = yearID, values_from = year)%>%
  select(-siteID)%>%
  as.matrix()%>%
  unname()


#make look up table for zones for each camsite
Zones <- sitecovs%>%mutate(X=X*attr(X, 'scaled:scale') + attr(X, 'scaled:center'), 
                           Y=Y*attr(Y, 'scaled:scale') + attr(Y, 'scaled:center'))%>%
  st_as_sf(., coords=c("X", "Y"), crs=3071)%>%
  st_join(.,st_transform(st_make_valid(get_spatial_data("bear_zones")), 3071))%>%
  rename(zone=bear_mgmt_zone_id)%>%st_drop_geometry()

zone <- Zones%>%select(yearID, siteID, zone)%>%
  mutate(zone = as.numeric(as.factor(zone)))%>%
  pivot_wider(names_from = yearID, values_from = zone)%>%
  select(-siteID)%>%
  as.matrix()%>%
  unname()

# make a siteID x year matrix of scaled X and Y coordinate
X <- sitecovs%>%select(year, siteID, X)%>%
  pivot_wider(names_from = year, values_from = X)%>%
  select(-siteID)%>%
  as.matrix()%>%
  unname()
Y <- sitecovs%>%select(year, siteID, Y)%>%
  pivot_wider(names_from = year, values_from = Y)%>%
  select(-siteID)%>%
  as.matrix()%>%
  unname()

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
write.csv(bearrangeknots2, "./bearrangeknots.csv")
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

sp.Z2 <- cbind.data.frame(coordsdf$cam_site_id, sp.Z)
colnames(sp.Z2)[1] <- "cam_site_id"
sp.Z3 <- left_join(camsiteskey, sp.Z2)
sp.Zarray <- array(NA, dim=c(maxsites,no.years,nrow(bearrangeknots2)))
for( t in 1:no.years ) {
  for( k in 1:nrow(bearrangeknots2)) {
    for( i in 1:nsite[t]){
      sp.Zarray[ i, t, k] <- as.numeric(sp.Z3[ c( sp.Z3$siteID == i & sp.Z3$yearID == t), 4+k])
    }
  }
}

#multi-scale state covariates
#Developed
Developed <- sitecovs%>%select(yearID, siteID, matches("Developed"))%>%
  pivot_longer(cols = matches("Developed"), names_pattern="(\\d+$)", names_to = "scale")

scales <- unique(Developed$scale)
Developed_array <- array(NA, dim=c(maxsites,no.years,5)) #sites, years, scales
for( t in 1:no.years ) {
  for( s in 1:5) {
    for( i in 1:nsite[t]){
      Developed_array[ i, t, s] <- Developed[ c( Developed$siteID == i & Developed$scale == scales[s] & Developed$yearID == t), "value"]$value
    }
  }
}
#reduce scales to 500m, 1000m, 5000m
Developed_array <- Developed_array[,,c(2,3,4,5)]

#Disturbance
Dist <- sitecovs%>%select(yearID, siteID, matches("Dist"))%>%
  pivot_longer(cols = matches("Dist"), names_pattern="(\\d+$)", names_to = "scale")

Dist_array <- array(NA, dim=c(maxsites,no.years,5)) #sites, years, scales
for( t in 1:no.years) {
  for( s in 1:5) {
    for( i in 1:nsite[t]){
      Dist_array[ i, t, s] <- Dist[ c( Dist$siteID == i & Dist$scale == scales[s] & Dist$yearID == t), "value"]$value
    }
  }
}
Dist_array <- Dist_array[,,c(2,3,4,5)]

#Proportion of Forest
Forest <- sitecovs%>%select(yearID, siteID, matches("Forest"))%>%
  pivot_longer(cols = matches("Forest"), names_pattern="(\\d+$)", names_to = "scale")

Forest_array <- array(NA, dim=c(maxsites,no.years,5)) #sites, years, scales
for( t in 1:no.years ) {
  for( s in 1:5) {
    for( i in 1:nsite[t]){
      Forest_array[ i, t, s] <- Forest[ c( Forest$siteID == i & Forest$scale == scales[s] & Forest$yearID == t), "value"]$value
    }
  }
}
Forest_array <- Forest_array[,,c(2,3,4,5)]

#Corn
Corn <- sitecovs%>%select(yearID, siteID, matches("Corn"))%>%
  pivot_longer(cols = matches("Corn"), names_pattern="(\\d+$)", names_to = "scale")

Corn_array <- array(NA, dim=c(maxsites,no.years,5)) #sites, years, scales
for( t in 1:no.years ) {
  for( s in 1:5) {
    for( i in 1:nsite[t]){
      Corn_array[ i, t, s] <- Corn[ c( Corn$siteID == i & Corn$scale == scales[s] & Corn$yearID == t), "value"]$value
    }
  }
}
Corn_array <- Corn_array[,,c(2,3,4,5)]

#detection covariates
detcovcols <- c("cam_site_id", "yearID", "year", "occ","camera_version", "meanEVI", "days_active")
detcovs <- ModelingDF2%>%select(all_of(detcovcols))%>%arrange(year, cam_site_id)%>%
  group_by(year)%>%mutate(siteID = as.numeric(factor(cam_site_id)))%>%ungroup()%>%relocate(siteID)


EVI <- readRDS("./EVIsiteyearoccSummer2019-2024.rds")%>%arrange(year, cam_site_id)%>%mutate(meanEVI=scale(meanEVI))
EVI2 <- left_join(Dethistlong2, EVI)

EVI_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      EVI_array[ i, j, t] <- EVI2[ c( EVI2$siteID == i & EVI2$occre == j & EVI2$yearID == t), "meanEVI"]$meanEVI
    }
  }
}
# knots:
EVI.num.knots = 5
EVI.knots<-quantile(unique(EVI2$meanEVI),
                seq(0,1,length=(EVI.num.knots+2))[-c(1,(EVI.num.knots+2))])

# construct design matrix Z:
EVI.Z_K <- (abs(outer(as.numeric(EVI2$meanEVI),EVI.knots,"-")))^3

EVI.OMEGA_all <- (abs(outer(EVI.knots,EVI.knots,"-")))^3

EVI.svd.OMEGA_all <- svd(EVI.OMEGA_all)

EVI.sqrt.OMEGA_all <- t(EVI.svd.OMEGA_all$v %*%
                      (t(EVI.svd.OMEGA_all$u)*sqrt(EVI.svd.OMEGA_all$d)))

EVI.Z <- t(solve(EVI.sqrt.OMEGA_all,t(EVI.Z_K)))
meanZ <- mean(EVI.Z)
sdZ <- sd(EVI.Z)
EVI.Z <- (EVI.Z-meanZ)/sdZ

EVI3 <- cbind(EVI2, EVI.Z)
EVIZ_array <- array(NA, dim=c(maxsites,no.occs,no.years, EVI.num.knots)) #sites, occasions, years
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      for(k in 1:EVI.num.knots){
        EVIZ_array[ i, j, t, k] <- EVI3[ c( EVI3$siteID == i & EVI3$occre == j & EVI3$yearID == t), 8+k]
      }
    }
  }
}




days_active <- detcovs%>%dplyr::select(yearID, siteID, occ, days_active)%>%
  left_join(Dethistlong2, .)

daysactive_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      daysactive_array[ i, j, t] <- days_active[ c( days_active$siteID == i & days_active$occre == j & days_active$yearID == t), "days_active"]$days_active
    }
  }
}


#Camera version  
unique(ModelingDF2$camera_version)
cam_version <- detcovs%>%select(yearID, siteID, occ, camera_version)%>%
  mutate(cam_versionID=as.numeric(as.factor(camera_version)),
         cam_versionID2=as.numeric(as.factor(if_else(camera_version == "V4", 2, 1))))%>%
  left_join(Dethistlong2, .)



camversion_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years. Do i need to remo
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      camversion_array[ i, j, t] <- cam_version[ c( cam_version$siteID == i & cam_version$occre == j & cam_version$yearID == t), "cam_versionID2"]$cam_versionID2
    }
  }
}

#what about scaling?
occ_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      occ_array[ i, j, t] <- Dethistlong2[ c( Dethistlong2$siteID == i & Dethistlong2$occre == j & Dethistlong2$yearID == t), "occ"]$occ
    }
  }
}
occscale_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years
Dethistlong3 <- Dethistlong2%>%mutate(occ = scale(occ))
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      occscale_array[ i, j, t] <- Dethistlong3[ c( Dethistlong3$siteID == i & Dethistlong3$occre == j & Dethistlong3$yearID == t), "occ"]$occ
    }
  }
}


############################################
############# R-N Model ####################
############################################
# Create constants for model from array.
constants <- list(
  nyear = dim(Bear_All)[[3]],
  nsite = nsite,
  camsites=camsites, #need this indexed by [i,t] and [i,t,k]? Or maybe not
  ncams=length(unique(ModelingDF2$cam_site_id)),
  nversions=length(unique(cam_version$cam_versionID2))-1,
  nsurveys=nsurveys2,
  camversion=camversion_array,
  sp.nknots=50,
  EVI.nknots=5
)
# Zone=zone,
# nZones=length(unique(Zones$bear_mgmt_unit_id))
# Bundle data (counts and covariates).
data <- list(
  y = Bear_All,
  Year = yr,
  #X=X,
  #Y=Y,
  Dev=Developed_array,
  Dist=Dist_array,
  Forest=Forest_array,
  Corn=Corn_array,
  EVI=EVI_array,
  daysactive=daysactive_array,
  catprobs = c(0.25, 0.25, 0.25, 0.25),
  sp.Z= sp.Zarray,
  EVI.Z= EVIZ_array
)

RNcode <- nimbleCode({
  
  # .............................................................
  # PRIORS
  # .............................................................
  
  ## Priors for SDM ##
  
  # CAR prior for spatial random effect SRE
  # s[1:ncell] ~ dcar_normal(adj[1:neigh], weights[1:neigh], num[1:ncell], tau)
  # precision of CAR prior
  # tau ~ dgamma(1, 1)
  
  # regression coefficient for ruffed grouse priority zone
  b_Yr ~ dnorm(0, 2)
  # regression coefficient for Human Footprint Index
  b_Dev ~ dnorm(0, 2)
  # regression coefficient for Dist
  b_Dist ~ dnorm(0, 2)
  # regression coefficient for Forest
  b_Forest ~ dnorm(0, 2)
  # regression coefficient for Corn
  b_Corn ~ dnorm(0, 2)
  
 
  
  
  for (k in 1:sp.nknots) {
    spat.spline.b[k] ~ dnorm(0,sigma.spat.spline.b)
  }
  for (k in 1:EVI.nknots) {
    b[k] ~ dnorm(0,sigma.EVI)
  }
  
  # spline random effect priors
  sigma.spat.spline.b~dunif(0,100)
  
  # spline random effect priors
  sigma.EVI~dunif(0,100)
  
  ## Priors for detection parameters ##
  # coefficient for occasion
  a_EVI ~ dlogis(0, 1)
  a_daysactive ~ dlogis(0,1)
  #a_occ ~ dlogis(0,1)
  for(v in 1:nversions){
    a_version[v] ~ dlogis(0,1)
  }
  
  # # camera site random effect for detection
  # for(c in 1:ncams){
  #   eps_p[c] ~ dnorm(0, sd_p)
  # }
  # # hyperprior for abundance random effect
  # sd_p ~ dgamma(1, 2)
  
  ## Priors for scales  ##
  abundance_scale[1] ~ dcat(catprobs[1:4])
  abundance_scale[2] ~ dcat(catprobs[1:4])
  abundance_scale[3] ~ dcat(catprobs[1:4])
  abundance_scale[4] ~ dcat(catprobs[1:4])
  
  
  for( t in 1:nyear ) { #loop over site then year?
    #state model
    # Loop through only the sites that are surveyed in a given year. 
    for( i in 1:nsite[t] ){
      N[i, t] ~ dpois( lambda[ i, t ] )
      log(lambda[i, t]) <-  b_Yr*Year[i, t] +  #make this a factor? b_ZoneYr[Zone[i,t]]*Year[i,t] +
        #scaled parameters
        b_Dev*Dev[i,t,abundance_scale[1]] + b_Dist*Dist[i,t,abundance_scale[2]] + b_Forest*Forest[i,t,abundance_scale[3]] + b_Corn*Corn[i,t,abundance_scale[4]] +
        #cam_site random effect
        s[i,t] #s[i,t] - spatial random effect
      
      s[i,t] <- inprod(spat.spline.b[1:sp.nknots],sp.Z[i,t,1:sp.nknots])
      
      #detection model  
      for(k in 1:nsurveys[i,t]){
        muy[i, k, t] <- 1 - pow(1-rho[i,k,t], N[i, t]) #
        logit(rho[i, k, t]) <- a_version[camversion[i,k,t]] + a_daysactive*daysactive[i,k,t] + 
          a_EVI * EVI[i,k,t] + EVI.spline[i,k,t] #+ eps_p[camsites[i,t]]
        y[i, k, t] ~ dbern(muy[i, k, t])
        
        EVI.spline[i,k, t] <- inprod(b[1:EVI.nknots],EVI.Z[i,k,t,1:EVI.nknots])
      }
    }
  }
  
  #PREDICTION
  

})




# function to provide random initial values for parameters
inits <- function() {
  base::list(N = matrix(data = rep(1, max(constants$nsite)*constants$nyear),
                        nrow = max(constants$nsite),
                        ncol = constants$nyear),
             b_Yr = runif(1, -1, 1),
             #b_ZoneYr = runif(constants$nZones, -1, 1),
             b_Dev = runif(1, -1, 1),
             b_Dist = runif(1, -1, 1),
             b_Forest= runif(1, -1, 1),
             b_Corn= runif(1, -1, 1),
             a_version = runif(constants$nversions, -1, 1),
             a_daysactive = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             spat.spline.b=rep(1,constants$sp.nknots),
             sigma.spat.spline.b=1,
             b=rep(1,constants$EVI.nknots),
             sigma.EVI=1,
             #a_occ = runif(1, -1, 1),
             #eps_N = rnorm(constants$ncams, 0, 2),
             # eps_p = rnorm(constants$ncams, 0, 2),
             # sd_p = runif(1, 0, 2),
             abundance_scale=rcat(4, c(0.25,0.25,0.25,0.25)),
             rho = array(data = runif(length(Bear_All), 0, 1),
                         dim=c(maxsites,no.occs,no.years))
  )
}

# parameters to monitor
keepers <- c('b_Dev', "b_Yr",
             "b_Dist", "b_Forest", "b_Corn",
             "spat.spline.b", "b",
             "a_version", "a_daysactive", "a_EVI", 
             "abundance_scale") #"b_X","b_Y",

# Will have to run chains for much longer (~40,000 iterations) to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 5000 # number of initial MCMC iterations to discard
ni <- 75000 # total number  of iterations

# .......................................................................
# RUN MODEL
# .......................................................................


# create model
model <- nimble::nimbleModel(code = RNcode, 
                             data = data, 
                             constants = constants, 
                             inits = inits())

# check to see if everything is initialized
model$initializeInfo()
model$lambda
which(is.na(model$lambda))
model$HFI
which(is.na(model$HFI))
which(is.na(model$Forest))
model$calculate()
model$calculate("abundance_scale[1]")

# compile the model
c_model <- nimble::compileNimble(model)

model_conf <- nimble::configureMCMC(model)# enableWAIC = TRUE

model_conf$addMonitors(keepers)

#reversible jump MCMC
# configureRJ(RN_code,
#             targetNodes = c("HFI", "Dist", "Forest"),
#             indicatorNodes = 'abundance_scale',
#             control = list(mean = 0, scale = .2)) #no idea what this does

#Stuber 2017 can i specify reversible jump among multiple scales and one of the scales will always be in the model instead of in and out indicator variable
#model_conf$printSamplers(c("abundance_scale[1]", "HFI"))

model_mcmc <- nimble::buildMCMC(model_conf)

c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)
c_model_mcmc$my_initializeModel


test <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = 0, 
                           niter = 100, 
                           nchains = nc, 
                           inits=inits())


samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc, 
                           thin= 5,
                           inits=inits())
samples2 <- append(samples, list("formula"= "log(lambda[i, t]) <-  Year  + Devsc + Distsc + Forestsc + Cornsc + spatialspline
                                 p <- version + daysactive + EVI + EVIspline"))
saveRDS(samples2, "./RNsamplesFullModelBearRange4.rds")

samples <- readRDS("./RNsamples50000Summer.rds")

MCMCsummary(samples,params = "b_Lat", round = 2)
MCMCsummary(samples,params = "b_HFI", round = 2)
MCMCsummary(samples,params = "abundance_scale", round = 2)
MCMCsummary(samples,params = "a_version", round = 2) 
MCMCsummary(samples[[1]],params = c("lambda[1000, 1]", "lambda[1255, 1]", "lambda[1262, 1]", "lambda[1066, 2]", "lambda[1067, 2]"), round = 2, ISB = FALSE)


PR <- rnorm(15000, 0, 2)
MCMCtrace(samplesSummer[1:3], 
          params = c( 'b_Dev', 'b_Forest', "b_Corn", "b_Dist", "b_Yr"),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c( "b_ZoneYr"),
          ISB = TRUE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c('spat.spline.b[1]', 'spat.spline.b[2]', 'spat.spline.b[3]', 'spat.spline.b[4]', 'spat.spline.b[5]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c('abundance_scale'),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c("a_version", "a_daysactive", "a_EVI"),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(test, 
          params = c("lambda[665, 1]", "lambda[719, 1]", "lambda[720, 1]", "lambda[769, 1]", "lambda[564, 2]"), #still has nodes for "lambda[1243, 1]" and such
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCsummary(samples, 
          params = c("spat.spline.b"), #still has nodes for "lambda[1243, 1]" and such
          ISB = TRUE,
          round=2)
lambdameans <- MCMCpstr( samples, params = c("lambda"), func=mean, type="chain")[[1]]

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

rhoplotdf <- cbind(pred.data2, rho.mean)
rhoplotdf$EVI <- rhoplotdf$EVI * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
ggplot(rhoplotdf, aes(x=EVI, y=rho.mean)) + geom_point() #+ geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5))



runCrossValidate
############################################################################################
####                           pop by zone                                           #######
############################################################################################
# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples, coda::mcmc))
samplesdf <- do.call(rbind, samples_mcmc)
lambdas <- samplesdf%>%as.data.frame()%>%select(starts_with("lambda"))%>%mutate("iter"=1:nrow(.))%>%relocate(iter)%>%
  pivot_longer(starts_with("lambda"), names_to = "param")%>%
  mutate("siteID"=str_extract(string = param, pattern = "^.*\\[(\\d+),", group=1), "yearID"=str_extract(string = param, pattern = ", (\\d{1}).*$", group=1))
View(head(lambdas, 50))
notthese <- paste0("lambda[", ,"]")
Zones$siteID <- as.character(Zones$siteID)
Zones$yearID <- as.character(Zones$yearID)
lambdaszone <- left_join(lambdas, Zones)
View(head(lambdaszone, 50))
popbyiterzone <- lambdaszone%>%group_by(yearID, bear_mgmt_zone_id, iter)%>%summarise(popbyzone=sum(value))
nsite2 <- data.frame("yearID"=as.character(1:6), "nsite"=nsite, nsiteprop=nsite/1243)
popbyyearzone <- popbyiterzone%>%group_by(yearID, bear_mgmt_zone_id)%>%left_join(., nsite2)%>%
  summarise(popbyzone2=mean(popbyzone)/unique(nsiteprop), sd=sd(popbyzone), lowci=quantile(popbyzone, .025)/unique(nsiteprop), highci=quantile(popbyzone, .975)/unique(nsiteprop), nsite=unique(nsite), nsiteprop2=unique(nsiteprop))

#
ggplot(popbyyearzone, aes(x=as.numeric(yearID), y=popbyzone2, color=bear_mgmt_zone_id)) + geom_point() +
  geom_line(aes(group = bear_mgmt_zone_id)) + geom_ribbon(aes(ymax = highci, ymin = lowci, fill = bear_mgmt_zone_id), alpha=0.3)


##############################  what were inits for scale selection #########################################
samples[["chain1"]][1,6:8]
samples[["chain2"]][1,6:8]
samples[["chain3"]][1,6:8]
############                      scrap           ################################################
badinits <- c("N[138, 1]", "N[665, 1]", "N[769, 1]", "N[564, 2]", "N[633, 3]", "N[734, 4]",
              "N[735, 4]", "N[736, 4]", "N[662, 5]", "N[715, 5]", "N[683, 6]", "N[756, 6]", "N[758, 6]",
               "N[817, 6]", "N[818, 6]")

badinits <- data.frame(siteID=str_extract(string = badinits, pattern = "(N\\[)(\\d+)", group = 2),
                       yearID=str_extract(string = badinits, pattern = "(N\\[\\d+, )(\\d+)", group = 2))
badinits2 <- camsiteskey%>%filter(paste0(siteID, yearID) %in% paste0(badinits$siteID, badinits$yearID))
camsitessf <- ModelingDF1%>%st_as_sf(., coords=c("X", "Y"), crs=3071)%>%distinct(cam_site_id, geometry)
#3/4 bad initis down in kenosha county, other one in superior, WI and has lambda of 0.58
badinits3 <- camsitessf%>%filter(cam_site_id %in% badinits2$cam_site_id)
notbadinits <- camsitessf%>%filter(!(cam_site_id %in% badinits2$cam_site_id))
badinitssitecovs <- sitecovs%>%filter(paste0(siteID, yearID) %in% paste0(badinits$siteID, badinits$yearID))
bearknotssf <- st_as_sf(bearrangeknots, coords=c("X", "Y"), crs=3071)
leaflet() %>% 
  # addProviderTiles("OpenStreetMap.Mapnik") %>%
  addTiles() %>%
  addCircleMarkers(data=st_transform(badinits3, 4326), popup=badinits3$cam_site_id, fillColor = "black", fillOpacity = 1,   stroke=F, radius=5)%>%
  addCircleMarkers(data=st_transform(bearknotssf, 4326), fillColor = "yellow", fillOpacity = 1,   stroke=F, radius=5)%>%
  addPolygons(data=st_transform(counties, crs=4326))
Zsbadinits <- matrix(NA, nrow=15, ncol=50)
for(i in 1:nrow(badinits)){
  print(as.numeric(badinits$siteID[i]),as.numeric(badinits$yearID[i]))
Zsbadinits[i,] <- sp.Zarray[as.numeric(badinits$siteID[i]),as.numeric(badinits$yearID[i]),]
}
spatsplinebadinits <- rowSums(Zsbadinits)
lambdabadinits <- spatsplinebadinits + intslist$b_Yr*badinitssitecovs$year + intslist$b_X*badinitssitecovs$X + intslist$b_Y*badinitssitecovs$Y +
  intslist$b_Dev*badinitssitecovs$Developed_5000 + intslist$b_Dist*badinitssitecovs$Dist_100 + intslist$b_Forest*badinitssitecovs$Forest_100 +
  intslist$b_Corn*badinitssitecovs$Corn100

ggplot(HFI, aes(x=value)) + facet_wrap(~scale) + geom_histogram()
ggplot(Forest, aes(x=value)) + facet_wrap(~scale) + geom_histogram()
ggplot(Dist, aes(x=value)) + facet_wrap(~scale) + geom_histogram()
ggplot(days_active, aes(x=days_active)) + geom_histogram()


NApopiterzone <- lambdaszone%>%filter(is.na(bear_mgmt_zone_id))


#make a list of specific occasions active for each camsite-year to avoid using NA surveys

nsurveys2 <- nsurveys%>%pivot_wider(id_cols=siteID,names_from=yearID, values_from=nsurveys)%>%select(-siteID)%>%
  as.matrix()%>%unname()
occasions <- ModelingDF2%>%group_by(year)%>%mutate(siteID=as.numeric(factor(cam_site_id)))%>%ungroup()%>%select(yearID, siteID, occ)%>%split(., .$yearID)
occasions2 <- lapply(occasions, function (x) split(x, x$siteID))
occasions3 <- lapply(occasions2, function (y) lapply(y, function (x) x$occ))


testing <- Dethistlong[177:187,]
testing%>%arrange(is.na(det), occ)
testing2 <- Dethistlong%>%group_by(yearID,siteID)%>%arrange(yearID, siteID, is.na(det), occ)%>%mutate(occre=row_number())%>%ungroup()
Test_All <- array(NA, dim=c(1243,11,6)) #sites, occasions, years

nsite <- sapply(dethist,nrow)

# In the "all" array (above), loop through each year of data (yearID) "t", each weekly
# count (rep) "j", each site (siteID) "i" (number of sites depends on year; addressed
# above and incorporated below), and assign counts.
for( t in 1:6 ) {
  for( j in 1:11) {
    for( i in 1:nsite[t]){
      Test_All[ i, j, t] <- testing2[ c( testing2$siteID == i & testing2$occre == j & testing2$yearID == t), "det"]$det
    }
  }
}

#what about scaling?
occ_array <- array(NA, dim=c(1243,11,6)) #sites, occasions, years
for( t in 1:6 ) {
  for( j in 1:11) {
    for( i in 1:nsite[t]){
      occscale_array[ i, j, t] <- testing2[ c( testing2$siteID == i & testing2$occre == j & testing2$yearID == t), "occ"]$occ
    }
  }
}
occscale_array <- array(NA, dim=c(1243,11,6)) #sites, occasions, years
testing2$occ <- scale(testing2$occ)
for( t in 1:6 ) {
  for( j in 1:11) {
    for( i in 1:nsite[t]){
      occscale_array[ i, j, t] <- testing2[ c( testing2$siteID == i & testing2$occre == j & testing2$yearID == t), "occ"]$occ
    }
  }
}

DAs <- BEARCPUE3[,c(1:3,8,11)]%>%st_drop_geometry()%>%rename(yearID=season)%>%
  group_by(year)%>%mutate(siteID=as.numeric(as.factor(cam_site_id)))%>%ungroup()%>%relocate(siteID)%>%
  arrange(year, cam_site_id)

EVI <- EVI2%>%mutate(yearID=year-2018)%>%group_by(year)%>%
  mutate(siteID = as.numeric(factor(cam_site_id)))#%>%
pivot_wider(names_from = occ, values_from = meanEVI)#%>%
select(-siteID)%>%
  as.matrix()%>%
  unname()

ifelse(apply(cam_version[,3:13], 1, unique))

#detcovs2 <- detcovs%>%group_by(yearID,siteID)%>%arrange(yearID, siteID, is.na(det), occ)%>%mutate(occre=row_number())%>%ungroup()


NANZmat <- coordsdf[which(is.na(Z[,1])),]
compare2knots <- cbind(NANZmat, as.data.frame(knots$design))


EVI <- readRDS("./EVIsiteyearoccSummer2019-2024.rds")%>%arrange(year, cam_site_id)%>%mutate(meanEVI=scale(meanEVI))
EVI2 <- left_join(Dethistlong2, EVI)
bsEVI <- bSpline(ModelingDF$meanEVI)

summary(fm1 <- mgcv::gam(BEAR_ADULT_AMT ~ s(meanEVI), data = ModelingDF, family = "poisson"))
## example of safe prediction
plot(ModelingDF[,c(14,8)], xlab = "meanEVI", ylab = "BEAR_ADULT_AMT")
EVI <- seq(2000, 8000, length.out = 100)
lines(EVI, predict(fm1, data.frame(meanEVI = EVI)), col="red")

inprod(beta[], X[i,])

plot(st_geometry(bearrange2))
plot(st_geometry(st_as_sf(test2, coords=c("X", "Y"), crs=3071)), col="red", add=TRUE)
plot(st_geometry(st_as_sf(bearrangeknots2, coords=c("X", "Y"), crs=3071)), col="blue", add=TRUE)


# make a pretty table with example of the input data to the model
Dethist2019 <- Dethistlong%>%mutate(occ=occ-3)%>%filter(year == 2019)%>%pivot_wider(names_from = occ, values_from = det)
Dethistdemo <- Dethist2019[c(45,49,50,51),1:11]
Dethistdemo2 <- Dethistdemo[, -c(3,4)]
colnames(Dethistdemo2) <- c("Camera Site", "Year", paste0("Week", 1:7))
library(gt)
Dethistdemo2%>%gt()%>%   tab_style(
  style = list(
    cell_text(weight = "bold", color = "white"),
    cell_fill(color = "#0072B2")
  ),
  locations = cells_column_labels(everything())
)%>%
  tab_style(
    style = list(
      cell_fill(color = "grey") # A light yellow color
    ),
    locations = cells_body(
      columns = c(cam_site_id, year)
    )
  )%>%
  cols_label(
    cam_site_id = "Camera Site",
    year = "Year",
    `1`= "Week1",
    `2`= "Week2",
    `3`= "Week3",
    `4`= "Week4",
    `5`= "Week5",
    `6`= "Week6",
    `7`= "Week7"
  )
