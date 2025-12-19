library(dplyr)
library(sf)
library(data.table)
library(tidyr)
library(nimble)
library(MCMCvis)
library(sswids)
library(stringr)

#make loo up table for zones for each camsite
Zones <- readRDS("./ModelingDFSummer.rds")%>%ungroup()%>%select(cam_site_id, season)%>%distinct()%>%st_transform(., 3071)%>%
  group_by(season)%>%mutate(siteID=row_number())%>%st_join(., st_transform(get_spatial_data("bear_zones"), 3071))%>%arrange(season, cam_site_id)%>%
  select(cam_site_id, yearID=season, siteID, bear_mgmt_zone_id)%>%mutate(cam_site_id_num = as.numeric(as.factor(cam_site_id)))%>%
  st_drop_geometry()


ModelingDF <- readRDS("./ModelingDFSummer.rds")%>%st_drop_geometry()%>%ungroup()%>%filter(occ > 3 & occ < 18) #5/21-8/26 this data frame doesn't have NAs for cam site-year-occs that dont have effort its just missing those rows
ModelingDF <- ModelingDF[-which(ModelingDF$camera_version %in% c("V2,V4","V2,V3")),] #remove occasions which have multiple camera versions
only1occ <- ModelingDF%>%group_by(cam_site_id, year)%>%summarise(N=n())%>%filter(N == 1)
ModelingDF <- ModelingDF%>%filter(
  !(paste0(cam_site_id, year) %in% paste0(only1occ$cam_site_id,only1occ$year)))
#split up data frame by year and reshape to wide format for detections histories for each year
dethist <- split(ModelingDF, ModelingDF$year)
dethist <- lapply(dethist, function (x) tidyr::pivot_wider(x, id_cols = c(cam_site_id, year), names_from = occ, values_from = BEAR_ADULT_AMT, names_sort = TRUE)%>%ungroup())
dethistall <- rbindlist(dethist)%>%group_by(year)%>%mutate(siteID=row_number())%>%ungroup()
Dethistlong <- dethistall%>%pivot_longer(cols=3:16, names_to="occ", values_to="det")%>%
  mutate(yearID=year-2018, det=ifelse(det > 0,1,0), occ=as.numeric(occ))%>%
  arrange(year, cam_site_id, occ)%>%ungroup()
Dethistlong2 <- Dethistlong%>%group_by(yearID,siteID)%>%arrange(yearID, siteID, is.na(det), occ)%>%mutate(occre=row_number())%>%ungroup()
nsite <- sapply(dethist,nrow)

maxsites <- max(nsite)
no.occs <- length(unique(ModelingDF$occ))
no.years <- length(unique(ModelingDF$year))
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
ModelingDF2 <- ModelingDF%>%mutate(across(c(10, 12, 14:35), scale))%>%rename(yearID=season)
#create vector of selection scale variables
scalevars <- grep(x = colnames(ModelingDF2), pattern = "_\\d+$", value = TRUE)
#create vector of site covariate variables
sitecovcols <- c("cam_site_id", "yearID", "year", "Lat", scalevars)
#make dataframe of site covariates and add site ID
sitecovs <- ModelingDF2%>%select(all_of(sitecovcols))%>%distinct()%>%arrange(year, cam_site_id)%>%
  group_by(year)%>%mutate(siteID = row_number())%>%ungroup()%>%relocate(siteID)
#make a siteID x year matrix of camera site IDs
camsites <- sitecovs%>%select(yearID, siteID, cam_site_id)%>%mutate(cam_site_id = as.numeric(as.factor(cam_site_id)))%>%
  pivot_wider(values_from = cam_site_id, names_from = yearID)%>%select(-siteID)%>%as.matrix()
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
# make a siteID x year matrix of scaled Latitude
Lat <- sitecovs%>%select(year, siteID, Lat)%>%
  pivot_wider(names_from = year, values_from = Lat)%>%
  select(-siteID)%>%
  as.matrix()%>%
  unname()

#spatial spline
coordsdf <- readRDS("./ModelingDF.rds")%>%cbind(., st_coordinates(.))%>%select("cam_site_id", "X", "Y")%>%distinct()
coordsmatrix <- coordsdf%>%st_drop_geometry()%>%select("X", "Y")%>%as.matrix()
knots <- read_csv("knots.csv")
# scale coordinates 
mean_x <- mean(coordsdf$X)
sd_x <- sd(coordsdf$X)
mean_y <- mean(coordsdf$Y)
sd_y <- sd(coordsdf$Y)

knots <- knots %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y) %>%
  dplyr::select(X.scale,Y.scale)

dat <- coordsdf %>%
  mutate(X.scale = (X-mean_x)/sd_x,
         Y.scale = (Y-mean_y)/sd_y)

# get matrix ready for spatial smoothing
knots.dist <- dist(knots,"euclidean",diag=T,upper=T)
omega_all = knots.dist^2*log(knots.dist) # basis
svd.omega_all <- svd(omega_all)
sqrt.omega_all <- t(svd.omega_all$v %*%
                      (t(svd.omega_all$u)*sqrt(svd.omega_all$d)))

# now for each spline
cov.dist_all = fields::rdist(x1=cbind(dat$X.scale,dat$Y.scale),x2=knots)
Z_K = cov.dist_all^2*log(cov.dist_all) # basis
Z <- t(solve(sqrt.omega_all,t(Z_K)))
Z <- (Z - meanZ)/sdZ #for prediction?


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

#detection covariates
detcovcols <- c("cam_site_id", "yearID", "year", "occ","camera_version", "meanEVI", "days_active")
detcovs <- ModelingDF2%>%select(all_of(detcovcols))%>%arrange(year, cam_site_id)%>%
  group_by(year)%>%mutate(siteID = as.numeric(factor(cam_site_id)))%>%ungroup()%>%relocate(siteID)


EVI <- readRDS("./EVIsiteyearoccSummer2019-2024.rds")%>%arrange(year, cam_site_id)%>%mutate(meanEVI=scale(meanEVI))
EVI2 <- left_join(Dethistlong2, EVI)
bsEVI <- bSpline(EVI$meanEVI)
inprod(beta[], X[i,])


EVI_array <- array(NA, dim=c(maxsites,no.occs,no.years)) #sites, occasions, years
for( t in 1:no.years ) {
  for( j in 1:no.occs) {
    for( i in 1:nsite[t]){
      EVI_array[ i, j, t] <- EVI2[ c( EVI2$siteID == i & EVI2$occre == j & EVI2$yearID == t), "meanEVI"]$meanEVI
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
  nknots=50
)

# Bundle data (counts and covariates).
data <- list(
  y = Bear_All,
  Year = yr,
  Latitude=Lat,
  HFI=HFI_array,
  Dist=Dist_array,
  Forest=Forest_array,
  EVI=EVI_array,
  EVI2=EVI2_array,
  daysactive=daysactive_array,
  occ=occscale_array,
  catprobs = c(0.2, 0.2, 0.2, 0.2, 0.2),
  Z=Z
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
  b_HFI ~ dnorm(0, 2)
  # regression coefficient for Lat
  b_Lat ~ dnorm(0, 2)
  # regression coefficient for Dist
  b_Dist ~ dnorm(0, 2)
  # regression coefficient for Forest
  b_Forest ~ dnorm(0, 2)
  
 
  
  
  for (k in 1:nknots) {
    spat.spline[k] ~ dnorm(0,sigma.spat.spline)
  }
  
  # spline random effect priors
  sigma.spat.spline~dunif(0,100)
  
  ## Priors for detection parameters ##
  # coefficient for occasion
  a_EVI ~ dlogis(0, 1)
  a_EVI2 ~ dlogis(0, 1)
  a_daysactive ~ dlogis(0,1)
  #a_occ ~ dlogis(0,1)
  for(v in 1:nversions){
    a_version[v] ~ dlogis(0,1)
  }
  
  
  ## Priors for scales  ##
  
  abundance_scale[1] ~ dcat(catprobs[1:5])
  abundance_scale[2] ~ dcat(catprobs[1:5])
  abundance_scale[3] ~ dcat(catprobs[1:5])
  
  
  for( t in 1:nyear ) { #loop over site then year?
    #state model
    # Loop through only the sites that are surveyed in a given year. 
    for( i in 1:nsite[t] ){
      N[i, t] ~ dpois( lambda[ i, t ] )
      log(lambda[i, t]) <-  b_Yr*Year[i, t] +#make this a factor? 
        b_Lat*Latitude[i, t] + 
        #scaled parameters
        b_HFI*HFI[i,t,abundance_scale[1]] + b_Dist*Dist[i,t,abundance_scale[2]] + b_Forest*Forest[i,t,abundance_scale[3]] + 
        #cam_site random effect
        s[camsites[i,t]] #s[i,t] - spatial random effect
      
      s[camsites[i,t]] <- inprod(spat.spline[1:nknots],Z[camsites[i,t],1:nknots])
      
      #detection model  
      for(k in 1:nsurveys[i,t]){
        muy[i, k, t] <- 1 - pow(1-rho[i,k,t], N[i, t]) #
        logit(rho[i, k, t]) <- a_version[camversion[i,k,t]] + a_daysactive*daysactive[i,k,t] + a_EVI*EVI[i,k,t] + a_EVI2*EVI[i,k,t]*EVI[i,k,t] #EVI and occ probably correlated + a_occ*occ[i,k,t] + eps_p[camsites[i,t]]
        y[i, k, t] ~ dbern(muy[i, k, t])
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
             b_Lat = runif(1, -1, 1),
             b_HFI = runif(1, -1, 1),
             b_Dist = runif(1, -1, 1),
             b_Forest= runif(1, -1, 1),
             a_version = runif(constants$nversions, -1, 1),
             a_daysactive = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             a_EVI2 = runif(1, -1, 1),
             #a_occ = runif(1, -1, 1),
             eps_N = rnorm(constants$ncams, 0, 2),
             #eps_p = rnorm(constants$ncams, 0, 2),
             abundance_scale=rcat(3, c(1/3,1/3,1/3)),
             rho = array(data = runif(length(Bear_All), 0, 1),
                         dim=c(maxsites,no.occs,no.years)),
             sd_n = runif(1, 0, 2)
             #sd_p = runif(1, 0, 2)
  )
}

# parameters to monitor
keepers <- c("lambda", 'b_HFI', "b_Lat", "b_Yr", "b_Dist", "b_Forest", 
             "a_version", "a_daysactive", "a_EVI", "a_EVI2",
             "abundance_scale", "muy") #"a_occ",

# Will have to run chains for much longer (~40,000 iterations) to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 5000 # number of initial MCMC iterations to discard
ni <- 30000 # total number  of iterations

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

model_conf <- nimble::configureMCMC(model, enableWAIC = TRUE)

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

samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc, 
                           thin= 5,
                           inits=inits(),
                           WAIC = TRUE)
samples2 <- append(samples, list("formula"= "log(lambda[i, t]) <-  Year + Lat + HFIsc + Distsc + Forestsc + 1|camsite
                                 p <- version + daysactive + EVI + EVI2"))
saveRDS(samples, "./RNsamples30000Summer.rds")

samples <- readRDS("./RNsamples50000Summer.rds")

MCMCsummary(samples,params = "b_Lat", round = 2)
MCMCsummary(samples,params = "b_HFI", round = 2)
MCMCsummary(samples,params = "abundance_scale", round = 2)
MCMCsummary(samples,params = "a_version", round = 2) #detection really low backtransformed to response scale, naive detection 0.07


PR <- rnorm(15000, 0, 2)
MCMCtrace(samples, 
          params = c('b_Lat', 'b_HFI', 'b_Forest', "b_Yr", "b_Dist"),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
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
MCMCtrace(samples[[1]], 
          params = c("a_version", "a_daysactive", "a_EVI", "a_EVI2"),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c("lambda[991, 1]", "lambda[992, 1]", "lambda[993, 1]", "lambda[994, 1]", "lambda[995, 1]"), #still has nodes for "lambda[1243, 1]" and such
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
lambdameans <- MCMCpstr( samples, params = c("lambda"), func=mean, type="chain")[[1]]

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
