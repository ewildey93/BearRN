library(ggplot2)
library(nimble)
library(MCMCvis)
# install.packages(c("here", "nimble", "coda", "stringr", "sf"))
library(here)
library(nimble)
library(coda)
library(stringr)
library(sf)
library(dplyr)
library(tidyr)
library(data.table)

ModelingDF <- readRDS("./ModelingDF.rds")%>%st_drop_geometry()%>%ungroup()
#ModelingDF <- RUGRCPUESaved%>%filter(cam_site_id != "cam_0246")%>%select(-24,-23) #wrong location for this cam in middle of lake cam_loc_seq_no == 
ModelingDF2 <- ModelingDF%>%mutate(across(14:58, scale))
ModelingDF2 <- ModelingDF2%>%mutate(across(c(2,3), as.factor, .names = "{col}_factor"))
ggplot(ModelingDF2, aes(x=CropRotation_1000,y=BEAR_ADULT_AMT)) + geom_point()
ggplot(RUGRCPUE3, aes(x=DaysOver20cm,y=GROUSE_AMT)) + geom_point()
ModelingDF2019 <- ModelingDF3%>%filter(year == 2019)%>%st_drop_geometry()


#separate Modeling DF into detection histories by year
dethist <- split(ModelingDF2, ModelingDF2$year)
dethist <- lapply(dethist, function (x) tidyr::pivot_wider(x, id_cols = c(cam_site_id, year), names_from = occ_factor, values_from = BEAR_ADULT_AMT, names_sort = TRUE)%>%ungroup())
dethist2 <- lapply(dethist, function(x) apply(x[,3:13],2, function(y) ifelse(y > 0, 1, 0)))
dethist2 <- Map(function(x, y, z) {rownames(x) <- paste0(y$cam_site_id, z); x}, x=dethist2, y=dethist, z=names(dethist))
dethistallyrs <- do.call(rbind, dethist2)
dethistcomplete <- dethistallyrs[complete.cases(dethistallyrs),]
completeguide <- which(complete.cases(dethistallyrs))
dethist2019 <- dethist2[["2019"]]


occNAs <- which(is.na(dethistallyrs))

#properly format site covariates to 1 row per camsite-year
sitecovs2019 <- ModelingDF2019%>%group_by(cam_site_id)%>%slice_head()
sitecovs <- ModelingDF2%>%group_by(cam_site_id, year)%>%select(14:57, 59)%>%slice_head()%>%arrange(year, cam_site_id)
sitecovscomplete <- sitecovs[completeguide,]
sitecovscomplete$HFI2020[is.na(sitecovscomplete$HFI2020)] <- 0


#properly format detection covariates to matrix of camsite-year x occasions
detcovs <- which(colnames(ModelingDF2) %in% c("meanEVI", "days_active"))#"occ_factor"
detcovs2 <- lapply(detcovs, y=ModelingDF2, 
                      function (x, y) pivot_wider(y, id_cols = c(cam_site_id, year), names_from = occ, 
                                                  values_from = all_of(x), values_fill = NA, names_sort = TRUE)%>%ungroup()%>%
                                                  arrange(year, cam_site_id))
names(detcovs2) <- c("days_active", "meanEVI")
detcovs3 <- lapply(detcovs2, function (x) as.matrix(x[,3:13]))
detcovs3 <- Map(function (x,y) {rownames(x) <- paste0(y$cam_site_id, y$year); x}, x=detcovs3, y=detcovs2)
detcovscomplete <- lapply(detcovs3, function(x) x[complete.cases(x),])
detcovscomplete2 <- lapply(detcovs3, function(x) x[completeguide,])




#put it together in data list for NIMBLE
data <- list(y=dethistallyrs, 
             HFI=as.vector(sitecovs$HFI2020), 
             Lat=as.vector(sitecovs$Lat), 
             Year=as.vector(as.numeric(sitecovs$year_factor)),
             effort=detcovs3[["days_active"]], 
             meanEVI=detcovs3[["meanEVI"]] )

which(is.na(sitecovscomplete$HFI2020))
hist(sitecovs$HFI2020)
HFIinits <- as.vector(ifelse(is.na(sitecovs$HFI2020), 0, NA))

#constants
#completeModelingDF <- ModelingDF%>%filter(paste0(cam_site_id, year) %in% rownames(dethistcomplete))
ncams <- length(unique(ModelingDF2$cam_site_id)) # number of unique cameras
cams <- str_extract(string = rownames(dethistallyrs), pattern = "cam_\\d{4}") #camera id for eachrow
camsites <- as.numeric(as.factor(str_extract(string = rownames(dethistallyrs), pattern = "cam_\\d{4}")))
ndeployments <- nrow(dethistallyrs)
nsurveys <- apply(dethistallyrs, MARGIN = 1, function (x) sum(!is.na(x)))
nyears <- 5
years <- sitecovs$year
constants <- list(ndeployments=ndeployments,
                  deployments=seq(1, ndeployments),
                  cams=cams, 
                  camsites=camsites,
                  ncams=ncams,
                  nsurveys=nsurveys) # this is wrong doesn't account for detection histories properly
                                     #   e.g. detection history: [0,0,0,0,NA,NA,NA,NA,0,0,0] this would only
                                     # count first 7 occasions




# .......................................................................
# MODEL CODE
# .......................................................................

camera <- nimble::nimbleCode( {
  
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
  
  
  
  ## Priors for detection parameters ##
  # coefficient for occasion
  a_EVI ~ dlogis(0, 1)
  a_effort ~ dlogis(0,1)
  
  
  # camera site random effect for both state level and detection level
  for(c in 1:ncams){
    eps_n[c] ~ dnorm(0, sd_n)
  }
  
  # hyperprior for detection and abundance random effect
  #sd_p ~ dgamma(1, 2)
  sd_n ~ dgamma(1, 2)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  #add year to this?
  # SDM - model for the latent state

    for(i in 1:ndeployments){
      log(lambda[i]) <-  eps_n[camsites[i]] + b_HFI*HFI[i] + b_Lat*Lat[i] + b_Yr*Year[i]#s[i] - spatial random effect
      n[i] ~ dpois(lambda[i])
    }
  
  # Camera submodel

  for(j in 1:ndeployments){
    for(k in 1:nsurveys[j]){
      muy[j, k] <- 1 - pow(rho[j, k], n[deployments[j]]) #
      logit(rho[j, k]) <-   a_effort*effort[j,k] + a_EVI*meanEVI[j,k] #no intercept in this model?
      y[j, k] ~ dbern(muy[j, k])
    }}
} )

# .......................................................................
# PREPARE MODEL TO RUN
# .......................................................................

# function to provide random initial values for parameters
inits <- function() {
  base::list(n = rep(1, constants$ndeployments),
             b_HFI = runif(1, -1, 1),
             b_Lat = runif(1, -1, 1),
             b_Yr = runif(1, -1, 1),
             a_effort = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             eps_n = rnorm(ncams, 0, 2),
             rho = matrix(data = runif(length(data$y), 0, 1),
                          nrow = nrow(data$y),
                          ncol = ncol(data$y)),
             sd_n = runif(1, 0, 2),
             HFI = HFIinits)
  }

# parameters to monitor
keepers <- c("lambda", 'b_HFI', "b_Lat", "b_Yr", "a_effort", "a_EVI")

# Will have to run chains for much longer (~40,000 iterations) to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 500 # number of initial MCMC iterations to discard
ni <- 5000 # total number  of iterations

# .......................................................................
# RUN MODEL
# .......................................................................


# create model
model <- nimble::nimbleModel(code = camera, 
                             data = data, 
                             constants = constants, 
                             inits = inits())

# check to see if everything is initialized
model$initializeInfo()
model$lambda
which(is.na(model$lambda))
model$HFI
which(is.na(model$HFI))

# compile the model
c_model <- nimble::compileNimble(model)

model_conf <- nimble::configureMCMC(model)

model_conf$addMonitors(keepers)

model_mcmc <- nimble::buildMCMC(model_conf)

c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)

# .......................................................................
# INSPECT RESULTS
# .......................................................................

# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples, coda::mcmc))

# Look at traceplots of the first five parameters
par(mfrow=c(3,2))
coda::traceplot(samples_mcmc[, 1:5])
# calculate Rhat convergence diagnostic for first five parameters
coda::gelman.diag(samples_mcmc[,1:5])

# extract mean and SD lambda of each grid cell
samplesdf <- do.call(rbind, samples_mcmc)
lambda <- samplesdf[, which(stringr::str_detect(string = colnames(samplesdf), pattern = 'lambda\\['))]
lambda_mean <- apply(lambda, 2, mean)
lambda_sd <- apply(lambda, 2, sd)

# map mean and standard deviation of expected abundance lambda
par(mfrow=c(1,1))
cov$lambda_mean <- lambda_mean
plot(cov["lambda_mean"], border = NA)
cov$lambda_sd <- lambda_sd
plot(cov["lambda_sd"], border = NA)

MCMCsummary(object = samples, round = 2, params=c('b_HFI', "b_Lat", "b_Yr", "a_effort", "a_EVI"))

samples %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = chain1[,"b_HFI"]), color = "white") + 
  labs(x = "beta HFI")

samples %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = chain1[,"b_Lat"]), color = "white") + 
  labs(x = "beta Lat")


MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "b_HFI")


##########################################################################################################
####                                  Long format detection non-detection data                        ####
##########################################################################################################
#Determine where we have data
has_data <- which(
  !is.na(dethistallyrs),
  arr.ind = TRUE
)

ob_cov_long <- matrix(
  1,
  nrow(has_data),
  ncol = 3
)

# Loop through the data and collect
#  the appropriate covariates.
for(i in 1:nrow(has_data)){
  ob_cov_long[i,2] <- detcovs3[["days_active"]][
    has_data[i,1],
    has_data[i,2]
  ]
  ob_cov_long[i,3] <- detcovs3[["meanEVI"]][
    has_data[i,1],
    has_data[i,2]
  ]
  
}

obcovsNAs <- cbind.data.frame(ob_cov_long, has_data)%>%rename(intercept =`1`, days.active=`2`, EVI=`3`)%>%
                              filter(is.na(EVI))

# We need to make the detection matrix long format as well.
#   We need this to have the same ordering as has_data.
#   This is actually really simple to do, we just remove
#   all the NA values and convert it to a vector like so.
y_long <- dethistallyrs[!is.na(dethistallyrs)]
deployment <- dethistallyrs[!is.na(dethistallyrs)]

# Finally, combine all the necessary data into a list for JAGS
# set up data list
data <- list(
  y = y_long,
  HFI=as.vector(sitecovs$HFI2020), 
  Lat=as.vector(sitecovs$Lat),
  Year=as.vector(as.numeric(sitecovs$year_factor)),
  effort=ob_cov_long[,2],
  meanEVI=ob_cov_long[,3]
)


ncams <- length(unique(ModelingDF2$cam_site_id)) # number of unique cameras
cams <- str_extract(string = rownames(dethistallyrs), pattern = "cam_\\d{4}") #camera id for eachrow
camsites <- as.numeric(as.factor(str_extract(string = rownames(dethistallyrs), pattern = "cam_\\d{4}")))
ndeployments <- nrow(dethistallyrs)
years <- sitecovs$year
nsurveys <- apply(dethistallyrs, MARGIN = 1, function (x) sum(!is.na(x)))
deploymentnames <- unlist(sapply(1:length(nsurveys), function (i) rep(rownames(dethistallyrs)[i], each=nsurveys[i])))
constants <- list(ndeployments=ndeployments,
                  deployments=as.numeric(as.factor(deploymentnames)),
                  cams=cams, 
                  camsites=camsites,
                  ncams=ncams,
                  nvisits=length(y_long))



camera <- nimble::nimbleCode( {
  
  
  # regression coefficient for ruffed grouse priority zone
  b_Yr ~ dnorm(0, 2)
  # regression coefficient for Human Footprint Index
  b_HFI ~ dnorm(0, 2)
  # regression coefficient for Lat
  b_Lat ~ dnorm(0, 2)
  
  ## Priors for detection parameters ##
  # coefficient for occasion
  a_EVI ~ dlogis(0, 1)
  a_effort ~ dlogis(0,1)
  alpha_0 ~ dunif(0,1)
  
  # camera site random intercept for both state level 
  for(c in 1:ncams){
    eps_n[c] ~ dnorm(0, sd_n)
  }
  
  # hyperprior for detection and abundance random effect
  sd_n ~ dgamma(1, 2)
  
  # .............................................................
  # LIKELIHOOD
  # .............................................................
  
  #add year to this?
  # SDM - model for the latent state
  #impute missing HFI values
  for(i in 1:ndeployments){
    HFI[i]~dnorm(mu_HFI, sigma_HFI)
  }
  mu_HFI ~ dnorm(0, 0.01)
  sigma_HFI ~ dunif(0,2) 
 
  
  for(i in 1:ndeployments){
    log(lambda[i]) <-  eps_n[camsites[i]] + b_HFI*HFI[i] + b_Lat*Lat[i] + b_Yr*Year[i]
    n[i] ~ dpois(lambda[i])
  }
  
  # Camera submodel
  
  for(j in 1:nvisits){
      muy[j] <- 1 - pow(1-rho[j], n[deployments[j]]) 
      logit(rho[j]) <-  alpha_0 + a_effort*effort[j] + a_EVI*meanEVI[j] 
      y[j] ~ dbern(muy[j])
    }
} )


# function to provide random initial values for parameters
HFIinits <- as.vector(ifelse(is.na(sitecovs$HFI2020), 0, NA))
inits <- function() {
  base::list(n = rep(1, constants$ndeployments),
             b_HFI = runif(1, -1, 1),
             b_Lat = runif(1, -1, 1),
             b_Yr = runif(1, -1, 1),
             a_effort = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             alpha_0 = runif(1, 0, 1),
             eps_n = rnorm(ncams, 0, 2),
             rho = runif(length(data$y), 0, 1),
             sd_n = runif(1, 0, 2),
             HFI = HFIinits,
             mu_HFI = runif(1, -1, 1),
             sigma_HFI = runif(1, 0, 2)
             )
}

# parameters to monitor
keepers <- c("lambda", 'b_HFI', "b_Lat", "b_Yr", "a_effort", "a_EVI")

# Will have to run chains for much longer (~40,000 iterations) to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 500 # number of initial MCMC iterations to discard
ni <- 5000 # total number  of iterations

# .......................................................................
# RUN MODEL
# .......................................................................


# create model
model <- nimble::nimbleModel(code = camera, 
                             data = data, 
                             constants = constants, 
                             inits = inits())

# check to see if everything is initialized
model$initializeInfo()
model$lambda
which(is.na(model$lambda))
model$HFI
which(is.na(model$HFI))

# compile the model
c_model <- nimble::compileNimble(model)

model_conf <- nimble::configureMCMC(model)

model_conf$addMonitors(keepers)

model_mcmc <- nimble::buildMCMC(model_conf)

c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)

samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc)

MCMCsummary(object = samples, round = 2, params=c('b_HFI', "b_Lat", "b_Yr", "a_effort", "a_EVI"))

samples %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = chain1[,"b_HFI"]), color = "white") + 
  labs(x = "beta HFI")

samples %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = chain1[,"b_Lat"]), color = "white") + 
  labs(x = "beta Lat")

samples %>%
  as_tibble() %>%
  ggplot() + 
  geom_histogram(aes(x = chain1[,"a_EVI"]), color = "white") + 
  labs(x = "beta EVI")


MCMCtrace(object = samples,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          params = "b_Lat")


# extract mean and SD lambda of each grid cell
samplesdf <- do.call(rbind, samples)
lambda <- samplesdf[, which(stringr::str_detect(string = colnames(samplesdf), pattern = 'lambda\\['))]
lambda_mean <- apply(lambda, 2, mean)
lambda_sd <- apply(lambda, 2, sd)


# map mean and standard deviation of expected abundance lambda
par(mfrow=c(1,1))
sitecovs$lambda_mean <- lambda_mean
plot(sitecovs["lambda_mean"], border = NA)
sitecovs$lambda_sd <- lambda_sd
plot(sitecovs["lambda_sd"], border = NA)

###########################################scrap###########################
cams2 <- as.data.frame(cams)
detcovsoff <- lapply(detcovscomplete2, function(x) rownames(x)[which(!complete.cases(x))])
wrongdetcovs<- ModelingDF%>%filter(paste0(cam_site_id, year) %in% detcovsoff[[1]])
