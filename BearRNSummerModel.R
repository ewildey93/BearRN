library(nimble)

load("path to file/BearRN.RData")

RNcode <- nimbleCode({
  
  # .............................................................
  # PRIORS
  # .............................................................
  
  
  # regression coefficient for bear zone
  for (z in 1:nZones){
    b_Zone[z] ~ dnorm(0, 2)
  }
  # regression coefficient for zone-year effect
  for (z in 1:nZones){
    b_ZoneYr[z] ~ dnorm(0, 2)
  }
  # regression coefficient for Developed
  b_Dev ~ dnorm(0, 2)
  # regression coefficient for Dist
  b_Dist ~ dnorm(0, 2)
  # regression coefficient for Forest
  b_Forest ~ dnorm(0, 2)
  # regression coefficient for Corn
  b_Corn ~ dnorm(0, 2)
  # regression coefficient for Wetland
  b_Wetland ~ dnorm(0, 2)
  
  
  
  #coefficients of spatial spline
  for (k in 1:sp.nknots) {
     spat.spline.b[k] ~ dnorm(0,sigma.spat.spline.b)
  }
  #coefficients of EVI spline
  for (k in 1:EVI.nknots) {
    b[k] ~ dnorm(0,sigma.EVI)
  }
  
  # spline random effect priors
  sigma.spat.spline.b~dunif(0,100)
  
  # spline random effect priors
  sigma.EVI~dunif(0,100)
  
  ## Priors for detection parameters ##
  # fixed effect coefficient for EVI
  a_EVI ~ dlogis(0, 1)
  # coefficient for effort
  a_daysactive ~ dlogis(0,1)
  # coefficient for camera version
  for(v in 1:nversions){
    a_version[v] ~ dlogis(0,1)
  }
  
  
  ## Priors for scales  ##
  abundance_scale[1] ~ dcat(catprobs[1:4])
  abundance_scale[2] ~ dcat(catprobs[1:4])
  abundance_scale[3] ~ dcat(catprobs[1:4])
  abundance_scale[4] ~ dcat(catprobs[1:4])
  abundance_scale[5] ~ dcat(catprobs[1:4])
  
  
  for( t in 1:nyear ) { 
    #state model
    # Loop through only the sites that are surveyed in a given year. 
    for( i in 1:nsite[t] ){
      N[i, t] ~ dpois( lambda[ i, t ] )
      log(lambda[i, t]) <-  b_Zone[Zone[i,t]] + b_ZoneYr[Zone[i,t]]*Year[i, t] +  
        #scaled parameters
        b_Corn*Corn[i,t,abundance_scale[1]] + b_Dev*Dev[i,t,abundance_scale[2]] +
        b_Dist*Dist[i,t,abundance_scale[3]] + b_Forest*Forest[i,t,abundance_scale[4]] + 
        b_Wetland*Wetland[i,t,abundance_scale[5]] +
        s[i,t] #s[i,t] - spatial random effect
      
      s[i,t] <- inprod(spat.spline.b[1:sp.nknots],sp.Z[i,t,1:sp.nknots])
      
      #detection model  
      for(k in 1:nsurveys[i,t]){
        muy[i, k, t] <- 1 - pow(1-rho[i,k,t], N[i, t]) #
        logit(rho[i, k, t]) <- a_version[camversion[i,k,t]] + a_daysactive*daysactive[i,k,t] + 
          a_EVI * EVI[i,k,t] + EVI.spline[i,k,t] 
        y[i, k, t] ~ dbern(muy[i, k, t])
        
        EVI.spline[i,k, t] <- inprod(b[1:EVI.nknots],EVI.Z[i,k,t,1:EVI.nknots])
        
      }
      
    }
  }
  
  #GOF
  for( t in 1:nyear ) {
    for( i in 1:nsite[t] ){
      for(k in 1:nsurveys[i,t]){
        ynew[i, k, t] ~ dbern(muy[i, k, t])
      }
      sum.y[i,t] <- sum(y[i,1:nsurveys[i,t],t])                    # Summation of observed detections per site
      sum.ynew[i,t] <- sum(ynew[i,1:nsurveys[i,t],t])              # Summation of replicated detection per site
      e.count[i,t] <- sum(muy[i,1:nsurveys[i,t],t])                # Expected detections per site
      
      # Chi-square discrepancy for the actual data.
      chi2.actual[i,t] <- pow((sum.y[i,t]-e.count[i,t]), 2) / (e.count[i,t] + e)   # e is a small constant to avoid any division by zero (Kery and Royle 2016)
      
      # Chi-square discrepancy for the simulated data
      chi2.sim[i,t] <- pow((sum.ynew[i,t]-e.count[i,t]), 2) / (e.count[i,t] + e)
    }
  }
  
  chifit.actual <- sum(chi2.actual[1:1071,1:6])
  chifit.sim <- sum(chi2.sim[1:1071,1:6])
  
  
})



# function to provide random initial values for parameters
inits <- function() {
  base::list(N = matrix(data = rep(1, max(constants$nsite)*constants$nyear),
                        nrow = max(constants$nsite),
                        ncol = constants$nyear),
             b_Zone = runif(constants$nZones, -1, 1),
             b_ZoneYr = runif(constants$nZones, -1, 1),
             b_Dev = runif(1, -1, 1),
             b_Dist = runif(1, -1, 1),
             b_Forest= runif(1, -1, 1),
             b_Corn= runif(1, -1, 1),
             b_Wetland= runif(1, -1, 1),
             a_version = runif(constants$nversions, -1, 1),
             a_daysactive = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             spat.spline.b=rep(1,constants$sp.nknots),
             sigma.spat.spline.b=1,
             b=rep(1,constants$EVI.nknots),
             sigma.EVI=1,
             abundance_scale=rcat(5, c(0.25,0.25,0.25,0.25)),
             rho = array(data = runif(length(Bear_All), 0, 1),
                         dim=c(maxsites,no.occs,no.years)),
             ynew=array(data = rep(1, length(Bear_All)),
                        dim=c(maxsites,no.occs,no.years)),
             sum.y=matrix(rep(5,max(constants$nsite)*constants$nyear),
                          nrow = max(constants$nsite),
                          ncol = constants$nyear),
             sum.ynew=matrix(rep(5,max(constants$nsite)*constants$nyear),
                             nrow = max(constants$nsite),
                             ncol = constants$nyear),
             e.count=matrix(rep(4.5,max(constants$nsite)*constants$nyear),
                            nrow = max(constants$nsite),
                            ncol = constants$nyear),
             chi2.actual=matrix(runif(max(constants$nsite)*constants$nyear, 1, 50),
                                nrow = max(constants$nsite),
                                ncol = constants$nyear),
             chi2.sim=matrix(runif(max(constants$nsite)*constants$nyear, 1, 50),
                             nrow = max(constants$nsite),
                             ncol = constants$nyear),
             chifit.actual=runif(1, 1, 100),
             chifit.sim=runif(1, 1, 100)
  )
}


# parameters to monitor
keepers <- c("b_Zone", "b_ZoneYr",'b_Dev',
             "b_Dist", "b_Forest", "b_Corn",
             "b_Wetland",
             "b", "spat.spline.b",
             "a_version", "a_daysactive", "a_EVI", 
             "abundance_scale", "chifit.actual",
             "chifit.sim") 

nc <- 3 # number of chains
nb <- 10000 # number of initial MCMC iterations to discard
ni <- 75000 # total number  of iterations

# .......................................................................
# RUN MODEL
# .......................................................................


# create model
model <- nimble::nimbleModel(code = RNcode, 
                             data = data, 
                             constants = constants, 
                             inits = inits())


# compile the model
c_model <- nimble::compileNimble(model)
model_conf <- nimble::configureMCMC(model)
model_conf$addMonitors(keepers)
model_mcmc <- nimble::buildMCMC(model_conf)
c_model_mcmc <- nimble::compileNimble(model_mcmc, project = model)
samples <- nimble::runMCMC(c_model_mcmc, 
                           nburnin = nb, 
                           niter = ni, 
                           nchains = nc, 
                           thin= 5,
                           inits=inits())