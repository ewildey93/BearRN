############################################
############# R-N Model ####################
############################################
# Create constants for model from array.
constants <- list(
  nyear = dim(Bear_All)[[3]],
  nsite = nsite,
  nversions=length(unique(cam_version$cam_versionID2))-1,
  nsurveys=nsurveys2,
  camversion=camversion_array,
  sp.nknots=50
)

# Bundle data (counts and covariates).
data <- list(
  y = Bear_All,
  Year = yr,
  X=X,
  Y=Y,
  sp.Z= sp.Zarray
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
  # regression coefficient for X coordinate
  b_X ~ dnorm(0, 2)
  # regression coefficient for Y coordinate
  b_Y ~ dnorm(0, 2)

  for(v in 1:nversions){
    a_version[v] ~ dlogis(0,1)
  }
  
  
  for (k in 1:sp.nknots) {
    spat.spline.b[k] ~ dnorm(0,sigma.spat.spline.b)
  }

  
  # spline random effect priors
  sigma.spat.spline.b~dunif(0,100)
  
  
  
  for( t in 1:nyear ) { #loop over site then year?
    #state model
    # Loop through only the sites that are surveyed in a given year. 
    for( i in 1:nsite[t] ){
      N[i, t] ~ dpois( lambda[ i, t ] )
      log(lambda[i, t]) <-  b_Yr*Year[i, t] + #make this a factor? 
        b_X*X[i, t] + b_Y*Y[i, t] +
        #cam_site random effect
        s[i,t] #s[i,t] - spatial random effect
      
      s[i,t] <- inprod(spat.spline.b[1:sp.nknots],sp.Z[i,t,1:sp.nknots])
      
      #detection model  
      for(k in 1:nsurveys[i,t]){
        muy[i, k, t] <- 1 - pow(1-rho[i,k,t], N[i, t]) #
        logit(rho[i, k, t]) <- a_version[camversion[i,k,t]]
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
             b_X = runif(1, -1, 1),
             b_Y = runif(1, -1, 1),
             a_version = runif(constants$nversions, -1, 1),
             spat.spline.b=rep(1,constants$sp.nknots),
             sigma.spat.spline.b=1,
             rho = array(data = runif(length(Bear_All), 0, 1),
                         dim=c(maxsites,no.occs,no.years))
  )
}

# parameters to monitor
keepers <- c("lambda",  "b_X","b_Y", "b_Yr", 
             "spat.spline.b", "muy", "a_version") 

# Will have to run chains for much longer (~40,000 iterations) to approach convergence
# running with 200 iterations took about 10 minutes on my laptop with 4 cores
# to speed things up, particularly for longer chains, you can run chains in parallel
# see: https://groups.google.com/g/nimble-users/c/RHH9Ybh7bSI
nc <- 3 # number of chains
nb <- 5000 # number of initial MCMC iterations to discard
ni <- 100 # total number  of iterations

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
mcmc.output <- nimble::runMCMC(c_model_mcmc, 
                              nburnin = 0, 
                              niter = 100, 
                              nchains = 1,
                              inits=inits(),
                              WAIC = FALSE)
PR <- rnorm(15000, 0, 2)
MCMCtrace(samples[[1]], 
          params = c('b_X', 'b_Y', "b_Yr"),
          ISB = FALSE,
          exact = TRUE,
          priors = PR,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples[[1]], 
          params = c('spat.spline.b[1]', 'spat.spline.b[2]', 'spat.spline.b[3]', 'spat.spline.b[4]', 'spat.spline.b[5]'),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCsummary(samples[[1]], 
          params = c('lambda[138, 1]', 'lambda[1255, 1]', 'lambda[1309, 1]', 'lambda[1107, 2]', 'lambda[1108, 2]',
                     'lambda[1209, 3]', 'lambda[1210, 3]', 'lambda[1303, 4]', 'lambda[1289, 5]', 'lambda[1290, 5]',
                     'lambda[1492, 6]'),
          ISB = FALSE,
          round=2)
MCMCtrace(mcmc.output, 
            params = c('lambda[138, 1]', 'lambda[1255, 1]', 'lambda[1309, 1]', 'lambda[1107, 2]', 'lambda[1108, 2]',
                       'lambda[1209, 3]', 'lambda[1210, 3]', 'lambda[1303, 4]', 'lambda[1289, 5]', 'lambda[1290, 5]',
                       'lambda[1492, 6]'),
            ISB = FALSE,
            pdf = FALSE,
            exact=TRUE)
MCMCsummary(mcmc.output, 
          params = c('lambda[138, 1]', 'lambda[1255, 1]', 'lambda[1309, 1]', 'lambda[1107, 2]', 'lambda[1108, 2]',
                     'lambda[1209, 3]', 'lambda[1210, 3]', 'lambda[1303, 4]', 'lambda[1289, 5]', 'lambda[1290, 5]',
                     'lambda[1492, 6]'),
          ISB = FALSE,
          round=2)
