

EVI.pred <- seq(from=range(EVI2$meanEVI)[1], to=range(EVI2$meanEVI)[2], length.out=100)
pred.data <- expand.grid(camversion=1, daysactive=0.2211643, EVI=EVI.pred)
# prediction dataset for spatial smoothing
pred.EVI.Z_K <- (abs(outer(as.numeric(pred.data$EVI),EVI.knots,"-")))^3

pred.EVI.Z <- t(solve(EVI.sqrt.OMEGA_all,t(pred.EVI.Z_K)))

# standardize for better performance
pred.EVI.Z <- (pred.EVI.Z-meanZ)/sdZ
pred.data2 <- cbind(pred.data, pred.EVI.Z)
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
  EVI.nknots=5
)

# Bundle data (counts and covariates).
data <- list(
  y = Bear_All,
  daysactive=daysactive_array,
  EVI=EVI_array,
  EVI.Z=EVIZ_array
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
  
  
  

  for (k in 1:EVI.nknots) {
    b[k] ~ dnorm(0,sigma.EVI)
  }
  
  # spline random effect priors
  sigma.EVI~dunif(0,100)
  
  ## Priors for detection parameters ##
  # coefficient for occasion
  #a_EVI ~ dlogis(0, 1)

for (t in 1:nyear){
  for (i in 1:nsite[t]){
  lambda[i,t] ~ dgamma(0.5, 0.5)
  }
}
 
  ## Priors for detection parameters ##
  # coefficient for occasion
  a_daysactive ~ dlogis(0,1)
  a_EVI ~ dlogis(0,1)
  #a_occ ~ dlogis(0,1)
  for(v in 1:nversions){
    a_version[v] ~ dlogis(0,1)
  }
  
  for( t in 1:nyear ) { #loop over site then year?
    #state model
    # Loop through only the sites that are surveyed in a given year. 
    for( i in 1:nsite[t] ){
      N[i, t] ~ dpois( lambda[ i, t ] )
      
      #detection model  
      for(k in 1:nsurveys[i,t]){
        muy[i, k, t] <- 1 - pow(1-rho[i,k,t], N[i, t]) #
        logit(rho[i, k, t]) <- a_version[camversion[i,k,t]] + a_daysactive*daysactive[i,k,t] + a_EVI*EVI[i,k,t] +
          EVI.spline[i,k,t] #EVI and occ probably correlated + a_occ*occ[i,k,t] + eps_p[camsites[i,t]]
        y[i, k, t] ~ dbern(muy[i, k, t])
        
        EVI.spline[i,k,t] <- inprod(b[1:EVI.nknots],EVI.Z[i,k,t,1:EVI.nknots])
        #EVI[i,k,t] <- b[1]*EVI.Z[i,k,t,1] + b[2]*EVI.Z[i,k,t,2] + b[3]*EVI.Z[i,k,t,3] + b[4]*EVI.Z[i,k,t,4] + b[5]*EVI.Z[i,k,t,5]
      }
    }
  }
  
  #PREDICTION
  # for (i in 1:npred){
  #   N.pred[i] ~ dpois(lambda.pred[i])
  #   
  #   
  #     muy.pred[i] <- 1 - pow(1-rho.pred[i], N.pred[i]) #
  #     logit(rho.pred[i]) <- a_version[camversion.pred[i]] + a_daysactive*daysactive.pred[i] + a_EVI*EVI.pred[i] + 
  #       pred.EVI.spline[i] #EVI and occ probably correlated + a_occ*occ[i,k,t] + eps_p[camsites[i,t]]
  #     y.pred[i] ~ dbern(muy.pred[i])
  #     
  #     pred.EVI.spline[i] <- inprod(b[1:EVI.nknots],pred.EVI.Z[i,1:EVI.nknots])
  #     #EVI[i,k,t] <- b[1]*EVI.Z[i,k,t,1] + b[2]*EVI.Z[i,k,t,2] + b[3]*EVI.Z[i,k,t,3] + b[4]*EVI.Z[i,k,t,4] + b[5]*EVI.Z[i,k,t,5]
  # 
  # }
  
})




# function to provide random initial values for parameters
inits <- function() {
  base::list(N = matrix(data = rep(1, max(constants$nsite)*constants$nyear),
                        nrow = max(constants$nsite),
                        ncol = constants$nyear),
             a_version = runif(constants$nversions, -1, 1),
             a_daysactive = runif(1, -1, 1),
             a_EVI = runif(1, -1, 1),
             b=rep(1,constants$EVI.nknots),
             sigma.EVI=1,
             lambda=matrix(data=rep(1,max(constants$nsite)*constants$nyear),
                           nrow = max(constants$nsite),
                           ncol = constants$nyear),
             rho = array(data = runif(length(Bear_All), 0, 1),
                         dim=c(maxsites,no.occs,no.years))
  )
}

# parameters to monitor
keepers <- c("lambda",
             "a_version", "a_daysactive", "a_EVI",
             "b",
             "muy[1,3,6]", "muy[1,9,6]")


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
model$calculate("rho.pred")

# compile the model
c_model <- nimble::compileNimble(model)

model_conf <- nimble::configureMCMC(model)

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
                           inits=inits())

MCMCtrace(samples, 
          params = c("a_version", "a_daysactive"),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c("muy[1217, 10, 5]"),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c("lambda[100, 2]", "lambda[100, 4]", "lambda[100, 5]"),
          ISB = FALSE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCtrace(samples, 
          params = c("b"),
          ISB = TRUE,
          exact = TRUE,
          pdf = FALSE,
          Rhat = TRUE,
          n.eff = TRUE)
MCMCsummary(samples,params = "muy[1217, 10, 5]", round = 3, ISB = FALSE)

muy.1.1 <- paste0("muy[1, ", 1:14, ", 1]")
muy1.1samples <- samples[["chain1"]][,which(colnames(samples[["chain1"]]) %in% muy.1.1)]
muy1.1means <- apply(muy1.1samples, 2, mean)
EVI1.1 <- EVI2$meanEVI[1:14] * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
muyEVI1.1plot <- cbind(muy1.1means, EVI1.1)
ggplot(muyEVI1.1plot, aes(x=EVI1.1, y=muy1.1means)) + geom_point() + geom_smooth()


muysamples <- samples[["chain1"]][,grep(pattern = "muy\\[\\d+, \\d+, 1\\]", x = colnames(samples[["chain1"]]))]
muymeans <- apply(muysamples, 2, mean)
muymeans <- muymeans[!is.na(muymeans)]
EVI1 <- EVI2$meanEVI[EVI2$yearID ==1 & !is.na(EVI2$det)] * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
muyEVI1plot <- cbind(muymeans, EVI1)
ggplot(muyEVI1plot, aes(x=EVI1, y=muymeans)) + geom_point() #+ geom_smooth()

#get mean response of spline
MCMCsummary(samples,params = "b", round = 3, ISB = T)
bmeans <- MCMCsummary(samples,params = "b", round = 3, ISB = T)
camversionmeans <- MCMCsummary(samples,params = "a_version", round = 3, ISB = T)
daysactivemean <- MCMCsummary(samples,params = "a_daysactive", round = 3, ISB = T)
daysactivemean <- MCMCsummary(samples,params = "a_EVI", round = 3, ISB = T)
rho.mean <- plogis(camversionmeans$mean[1] + daysactivemean$mean*0.2211643 + pred.data$EVI +bmeans$mean[1]*pred.EVI.Z[,1] + bmeans$mean[2]*pred.EVI.Z[,2] + bmeans$mean[3]*pred.EVI.Z[,3] +
                                                bmeans$mean[4]*pred.EVI.Z[,4] + bmeans$mean[5]*pred.EVI.Z[,5])
#get CRIs for spline
#combine MCMC chains into one
allchains <- MCMCchains(samples, params =c("a_version[1]", "a_daysactive", "a_EVI", "b"))
#loop through estimated parameter at each iteration
rho.preds <- array(dim = c(length(pred.data$EVI), length(out1N[,"alpha0"])))
for(j in 1:length(out1N[,"alpha0"])){
  rho.preds[,j] <- plogis(allchains[,"a_version[1]"][j] + allchains[,"a_daysactive"][j]*0.2211643 + allchains[,"a_EVI"][j]*pred.data$EVI + 
                            allchains[,"b[1]"][j]*pred.EVI.Z[,1] + allchains[,"b[2]"][j]*pred.EVI.Z[,2] + allchains[,"b[3]"][j]*pred.EVI.Z[,3] +
                            allchains[,"b[4]"][j]*pred.EVI.Z[,4] + allchains[,"b[5]"][j]*pred.EVI.Z[,5])
}
#calculate interval
CL <- apply(rho.preds, 1, function(x){quantile(x, prob = c(0.025, 0.975))})

rhoplotdf <- cbind(pred.data2, rho.mean, t(CL))
rhoplotdf$EVI <- rhoplotdf$EVI * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
ggplot(rhoplotdf, aes(x=EVI, y=rho.pred)) + geom_point() + geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5))

rhoplotdf <- cbind(pred.data2, rho.mean)
rhoplotdf$EVI <- rhoplotdf$EVI * attr(EVI2$meanEVI, "scaled:scale") + attr(EVI2$meanEVI, "scaled:center")
ggplot(rhoplotdf, aes(x=EVI, y=rho.mean)) + geom_point() #+ geom_ribbon(aes(ymin = CL2.5, ymax = CL97.5))


##########################################scrap############################################
which(grep(pattern = "muy\\[\\d+, \\d+, 1\\]", x = colnames(samples[["chain1"]]), value = T)=="muy[2, 12, 7]")
muynamesyear1 <- as.data.frame(colnames(muysamples))
muynamesyear1$site <- str_extract(string = muynamesyear1$`colnames(muysamples)`, pattern = "(muy\\[)(\\d+)", group = 2)
muynamesyear1$occ <- str_extract(string = muynamesyear1$`colnames(muysamples)`, pattern = "(muy\\[\\d+, )(\\d+)", group = 2)
muynamesyear1$yr <- str_extract(string = muynamesyear1$`colnames(muysamples)`, pattern = "(muy\\[\\d+, \\d+, )(\\d+)", group = 2)
muynamesyear1$means <- muymeans


#for predicting within MCMC that didn't work
,
N.pred = rep(1,20),
rho.pred = runif(20, 0 ,1),
lambda.pred=rep(1,20)
# parameters to monitor
 #maybe don't save muys in the future