library(dplyr)
library(sf)
library(data.table)
library(tidyr)
library(nimble)
library(MCMCvis)
library(sswids)
library(stringr)
library(leaflet)
library(purrr)
library(units)

#make loo up table for zones for each camsite
Zones <- readRDS("./ModelingDF.rds")%>%ungroup()%>%select(cam_site_id, season)%>%distinct()%>%
  group_by(season)%>%mutate(siteID=row_number())%>%st_join(., st_transform(get_spatial_data("bear_zones"), 3071)%>%st_make_valid())%>%
  arrange(season, cam_site_id)%>%select(cam_site_id, yearID=season, siteID, bear_mgmt_zone_id)%>%
  mutate(cam_site_id_num = as.numeric(as.factor(cam_site_id)))%>%
  st_drop_geometry()

nsiteszoneyr <- Zones%>%group_by(yearID, bear_mgmt_zone_id)%>%summarise(N=n())%>%mutate(nyear=sum(N))

samples <- readRDS("./RNsamples.rds")

MCMCsummary(samples,params = "b_Lat", round = 2)
MCMCsummary(samples,params = "b_HFI", round = 2)
MCMCsummary(samples,params = "abundance_scale", round = 2)

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
MCMCtrace(samples, 
          params = c("a_version", "a_daysactive", "a_EVI", "a_occ"),
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
lambdameans <- MCMCpstr( samples, params = c("lambda"), func=mean, type="summary")[[1]]%>%
  as_tibble(rownames = "siteID") %>% pivot_longer(starts_with("V"), names_to="yearID") %>% mutate(yearID = parse_number(yearID), siteID=as.numeric(siteID))%>%
  drop_na(value)%>%left_join(., Zones)
lambdachains <- MCMCpstr( samples, params = c("lambda"), type="chains")[[1]]#%>%
lambdachains2 <- apply(lambdachains, 3, c)
lambdasd <- MCMCpstr( samples, params = c("lambda"), func=sd, type="summary")


############################################################################################
####                           pop by zone                                           #######
############################################################################################
# convert to mcmc object for inspection via coda package
samples_mcmc <- coda::as.mcmc.list(lapply(samples, coda::mcmc))
samplesdf <- do.call(rbind, samples_mcmc)
lambdas <- samplesdf%>%as.data.frame()%>%select(starts_with("lambda"))%>%mutate("iter"=1:nrow(.))%>%relocate(iter)%>%
  pivot_longer(starts_with("lambda"), names_to = "param")%>%drop_na(.,value)%>%
  mutate("siteID"=str_extract(string = param, pattern = "^.*\\[(\\d+),", group=1), "yearID"=str_extract(string = param, pattern = ", (\\d{1}).*$", group=1))
View(head(lambdas, 50))
lambdas%>%group_by(yearID)%>%summarise(N=n_distinct(siteID))
Zones$siteID <- as.character(Zones$siteID)
Zones$yearID <- as.character(Zones$yearID)
lambdaszone <- left_join(lambdas, Zones)
View(head(lambdaszone, 50))
lambdaszone2 <- lapply(1:6, function(x) filter(lambdaszone, as.numeric(yearID) == x & as.numeric(siteID) <= nsite[x]))%>%list_rbind()
lambdaszone2%>%group_by(yearID)%>%summarise(N=n_distinct(siteID))
#get average total pop by zone or average camsite lambda first?
popbyiterzone <- lambdaszone2%>%group_by(yearID, bear_mgmt_zone_id, iter)%>%summarise(sumsites=sum(value), sites=n_distinct(siteID), areasurveyed=sites*15, dens=sumsites/areasurveyed)
bearzones <- get_spatial_data("bear_zones")%>%st_make_valid()%>%mutate(area=set_units(st_area(.), "km^2"), aream2=st_area(.))
popbyyearzone <- popbyiterzone%>%group_by(yearID, bear_mgmt_zone_id)%>%left_join(., st_drop_geometry(bearzones[,c(1,6)]))%>%
  summarise(densbyzone=mean(dens), lowdens=quantile(dens, .025), highdens=quantile(dens, .975),
            pop=mean(dens)*unique(area), lowci=quantile(dens, .025)*unique(area), highci=quantile(dens, .975)*unique(area), nsite=unique(sites))

#pop by zone
ggplot(popbyyearzone, aes(x=as.numeric(yearID), y=pop, color=bear_mgmt_zone_id)) + geom_point() +
  geom_line(aes(group = bear_mgmt_zone_id)) + geom_ribbon(aes(ymax = highci, ymin = lowci, fill = bear_mgmt_zone_id), alpha=0.3)
#dens by zone
ggplot(popbyyearzone, aes(x=as.numeric(yearID), y=densbyzone, color=bear_mgmt_zone_id)) + geom_point() +
  geom_line(aes(group = bear_mgmt_zone_id)) + geom_ribbon(aes(ymax = highdens, ymin = lowdens, fill = bear_mgmt_zone_id), alpha=0.3)
popbyyearzone%>%group_by(yearID)%>%summarise(totalpop=sum(popbyzone2))

################ scrap ############################
cam2927 <- readRDS("./ModelingDF.rds")%>%ungroup()%>%select(cam_site_id, season)%>%distinct()%>%filter(cam_site_id == "cam_2927")
bearzones <- get_spatial_data("bear_zones")%>%st_make_valid()%>%mutate(area=set_units(st_area(.), "km^2"), aream2=st_area(.))
cam2927lambda <- lambdaszone2%>%filter(cam_site_id == "cam_2927")
unique(cam2927lambda$param)

leaflet() %>% 
  addProviderTiles('Esri.WorldImagery') %>%
  addPolygons(data=st_transform(bearzones, 4326),  fillOpacity = 0.75, label = bearzones$bear_mgmt_zone_id)%>%
  addCircleMarkers(data=st_transform(cam2927, 4326), fillOpacity = 0.75, fillColor = "black",   stroke=F, radius=5,
                   label = cam2927$cam_site_id)
st_join(cam2927, bearzones)
