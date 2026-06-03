library(ggmap)
library(sf)
library(CropScapeR)
library(landscapemetrics)
library(terra)
library(purrr)
library(tigris)
library(cowplot)
library(ggplot2)
library(sswids)
library(shadowtext)
library(units)

location <- data.frame("X"=-91.312287, "Y"=45.489776)
locationSF <- st_as_sf(location, coords=c("X", "Y"), crs=4326)
locationSF3071 <- st_transform(locationSF, 3071)
locationSFbuffers <- lapply(c(500, 1000, 2500, 5000), function(x) st_buffer(locationSF3071, dist = x, ))
locationSFbuffers4326 <- lapply(locationSFbuffers, function (x) st_transform(x, crs=4326))
locationSFbuffers4326 <- do.call(rbind, locationSFbuffers4326)
locationSFbuffers4326$BufferSize <- as.character(c(0.5, 1, 2.5, 5))

CropRasts <- list.files("C:/Users/wildeefb/Documents/GeoSpatial/BearCrops/", full.names = TRUE)
CropRasts2 <- lapply(CropRasts, function (x) rast(x))
data("linkdata")
locationSF2 <- st_transform(locationSF, crs=crs(CropRasts2[[1]]))


corndemo <- c(500, 1000, 2500, 5000) %>% 
  set_names() %>% 
  # produce a dataframe after this is all done
  map_dfr( 
  ~sample_lsm(
            # raster layer
            landscape = CropRasts2[[1]],
            # camera locations
            y = locationSF2,
            # get landcover class level metrics
            level = "class",
            # return NA values for classes not in buffer
            # all_classes = TRUE, 
            # can do multiple metrics at once
            what = 'lsm_c_pland',
            # buffer sizes to use
            size = .x, 
            # default is square buffer
            shape = "circle", 
            # turn warnings on or off
            verbose = FALSE 
          ), # get buffer size column in the output
  .id = "buffer_size"
  )

corndemo2 <- left_join(corndemo, linkdata, by=join_by(class==MasterCat))
corndemo2 <- corndemo2[grep(pattern = "Corn", x = corndemo2$Crop, ignore.case = TRUE),]
corndemo2 <- corndemo2%>%group_by(buffer_size)%>%summarise(CornProp=sum(value))
corndemo2$Abundance <- 3
corndemo2$buffer_size <- as.numeric(corndemo2$buffer_size)
noise <- data.frame("buffer_size"=sample(corndemo2$buffer_size, size = 30, replace = T),
                    "Abundance"=runif(30, 0 , 2.75))
us_states = states(cb = FALSE, class = "sf")
us_states = subset(us_states, 
                   !NAME %in% c(
                     "United States Virgin Islands",
                     "Commonwealth of the Northern Mariana Islands",
                     "Guam",
                     "American Samoa",
                     "Puerto Rico",
                     "Alaska",
                     "Hawaii"
                   ))
us_states3071 <- st_transform(us_states, crs=3071)
WisconsinBB <- st_as_sfc(st_bbox(Wisconsin))
WisconsinBBBuffer <- st_buffer(WisconsinBB, dist=60000)
inset <- ggplot() + 
  geom_sf(data = us_states3071, fill = "white", size = 0.2) + 
  geom_sf(data = WisconsinBBBuffer, fill = NA, color = "red", size = 1.2) +
  theme_void()

#download basemap
register_google()
map2 <- get_map(location = st_coordinates(locationSF),
                maptype="satellite", source="google", zoom=12)
cbPalette <- c("#F0E442", "#56B4E9", "#009E73", "#E69F00" , "#0072B2", "#D55E00")

ggmap(map2, baselayer=Wisconsin4326, maprange=TRUE) +
  geom_sf(data=locationSFbuffers4326, inherit.aes=FALSE,aes(color=BufferSize), fill=NA,  linewidth=2) +
  scale_color_brewer(palette = "Set2",
                    name="Proportion of Corn",
                    labels = c("0.5" = "51.1%",
                               "1" = "35.7%",
                               "2.5" = "13.1%",
                               "5" = "11.7%")) +
  theme(axis.title = element_blank(), # Removes both x and y axis titles
        axis.text = element_blank(),  # Removes both x and y axis tick labels
        axis.ticks = element_blank(),
        legend.position=c(.151,.852)
        )



ggplot(corndemo2, aes(x=CornProp, y=Abundance, color=as.factor(buffer_size))) + geom_point(size=4) +
  #geom_point(data = noise, aes(x=buffer_size, y=Abundance), color="black", size=4) +
  scale_color_manual(values=cbPalette, name="Buffer Size") +
  ylim(0, 4) +
  theme(panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        panel.background = element_blank(), # Ensure the panel background is blank/transparent
        axis.line = element_line(colour = "black"), # Explicitly draw the bottom/left axis lines
        panel.border = element_blank()) +
  xlab("Corn Land Cover") +
  ylab("Site Abundance")
  

bearzones <- get_spatial_data("bear_zones")%>%st_transform(., 3071)
bearrangeU <- st_sf(data.frame("BearRange"="1") ,geometry = st_union(bearrange2))
coordsdf$label <- "Unique Camera Site"
bearea <- st_area(st_union(bearrange2))
units(bearea) <- "km^2"

bearzonesrange <- ggplot() + 
   geom_sf(data=bearrangeU, aes(fill=BearRange)) +
  geom_sf(data=st_as_sf(coordsdf, coords=c("X", "Y"), crs=3071),aes(color=label)) +
  geom_sf(data=bearzones, fill=NA, color="black") + 
  geom_shadowtext(data = bearzones, 
                  aes(label = bear_mgmt_zone_id, 
                      geometry = geoms), # specify geometry for sf data
                  stat = "sf_coordinates",
                  size = 9,
                  color = "black",
                  bg.color = "white",      # color of the outline
                  bg.r = 0.05,              # thickness of the outline
                  show.legend = FALSE) +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        panel.background = element_blank(),
        axis.text = element_blank(),  # Removes both x and y axis tick labels
        axis.ticks = element_blank(),# Ensure the panel background is blank/transparent
        panel.border = element_blank(),
        legend.title = element_blank()) +
  scale_color_manual(values="black", label="Unique\n Camera\n Site") +
  scale_fill_manual(values="#0072B2", label="Bear Range") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme(legend.spacing = unit(0.05, "cm"))

bearzones <- ggplot() + 
  #geom_sf(data=bearrangeU, aes(fill=BearRange)) +
  geom_sf(data=bearzones, fill=NA, color="black") + geom_sf_text(data = bearzones, aes(label = bear_mgmt_zone_id), size=9, color="black") + 
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        panel.background = element_blank(),
        axis.text = element_blank(),  # Removes both x and y axis tick labels
        axis.ticks = element_blank(),# Ensure the panel background is blank/transparent
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") #+
  scale_fill_manual(values="#0072B2", label="Bear Range")

gg_inset_map1 = ggdraw() +
  draw_plot(bearzonesrange) +
  draw_plot(inset, x = 0.6, y = 0.6, width = 0.4, height = 0.4)
  
gg_inset_map1
ggplot(st_union(bearrange2)) + geom_sf(fill="#0072B2", color="white") +  geom_sf(data=Wisconsin2, inherit.aes = FALSE, fill=NA) +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        panel.background = element_blank(),
        axis.text = element_blank(),  # Removes both x and y axis tick labels
        axis.ticks = element_blank(),# Ensure the panel background is blank/transparent
        panel.border = element_blank())

ggplot(bearrange2) +   geom_sf(Wisconsin2, inherit.aes = FALSE, aes(fill=NA)) +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(), # Remove major gridlines
        panel.grid.minor = element_blank(), # Remove minor gridlines
        panel.background = element_blank(),
        axis.text = element_blank(),  # Removes both x and y axis tick labels
        axis.ticks = element_blank(),# Ensure the panel background is blank/transparent
        panel.border = element_blank())

geom_label(data=siteSF,inherit.aes = FALSE, aes(x= lon, y=lat, label=1:nrow(siteSF)), # conflict with 1 and 2
             show.legend=FALSE, color="black", size=4) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

sitemap