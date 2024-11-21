library(sf)
library(terra)
library(tidyverse)

# ---- load dataset x-y locations ----
datapoints <- read.csv('AgeSt_mapping/dataset_locations_updated.csv')


# specify a crs
reg.proj <- st_crs(4326)

## convert df to points geom

# # clim points
# clm.point <- climpoints %>%
#   st_as_sf(., coords = c('long', 'lat'), crs = reg.proj)

# # jitter coords for overlapping points
# jtr.point <- datapoints %>% 
#   filter(Notes == "Locations are close or in distinguishable from each other") %>%
#   st_as_sf(., coords = c('long', 'lat'), crs = reg.proj) 

## 3 types of jitter/vs non-jitter needed

# jitter coords for Barger datasets to a small extent
lo.jtr.point <- datapoints %>% 
  filter(Notes == "Locations are close") %>% 
  distinct(.,mapping_loc, .keep_all = TRUE) %>% # remove duplicated sites
  st_as_sf(., coords = c('long', 'lat'), crs = reg.proj) %>% 
  st_jitter(., .15)

# no jitter for non overlapping datasets
non.jtr.point <- datapoints %>% 
  distinct(.,mapping_loc, .keep_all = TRUE) %>% # remove duplicated sites
  filter(Notes == "") %>%
  st_as_sf(., coords = c('long', 'lat'), crs = reg.proj)

# manually Jittered datasets
amt <- .17

lat.jitter <- c(0,amt,0,0,amt,-amt*1.2,amt,-amt,0,-amt,amt*1.5)
long.jitter <- c(0,0,amt,0,0,amt,amt,0,-amt,-amt,-amt*.5)

jtr.point <- datapoints %>% 
  distinct(.,mapping_loc, .keep_all = TRUE) %>% # remove duplicated sites
  filter(Notes == "Locations are close or in distinguishable from each other") %>% 
  mutate(lat = lat + lat.jitter) %>% 
  mutate(long = long + long.jitter) %>% 
  st_as_sf(., coords = c('long', 'lat'), crs = reg.proj) %>%
  st_jitter(., .05)

#
# ---- load state boundaries ----
states <- st_read('AgeSt_mapping/states_outline/tl_2012_us_state.shp') %>% filter(STUSPS=='CA'|STUSPS=='OR'|STUSPS=='NV'|STUSPS=='UT'|STUSPS=='CO'|STUSPS=='NM'|STUSPS=='AZ'|STUSPS=='ID'|STUSPS=='WY') %>% 
  st_transform(., crs = 4326)

# create area boundary to clip PJ distribution to
st_simp <- st_union(states, by_feature = FALSE) %>% st_sf()

# ---- load PJ species distributions ----
pied <- st_read('AgeSt_mapping/PiEd_distribution/pinuedul.shp') %>% 
  st_transform(., crs = 4326) %>% 
  st_intersection(.,st_simp) %>%
  filter(CODE == 1)

pimo <- st_read('AgeSt_mapping/PiMo_distribution/pinumono.shp') %>% 
  st_transform(., crs = 4326) %>% 
  st_intersection(.,st_simp) %>%
  filter(CODE == 1)

# Not including Rocky mountain juniper
# jusc <- st_read('JuSc_distribution/juniscop.shp') %>% 
#   st_transform(., crs = 4326) %>% 
#   st_intersection(.,st_simp) %>%
#   filter(CODE == 1)

juoc <- st_read('AgeSt_mapping/JuOc_distribution/juniocci.shp') %>% 
  st_transform(., crs = 4326) %>% 
  st_intersection(.,st_simp) %>%
  filter(CODE == 1)

juos <- st_read('AgeSt_mapping/JuOs_distribution/junioste.shp') %>% 
  st_transform(., crs = 4326) %>% 
  st_intersection(.,st_simp) %>%
  filter(CODE == 1)

jumo <- st_read('AgeSt_mapping/JuMo_distribution/junimono.shp') 
st_crs(jumo) <- 4267
jumo <- jumo %>% st_transform(., crs = 4326) %>% 
  st_intersection(.,st_simp) %>%
  filter(CODE == 1)

# ---- Plot map of study locations on top of PJ distribution ----

col1 = c('PJ Distribution' = '#66c2a4')
col2 = c('Study Locations' = 'black')

map<-ggplot(states) + 
  geom_sf(fill = '#f0f0f0', col = 'grey', size = .3) + 
  
  # pj distribution
  geom_sf(data = pied, aes(fill = 'PJ Distribution'), col = 'NA') +
  geom_sf(data = pimo, aes(fill = 'PJ Distribution'), col = 'NA') +
  geom_sf(data = juoc, aes(fill = 'PJ Distribution'), col = 'NA') +
  geom_sf(data = juos, aes(fill = 'PJ Distribution'), col = 'NA') +
  geom_sf(data = jumo, aes(fill = 'PJ Distribution'), col = 'NA') +
  
  # # clim data
  # geom_sf(data = clm.point, aes(shape = 'Temperature Data'), col = 'red', fill = 'red', alpha = 0.6, size = .7) +
  # 
  # study locs
  geom_sf(data = jtr.point, aes(col = 'Study Locations'), pch = 1, size = 1) +
  geom_sf(data = lo.jtr.point, aes(col = 'Study Locations'), pch = 1, size = 1) +
  geom_sf(data = non.jtr.point, aes(col = 'Study Locations'), pch = 1, size = 1) +
  
  #scale_shape_manual(name = '', values = 23, labels = 'Temperature Reconstructions') +
  scale_color_manual(name = '', values = col2) +
  scale_fill_manual(name = '', values = col1) +
  
  theme(
    text = element_text(size = 7.5),
    legend.justification = c("right", "top"),
    legend.position = c(.99, .99),
    legend.title=element_blank(),
    legend.spacing.y = unit(-.2, 'cm'),
    legend.box.margin = margin(4,4,4,4),
    legend.key=element_blank(),
    legend.key.size = unit(.3, 'cm'),
    legend.text=element_text(size=5),
    panel.background = element_rect(linetype = "solid",fill = NA),
    panel.border = element_rect(linetype = "solid", fill = NA),
    panel.grid.major = element_line(colour = "grey", size = .3)
        )+geom_text(x = 125, y = 45
                    
                    , label = "C"
                    , color = "black"
                    , size=8)
map

#ggsave("study_locations_mapped_v3.png", width = 4, height = 5, dpi = 500)

## misc
# plot each different species after different color to check extent
ggplot(states) + 
  geom_sf(size = .3) + 
  geom_sf(data = pied, fill = 'blue') +
  geom_sf(data = pimo, fill = 'red') +
  geom_sf(data = juoc, fill = 'green') +
  geom_sf(data = juos, fill = 'purple') +
  geom_sf(data = jumo, fill = 'orange')
