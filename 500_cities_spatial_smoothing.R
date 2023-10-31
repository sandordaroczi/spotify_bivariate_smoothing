library(rgdal)          # For readOGR()
library(rgeos)          # For unionSpatialPolygons(), gIntersection(), gBuffer()
library(maptools)       # For unionSpatialPolygons()
library(ggplot2)        # For fortify(), ggplot()
library(readxl)         # For read_excel()
library(magrittr)       # For the pipe operator %>%
library(scales)         # For rescale()
library(dplyr)          # For inner_join(), bind_rows(),                   ### between(), mutate()
library(gridExtra)      # For grid.arrange()
library(tidyr)          # For gather()
library(devtools)       # For guide_colourbar()
library(gstat)          # For kriging()
library(sp)
library(spdep)
library(mgcv)           # For Markov Random Fields (MRF) smoothing

## Sources:
## https://rpubs.com/Earlien/Publicly_Available_Spatial_Data_Sets
## https://fromthebottomoftheheap.net/2017/10/19/first-steps-with-mrf-smooths/


####################################################################
# 1. Fill colours
####################################################################

# Fill colours for the spatial covariate(s)
Fill.colours.style1 <- c("#4dac26", "#b8e186", "#f7f7f7", "#f1b6da", "#d01c8b")

# Fill colours for the response
Fill.colours.style2 <- c("#2c7bb6", "#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c", "#d7191c")






####################################################################
# 2. Read in dataset
####################################################################

setwd("/Users/sandordaroczi/R/regression_seminar")
getwd()

city_name = "Atlanta"

cities_data <- read.csv("500_Cities__Census_Tract-level_Data__GIS_Friendly_Format___2019_release.csv", header=TRUE)
head(cities_data)

df <- subset(cities_data, PlaceName == city_name)
df$response <- df$STROKE_CrudePrev
df$x <- df$CSMOKING_CrudePrev

head(df)
n = nrow(df)
cols = colnames(df)

df$Population2010 = as.numeric(sub(",", "", df$Population2010, fixed = TRUE))



####################################################################
# 3. Import the shapefile
####################################################################

map_full <- readOGR("shapefiles/500Cities_Tracts_Clip.shp", verbose = FALSE)

# Subset the map for Boston only
map <- map_full[map_full$PlaceName == city_name,]

# Create border map from full map
N <- length(map)
N_original <- length(map)
map.border <- unionSpatialPolygons(map, IDs = rep(1, N))
map.border <- fortify(map.border)

# Fortify the map into a dataframe which can be more easily handled
map.df <- fortify(map)

# Quick plot
plot(map)
head(map)

# Calculating id's of missing records
missing.record <- which(!(unique(map$tract2010) %in% df$TractFIPS))

# Make shapefile dataframe IDs numeric
map.df$id <- as.numeric(map.df$id)
if(!all(range(map.df$id) == c(1, N))){
  map.df$id <- as.numeric(map.df$id) - min(map.df$id) + 1  # Need to add 1 so range is 1:N
  print(range(map.df$id))
}

# Transforming IDs into the range 1:N
unique_ids = unique(map.df$id)
for (i in 1:length(unique_ids)){
  map.df$id[map.df$id == unique_ids[i]] = i
}

# Filtering for missing records
missing.records <- which(!(unique(map$tract2010) %in% df$TractFIPS))
map.df <- filter(map.df, !(id %in% missing.records))

# Transforming IDs once again, could not figure out a more elegant way of handling this
unique_ids = unique(map.df$id)
for (i in 1:length(unique_ids)){
  map.df$id[map.df$id == unique_ids[i]] = i
}

N <- length(unique(map.df$id))
stopifnot(all(unique(map.df$id) == 1:N))


# Calculating average lat and long values
map.df$lat.mean <- ave(map.df$lat, map.df$id)
map.df$long.mean <- ave(map.df$long, map.df$id)




####################################################################
# 4. Creating subset of data
####################################################################

df_subset = df[,c("TractFIPS", "Population2010", "response", "x")]
head(df_subset)
n = nrow(df_subset)

colnames(df_subset) <- c("tract.fips", "population", "response", "x")

# Dataframe of the variables we need to append to the shapefile
Append <- data.frame(
  id = unique(map.df$id),
  response = df_subset$response,
  x = df_subset$x,
  population = df_subset$population,
  tract.fips = df_subset$tract.fips
)

# Merge values with shapefile dataframe
df_joined <- inner_join(map.df, Append, by = "id")
df_joined$id <- as.factor(df_joined$id)



####################################################################
# 5. EDA + Plots
####################################################################

glimpse(map)
summary(df_joined)


cor(df_joined$x, df_joined$response)
plot(df_joined$x, df_joined$response)


gg.locations <- df_joined %>% as.data.frame %>%
  ggplot(aes(long.mean, lat.mean)) + geom_point(size=1) + coord_equal() +
  ggtitle("Data locations")

# Set up basic components of the plots
gg.base <- ggplot(data = NULL, aes(x = long, y = lat, group = group)) +
  geom_polygon(data = map.border, fill = "grey30", color = "black", size = 0.8) +
  scale_x_continuous(expand = c(0.05, 0.05)) +
  scale_y_continuous(expand = c(0.05, 0.05)) +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_fill_continuous(guide = guide_colourbar()) +
  coord_fixed() # Keeps x-y ratio consistent

# Plotting the response variable
gg.response <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = response), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,max(df_joined$response))) +
  ggtitle("Prevalence of Stroke")

# Plotting the x variable
gg.x <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = x), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,max(df_joined$x))) +
  ggtitle("Prevalence of Smoking")

grid.arrange(
  gg.locations, gg.response, gg.x,
  nrow = 1
)

####################################################################
# 6. Defining neighbour structure
####################################################################


df_joined_aggr <- aggregate(cbind(long, lat, response, x, population, tract.fips) ~ id + group, data = df_joined, FUN = mean)

nb <- poly2nb(map, row.names = 1:N_original)
coords <- coordinates(map)
plot(nb, coords=coords, col="blue")
plot(map, add=TRUE)

nb <- subset(nb, !(1:length(nb) %in% missing.records))
names(nb) <- 1:length(unique(df_joined_aggr$id))



####################################################################
# 7. Running MRF smoother
####################################################################



df_joined_aggr$id <- as.factor(df_joined_aggr$id)

mrf.m1 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), sp = 1), # define MRF smooth
          data = df_joined_aggr,
          method = 'REML')

mrf.m2 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), sp = 10), # define MRF smooth
          data = df_joined_aggr,
          method = 'REML')

mrf.m3 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), sp = 100), # define MRF smooth
          data = df_joined_aggr,
          method = 'REML')

mrf.m4 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), k = 80), # define MRF smooth
              data = df_joined_aggr,
              method = 'REML')

mrf.m5 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), k = 30), # define MRF smooth
              data = df_joined_aggr,
              method = 'REML')

mrf.m6 <- gam(response ~ s(id, bs = 'mrf', xt = list(nb = nb), k = 10), # define MRF smooth
              data = df_joined_aggr,
              method = 'REML')

df_joined_aggr <- transform(df_joined_aggr, mrf.m1 = predict(mrf.m1, type = 'response'))
df_joined_aggr <- transform(df_joined_aggr, mrf.m2 = predict(mrf.m2, type = 'response'))
df_joined_aggr <- transform(df_joined_aggr, mrf.m3 = predict(mrf.m3, type = 'response'))
df_joined_aggr <- transform(df_joined_aggr, mrf.m4 = predict(mrf.m4, type = 'response'))
df_joined_aggr <- transform(df_joined_aggr, mrf.m5 = predict(mrf.m5, type = 'response'))
df_joined_aggr <- transform(df_joined_aggr, mrf.m6 = predict(mrf.m6, type = 'response'))

summary(mrf.m1)
summary(mrf.m2)
summary(mrf.m3)
summary(mrf.m4)
summary(mrf.m5)
summary(mrf.m6)

df_joined <- inner_join(df_joined, df_joined_aggr[,c("id", "mrf.m1", "mrf.m2", "mrf.m3", "mrf.m4", "mrf.m5", "mrf.m6")], by = "id")




####################################################################
# 8. Plotting results
####################################################################


gg.mrf.m1 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m1), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with sp=1")

gg.mrf.m2 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m2), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with sp=10")

gg.mrf.m3 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m3), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with sp=100")



gg.mrf.m4 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m4), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with k=80")

gg.mrf.m5 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m5), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with k=30")

gg.mrf.m6 <- gg.base + 
  geom_polygon(data = df_joined, aes(fill = mrf.m6), color = NA) +
  scale_fill_gradientn(colours = Fill.colours.style2, limits = c(0,9)) +
  ggtitle("MRF with k=10")



grid.arrange(
  gg.response, gg.mrf.m1, gg.mrf.m2, gg.mrf.m3, gg.response, gg.mrf.m4, gg.mrf.m5, gg.mrf.m6,
  nrow = 2
)













