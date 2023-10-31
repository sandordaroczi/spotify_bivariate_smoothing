library(tidyverse)
library(skimr)
library(ggplot2)
library(tidyr)
library(gridExtra)     # For grid.arrange()
library(caTools)       # For sample.split
library(caret)         # For RMSE
library(mgcv)          # For gam() function
library(npreg)         # For sm() function
library(cowplot)

####################################################################
# 1. Read in dataset
####################################################################

setwd("/Users/sandordaroczi/R/regression_seminar")
getwd()

df <- read.csv("songs_new.csv")

head(df)
n = nrow(df)
cols = colnames(df)


####################################################################
# 2. Data preprocessing
####################################################################

genre_names = cols[18:28]
genre <- rep(0,n)

for (j in genre_names){
  for (i in 1:n){
    if (df[i,j] == 1){
      genre[i] = j
    }
  }
}

df["genre"] = genre
df_new = df[,-18:-28]
df_new["genre"] = factor(genre)
df_new["track.mode"] = factor(df$track.mode)
df_new["track.explicit"] = factor(df$track.explicit)
df_new["track.key"] = factor(df$track.key)
df_new["time_signature"] = factor(df$time_signature)
new_order = c(17,1,2,4,6:11,14,16,18,3,5,12,15,19)
df_new = df_new[,new_order]
colnames(df_new) = c("popularity", "danceability", "energy", "loudness",
                 "speechiness", "acousticness", "instrumentalness",
                 "liveness", "valence", "tempo", "duration_ms", "number",
                 "num_available_markets", "key", "mode", "time_signature",
                 "explicit", "genre")
df <- df_new


# Removing outliers: in total, 12 observations are removed
dim(df)
df <- df[df$loudness > min(df$loudness),]
df <- df[df$speechiness > min(df$speechiness),]
df <- df[df$danceability > min(df$danceability),]
df <- df[df$tempo > min(df$tempo),]
df <- df[df$tempo < max(df$tempo),]
dim(df)

df <- df[df$valence > max(head(sort(df$valence), 3)),]
dim(df)

df <- df[df$duration_ms < min(tail(sort(df$duration_ms), 4)),]
dim(df)



####################################################################
# 3. Plotting scatter plots agains continuous covariates
####################################################################

df[1:5] %>%
  gather(-popularity, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = popularity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

df[c(1,6:9)] %>%
  gather(-popularity, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = popularity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

df[c(1,10:13)] %>%
  gather(-popularity, key = "var", value = "value") %>% 
  ggplot(aes(x = value, y = popularity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~ var, scales = "free") +
  theme_bw()

####################################################################
# 4. Full model
####################################################################

# Full model, with credits to Brian Zeiser Dietrich
full_model <- lm(popularity ~ . 
                 + poly(loudness,2) + poly(energy,2) 
                 + poly(duration_ms,2) + poly(num_available_markets,2)
                 + mode * key + explicit
                 + explicit * key
                 - loudness - energy - duration_ms - num_available_markets, 
                 data = df)
full <- summary(full_model)
full



####################################################################
# 5. Choosing covariates
####################################################################

# Possible candidates for x and y variables with nonlinear effect on response
# Danceability
# Acousticness
# Energy
# Loudness
# Duration_ms

# Choosing x and y regressor variables
df$x <- df$danceability
df$y <- df$energy

# Scale values
standardize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

df$x <- standardize(df$x)
df$y <- standardize(df$y)

# Choosing z response variable
df$z <- df$popularity
df <- df[,c("x", "y", "z")]




####################################################################
# 5. Linear models
####################################################################


# Small model, with two covariates
bivariate_model <- lm(z ~ x + y, data = df)
bivariate <- summary(bivariate_model)
bivariate

# Small model, with two covariates, including interaction
bivariate_model.inter <- lm(z ~ x * y, data = df)
bivariate.inter <- summary(bivariate_model.inter)
bivariate.inter


####################################################################
# 6. Fitting bivariate smoothing models
####################################################################



# Fitting with gam using te (tensor product smooths)
model.te1 <- gam(z ~ te(x, y, bs = 'ps', sp = c(0.05, 0.05)),
                 family=gaussian,
                 method = "REML",
                 data=df)

model.te2 <- gam(z ~ te(x, y, bs = 'ps', sp = c(1, 1)),
                family=gaussian,
                method = "REML",
                data=df)

model.te3 <- gam(z ~ te(x, y, bs = 'ps', sp = c(10, 10)),
                 family=gaussian,
                 method = "REML",
                 data=df)

df$pred.te1 <- model.te1 %>%
  predict(df)
df$pred.te2 <- model.te2 %>%
  predict(df)
df$pred.te3 <- model.te3 %>%
  predict(df)


####################################################################
# 7. Predicting on the unit grid
####################################################################

newd <- crossing(x = seq(0,1, length = 100), y = seq(0,1, length = 100))
newd <- mutate(
  newd,
  pred.te1      = predict(model.te1, newd),
  pred.te2      = predict(model.te2, newd),
  pred.te3      = predict(model.te3, newd),
  pred.lm       = predict(bivariate_model, newdata = newd),
  pred.lm.inter = predict(bivariate_model.inter, newdata = newd)
  )


####################################################################
# 8. Plotting results
####################################################################

alpha <- 1
size <- 1

gg.z <- ggplot(df, aes(x, y, color=z)) +
  geom_point(alpha=alpha, size=size) +
  scale_color_gradient(low = "white", high = "red") +
  theme_bw() +
  labs(title = "popularity")

gg.pred.te1 <- ggplot(newd, aes(x = x, y = y, fill = pred.te1)) +
  geom_raster() +
  geom_contour(aes(z = pred.te1)) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  labs(title = "pred.te1", subtitle = "sp = 0.05")

gg.pred.te2 <- ggplot(newd, aes(x = x, y = y, fill = pred.te2)) +
  geom_raster() +
  geom_contour(aes(z = pred.te2)) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  labs(title = "pred.te2", subtitle = "sp = 1")

gg.pred.te3 <- ggplot(newd, aes(x = x, y = y, fill = pred.te3)) +
  geom_raster() +
  geom_contour(aes(z = pred.te3)) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  labs(title = "pred.te3", subtitle = "sp = 10")

gg.pred.lm <- ggplot(newd, aes(x = x, y = y, fill = pred.lm)) +
  geom_raster() +
  geom_contour(aes(z = pred.lm)) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  labs(title = "pred.lm")

gg.pred.lm.inter <- ggplot(newd, aes(x = x, y = y, fill = pred.lm.inter)) +
  geom_raster() +
  geom_contour(aes(z = pred.lm.inter)) +
  scale_fill_distiller(palette = "RdBu", type = "div") +
  labs(title = "pred.lm.inter")

grid.arrange(
  gg.z, gg.pred.lm, gg.pred.lm.inter, gg.pred.te1, gg.pred.te2, gg.pred.te3,
  nrow = 2
)

grid.arrange(
  gg.z, gg.pred.lm.inter,
  nrow = 1
)

grid.arrange(
  gg.z, gg.pred.te1,
  nrow = 1
)

grid.arrange(
  gg.z, gg.pred.te2,
  nrow = 1
)

grid.arrange(
  gg.z, gg.pred.te3,
  nrow = 1
)

plot_ly(z = ~xtabs(pred.lm ~ x + y, data = newd)) %>%
  add_surface()  %>%
  layout(
    title = "pred.lm",
    scene = list(zaxis = list(title = "pred.lm")))

plot_ly(z = ~xtabs(pred.lm.inter ~ x + y, data = newd)) %>%
  add_surface()  %>%
  layout(
    title = "pred.lm.inter",
    scene = list(zaxis = list(title = "pred.lm.inter")))

plot_ly(z = ~xtabs(pred.te1 ~ x + y, data = newd)) %>%
  add_surface()  %>%
  layout(
    title = "pred.te1",
    scene = list(zaxis = list(title = "pred.te1")))

plot_ly(z = ~xtabs(pred.te2 ~ x + y, data = newd)) %>%
  add_surface()  %>%
  layout(
    title = "pred.te2",
    scene = list(zaxis = list(title = "pred.te2")))

plot_ly(z = ~xtabs(pred.te3 ~ x + y, data = newd)) %>%
  add_surface()  %>%
  layout(
    title = "pred.te3 ~ x + y, data = newd",
    scene = list(zaxis = list(title = "pred.te3 ~ x + y, data = newd")))

