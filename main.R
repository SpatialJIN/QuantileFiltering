require(tidyverse)
require(dplyr)
require(CVXR)
require(sf)
require(spdep)
require(rgdal)
require(SpatialEpi)
require(ggplot2)
require(CVXR)

# Lung cancer in Ohio
dohio <- read.csv("./Data/dataohiocomplete.csv")
head(dohio)

map <- readOGR("./Data/fe_2007_39_county/fe_2007_39_county.shp", verbose = FALSE)
plot(map)

d <- aggregate(
  x = dohio$y,
  by = list(county = dohio$NAME, year = dohio$year),
  FUN = sum
)
names(d) <- c("county", "year", "Y")
# head(d)

dohio <- dohio[order(
  dohio$county,
  dohio$year,
  dohio$gender,
  dohio$race
), ]

n.strata <- 4
E <- expected(
  population = dohio$n,
  cases = dohio$y,
  n.strata = n.strata
)

nyears <- length(unique(dohio$year))
countiesE <- rep(unique(dohio$NAME),
                 each = nyears)

ncounties <- length(unique(dohio$NAME))
yearsE <- rep(unique(dohio$year),
              times = ncounties)

dE <- data.frame(county = countiesE, year = yearsE, E = E)

d <- d %>%
  left_join(dE, by = c("county", "year")) %>%
  mutate("SIR" = .$Y / .$E)

dw <- reshape(d,
              timevar = "year",
              idvar = "county",
              direction = "wide"
)

map <- merge(map, dw, by.x = "NAME", by.y = "county")
map_sf <- st_as_sf(map)
map_sf <- gather(map_sf, year, SIR, paste0("SIR.", 1968:1988))
map_sf$year <- as.integer(substring(map_sf$year, 5, 8))

ggplot(map_sf) + geom_sf(aes(fill = SIR)) +
  facet_wrap(~year, dir = "h", ncol = 7) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )


#===============================================#

# nb <- poly2nb(map)
# head(nb)
# lapply(nb, write, "./Data/nb.txt", append=TRUE, ncolumns=1000) # save nb

W <- data.frame(matrix(ncol=88,nrow=88))
colnames(W) <- map$NAME
rownames(W) <- map$NAME

for (i in 1:88) {
  for (j in nb[[i]]) {
    W[i,j] <- 1
  }
}
W[is.na(W)] <- 0
W <- as.matrix(W)
iota <- c(rep(1,88)) %>% as.matrix()
D <- W%*%iota %>%
  as.numeric() %>%
  diag()
L <- D-W
# write.csv(D, "./Data/AM.csv")

# quantile_filter

quantile_filter=function(y,Omega,lambda,tau){
  n=nrow(y)
  x=Variable(n)
  p=Problem( Minimize( 0.5*p_norm(y-x,1)+(tau-0.5)*sum(y-x)+lambda*quad_form(x,Omega) ) )
  result=solve(p)
  xhat=result$getValue(x)
  return(xhat)
}

SIR <- data.frame()
SIR <-select(map@data, starts_with("SIR."))
rownames(SIR) <- map$NAME

xhat_1 <- map$NAME
xhat_5 <- map$NAME
xhat_9 <- map$NAME
for (i in 1:21) {
  xhat_1 <- cbind(xhat_1, quantile_filter(SIR[i], L, 1, 0.1))
  xhat_5 <- cbind(xhat_5, quantile_filter(SIR[i], L, 1, 0.5))
  xhat_9 <- cbind(xhat_9, quantile_filter(SIR[i], L, 1, 0.9))
}
colnames(xhat_1) <- c("county","xhat.1.SIR.1968","xhat.1.SIR.1969","xhat.1.SIR.1970","xhat.1.SIR.1971","xhat.1.SIR.1972","xhat.1.SIR.1973",
                      "xhat.1.SIR.1974","xhat.1.SIR.1975","xhat.1.SIR.1976","xhat.1.SIR.1977","xhat.1.SIR.1978","xhat.1.SIR.1979",
                      "xhat.1.SIR.1980","xhat.1.SIR.1981","xhat.1.SIR.1982","xhat.1.SIR.1983","xhat.1.SIR.1984","xhat.1.SIR.1985",
                      "xhat.1.SIR.1986","xhat.1.SIR.1987","xhat.1.SIR.1988")
colnames(xhat_5) <- c("county","xhat.5.SIR.1968","xhat.5.SIR.1969","xhat.5.SIR.1970","xhat.5.SIR.1971","xhat.5.SIR.1972","xhat.5.SIR.1973",
                      "xhat.5.SIR.1974","xhat.5.SIR.1975","xhat.5.SIR.1976","xhat.5.SIR.1977","xhat.5.SIR.1978","xhat.5.SIR.1979",
                      "xhat.5.SIR.1980","xhat.5.SIR.1981","xhat.5.SIR.1982","xhat.5.SIR.1983","xhat.5.SIR.1984","xhat.5.SIR.1985",
                      "xhat.5.SIR.1986","xhat.5.SIR.1987","xhat.5.SIR.1988")
colnames(xhat_9) <- c("county","xhat.9.SIR.1968","xhat.9.SIR.1969","xhat.9.SIR.1970","xhat.9.SIR.1971","xhat.9.SIR.1972","xhat.9.SIR.1973",
                      "xhat.9.SIR.1974","xhat.9.SIR.1975","xhat.9.SIR.1976","xhat.9.SIR.1977","xhat.9.SIR.1978","xhat.9.SIR.1979",
                      "xhat.9.SIR.1980","xhat.9.SIR.1981","xhat.9.SIR.1982","xhat.9.SIR.1983","xhat.9.SIR.1984","xhat.9.SIR.1985",
                      "xhat.9.SIR.1986","xhat.9.SIR.1987","xhat.9.SIR.1988")
xhat_1 <- xhat_1 %>%
  as.data.frame()
xhat_5 <- xhat_5 %>%
  as.data.frame()
xhat_9 <- xhat_9 %>%
  as.data.frame()


map_qf <- readOGR("./Data/fe_2007_39_county/fe_2007_39_county.shp", verbose = FALSE)

# xhat_1
map_xhat_1 <- merge(map_qf, xhat_1, by.x = "NAME", by.y = "county")
map_xhat_1_sf <- st_as_sf(map_xhat_1)
map_xhat_1_sf <- gather(map_xhat_1_sf, year, xhat.1.SIR, paste0("xhat.1.SIR.", 1968:1988))
map_xhat_1_sf$year <- as.integer(substring(map_xhat_1_sf$year, 12, 15))
map_xhat_1_sf$xhat.1.SIR <- map_xhat_1_sf$xhat.1.SIR %>%
  as.numeric()

ggplot(map_xhat_1_sf) + geom_sf(aes(fill = xhat.1.SIR)) +
  facet_wrap(~year, dir = "h", ncol = 7) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )

# xhat_5

map_xhat_5 <- merge(map_qf, xhat_5, by.x = "NAME", by.y = "county")
map_xhat_5_sf <- st_as_sf(map_xhat_5)
map_xhat_5_sf <- gather(map_xhat_5_sf, year, xhat.5.SIR, paste0("xhat.5.SIR.", 1968:1988))
map_xhat_5_sf$year <- as.integer(substring(map_xhat_5_sf$year, 12, 15))
map_xhat_5_sf$xhat.5.SIR <- map_xhat_5_sf$xhat.5.SIR %>%
  as.numeric()

ggplot(map_xhat_5_sf) + geom_sf(aes(fill = xhat.5.SIR)) +
  facet_wrap(~year, dir = "h", ncol = 7) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )

# xhat_9

map_xhat_9 <- merge(map_qf, xhat_9, by.x = "NAME", by.y = "county")
map_xhat_9_sf <- st_as_sf(map_xhat_9)
map_xhat_9_sf <- gather(map_xhat_9_sf, year, xhat.9.SIR, paste0("xhat.9.SIR.", 1968:1988))
map_xhat_9_sf$year <- as.integer(substring(map_xhat_9_sf$year, 12, 15))
map_xhat_9_sf$xhat.9.SIR <- map_xhat_9_sf$xhat.9.SIR %>%
  as.numeric()

ggplot(map_xhat_9_sf) + geom_sf(aes(fill = xhat.9.SIR)) +
  facet_wrap(~year, dir = "h", ncol = 7) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red"
  )





# check
if(FALSE){
df <- data.frame(county = map$NAME, neigh = rep(0, length(map)))
rownames(df) <- map$NAME
map <- SpatialPolygonsDataFrame(map, df, match.ID = FALSE)
map$neigh[nb[[2]]] <- 1
map$neigh[nb[[44]]] <- 1
map$neigh[nb[[58]]] <- 1
coord <- coordinates(map)
map$long <- coord[, 1]
map$lat <- coord[, 2]
map$ID <- 1:dim(map@data)[1]
mapsf <- st_as_sf(map)
ggplot(mapsf) + geom_sf(aes(fill = as.factor(neigh))) +
  geom_text(aes(long, lat, label = mapsf$county), color = "white") +
  theme_bw() + guides(fill = FALSE)
}

# 1988's data only

map88_1_sf <- map_xhat_1_sf %>%
  filter(year == 1988) %>%
  rename("SIR" = "xhat.1.SIR") %>%
  mutate(tau = 0.1)
map88_5_sf <- map_xhat_5_sf %>%
  filter(year == 1988) %>%
  rename("SIR" = "xhat.5.SIR") %>%
  mutate(tau = 0.5)
map88_9_sf <- map_xhat_9_sf %>%
  filter(year == 1988) %>%
  rename("SIR" = "xhat.9.SIR") %>%
  mutate(tau = 0.9)

map88_sf <- rbind(map88_1_sf, map88_5_sf, map88_9_sf)
Att.labs <- c("tau = 0.1",
              "tau = 0.5",
              "tau = 0.9")
names(Att.labs) <- c("0.1",
                     "0.5",
                     "0.9")

ggplot(map88_sf) + geom_sf(aes(fill = SIR)) +
  facet_wrap(~tau, dir = "h", ncol = 3,
             labeller = labeller(tau = Att.labs)) +
  ggtitle("1988's SIR of Lung cancer in Ohio") + theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "chartreuse4", mid = "white", high = "firebrick1"
  )

# normal graph of 1988
map88_lc_sf <- map_sf %>%
  filter(year == 1988)
ggplot(map_sf) + geom_sf(aes(fill = SIR)) +
  ggtitle("1988's SIR of Lung cancer in Ohio") + theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "chartreuse4", mid = "white", high = "firebrick1"
  )

# residual error
map88_1_re <- map88_1_sf %>%
  mutate(ita = map88_lc_sf$SIR - .$SIR)
map88_5_re <- map88_5_sf %>%
  mutate(ita = map88_lc_sf$SIR - .$SIR)
map88_9_re <- map88_9_sf %>%
  mutate(ita = map88_lc_sf$SIR - .$SIR)
for (i in 1:88) {
  if (map88_1_re$ita[i]<0) {
    map88_1_re$ita[i] = -1
  }else if (map88_1_re$ita[i]>0){
    map88_1_re$ita[i] = 1
  }
}
for (i in 1:88) {
  if (map88_5_re$ita[i]<0) {
    map88_5_re$ita[i] = -1
  }else if (map88_5_re$ita[i]>0){
    map88_5_re$ita[i] = 1
  }
}
for (i in 1:88) {
  if (map88_9_re$ita[i]<0) {
    map88_9_re$ita[i] = -1
  }else if (map88_9_re$ita[i]>0){
    map88_9_re$ita[i] = 1
  }
}
map88_re <- rbind(map88_1_re, map88_5_re, map88_9_re)
Att.labs <- c("tau = 0.1",
              "tau = 0.5",
              "tau = 0.9")
names(Att.labs) <- c("0.1",
                     "0.5",
                     "0.9")

ggplot(map88_re) + geom_sf(aes(fill = ita)) +
  facet_wrap(~tau, dir = "h", ncol = 3,
             labeller = labeller(tau = Att.labs)) +
  ggtitle("Ita") + theme_void() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 0, low = "lightskyblue1", mid = "cadetblue", high = "royalblue1"
  )


# average
summary(map88_1_sf$SIR)
summary(map88_5_sf$SIR)
summary(map88_9_sf$SIR)
