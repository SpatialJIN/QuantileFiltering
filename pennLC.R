require(tidyverse)
require(dplyr)
require(CVXR)
require(SpatialEpi)
require(spdep)
require(sf)
require(ggplot2)
require(tibble)
require(leaflet)
require(viridis)

map <- pennLC$spatial.polygon
# plot(map)

nb <- poly2nb(map)

coord <- coordinates(map)
map$long <- coord[, 1]
map$lat <- coord[, 2]

#===============================================#

mapsf <- st_as_sf(map) %>% 
  rownames_to_column("county")

W <- data.frame(matrix(ncol=67,nrow=67))
colnames(W) <- mapsf$county
rownames(W) <- mapsf$county

for (i in 1:67) {
  for (j in nb[[i]]) {
    W[i,j] <- 1
  }
}
W[is.na(W)] <- 0
W <- as.matrix(W)
iota <- c(rep(1,67)) %>% as.matrix()
D <- W%*%iota %>%
  as.numeric() %>%
  diag()
L <- D - W

# quantile_filter

quantile_filter=function(y,Omega,lambda,tau){
  n=nrow(y)
  x=Variable(n)
  p=Problem( Minimize( 0.5*p_norm(y-x,1)+(tau-0.5)*sum(y-x)+lambda*quad_form(x,Omega) ) )
  result=solve(p)
  xhat=result$getValue(x)
  return(xhat)
}

# SIR

d <- group_by(pennLC$data, county) %>% summarize(Y = sum(cases))
d <- aggregate(
  x = pennLC$data$cases,
  by = list(county = pennLC$data$county),
  FUN = sum
)
names(d) <- c("county", "Y")
pennLC$data <- pennLC$data[order(
  pennLC$data$county,
  pennLC$data$race,
  pennLC$data$gender,
  pennLC$data$age
), ]
E <- expected(
  population = pennLC$data$population,
  cases = pennLC$data$cases, n.strata = 16
)
d$E <- E[match(d$county, unique(pennLC$data$county))]
d$SIR <- d$Y / d$E
mapsf <- left_join(mapsf, d, by = "county")

# quantile filtering for SIRs

sir <- mapsf$SIR %>% 
  as.data.frame()
rownames(sir) <- mapsf$county

xhat_1 <- mapsf$county
xhat_5 <- mapsf$county
xhat_9 <- mapsf$county

xhat_1 <- cbind(xhat_1, quantile_filter(sir, L, 4.8, 0.1))
xhat_5 <- cbind(xhat_5, quantile_filter(sir, L, 4.8, 0.5))
xhat_9 <- cbind(xhat_9, quantile_filter(sir, L, 4.8, 0.9))


colnames(xhat_1) <- c("county", "xhat.1.sir")
colnames(xhat_5) <- c("county", "xhat.5.sir")
colnames(xhat_9) <- c("county", "xhat.9.sir")

xhat_1 <- xhat_1 %>%
  as.data.frame()
xhat_5 <- xhat_5 %>%
  as.data.frame()
xhat_9 <- xhat_9 %>%
  as.data.frame()

mapqf_1 <- merge(mapsf, xhat_1, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.1.sir") %>%
  mutate(tau = 0.1)
mapqf_1$SIR_QF <- mapqf_1$SIR_QF %>% 
  as.numeric()
mapqf_5 <- merge(mapsf, xhat_5, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.5.sir") %>%
  mutate(tau = 0.5)
mapqf_5$SIR_QF <- mapqf_5$SIR_QF %>% 
  as.numeric()
mapqf_9 <- merge(mapsf, xhat_9, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.9.sir") %>%
  mutate(tau = 0.9)
mapqf_9$SIR_QF <- mapqf_9$SIR_QF %>% 
  as.numeric()
mapQF <- rbind(mapqf_1, mapqf_5, mapqf_9)



# Visualization
ggplot(mapsf) + geom_sf(aes(fill = SIR)) +
  ggtitle("The lung cancer SIRs in Pennsylvania") +
  theme_gray() +
  scale_fill_gradient2(
    midpoint = 1, low = "chartreuse4", mid = "white", high = "firebrick1"
  )

#===============================================#

ggplot(mapQF) + geom_sf(aes(fill = SIR_QF)) +
  facet_wrap(~tau, dir = "h", ncol = 3,
             labeller = label_bquote(col = tau == .(tau)))+
  ggtitle("Results with quantile filtering") +
  theme_gray() +
  scale_fill_gradient2(
    midpoint = 1, low = "chartreuse4", mid = "white", high = "firebrick1"
  )

# residual error
mapqf_1_re <- mapqf_1 %>%
  mutate(eta = .$SIR - .$SIR_QF)
mapqf_5_re <- mapqf_5 %>%
  mutate(eta = .$SIR - .$SIR_QF)
mapqf_9_re <- mapqf_9 %>%
  mutate(eta = .$SIR - .$SIR_QF)
for (i in 1:67) {
  if (mapqf_1_re$eta[i] < -1e-09) {
    mapqf_1_re$eta[i] = -1
  }else if (mapqf_1_re$eta[i] > 1e-09){
    mapqf_1_re$eta[i] = 1
  }
}
for (i in 1:67) {
  if (mapqf_5_re$eta[i] < -1e-09) {
    mapqf_5_re$eta[i] = -1
  }else if (mapqf_5_re$eta[i] > 1e-09){
    mapqf_5_re$eta[i] = 1
  }
}
for (i in 1:67) {
  if (mapqf_9_re$eta[i]< -1e-09) {
    mapqf_9_re$eta[i] = -1
  }else if (mapqf_9_re$eta[i] > 1e-09){
    mapqf_9_re$eta[i] = 1
  }
}
mapqf_re <- rbind(mapqf_1_re, mapqf_5_re, mapqf_9_re)

ggplot(mapqf_re) + geom_sf(aes(fill = eta)) +
  facet_wrap(~tau, dir = "h", ncol = 3,
             labeller = label_bquote(col = tau == .(tau)))+
  ggtitle("Graphs of eta") + theme_gray() +
  labs(fill = expression(eta)) +
  scale_fill_gradient2(
    midpoint = 0, low = "lightskyblue1", mid = "cadetblue", high = "royalblue1"
  )

# average
summary(mapqf_1$SIR_QF)
summary(mapqf_5$SIR_QF)
summary(mapqf_9$SIR_QF)


# lambda = 0.1, 1, 100

xhat_1_01 <- mapsf$county
xhat_5_01 <- mapsf$county
xhat_9_01 <- mapsf$county
xhat_1_1 <- mapsf$county
xhat_5_1 <- mapsf$county
xhat_9_1 <- mapsf$county
xhat_1_100 <- mapsf$county
xhat_5_100 <- mapsf$county
xhat_9_100 <- mapsf$county

xhat_1_01 <- cbind(xhat_1_01, quantile_filter(sir, L, 0.1, 0.1))
xhat_5_01 <- cbind(xhat_5_01, quantile_filter(sir, L, 0.1, 0.5))
xhat_9_01 <- cbind(xhat_9_01, quantile_filter(sir, L, 0.1, 0.9))
xhat_1_1 <- cbind(xhat_1_1, quantile_filter(sir, L, 1, 0.1))
xhat_5_1 <- cbind(xhat_5_1, quantile_filter(sir, L, 1, 0.5))
xhat_9_1 <- cbind(xhat_9_1, quantile_filter(sir, L, 1, 0.9))
xhat_1_100 <- cbind(xhat_1_100, quantile_filter(sir, L, 100, 0.1))
xhat_5_100 <- cbind(xhat_5_100, quantile_filter(sir, L, 100, 0.5))
xhat_9_100 <- cbind(xhat_9_100, quantile_filter(sir, L, 100, 0.9))

colnames(xhat_1_01) <- c("county", "xhat.1.01.sir")
colnames(xhat_5_01) <- c("county", "xhat.5.01.sir")
colnames(xhat_9_01) <- c("county", "xhat.9.01.sir")
colnames(xhat_1_1) <- c("county", "xhat.1.1.sir")
colnames(xhat_5_1) <- c("county", "xhat.5.1.sir")
colnames(xhat_9_1) <- c("county", "xhat.9.1.sir")
colnames(xhat_1_100) <- c("county", "xhat.1.100.sir")
colnames(xhat_5_100) <- c("county", "xhat.5.100.sir")
colnames(xhat_9_100) <- c("county", "xhat.9.100.sir")

xhat_1_01 <- xhat_1_01 %>%
  as.data.frame()
xhat_5_01 <- xhat_5_01 %>%
  as.data.frame()
xhat_9_01 <- xhat_9_01 %>%
  as.data.frame()
xhat_1_1 <- xhat_1_1 %>%
  as.data.frame()
xhat_5_1 <- xhat_5_1 %>%
  as.data.frame()
xhat_9_1 <- xhat_9_1 %>%
  as.data.frame()
xhat_1_100 <- xhat_1_100 %>%
  as.data.frame()
xhat_5_100 <- xhat_5_100 %>%
  as.data.frame()
xhat_9_100 <- xhat_9_100 %>%
  as.data.frame()

mapqf_1_01 <- merge(mapsf, xhat_1_01, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.1.01.sir") %>%
  mutate(tau = 0.1, lambda = 0.1)
mapqf_1_01$SIR_QF <- mapqf_1_01$SIR_QF %>% 
  as.numeric()
mapqf_5_01 <- merge(mapsf, xhat_5_01, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.5.01.sir") %>%
  mutate(tau = 0.5, lambda = 0.1)
mapqf_5_01$SIR_QF <- mapqf_5_01$SIR_QF %>% 
  as.numeric()
mapqf_9_01 <- merge(mapsf, xhat_9_01, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.9.01.sir") %>%
  mutate(tau = 0.9, lambda = 0.1)
mapqf_9_01$SIR_QF <- mapqf_9_01$SIR_QF %>% 
  as.numeric()
mapqf_1_1 <- merge(mapsf, xhat_1_1, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.1.1.sir") %>%
  mutate(tau = 0.1, lambda = 1)
mapqf_1_1$SIR_QF <- mapqf_1_1$SIR_QF %>% 
  as.numeric()
mapqf_5_1 <- merge(mapsf, xhat_5_1, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.5.1.sir") %>%
  mutate(tau = 0.5, lambda = 1)
mapqf_5_1$SIR_QF <- mapqf_5_1$SIR_QF %>% 
  as.numeric()
mapqf_9_1 <- merge(mapsf, xhat_9_1, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.9.1.sir") %>%
  mutate(tau = 0.9, lambda = 1)
mapqf_9_1$SIR_QF <- mapqf_9_1$SIR_QF %>% 
  as.numeric
mapqf_1_100 <- merge(mapsf, xhat_1_100, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.1.100.sir") %>%
  mutate(tau = 0.1, lambda = 100)
mapqf_1_100$SIR_QF <- mapqf_1_100$SIR_QF %>% 
  as.numeric()
mapqf_5_100 <- merge(mapsf, xhat_5_100, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.5.100.sir") %>%
  mutate(tau = 0.5, lambda = 100)
mapqf_5_100$SIR_QF <- mapqf_5_100$SIR_QF %>% 
  as.numeric()
mapqf_9_100 <- merge(mapsf, xhat_9_100, by.x = "county", by.y = "county") %>% 
  rename("SIR_QF" = "xhat.9.100.sir") %>%
  mutate(tau = 0.9, lambda = 100)
mapqf_9_100$SIR_QF <- mapqf_9_100$SIR_QF %>% 
  as.numeric()

mapQF_l <- rbind(mapqf_1_01, mapqf_5_01, mapqf_9_01,
                 mapqf_1_1, mapqf_5_1, mapqf_9_1,
                 mapqf_1_100, mapqf_5_100, mapqf_9_100)

mapQF_l$tau <- mapQF_l$tau %>% 
  as.character() %>% 
  paste("tau ==",.)

mapQF_l$lambda <- mapQF_l$lambda %>% 
  as.character() %>% 
  paste("lambda ==", .)


ggplot(mapQF_l) + geom_sf(aes(fill = SIR_QF)) +
  facet_wrap(~tau+lambda, dir = "h", ncol = 3,
             labeller = label_parsed) +
  ggtitle("Results with quantile filtering") +
  theme_gray() +
  scale_fill_gradient2(
    midpoint = 1, low = "chartreuse4", mid = "white", high = "firebrick1"
  )

# Empirical eumulative distribution function plots of residuals
xlim=c(-0.5,0.5)
ylim=c(0,1)
cols <- c("blue", "green", "red")
ltys <- c(1, 2, 3)
plot(ecdf(mapqf_1$SIR-mapqf_1$SIR_QF),verticals=T,do.points=F,col=cols[1],lty = ltys[1],xlim=xlim,ylim=ylim,xlab="",ylab="",main="");par(new=T)
plot(ecdf(mapqf_5$SIR-mapqf_5$SIR_QF),verticals=T,do.points=F,col=cols[2],lty = ltys[2],xlim=xlim,ylim=ylim,xlab="",ylab="",main="");par(new=T)
plot(ecdf(mapqf_9$SIR-mapqf_9$SIR_QF),verticals=T,do.points=F,col=cols[3],lty = ltys[3],xlim=xlim,ylim=ylim,xlab="",ylab="",main="");par(new=T)
plot(c(0,0,0),c(0.1,0.5,0.9),type="p",pch=19,cex=1,xlim=xlim,ylim=ylim,xlab="",ylab="")
abline(v=0,lty="dashed",col="gray70")
abline(h=c(0.1,0.5,0.9),lty="dashed",col="gray70")
legend("topleft", col = cols, lty = ltys, legend = c(expression(tau == 0.1),expression(tau == 0.5),expression(tau == 0.9)))
title("Empirical eumulative distribution function plots of residuals")


# Leaflet
pal4sir <- colorNumeric("viridis", domain = mapsf$SIR)
pal4tau01 <- colorNumeric("viridis", domain = mapqf_1$SIR_QF)
pal4tau05 <- colorNumeric("viridis", domain = mapqf_5$SIR_QF)
pal4tau09 <- colorNumeric("viridis", domain = mapqf_9$SIR_QF)

labels4sir <- sprintf("<strong> %s <strong/> <br/>
                      Longitude: %s <br/>
                      Latitude: %s <br/>
                      Observed: %s <br/>
                      Expected: %s <br/>
                      SIR: %s <br/>",
                      mapsf$county, mapsf$long, mapsf$lat,
                      mapsf$Y, mapsf$E, mapsf$SIR) %>%
  lapply(htmltools::HTML)

labels4tau01 <- sprintf("<strong> %s <strong/> <br/>
                          Longitude: %s <br/>
                          Latitude: %s <br/>
                          tau: %s <br/>
                          SIR_QF: %s <br/>",
                           mapqf_1$county, mapqf_1$long, mapqf_1$lat,
                           mapqf_1$tau, mapqf_1$SIR_QF) %>%
  lapply(htmltools::HTML)

labels4tau05 <- sprintf("<strong> %s <strong/> <br/>
                          Longitude: %s <br/>
                          Latitude: %s <br/>
                          tau: %s <br/>
                          SIR_QF: %s <br/>",
                           mapqf_5$county, mapqf_5$long, mapqf_5$lat,
                           mapqf_5$tau, mapqf_5$SIR_QF) %>%
  lapply(htmltools::HTML)

labels4tau09 <- sprintf("<strong> %s <strong/> <br/>
                          Longitude: %s <br/>
                          Latitude: %s <br/>
                          tau: %s <br/>
                          SIR_QF: %s <br/>",
                           mapqf_9$county, mapqf_9$long, mapqf_9$lat,
                           mapqf_9$tau, mapqf_9$SIR_QF) %>%
  lapply(htmltools::HTML)

leaflet(data = c(mapsf, mapqf_1, mapqf_5, mapqf_9)) %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles(providers$CartoDB.Positron, group = "CartoDB") %>%
  addPolygons(data = mapsf, color = "grey", weight = 1, fillColor = ~pal4sir(mapsf$SIR),
              fillOpacity = 0.25, label = labels4sir, group = "SIR") %>%
  addLegend("bottomright",
            pal = pal4sir, values = ~mapsf$SIR,
            title = "SIR", group = "SIR", position = "topleft") %>% 
  addPolygons(data = mapqf_1, color = "grey", weight = 1, fillColor = ~pal4tau01(mapqf_1$SIR_QF),
              fillOpacity = 0.25, label = labels4tau01, group = "tau01") %>%
  addLegend("bottomright",
            pal = pal4tau01, values = ~mapqf_1$SIR_QF,
            title = "SIR_QF(tau=0.1)", group = "tau01", position = "topleft") %>%
  addPolygons(data = mapqf_5, color = "grey", weight = 1, fillColor = ~pal4tau05(mapqf_5$SIR_QF),
              fillOpacity = 0.25, label = labels4tau05, group = "tau05") %>%
  addLegend("bottomright",
            pal = pal4tau05, values = ~mapqf_5$SIR_QF,
            title = "SIR_QF(tau=0.5)", group = "tau05", position = "topleft") %>%
  addPolygons(data = mapqf_9, color = "grey", weight = 1, fillColor = ~pal4tau09(mapqf_9$SIR_QF),
              fillOpacity = 0.25, label = labels4tau09, group = "tau09") %>%
  addLegend("bottomright",
            pal = pal4tau09, values = ~mapqf_9$SIR_QF,
            title = "SIR_QF(tau=0.9)", group = "tau09", position = "topleft") %>%
  addScaleBar(position = c("bottomleft")) %>%
  addLayersControl(
    baseGroups = c("OpenStreetMap", "CartoDB"),
    overlayGroups = c("SIR", "tau01", "tau05", "tau09")
  )
