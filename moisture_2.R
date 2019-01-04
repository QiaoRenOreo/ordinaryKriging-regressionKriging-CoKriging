library(gstat)
library(maptools)
library(rgdal)
library(sp)
library(e1071)

# step 1. read csv file
getwd()
setwd("D:/Study/Module12/Project/BBdata/4ManderDataset")
d <- import("mander_MD30zeroRemoved.csv")
# d <- read.csv(file="ManderMHNew.csv", header=TRUE, sep=",")
coordinates(d) <- ~X + Y
head(d)

# step 2. check do we need to do log transformation on MHW? do log 3 times on MHW
#proj4string(d) <- CRS("+init=epsg:3035 +units=km")
plot(d$MHW, axes=T)
hist(d$MHW)
lMHW <- log(d$MHW) # firt log
hist(lMHW)
skewness(lMHW)
l2MHW <- log(lMHW) # second log
skewness(l2MHW)
hist(l2MHW)
l3MHW <- log(l2MHW) #third log
skewness(l3MHW)
hist(l3MHW)

# step 3. check the correlation between MHW and MD30
x <- cor(d$MHW, d$MoistDefAvg30)

# step 4. define the subsample. 
# subsample will be used in validation. they don't participate prediction
subsamp <- c(1, 6, 10, 11, 19, 33, 36, 37, 41, 42, 47, 49, 53, 54,
             57, 66, 67, 68, 72, 86, 89, 98, 100, 127, 129, 133, 144, 148, 150, 155)
d.aa <- d[subsamp,]
names(d.aa)
names(d.aa)[2] <- "MD30"
d.c <- d[-subsamp,] # d.c is original sample exclude subsample
class(d)
###########################################################

# step 5. create variogram of MD30
MD30.ev <- variogram(log(MoistDefAvg30)~ 1, data=d, cutoff = 1200, width=70)
plot(MD30.ev, main="semivariance of MD30")
#MD30.ev


# step 5. create variogram on 3rd-log transformed MHW

# MHW1: default variogram
MHW1.ev <- variogram(log(log(log(MHW)))~ 1, data=d.c) 
MHW1.ev
head(d.c)

# MHW2: self-defined cutoff and width
MHW2.ev <- variogram(log(log(log(MHW)))~ 1, data=d.c, cutoff = 1200, width=70)
plot(MHW2.ev)
MHW2.ev

####################################################################
# step 6. use a model to fit MD30 variogram, then check SSE
MD30.mv <- fit.variogram(MD30.ev, model=vgm(0.5,"Exp", 500, 0.5))
X11()
plot(MD30.ev, model=MD30.mv, main="variogram of MD30 with a fitting model")
str(MD30.mv)
sse_MD30=attr(MD30.mv, "SSErr")
sse_MD30

# step 6. use a model to fit MHW1 and MHW2 variogram, then check SSE


# MHW1: default variogram
plot(MHW1.ev)
MHW1.mv <- fit.variogram(MHW1.ev, model=vgm(0.5,"Exp", 500, 0.5))
plot(MHW1.ev, model=MHW1.mv, main="variogram of MHW1 with a fitting model")
str(MHW1.mv)
sse_MHW1=attr(MHW1.mv, "SSErr")
sse_MHW1 # =1.909377e-08

# MHW2: variogram with self-defined cutoff, width
plot(MHW2.ev)
MHW2.mv <- fit.variogram(MHW2.ev, model=vgm(0.5,"Exp", 500, 0.5))
plot(MHW2.ev, model=MHW2.mv, main="variogram of MHW2 with a fitting model")
str(MHW2.mv)
sse_MHW2=attr(MHW2.mv, "SSErr")
sse_MHW2 # =1.909377e-08

#########################################
#step 7. cokriging
g <- gstat(NULL, "ln.MD30", log(MoistDefAvg30)~ 1, d.c)
g <- gstat(g, "ln.MHW", log(log(log(MHW)))~ 1, d.c)
v <- variogram(g,cutoff = 1000, width=70)
plot(v)
v
str(v)
g <- gstat(g, model=vgm(1, "Exp", 500, 1), fill.all=TRUE)
g
g.fit <- fit.lmc(v, g)
g <- g.fit
plot(v, g.fit)
str(g)
rm(g)
##########################################
# step 8. create grid
xy <- expand.grid(x=seq(249000, 253500, by=50), y=seq(493500,496500, by=50))
xys <- SpatialPoints(xy)
gridded(xys) <- TRUE
#proj4string(xys) <- CRS("+init=epsg:3035 +units=km")
points(as.data.frame(d)[,1:2], col=4, pch=19)
plot(xys)

str(xys)

# step 9. perform co-kriging
str(d)
class(ln.MHW)
x <- predict(g, xys)

str(g)

# step 10. MD30 ordinary kriging: prediction and variance
aa2 <-krige(log(MoistDefAvg30)~1,d.c,newdata=xys, model=MD30.mv)
X11()
spplot(aa2, "var1.pred",sp.layout=list("sp.points", pch="+", d.c), 
       scales=list(draw=TRUE),main ="Kriging Predictions of MD30")
X11()
spplot(aa2, "var1.var",sp.layout=list("sp.points", pch="+", d.c), 
       scales=list(draw=TRUE),main ="Kriging Variance of MD30")


# step 11. cross validation on MD30
aa2.cv <- krige.cv((MoistDefAvg30)~1, d, model=MD30.mv)  
str(aa2.cv)
# head(aa2)
me_ok <- sum(aa2.cv$residual) / length(aa2.cv$residual)
mse_ok <- sum(aa2.cv$residual^2) / length(aa2.cv$residual)
rmse_ok <- sqrt(mse_ok)

me_ok
mse_ok
rmse_ok


# step 10. MHW ordinary kriging: prediction and variance
aa3 <-krige(log(log(log(MHW)))~ 1,d,newdata=xys, model=MHW2.mv)
spplot(aa3, "var1.pred",main ="Kriged predictions of MHW2")
spplot(aa3, "var1.var",main ="Kriging variance of MHW2")

# step 11. cross validation on log transformed MHW
aa3.cv <- krige.cv(log(log(log(MHW)))~ 1, d, model=MHW2.mv)
str(aa3.cv)
me_ok2 <- sum(aa3.cv$residual) / length(aa3.cv$residual)
mse_ok2 <- sum(aa3.cv$residual^2) / length(aa3.cv$residual)
rmse_ok2 <- sqrt(mse_ok1)

me_ok1
me_ok2

mse_ok1
mse_ok2

rmse_ok1
rmse_ok2


#########################################
# step 12. ordinary cokriging and accuracy assesment
aa1  <- predict(g, d.aa)
MD30.err1 <-log(d.aa$MD30)-aa1$ln.MD30.pred
ME<-sum(MD30.err1)/length(MD30.err1)
RMSE<-sum(MD30.err1^2)/length(MD30.err1)
#aa2 <-krige(d$MoistDefAvg30~1,d,newdata=xys, model=MD30.mv)
#aa3 <-krige(d$MoistDefAvg30~d$MHW, data=d, newdata=xys, model=MD30.mv)
#spplot(aa2)


# step 13. cross validation for cokring 
v_ck <-gstat.cv(g)
length(d$MoistDefAvg30)
length(d$MHW)
length(d.c)
summary(v_ck)
x_RMSE=sqrt(sum(v_ck$residual^2)/(length(v_ck$residual))) 
rm(x)


nrow(xys)