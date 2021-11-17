# library(OpenML,farff)
# freMTPL2freq <- getOMLDataSet(data.id = 41214)$data
library (CASdatasets)
library(gamlss)
library(fitdistrplus)
library(nnet)
library(mboost)
#library(ModTools)
library(xgboost)
#install.packages("remotes")
#remotes::install_github("f1kidd/fmlogit")
#library(fmlogit)
data (freMTPL2freq)
dat <- freMTPL2freq
dat$VehGas <- factor (dat$VehGas)
dat$ClaimNb <- as.numeric(dat$ClaimNb)
data (freMTPL2sev)
sev <- freMTPL2sev
sev$ClaimNb0 <- 1
dat0 <-
  aggregate (sev , by = list (IDpol = sev$IDpol), FUN = sum)[c (1 , 3:4)]
names (dat0)[2] <- "ClaimTotal"
dat <- merge (x = dat ,
              y = dat0 ,
              by = "IDpol" ,
              all.x = TRUE)
dat [is.na(dat)] <- 0
dat <- dat [which (dat$ClaimNb0 <= 5) , ]
cor(dat$ClaimNb, dat$ClaimNb0)
dat$ClaimNb <- dat$ClaimNb0
dat <- dat[, -14]
dat$Exposure <- pmin (dat$Exposure , 1)
sev <- sev [which (sev$IDpol %in% dat$IDpol), c (1 , 2)]
dat$VehBrand <- factor (
  dat$VehBrand,
  levels = c (
    "B1" ,
    "B2" ,
    "B3" ,
    "B4" ,
    "B5" ,
    "B6",
    "B10" ,
    "B11" ,
    "B12" ,
    "B13" ,
    "B14"
  )
)
####
dat$AreaGLM <- as.integer (dat$Area)
dat$VehPowerGLM <- as.factor (pmin (dat$VehPower , 9))
dat$VehAgeGLM <- as.factor (cut (
  dat$VehAge ,
  c (0 , 5 , 12 , 101) ,
  labels = c ("0 -5" , "6 -12" , "12+") ,
  include.lowest = TRUE
))
dat$DrivAgeGLM <-
  as.factor (cut (
    dat$DrivAge ,
    c (18 , 20 , 25 , 30 , 40 , 50 , 70 , 101) ,
    labels = c ("18 -20" , "21 -25" , "26 -30" , "31 -40" , "41 -50" , "51 -70" , "71+") ,
    include.lowest = TRUE
  ))
dat$BonusMalusGLM <- pmin (dat$BonusMalus , 150)
dat$DensityGLM <- log (dat$Density)
dat$VehGasGLM <- as.integer(dat$VehGas)
####

# str(dat)
# str(sev)
dat<-dat[which(dat$ClaimNb>0),]
dat$Sev<-dat$ClaimTotal/dat$ClaimNb

VehBrandOrd<-aggregate(dat$Sev,by=list(dat$VehBrand),mean)
VehBrandOrd<-VehBrandOrd[order(VehBrandOrd$x),]
VehBrandOrd$order<-1:nrow(VehBrandOrd)
dat$VehBrandOrd<-NA
for (k in 1:nrow(VehBrandOrd)){
  dat$VehBrandOrd[dat$VehBrand==VehBrandOrd$Group.1[k]]<-VehBrandOrd$order[k]
}
aggregate(dat$Sev,by=list(dat$VehBrandOrd),mean)

RegionOrd<-aggregate(dat$Sev,by=list(dat$Region),mean)
RegionOrd<-RegionOrd[order(RegionOrd$x),]
RegionOrd$order<-1:nrow(RegionOrd)
dat$RegionOrd<-NA
for (k in 1:nrow(RegionOrd)){
  dat$RegionOrd[dat$Region==RegionOrd$Group.1[k]]<-RegionOrd$order[k]
}
aggregate(dat$Sev,by=list(dat$RegionOrd),mean)
