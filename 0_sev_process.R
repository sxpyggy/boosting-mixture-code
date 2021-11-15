# library(OpenML,farff)
# freMTPL2freq <- getOMLDataSet(data.id = 41214)$data
rm(list=ls())
library ( CASdatasets )
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
####

str(dat)
str(sev)

gamma_mle<-function(y,z){
  mu_mle<-sum(z*y)/sum(z)
  n<-length(y)
  logL_shape<-function(shape_mle){
    sum(z)*shape_mle*log(shape_mle)-
      sum(z)*shape_mle*log(mu_mle)-
      sum(z)*lgamma(shape_mle)+
      (shape_mle-1)*sum(z*log(y))-
      shape_mle/mu_mle*sum(z*y)
  }
  shape_mle<-optimise(logL_shape,maximum = T,
                      interval = c(0,1000))$maximum
  list(mu=mu_mle,shape=shape_mle,
       scale=mu_mle/shape_mle,rate=shape_mle/mu_mle)
}

mean_excess<-function(y){
  yy<-sort(y,decreasing = T)
  yy_sum<-cumsum(yy)
  nn_sum<-cumsum(rep(1,length(y)))
  list(x=yy,mean_ex=yy_sum/nn_sum)
}

hill_alpha<-function(y){
  n<-length(y)
  yy<-sort(y,decreasing = T)
  log_yy<-log(yy)
  csum<-cumsum(log_yy)
  kk<-(n):1
  alpha<-((csum-log_yy*(n-kk+1))/(n-kk+1))^-1
  alpha_se<-sqrt(alpha^2*(n-kk+1)^2/(n-kk)^2/(n-kk-1))
  list(alpha=alpha,alpha_se=alpha_se,n=3:n)
}

dpareto<-function(y,theta,alpha){
  ifelse(y<theta,0,alpha/theta*(y/theta)^(-alpha-1))
}

alpha_wmle<-function(y,z,theta){
  z[which(y<theta)]<-0
  (sum(z*log(y))/sum(z)-log(theta))^(-1)
}

## mixed density function

f_mix_cp<-function(y_v,p_hat,par_mat,alpha,theta){
  # y_v: data vector
  # p_hat: mixing probs
  # par_mat: gamma parameters
  # alpha: tail index
  # theta: pareto threshold
  fy<-rep(NA,length(y_v))
  for (i in 1:length(y_v))
  {
    y<-y_v[i]
    dgam<-rep(NA,nrow(par_mat))
    for (k in 1:nrow(par_mat)){
      dgam[k]<-dgamma(y,shape=par_mat$shape[k],scale=par_mat$scale[k])
    }
    dpar<-dpareto(y,theta,alpha)
    fy[i]<-sum(p_hat*c(dgam,dpar))
  }
  fy
}

f_mix_vp<-function(y_v,p_hat_mat,par_mat,alpha,theta){
  # y_v: data vector
  # p_hat: mixing probs
  # par_mat: gamma parameters
  # alpha: tail index
  # theta: pareto threshold
  fy<-rep(NA,length(y_v))
  for (i in 1:length(y_v))
  {
    y<-y_v[i]
    dgam<-rep(NA,nrow(par_mat))
    for (k in 1:nrow(par_mat)){
      dgam[k]<-dgamma(y,shape=par_mat$shape[k],scale=par_mat$scale[k])
    }
    dpar<-dpareto(y,theta,alpha)
    fy[i]<-sum(p_hat_mat[i,]*c(dgam,dpar))
  }
  fy
}

binomial_p <- function(){
  Family(
    loss = function(y, f) {
      y*f-log(1+exp(f))
    },
    ngradient = function(y, f, w=1){
      y - (exp(f)/(1+exp(f)))
    },
    offset = function(y, w) {
      p = weighted.mean(y, w)
      log(p/(1-p))
    },
    nuisance = function() return(NA),
    response = function(f) exp(f)/(1+exp(f)),
    rclass = function(f) NA,
    name = "binomial_p")
}

p_linear<-function(x){
  log(x/(1-x))
}
