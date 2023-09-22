gamma_mle<-function(y,z){
  # the MLE of gamma with response y and weights z
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
                      interval = c(0,5000))$maximum
  list(mu=mu_mle,shape=shape_mle,
       scale=mu_mle/shape_mle,rate=shape_mle/mu_mle)
}

mean_excess<-function(y){
  # the threshold in pareto distribution
  yy<-sort(y,decreasing = T)
  yy_sum<-cumsum(yy)
  nn_sum<-cumsum(rep(1,length(y)))
  list(x=yy,mean_ex=yy_sum/nn_sum)
}

hill_alpha<-function(y){
  # the tail index estimate in pareto distribution
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
  # density of pareto distribution
  ifelse(y<theta,0,alpha/theta*(y/theta)^(-alpha-1))
}

alpha_wmle<-function(y,z,theta){
  # scale MLE in the gamma distribution
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

logL_gamma<-function(weights, Y, Y_hat, shape){
  sum(weights) * shape * log(shape) -
    sum(weights * log(Y_hat)) * shape -
    sum(weights) * lgamma(shape) +
    (shape - 1) * sum(weights * log(Y)) -
    shape * sum(weights * Y / Y_hat)
}

# binomial_p <- function(){
#   Family(
#     loss = function(y, f) {
#       y*f-log(1+exp(f))
#     },
#     ngradient = function(y, f, w=1){
#       y - (exp(f)/(1+exp(f)))
#     },
#     offset = function(y, w) {
#       p = weighted.mean(y, w)
#       log(p/(1-p))
#     },
#     nuisance = function() return(NA),
#     response = function(f) exp(f)/(1+exp(f)),
#     rclass = function(f) NA,
#     name = "binomial_p")
# }
#
# p_linear<-function(x){
#   # log odds
#   log(x/(1-x))
# }

