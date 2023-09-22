library("gbm")
library("splines")
neg_ll3 <- function(y, mu, sigma, p) {
  ## y has dimension of n; mu, sigma, p have dimension of K*n ##
  return(-sum(log(
    p[, 1] * dnorm(y, mean = mu[, 1], sd = sigma[, 1]) +
      p[, 2] * dnorm(y, mean = mu[, 2], sd = sigma[, 2]) +
      p[, 3] * dnorm(y, mean = mu[, 3], sd = sigma[, 3])
  )) / length(y))
}
neg_ll5 <- function(y, mu, sigma, p) {
  ## y has dimension of n; mu, sigma, p have dimension of K*n ##
  return(-sum(log(
      p[, 1] * dnorm(y, mean = mu[, 1], sd = sigma[, 1]) +
      p[, 2] * dnorm(y, mean = mu[, 2], sd = sigma[, 2]) +
      p[, 3] * dnorm(y, mean = mu[, 3], sd = sigma[, 3]) +
      p[, 4] * dnorm(y, mean = mu[, 4], sd = sigma[, 4]) +
      p[, 5] * dnorm(y, mean = mu[, 5], sd = sigma[, 5])
  )) / length(y))
}

sim_gaussian <- function(n, seed) {
  # mixing probs are related to covariates, means are not related to covariates
  set.seed(seed)

  X1 = rnorm(n)
  X2 = rnorm(n)
  X3 = rnorm(n)
  X4 = rbinom(n,1,0.5)
  X5 = X1+rnorm(n)
  X6 = rnorm(n)

  F1 = tanh(X1^2+X2)
  F2 = tanh(X1+X2^2)
  F3 = tanh(2*X4+X3+X3^2+2*X3*X4)

  # pairs(cbind(F1,F2,F3,X1,X2,X3,X4))

  P <- FtoP(FF = cbind(F1, F2, F3))

  # pairs(cbind(P,X1,X2,X3,X4))
  # boxplot(P)

  MU1 <- rep(-5,n)
  MU2 <- rep(0,n)
  MU3 <- rep(5,n)
  # MU1 <- 7 + 2 * X1 + exp(0.3 * X2)
  # MU2 <- 1 - 0.1 * X1 ^ 2 + 5 * X3 - X1 * X2
  # MU3 <- -5 - 2 * X1 + 4 * sin(X2) + 0.2 * log(X4)
  dat <-
    data.frame(
      X1,
      X2,
      X3,
      X4,
      X5,
      X6,
      F1,
      F2,
      F3,
      P1 = P[, 1],
      P2 = P[, 2],
      P3 = P[, 3],
      MU1,
      MU2,
      MU3,
      SG1 = 0.5,
      SG2 = 1,
      SG3 = 2,
      Z1 = NA,
      Z2 = NA,
      Z3 = NA,
      Y1 = NA,
      Y2 = NA,
      Y3 = NA,
      Y = NA
    )
  for (i in 1:nrow(dat)) {
    dat[i, c("Z1", "Z2", "Z3")] <-
      t(rmultinom(1, size = 1, dat[i, c("P1", "P2", "P3")]))
    dat[i, c("Y1", "Y2", "Y3")] <-
      t(rnorm(3, mean = as.numeric(dat[i, c("MU1", "MU2", "MU3")]), sd =
                as.numeric(dat[i, c("SG1", "SG2", "SG3")])))
    dat$Y[i] <-
      sum(dat[i, c("Z1", "Z2", "Z3")] * dat[i, c("Y1", "Y2", "Y3")])
  }
  dat[,c("F1","F2","F3")]<-PtoF(dat[,c("P1","P2","P3")])
  dat
}

dat_Norm<-function(K,n){
  # put into a suitable format for mixture of Gaussian examples
  # z is the indictor of component; f is pdf of components; mu is the mean of component
  # p is the mixing probabilities; plinear is  the linear predictor in mixing probs
  # sigma is the variance; weight is the weight for mu regression
  dat_bst<-data.frame(BST=matrix(NA,ncol=7*K,nrow=n))
  bst_names<-NULL
  for (k in 1:K){
    vnames<-c(paste("z",k,sep=""),paste("f",k,sep=""),paste("mu",k,sep=""),
              paste("p",k,sep=""),paste("plinear",k,sep=""),
              paste("sigma",k,sep=""),paste("weight",k,sep=""))
    bst_names<-c(bst_names,vnames)
  }
  names(dat_bst)<-bst_names
  dat_bst
}

EM_gaussian <-
  function(X, Y, Xtest, Ytest, M0, model, structure, trace, patience) {

    # model = "glm" or "null"; structure = "both", "mu", "p"

    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_test <- dat_Norm(K = 3, n = nrow(Xtest))

    #initialization

    par_mat$p1 <- 1/3
    par_mat$p2 <- 1/3
    par_mat$p3 <- 1/3
    par_mat$mu1 <- 1
    par_mat$mu2 <- 0
    par_mat$mu3 <- -1
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    learn_loss <- NULL
    test_loss <- NULL
    learn_impT <- NULL
    learn_impT[1]<-T

    for (m in 1:M0) {
      #expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)

      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz
      # the expectation of z

      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2

      # null model for mu
      mu1_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight3)

      # null model for z
      p_glm <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ 1,
                 data = X,
                 trace = F)


      if (model == "glm" &
          (structure == "mu" | structure == "both")) {
        mu1_glm <- lm(Y ~ ., data = X, weights = par_mat$weight1)
        mu2_glm <- lm(Y ~ ., data = X, weights = par_mat$weight2)
        mu3_glm <- lm(Y ~ ., data = X, weights = par_mat$weight3)
      }

      if (model == "glm" &
          (structure == "p" | structure == "both")) {
        p_glm <-
          multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ .,
                   data = X,
                   trace = F)
      }

      par_mat$mu1 <- predict(mu1_glm, newdata = X)
      par_mat$mu2 <- predict(mu2_glm, newdata = X)
      par_mat$mu3 <- predict(mu3_glm, newdata = X)

      par_mat[, c("p1", "p2", "p3")] <-
        predict(p_glm, newdata = X, type = "probs")

      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))

      par_mat_test$mu1 <- predict(mu1_glm, newdata = Xtest)
      par_mat_test$mu2 <- predict(mu2_glm, newdata = Xtest)
      par_mat_test$mu3 <- predict(mu3_glm, newdata = Xtest)

      par_mat_test[, c("p1", "p2", "p3")] <-
        predict(p_glm, newdata = Xtest, type = "probs")

      par_mat_test$sigma1 <- unique(par_mat$sigma1)
      par_mat_test$sigma2 <- unique(par_mat$sigma2)
      par_mat_test$sigma3 <- unique(par_mat$sigma3)

      learn_loss[m] <-
        neg_ll3(Y, par_mat[, c("mu1", "mu2", "mu3")], par_mat[, c("sigma1", "sigma2", "sigma3")], par_mat[, c("p1", "p2", "p3")])
      test_loss[m] <-
        neg_ll3(Ytest, par_mat_test[, c("mu1", "mu2", "mu3")], par_mat_test[, c("sigma1", "sigma2", "sigma3")], par_mat_test[, c("p1", "p2", "p3")])
      if (m>1){
        learn_impT[m] <- (learn_loss[m-1]-learn_loss[m])>10^-4
      }
       if (trace == T) {
        print(paste(
          "iteration:",
          m,
          "; ",
          "learn_loss:",
          round(learn_loss[m], 4),
          "; ",
          "test_loss:",
          round(test_loss[m], 4),
          "; ",
          learn_impT[m],
          sep=""
        ))
       }
      if (m>patience){
        if (sum(learn_impT[(m-patience+1):m])==0){
          break
      }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat[, c("p1", "p2", "p3")])
    par_mat_test[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat_test[, c("p1", "p2", "p3")])

    list(
      learn_loss = learn_loss,
      test_loss = test_loss,
      fitted = par_mat,
      test = par_mat_test,
      mu_models = list(mu1_glm, mu2_glm, mu3_glm),
      p_model = p_glm,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3)
      )
    )
  }

EM_gaussian_r <-
  function(X, Y, Xtest, Ytest, M0, model, structure, trace, patience) {

    # model = "glm" or "null"; structure = "both", "mu", "p"

    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_test <- dat_Norm(K = 3, n = nrow(Xtest))

    #initialization

    par_mat$p1 <- 1/3
    par_mat$p2 <- 1/3
    par_mat$p3 <- 1/3
    par_mat$mu1 <- 1
    par_mat$mu2 <- 0
    par_mat$mu3 <- -1
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    learn_loss <- NULL
    test_loss <- NULL
    learn_impT <- NULL
    learn_impT[1]<-T

    for (m in 1:M0) {
      #expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)

      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz
      # the expectation of z

      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2

      # null model for mu
      mu1_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight3)

      # null model for z
      p_glm <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ 1,
                 data = X,
                 trace = F)


      if (model == "glm" &
          (structure == "mu" | structure == "both")) {
        mu1_glm <- lm(Y ~ ., data = X, weights = par_mat$weight1)
      }

      if (model == "glm" &
          (structure == "p" | structure == "both")) {
        p_glm <-
          multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ .,
                   data = X,
                   trace = F)
      }

      par_mat$mu1 <- predict(mu1_glm, newdata = X)
      par_mat$mu2 <- predict(mu2_glm, newdata = X)
      par_mat$mu3 <- predict(mu3_glm, newdata = X)

      par_mat[, c("p1", "p2", "p3")] <-
        predict(p_glm, newdata = X, type = "probs")

      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))

      par_mat_test$mu1 <- predict(mu1_glm, newdata = Xtest)
      par_mat_test$mu2 <- predict(mu2_glm, newdata = Xtest)
      par_mat_test$mu3 <- predict(mu3_glm, newdata = Xtest)

      par_mat_test[, c("p1", "p2", "p3")] <-
        predict(p_glm, newdata = Xtest, type = "probs")

      par_mat_test$sigma1 <- unique(par_mat$sigma1)
      par_mat_test$sigma2 <- unique(par_mat$sigma2)
      par_mat_test$sigma3 <- unique(par_mat$sigma3)

      learn_loss[m] <-
        neg_ll3(Y, par_mat[, c("mu1", "mu2", "mu3")], par_mat[, c("sigma1", "sigma2", "sigma3")], par_mat[, c("p1", "p2", "p3")])
      test_loss[m] <-
        neg_ll3(Ytest, par_mat_test[, c("mu1", "mu2", "mu3")], par_mat_test[, c("sigma1", "sigma2", "sigma3")], par_mat_test[, c("p1", "p2", "p3")])
      if (m>1){
        learn_impT[m] <- (learn_loss[m-1]-learn_loss[m])>10^-4
      }
      if (trace == T) {
        print(paste(
          "iteration:",
          m,
          "; ",
          "learn_loss:",
          round(learn_loss[m], 4),
          "; ",
          "test_loss:",
          round(test_loss[m], 4),
          "; ",
          learn_impT[m],
          sep=""
        ))
      }
      if (m>patience){
        if (sum(learn_impT[(m-patience+1):m])==0){
          break
        }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat[, c("p1", "p2", "p3")])
    par_mat_test[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat_test[, c("p1", "p2", "p3")])

    list(
      learn_loss = learn_loss,
      test_loss = test_loss,
      fitted = par_mat,
      test = par_mat_test,
      mu_models = list(mu1_glm, mu2_glm, mu3_glm),
      p_model = p_glm,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3)
      )
    )
  }

EB_gaussian <-

  function(X, Y, Xval, Yval, M0, n_tree_mu, n_tree_p, structure, trace,
           patience, parallel,int_mu1,int_mu2,int_mu3) {

    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_val <- NULL

    #initialization
    par_mat$p1 <- 1/3
    par_mat$p2 <- 1/3
    par_mat$p3 <- 1/3
    par_mat$mu1 <- -1
    par_mat$mu2 <- 0
    par_mat$mu3 <- 1
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    if (is.null(Yval)!=TRUE){
      par_mat_val <- dat_Norm(K = 3, n = nrow(Xval))
      par_mat_val$p1 <- 1/3
      par_mat_val$p2 <- 1/3
      par_mat_val$p3 <- 1/3
      par_mat_val$mu1 <- -1
      par_mat_val$mu2 <- 0
      par_mat_val$mu3 <- 1
      par_mat_val$sigma1 <- 0.5
      par_mat_val$sigma2 <- 0.5
      par_mat_val$sigma3 <- 0.5
    }

    train_loss <- NULL
    valid_loss <- NULL
    mu_iter<-NULL
    train_impT <- NULL
    train_impT[1]<-T
    train_imp <- 0

    for (m in 1:M0) {

      # expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)

      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz

      par_mat$weight1 = par_mat$z1 #/ par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 #/ par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 #/ par_mat$sigma3 ^ 2

      if (is.null(Yval)!=TRUE){

      par_mat_val$f1 = dnorm(Yval, par_mat_val$mu1, par_mat_val$sigma1)
      par_mat_val$f2 = dnorm(Yval, par_mat_val$mu2, par_mat_val$sigma2)
      par_mat_val$f3 = dnorm(Yval, par_mat_val$mu3, par_mat_val$sigma3)

      sum_fz_val <-
        par_mat_val$p1 * par_mat_val$f1 +
        par_mat_val$p2 * par_mat_val$f2 +
        par_mat_val$p3 * par_mat_val$f3

      par_mat_val$z1 = par_mat_val$p1 * par_mat_val$f1 / sum_fz_val
      par_mat_val$z2 = par_mat_val$p2 * par_mat_val$f2 / sum_fz_val
      par_mat_val$z3 = par_mat_val$p3 * par_mat_val$f3 / sum_fz_val

      par_mat_val$weight1 = par_mat_val$z1 #/ par_mat_val$sigma1 ^ 2
      par_mat_val$weight2 = par_mat_val$z2 #/ par_mat_val$sigma2 ^ 2
      par_mat_val$weight3 = par_mat_val$z3 #/ par_mat_val$sigma3 ^ 2
      }

      # null model for mu and p

      mu1_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
      par_mat$mu1 <- predict(mu1_bst, newdata = X)
      par_mat$mu2 <- predict(mu2_bst, newdata = X)
      par_mat$mu3 <- predict(mu3_bst, newdata = X)
      if (is.null(Yval)!=TRUE){
        par_mat_val$mu1 <- predict(mu1_bst, newdata = Xval)
        par_mat_val$mu2 <- predict(mu2_bst, newdata = Xval)
        par_mat_val$mu3 <- predict(mu3_bst, newdata = Xval)
      }

      p_bst <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ 1,
                 data = X,
                 trace = F)
      par_mat[, c("p1", "p2", "p3")] <-
        predict(p_bst, newdata = X, type = "probs")
      if (is.null(Yval)!=TRUE){
      par_mat_val[, c("p1", "p2", "p3")] <-
        predict(p_bst, newdata = Xval, type = "probs")
      }

      # boosting for mu

      if (structure == "mu"| structure == "both") {
        param<-list(max_depth=4, eta =0.5, objective="reg:squarederror")
        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight1)
        setinfo(dtrain, "base_margin", rep(int_mu1,nrow(X)))
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval,weight=par_mat_val$weight1)
        setinfo(dvalid, "base_margin", rep(int_mu1,nrow(Xval)))
        watchlist=list(train=dtrain, eval= dvalid)
        mu1_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu1:", mu1_bst$niter))
        par_mat$mu1 <- predict(mu1_bst, newdata = dtrain)
        if (is.null(Yval)!=TRUE) {
          par_mat_val$mu1 <- predict(mu1_bst, newdata = dvalid)
        }

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y, weight=par_mat$weight2)
        setinfo(dtrain, "base_margin", rep(int_mu2,nrow(X)))
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval,weight=par_mat_val$weight2)
        setinfo(dvalid, "base_margin", rep(int_mu2,nrow(Xval)))
        watchlist=list(train=dtrain, eval= dvalid)
        mu2_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu2:", mu2_bst$niter))
        par_mat$mu2 <- predict(mu2_bst, newdata = dtrain)
        if (is.null(Yval)!=TRUE){
          par_mat_val$mu2 <- predict(mu2_bst, newdata = dvalid)
        }

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight3)
        setinfo(dtrain, "base_margin", rep(int_mu3,nrow(X)))
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval,weight=par_mat_val$weight3)
        setinfo(dvalid, "base_margin", rep(int_mu3,nrow(Xval)))
        watchlist=list(train=dtrain, eval= dvalid)
        mu3_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu3:", mu3_bst$niter))
        par_mat$mu3 <- predict(mu3_bst, newdata = dtrain)
        if (is.null(Yval)!=TRUE){
          par_mat_val$mu3 <- predict(mu3_bst, newdata = dvalid)
        }
      }

      # boosting for p

      if (structure == "p" |structure == "both") {
        p_bst <-
          BST(
            X = X,
            Y = par_mat[,c("z1","z2","z3")],
            Pinit = NULL,
            Xval = Xval,
            Yval = par_mat_val[, c("z1", "z2", "z3")],
            Pvalinit = NULL,
            M = n_tree_p,
            cp = 0.001,
            maxdepth = 4,
            lr = 0.2,
            trace =T,
            patience=2,
            parallel = parallel
          )
        par_mat[, c("p1", "p2", "p3")] <-
          p_bst$fitted[,c("BST_p1","BST_p2","BST_p3")]
        if (is.null(Yval)!=TRUE){
        par_mat_val[, c("p1", "p2", "p3")] <-
          p_bst$valid[,c("BST_p1","BST_p2","BST_p3")]
        }
      }

      # sigma
      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))

      if (is.null(Yval)!=TRUE){
      par_mat_val$sigma1 <- unique(par_mat$sigma1)
      par_mat_val$sigma2 <- unique(par_mat$sigma2)
      par_mat_val$sigma3 <- unique(par_mat$sigma3)
      }

      # outer train loss and validation loss

      train_loss[m] <-
        neg_ll3(Y, par_mat[, c("mu1", "mu2", "mu3")], par_mat[, c("sigma1", "sigma2", "sigma3")], par_mat[, c("p1", "p2", "p3")])
      if (is.null(Yval)!=TRUE){
      valid_loss[m] <-
        neg_ll3(Yval, par_mat_val[, c("mu1", "mu2", "mu3")], par_mat_val[, c("sigma1", "sigma2", "sigma3")], par_mat_val[, c("p1", "p2", "p3")])
      }
      if (m>1){
        train_imp <- train_loss[m-1]-train_loss[m]
        train_impT[m] <- (train_imp>10^-4)
      }
      if (trace == T) {
        print(paste(
          "EB-iteration:",
          m,
          "; ",
          "train_loss:",
          round(train_loss[m], 4),
          "; ",
          "valid_loss:",
          round(valid_loss[m], 4),
          "; ",
          "train_loss_improve",
          round(train_imp, 4),
          "; ",
          train_impT[m],
          sep=""
        ))
      }
      if (m>patience){
        # minimize the training loss rather than early stop
        if (sum(train_impT[(m-patience+1):m])==0){
          break
        }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat[, c("p1", "p2", "p3")])
    if (is.null(Yval)!=TRUE){
    par_mat_val[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat_val[, c("p1", "p2", "p3")])
    }

    list(
      train_loss = train_loss,
      valid_loss = valid_loss,
      fitted = par_mat,
      valid = par_mat_val,
      mu_models = list(mu1_bst, mu2_bst, mu3_bst),
      mu_iter = c(mu1_bst$niter,mu2_bst$niter,mu3_bst$niter),
      p_model = p_bst,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3)
      )
    )
  }

EB_gaussian_r <-

  function(X, Y, Xval, Yval, M0, n_tree_mu, n_tree_p, structure, trace, patience,
           parallel,int_mu1,int_mu2,int_mu3) {

    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_val <- NULL

    #initialization
    par_mat$p1 <- 1/3
    par_mat$p2 <- 1/3
    par_mat$p3 <- 1/3
    par_mat$mu1 <- -1
    par_mat$mu2 <- 0
    par_mat$mu3 <- 1
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    if (is.null(Yval)!=TRUE){
      par_mat_val <- dat_Norm(K = 3, n = nrow(Xval))
      par_mat_val$p1 <- 1/3
      par_mat_val$p2 <- 1/3
      par_mat_val$p3 <- 1/3
      par_mat_val$mu1 <- -1
      par_mat_val$mu2 <- 0
      par_mat_val$mu3 <- 1
      par_mat_val$sigma1 <- 0.5
      par_mat_val$sigma2 <- 0.5
      par_mat_val$sigma3 <- 0.5
    }

    train_loss <- NULL
    valid_loss <- NULL
    mu_iter<-NULL
    train_impT <- NULL
    train_impT[1]<-T
    train_imp <- 0

    for (m in 1:M0) {

      # expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)

      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz

      par_mat$weight1 = par_mat$z1 #/ par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 #/ par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 #/ par_mat$sigma3 ^ 2

      if (is.null(Yval)!=TRUE){

        par_mat_val$f1 = dnorm(Yval, par_mat_val$mu1, par_mat_val$sigma1)
        par_mat_val$f2 = dnorm(Yval, par_mat_val$mu2, par_mat_val$sigma2)
        par_mat_val$f3 = dnorm(Yval, par_mat_val$mu3, par_mat_val$sigma3)

        sum_fz_val <-
          par_mat_val$p1 * par_mat_val$f1 +
          par_mat_val$p2 * par_mat_val$f2 +
          par_mat_val$p3 * par_mat_val$f3

        par_mat_val$z1 = par_mat_val$p1 * par_mat_val$f1 / sum_fz_val
        par_mat_val$z2 = par_mat_val$p2 * par_mat_val$f2 / sum_fz_val
        par_mat_val$z3 = par_mat_val$p3 * par_mat_val$f3 / sum_fz_val

        par_mat_val$weight1 = par_mat_val$z1 #/ par_mat_val$sigma1 ^ 2
        par_mat_val$weight2 = par_mat_val$z2 #/ par_mat_val$sigma2 ^ 2
        par_mat_val$weight3 = par_mat_val$z3 #/ par_mat_val$sigma3 ^ 2
      }

      # null model for mu and p

      mu1_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
      par_mat$mu1 <- predict(mu1_bst, newdata = X)
      par_mat$mu2 <- predict(mu2_bst, newdata = X)
      par_mat$mu3 <- predict(mu3_bst, newdata = X)
      if (is.null(Yval)!=TRUE){
        par_mat_val$mu1 <- predict(mu1_bst, newdata = Xval)
        par_mat_val$mu2 <- predict(mu2_bst, newdata = Xval)
        par_mat_val$mu3 <- predict(mu3_bst, newdata = Xval)
      }

      p_bst <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ 1,
                 data = X,
                 trace = F)
      par_mat[, c("p1", "p2", "p3")] <-
        predict(p_bst, newdata = X, type = "probs")
      if (is.null(Yval)!=TRUE){
        par_mat_val[, c("p1", "p2", "p3")] <-
          predict(p_bst, newdata = Xval, type = "probs")
      }

      # boosting for mu

      if (structure == "mu"| structure == "both") {
        param<-list(max_depth=4, eta =0.5, objective="reg:squarederror")
        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight1)
        setinfo(dtrain, "base_margin", rep(int_mu1,nrow(X)))
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval,weight=par_mat_val$weight1)
        setinfo(dvalid, "base_margin", rep(int_mu1,nrow(Xval)))
        watchlist=list(train=dtrain, eval= dvalid)
        mu1_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu1:", mu1_bst$niter))
        par_mat$mu1 <- predict(mu1_bst, newdata = dtrain)
        if (is.null(Yval)!=TRUE) {
          par_mat_val$mu1 <- predict(mu1_bst, newdata = dvalid)
        }

        # dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight3)
        # setinfo(dtrain, "base_margin", rep(int_mu3,nrow(X)))
        # dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval,weight=par_mat_val$weight3)
        # setinfo(dvalid, "base_margin", rep(int_mu3,nrow(Xval)))
        # watchlist=list(train=dtrain, eval= dvalid)
        # mu3_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        # print(paste("best iteration for mu3:", mu3_bst$niter))
        # par_mat$mu3 <- predict(mu3_bst, newdata = dtrain)
        # if (is.null(Yval)!=TRUE){
        #   par_mat_val$mu3 <- predict(mu3_bst, newdata = dvalid)
        # }
      }


      # boosting for p

      if (structure == "p" |structure == "both") {
        p_bst <-
          BST(
            X = X,
            Y = par_mat[,c("z1","z2","z3")],
            Pinit = NULL,
            Xval = Xval,
            Yval = par_mat_val[, c("z1", "z2", "z3")],
            Pvalinit = NULL,
            M = n_tree_p,
            cp = 0.001,
            maxdepth = 4,
            lr = 0.2,
            trace =T,
            patience=2,
            parallel = parallel
          )
        par_mat[, c("p1", "p2", "p3")] <-
          p_bst$fitted[,c("BST_p1","BST_p2","BST_p3")]
        if (is.null(Yval)!=TRUE){
          par_mat_val[, c("p1", "p2", "p3")] <-
            p_bst$valid[,c("BST_p1","BST_p2","BST_p3")]
        }
      }

      # sigma
      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))

      if (is.null(Yval)!=TRUE){
        par_mat_val$sigma1 <- unique(par_mat$sigma1)
        par_mat_val$sigma2 <- unique(par_mat$sigma2)
        par_mat_val$sigma3 <- unique(par_mat$sigma3)
      }

      # outer train loss and validation loss

      train_loss[m] <-
        neg_ll3(Y, par_mat[, c("mu1", "mu2", "mu3")], par_mat[, c("sigma1", "sigma2", "sigma3")], par_mat[, c("p1", "p2", "p3")])
      if (is.null(Yval)!=TRUE){
        valid_loss[m] <-
          neg_ll3(Yval, par_mat_val[, c("mu1", "mu2", "mu3")], par_mat_val[, c("sigma1", "sigma2", "sigma3")], par_mat_val[, c("p1", "p2", "p3")])
      }
      if (m>1){
        train_imp <- train_loss[m-1]-train_loss[m]
        train_impT[m] <- (train_imp>10^-4)
      }
      if (trace == T) {
        print(paste(
          "EB-iteration:",
          m,
          "; ",
          "train_loss:",
          round(train_loss[m], 4),
          "; ",
          "valid_loss:",
          round(valid_loss[m], 4),
          "; ",
          "train_loss_improve",
          round(train_imp, 4),
          "; ",
          train_impT[m],
          sep=""
        ))
      }
      if (m>patience){
        # minimize the training loss rather than early stop
        if (sum(train_impT[(m-patience+1):m])==0){
          break
        }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat[, c("p1", "p2", "p3")])
    if (is.null(Yval)!=TRUE){
      par_mat_val[, c("plinear1", "plinear2", "plinear3")] <-
        PtoF(par_mat_val[, c("p1", "p2", "p3")])
    }

    list(
      train_loss = train_loss,
      valid_loss = valid_loss,
      fitted = par_mat,
      valid = par_mat_val,
      mu_models = list(mu1_bst, mu2_bst, mu3_bst),
      mu_iter = c(mu1_bst$niter, mu2_bst$niter, mu3_bst$niter),
      p_model = p_bst,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3)
      )
    )
  }

EM_gaussian5 <-
  function(X, Y, Xtest, Ytest, M0, model, structure, trace, patience) {

    # model = "glm" or "null"; structure = "both", "mu", "p"

    par_mat <- dat_Norm(K = 5, n = nrow(X))
    par_mat_test <- dat_Norm(K = 5, n = nrow(Xtest))

    #initialization

    par_mat$p1 <- 1/5
    par_mat$p2 <- 1/5
    par_mat$p3 <- 1/5
    par_mat$p4 <- 1/5
    par_mat$p5 <- 1/5

    par_mat$mu1 <- -2
    par_mat$mu2 <- -1
    par_mat$mu3 <- 0
    par_mat$mu4 <- 1
    par_mat$mu5 <- 2

    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5
    par_mat$sigma4 <- 0.5
    par_mat$sigma5 <- 0.5

    learn_loss <- NULL
    test_loss <- NULL
    learn_impT <- NULL
    learn_impT[1]<-T

    for (m in 1:M0) {
      #expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)
      par_mat$f4 = dnorm(Y, par_mat$mu4, par_mat$sigma4)
      par_mat$f5 = dnorm(Y, par_mat$mu5, par_mat$sigma5)

      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3 +
        par_mat$p4 * par_mat$f4 + par_mat$p5 * par_mat$f5

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz
      par_mat$z4 = par_mat$p4 * par_mat$f4 / sum_fz
      par_mat$z5 = par_mat$p5 * par_mat$f5 / sum_fz
      # the expectation of z

      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2
      par_mat$weight4 = par_mat$z4 / par_mat$sigma4 ^ 2
      par_mat$weight5 = par_mat$z5 / par_mat$sigma5 ^ 2

      # null model for mu
      mu1_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
      mu4_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight4)
      mu5_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight5)

      # null model for z
      p_glm <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3", "z4", "z5")]) ~ 1,
                 data = X,
                 trace = F)


      if (model == "glm" &
          (structure == "mu" | structure == "both")) {
        mu1_glm <- lm(Y ~ ., data = X, weights = par_mat$weight1)
        mu2_glm <- lm(Y ~ ., data = X, weights = par_mat$weight2)
        mu3_glm <- lm(Y ~ ., data = X, weights = par_mat$weight3)
        mu4_glm <- lm(Y ~ ., data = X, weights = par_mat$weight4)
        mu5_glm <- lm(Y ~ ., data = X, weights = par_mat$weight5)
      }

      if (model == "glm" &
          (structure == "p" | structure == "both")) {
        p_glm <-
          multinom(as.matrix(par_mat[, c("z1", "z2", "z3", "z4", "z5")]) ~.,
                   data = X,
                   trace = F)
      }

      par_mat$mu1 <- predict(mu1_glm, newdata = X)
      par_mat$mu2 <- predict(mu2_glm, newdata = X)
      par_mat$mu3 <- predict(mu3_glm, newdata = X)
      par_mat$mu4 <- predict(mu4_glm, newdata = X)
      par_mat$mu5 <- predict(mu5_glm, newdata = X)

      par_mat[, c("p1", "p2", "p3", "p4", "p5")] <-
        predict(p_glm, newdata = X, type = "probs")

      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))
      par_mat$sigma4 <-
        sqrt(sum(par_mat$z4 * (Y - par_mat$mu4) ^ 2) / sum(par_mat$z4))
      par_mat$sigma5 <-
        sqrt(sum(par_mat$z5 * (Y - par_mat$mu5) ^ 2) / sum(par_mat$z5))

      par_mat_test$mu1 <- predict(mu1_glm, newdata = Xtest)
      par_mat_test$mu2 <- predict(mu2_glm, newdata = Xtest)
      par_mat_test$mu3 <- predict(mu3_glm, newdata = Xtest)
      par_mat_test$mu4 <- predict(mu4_glm, newdata = Xtest)
      par_mat_test$mu5 <- predict(mu5_glm, newdata = Xtest)


      par_mat_test[, c("p1", "p2", "p3", "p4", "p5")] <-
        predict(p_glm, newdata = Xtest, type = "probs")

      par_mat_test$sigma1 <- unique(par_mat$sigma1)
      par_mat_test$sigma2 <- unique(par_mat$sigma2)
      par_mat_test$sigma3 <- unique(par_mat$sigma3)
      par_mat_test$sigma4 <- unique(par_mat$sigma4)
      par_mat_test$sigma5 <- unique(par_mat$sigma5)

      learn_loss[m] <-
        neg_ll5(Y, par_mat[, c("mu1", "mu2", "mu3","mu4","mu5")],
                par_mat[, c("sigma1", "sigma2", "sigma3","sigma4","sigma5")],
                par_mat[, c("p1", "p2", "p3", "p4", "p5")])
      test_loss[m] <-
        neg_ll5(Ytest, par_mat_test[, c("mu1", "mu2", "mu3", "mu4", "mu5")],
                par_mat_test[, c("sigma1", "sigma2", "sigma3", "sigma4", "sigma5")],
                par_mat_test[, c("p1", "p2", "p3", "p4", "p5")])
      if (m>1){
        learn_impT[m] <- (learn_loss[m-1]-learn_loss[m])>10^-4
      }
      if (trace == T) {
        print(paste(
          "iteration:",
          m,
          "; ",
          "learn_loss:",
          round(learn_loss[m], 4),
          "; ",
          "test_loss:",
          round(test_loss[m], 4),
          "; ",
          learn_impT[m],
          sep=""
        ))
      }
      if (m>patience){
        if (sum(learn_impT[(m-patience+1):m])==0){
          break
        }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3", "plinear4", "plinear5")] <-
      PtoF(par_mat[, c("p1", "p2", "p3", "p4", "p5")])
    par_mat_test[, c("plinear1", "plinear2", "plinear3", "plinear4", "plinear5")] <-
      PtoF(par_mat_test[, c("p1", "p2", "p3", "p4", "p5")])

    list(
      learn_loss = learn_loss,
      test_loss = test_loss,
      fitted = par_mat,
      test = par_mat_test,
      mu_models = list(mu1_glm, mu2_glm, mu3_glm, mu4_glm, mu5_glm),
      p_model = p_glm,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3),
        unique(par_mat$sigma4),
        unique(par_mat$sigma5)
      )
    )
  }

EB_gaussian5 <-

  function(X, Y, Xval, Yval, M0, n_tree_mu, n_tree_p, structure, trace, patience, parallel) {

    par_mat <- dat_Norm(K = 5, n = nrow(X))
    par_mat_val <- NULL

    #initialization
    par_mat$p1 <- 1/5
    par_mat$p2 <- 1/5
    par_mat$p3 <- 1/5
    par_mat$p4 <- 1/5
    par_mat$p5 <- 1/5

    par_mat$mu1 <- -2
    par_mat$mu2 <- -1
    par_mat$mu3 <- 0
    par_mat$mu4 <- 1
    par_mat$mu5 <- 2

    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5
    par_mat$sigma4 <- 0.5
    par_mat$sigma5 <- 0.5


    if (is.null(Yval)!=TRUE){
      par_mat_val <- dat_Norm(K = 5, n = nrow(Xval))
      par_mat_val$p1 <- 1/5
      par_mat_val$p2 <- 1/5
      par_mat_val$p3 <- 1/5
      par_mat_val$p4 <- 1/5
      par_mat_val$p5 <- 1/5

      par_mat_val$mu1 <- -2
      par_mat_val$mu2 <- -1
      par_mat_val$mu3 <- 0
      par_mat_val$mu4 <- 1
      par_mat_val$mu5 <- 2

      par_mat_val$sigma1 <- 0.5
      par_mat_val$sigma2 <- 0.5
      par_mat_val$sigma3 <- 0.5
      par_mat_val$sigma4 <- 0.5
      par_mat_val$sigma5 <- 0.5
    }

    train_loss <- NULL
    valid_loss <- NULL
    mu_iter<-NULL
    train_impT <- NULL
    train_impT[1]<-T
    train_imp <- 0

    for (m in 1:M0) {

      # expectation of latent variable

      par_mat$f1 = dnorm(Y, par_mat$mu1, par_mat$sigma1)
      par_mat$f2 = dnorm(Y, par_mat$mu2, par_mat$sigma2)
      par_mat$f3 = dnorm(Y, par_mat$mu3, par_mat$sigma3)
      par_mat$f4 = dnorm(Y, par_mat$mu4, par_mat$sigma4)
      par_mat$f5 = dnorm(Y, par_mat$mu5, par_mat$sigma5)


      sum_fz <-
        par_mat$p1 * par_mat$f1 + par_mat$p2 * par_mat$f2 + par_mat$p3 * par_mat$f3 +
        par_mat$p4 * par_mat$f4 + par_mat$p5 * par_mat$f5

      par_mat$z1 = par_mat$p1 * par_mat$f1 / sum_fz
      par_mat$z2 = par_mat$p2 * par_mat$f2 / sum_fz
      par_mat$z3 = par_mat$p3 * par_mat$f3 / sum_fz
      par_mat$z4 = par_mat$p4 * par_mat$f4 / sum_fz
      par_mat$z5 = par_mat$p5 * par_mat$f5 / sum_fz


      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2
      par_mat$weight4 = par_mat$z4 / par_mat$sigma4 ^ 2
      par_mat$weight5 = par_mat$z5 / par_mat$sigma5 ^ 2

      if (is.null(Yval)!=TRUE){

        par_mat_val$f1 = dnorm(Yval, par_mat_val$mu1, par_mat_val$sigma1)
        par_mat_val$f2 = dnorm(Yval, par_mat_val$mu2, par_mat_val$sigma2)
        par_mat_val$f3 = dnorm(Yval, par_mat_val$mu3, par_mat_val$sigma3)
        par_mat_val$f4 = dnorm(Yval, par_mat_val$mu4, par_mat_val$sigma4)
        par_mat_val$f5 = dnorm(Yval, par_mat_val$mu5, par_mat_val$sigma5)


        sum_fz_val <-
          par_mat_val$p1 * par_mat_val$f1 +
          par_mat_val$p2 * par_mat_val$f2 +
          par_mat_val$p3 * par_mat_val$f3 +
          par_mat_val$p4 * par_mat_val$f4 +
          par_mat_val$p5 * par_mat_val$f5

        par_mat_val$z1 = par_mat_val$p1 * par_mat_val$f1 / sum_fz_val
        par_mat_val$z2 = par_mat_val$p2 * par_mat_val$f2 / sum_fz_val
        par_mat_val$z3 = par_mat_val$p3 * par_mat_val$f3 / sum_fz_val
        par_mat_val$z4 = par_mat_val$p4 * par_mat_val$f4 / sum_fz_val
        par_mat_val$z5 = par_mat_val$p5 * par_mat_val$f5 / sum_fz_val

        par_mat_val$weight1 = par_mat_val$z1 / par_mat_val$sigma1 ^ 2
        par_mat_val$weight2 = par_mat_val$z2 / par_mat_val$sigma2 ^ 2
        par_mat_val$weight3 = par_mat_val$z3 / par_mat_val$sigma3 ^ 2
        par_mat_val$weight4 = par_mat_val$z4 / par_mat_val$sigma4 ^ 2
        par_mat_val$weight5 = par_mat_val$z5 / par_mat_val$sigma5 ^ 2
      }

      # null model for mu and p

      mu1_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
      mu4_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight4)
      mu5_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight5)
      par_mat$mu1 <- predict(mu1_bst, newdata = X)
      par_mat$mu2 <- predict(mu2_bst, newdata = X)
      par_mat$mu3 <- predict(mu3_bst, newdata = X)
      par_mat$mu4 <- predict(mu4_bst, newdata = X)
      par_mat$mu5 <- predict(mu5_bst, newdata = X)
      if (is.null(Yval)!=TRUE){
        par_mat_val$mu1 <- predict(mu1_bst, newdata = Xval)
        par_mat_val$mu2 <- predict(mu2_bst, newdata = Xval)
        par_mat_val$mu3 <- predict(mu3_bst, newdata = Xval)
        par_mat_val$mu4 <- predict(mu4_bst, newdata = Xval)
        par_mat_val$mu5 <- predict(mu5_bst, newdata = Xval)
      }

      p_bst <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3", "z4", "z5")]) ~ 1,
                 data = X,
                 trace = F)
      par_mat[, c("p1", "p2", "p3", "p4", "p5")] <-
        predict(p_bst, newdata = X, type = "probs")
      if (is.null(Yval)!=TRUE){
        par_mat_val[, c("p1", "p2", "p3", "p4", "p5")] <-
          predict(p_bst, newdata = Xval, type = "probs")
      }

      # boosting for mu

      if (structure == "mu"| structure == "both") {
        param<-list(max_depth=4, eta =0.5, objective="reg:squarederror")
        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight1)
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval, weight=par_mat_val$weight1)
        watchlist=list(train=dtrain, eval= dvalid)
        mu1_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu1:", mu1_bst$niter))

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight2)
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval, weight=par_mat_val$weight2)
        watchlist=list(train=dtrain, eval= dvalid)
        mu2_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu2:", mu2_bst$niter))

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight3)
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval, weight=par_mat_val$weight3)
        watchlist=list(train=dtrain, eval= dvalid)
        mu3_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu3:", mu3_bst$niter))

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight4)
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval, weight=par_mat_val$weight4)
        watchlist=list(train=dtrain, eval= dvalid)
        mu4_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu4:", mu4_bst$niter))

        dtrain<-xgb.DMatrix(data=as.matrix(X), label = Y,weight=par_mat$weight5)
        dvalid<-xgb.DMatrix(data=as.matrix(Xval), label= Yval, weight=par_mat_val$weight5)
        watchlist=list(train=dtrain, eval= dvalid)
        mu5_bst <-xgb.train(param, dtrain, nrounds=100,verbose = 0,watchlist,early_stopping_rounds = 2)
        print(paste("best iteration for mu5:", mu5_bst$niter))

        par_mat$mu1 <- predict(mu1_bst, newdata = dtrain)
        par_mat$mu2 <- predict(mu2_bst, newdata = dtrain)
        par_mat$mu3 <- predict(mu3_bst, newdata = dtrain)
        par_mat$mu4 <- predict(mu4_bst, newdata = dtrain)
        par_mat$mu5 <- predict(mu5_bst, newdata = dtrain)
        if (is.null(Yval)!=TRUE){
          par_mat_val$mu1 <- predict(mu1_bst, newdata = dvalid)
          par_mat_val$mu2 <- predict(mu2_bst, newdata = dvalid)
          par_mat_val$mu3 <- predict(mu3_bst, newdata = dvalid)
          par_mat_val$mu4 <- predict(mu4_bst, newdata = dvalid)
          par_mat_val$mu5 <- predict(mu5_bst, newdata = dvalid)
        }
      }

      # boosting for p

      if (structure == "p" |structure == "both") {
        p_bst <-
          BST(
            X = X,
            Y = par_mat[,c("z1","z2","z3","z4","z5")],
            Pinit = NULL,
            Xval = Xval,
            Yval = par_mat_val[, c("z1", "z2", "z3", "z4", "z5")],
            Pvalinit = NULL,
            M = n_tree_p,
            cp = 0.001,
            maxdepth = 4,
            lr = 0.2,
            trace =T,
            patience=2,
            parallel = parallel
          )
        par_mat[, c("p1", "p2", "p3", "p4", "p5")] <-
          p_bst$fitted[,c("BST_p1","BST_p2","BST_p3", "BST_p4", "BST_p5")]
        if (is.null(Yval)!=TRUE){
          par_mat_val[, c("p1", "p2", "p3","p4","p5")] <-
            p_bst$valid[,c("BST_p1","BST_p2","BST_p3", "BST_p4", "BST_p5")]
        }
      }

      # sigma estimation

      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))
      par_mat$sigma4 <-
        sqrt(sum(par_mat$z4 * (Y - par_mat$mu4) ^ 2) / sum(par_mat$z4))
      par_mat$sigma5 <-
        sqrt(sum(par_mat$z5 * (Y - par_mat$mu5) ^ 2) / sum(par_mat$z5))

      if (is.null(Yval)!=TRUE){
        par_mat_val$sigma1 <- unique(par_mat$sigma1)
        par_mat_val$sigma2 <- unique(par_mat$sigma2)
        par_mat_val$sigma3 <- unique(par_mat$sigma3)
        par_mat_val$sigma4 <- unique(par_mat$sigma4)
        par_mat_val$sigma5 <- unique(par_mat$sigma5)
      }

      # outer train loss and validation loss

      train_loss[m] <-
        neg_ll5(Y, par_mat[, c("mu1", "mu2", "mu3","mu4","mu5")],
                par_mat[, c("sigma1", "sigma2", "sigma3","sigma4","sigma5")],
                par_mat[, c("p1", "p2", "p3","p4","p5")])
      if (is.null(Yval)!=TRUE){
        valid_loss[m] <-
          neg_ll5(Yval, par_mat_val[, c("mu1", "mu2", "mu3","mu4","mu5")],
                  par_mat_val[, c("sigma1", "sigma2", "sigma3","sigma4","sigma5")],
                  par_mat_val[, c("p1", "p2", "p3","p4","p5")])
      }
      if (m>1){
        train_imp <- train_loss[m-1]-train_loss[m]
        train_impT[m] <- (train_imp>10^-4)
      }
      if (trace == T) {
        print(paste(
          "EB-iteration:",
          m,
          "; ",
          "train_loss:",
          round(train_loss[m], 4),
          "; ",
          "valid_loss:",
          round(valid_loss[m], 4),
          "; ",
          "train_loss_improve",
          round(train_imp, 4),
          "; ",
          train_impT[m],
          sep=""
        ))
      }
      if (m>patience){
        # minimize the training loss rather than early stop
        if (sum(train_impT[(m-patience+1):m])==0){
          break
        }
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3","plinear4","plinear5")] <-
      PtoF(par_mat[, c("p1", "p2", "p3","p4","p5")])
    if (is.null(Yval)!=TRUE){
      par_mat_val[, c("plinear1", "plinear2", "plinear3","plinear4","plinear5")] <-
        PtoF(par_mat_val[, c("p1", "p2", "p3","p4","p5")])
    }

    list(
      train_loss = train_loss,
      valid_loss = valid_loss,
      fitted = par_mat,
      valid = par_mat_val,
      mu_models = list(mu1_bst, mu2_bst, mu3_bst, mu4_bst, mu5_bst),
      mu_iter = c(mu1_bst$niter,mu2_bst$niter,mu3_bst$niter,mu4_bst$niter,mu5_bst$niter),
      p_model = p_bst,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3),
        unique(par_mat$sigma4),
        unique(par_mat$sigma5)
      )
    )
  }




