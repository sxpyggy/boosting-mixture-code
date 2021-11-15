library("gbm")
neg_ll3 <- function(y, mu, sigma, p) {
  ## y has dimension of n; mu, sigma, p have dimension of K*n ##
  return(-sum(log(
    p[, 1] * dnorm(y, mean = mu[, 1], sd = sigma[, 1]) +
      p[, 2] * dnorm(y, mean = mu[, 2], sd = sigma[, 2]) +
      p[, 3] * dnorm(y, mean = mu[, 3], sd = sigma[, 3])
  )) / length(y))
}

sim_gaussian <- function(n, seed) {
  set.seed(seed)
  X1 = rnorm(n, 2, 1)
  X2 = rexp(n, 2)
  X3 = rbinom(n, 1, 0.5)
  X4 = rgamma(n,0.5,0.5)
  F1 = X1+log(X2)
  F2 = 1-X1+2*X3-X1*X2
  F3 = log(X2)+X3*X1+2*sin(X1)
  P <- FtoP(FF = cbind(F1, F2, F3))
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
      F1,
      F2,
      F3,
      P1 = P[, 1],
      P2 = P[, 2],
      P3 = P[, 3],
      MU1,
      MU2,
      MU3,
      SG1 = 1,
      SG2 = 1,
      SG3 = 1,
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

EM_gaussian <-
  function(X, Y, Xtest, Ytest, M0, model, structure, trace) {
    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_test <- dat_Norm(K = 3, n = nrow(Xtest))

    #initialization

    par_mat$p1 <- 0.3
    par_mat$p2 <- 0.3
    par_mat$p3 <- 0.4
    par_mat$mu1 <- -10
    par_mat$mu2 <- 0
    par_mat$mu3 <- 15
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    learn_loss <- rep(0, M0)
    test_loss <- rep(0, M0)

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

      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2

      mu1_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_glm <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
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
      if (trace == T) {
        print(paste(
          "iteration:",
          m,
          ";",
          "learn_loss:",
          round(learn_loss[m], 4),
          ";",
          "test_loss:",
          round(test_loss[m], 4)
        ))
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

  function(X, Y, Xval, Yval, MUinit, Pinit, MUVinit, PVinit,
           style, M0, structure, trace) {

    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_val <- NULL
    if (is.null(Yval)!=TRUE){
      par_mat_val <- dat_Norm(K = 3, n = nrow(Xval))
    }

    #initialization

    par_mat[]

    par_mat$p1 <- 0.3
    par_mat$p2 <- 0.3
    par_mat$p3 <- 0.4
    par_mat$mu1 <- -10
    par_mat$mu2 <- 0
    par_mat$mu3 <- 15
    par_mat$sigma1 <- 0.5
    par_mat$sigma2 <- 0.5
    par_mat$sigma3 <- 0.5

    train_loss <- rep(0, M0)
    valid_loss <- rep(0, M0)

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

      par_mat$weight1 = par_mat$z1 / par_mat$sigma1 ^ 2
      par_mat$weight2 = par_mat$z2 / par_mat$sigma2 ^ 2
      par_mat$weight3 = par_mat$z3 / par_mat$sigma3 ^ 2

      mu1_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight1)
      mu2_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight2)
      mu3_bst <- lm(Y ~ 1, data = X, weights = par_mat$weight3)
      p_bst <-
        multinom(as.matrix(par_mat[, c("z1", "z2", "z3")]) ~ 1,
                 data = X,
                 trace = F)
      par_mat[, c("p1", "p2", "p3")] <-
        predict(p_bst, newdata = X, type = "probs")
      par_mat_val[, c("p1", "p2", "p3")] <-
        predict(p_bst, newdata = Xval, type = "probs")

      if (structure == "mu"| structure == "both") {
        mu1_bst <- gbm(Y ~ ., data = X, weights = par_mat$weight1,
                       distribution = "gaussian", n.trees = 50,
                       interaction.depth = 4, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 0.8, verbose = F)
        mu2_bst <- gbm(Y ~ ., data = X, weights = par_mat$weight2,
                       distribution = "gaussian", n.trees = 50,
                       interaction.depth = 4, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 0.8, verbose = F)
        mu3_bst <- gbm(Y ~ ., data = X, weights = par_mat$weight3,
                       distribution = "gaussian", n.trees = 50,
                       interaction.depth = 4, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 0.8, verbose = F)
      }

      if (structure == "p" |structure == "both") {
        p_bst <-
          BST(
            X = X,
            Y = par_mat[,c("z1","z2","z3")],
            Xval = Xval,
            Yval = par_mat_val[, c("z1", "z2", "z3")],
            M = 15,
            cp = 0.01,
            maxdepth = 3,
            lr = 0.2,
            trace=T
          )
        par_mat[, c("p1", "p2", "p3")] <-
          p_bst$fitted[,c("BST_p1","BST_p2","BST_p3")]
        par_mat_val[, c("p1", "p2", "p3")] <-
          p_bst$valid[,c("BST_p1","BST_p2","BST_p3")]
      }

      par_mat$mu1 <- predict(mu1_bst, newdata = X)
      par_mat$mu2 <- predict(mu2_bst, newdata = X)
      par_mat$mu3 <- predict(mu3_bst, newdata = X)

      par_mat$sigma1 <-
        sqrt(sum(par_mat$z1 * (Y - par_mat$mu1) ^ 2) / sum(par_mat$z1))
      par_mat$sigma2 <-
        sqrt(sum(par_mat$z2 * (Y - par_mat$mu2) ^ 2) / sum(par_mat$z2))
      par_mat$sigma3 <-
        sqrt(sum(par_mat$z3 * (Y - par_mat$mu3) ^ 2) / sum(par_mat$z3))

      par_mat_val$mu1 <- predict(mu1_bst, newdata = Xval)
      par_mat_val$mu2 <- predict(mu2_bst, newdata = Xval)
      par_mat_val$mu3 <- predict(mu3_bst, newdata = Xval)

      par_mat_val$sigma1 <- unique(par_mat$sigma1)
      par_mat_val$sigma2 <- unique(par_mat$sigma2)
      par_mat_val$sigma3 <- unique(par_mat$sigma3)

      train_loss[m] <-
        neg_ll3(Y, par_mat[, c("mu1", "mu2", "mu3")], par_mat[, c("sigma1", "sigma2", "sigma3")], par_mat[, c("p1", "p2", "p3")])
      valid_loss[m] <-
        neg_ll3(Yval, par_mat_val[, c("mu1", "mu2", "mu3")], par_mat_val[, c("sigma1", "sigma2", "sigma3")], par_mat_val[, c("p1", "p2", "p3")])
      if (trace == T) {
        print(paste(
          "iteration:",
          m,
          ";",
          "train_loss:",
          round(train_loss[m], 4),
          ";",
          "valid_loss:",
          round(valid_loss[m], 4)
        ))
      }
    }
    par_mat[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat[, c("p1", "p2", "p3")])
    par_mat_val[, c("plinear1", "plinear2", "plinear3")] <-
      PtoF(par_mat_val[, c("p1", "p2", "p3")])

    list(
      train_loss = train_loss,
      valid_loss = valid_loss,
      fitted = par_mat,
      valid = par_mat_val,
      mu_models = list(mu1_bst, mu2_bst, mu3_bst),
      p_model = p_bst,
      sigma = c(
        unique(par_mat$sigma1),
        unique(par_mat$sigma2),
        unique(par_mat$sigma3)
      )
    )
  }



