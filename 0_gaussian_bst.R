library("gbm")
neg_ll3 <- function(y, mu, sigma, p) {
  ## y has dimension of n; mu, sigma, p have dimension of K*n ##
  return(-sum(log(
    p[, 1] * dnorm(y, mean = mu[, 1], sd = sigma[, 1]) +
      p[, 2] * dnorm(y, mean = mu[, 2], sd = sigma[, 2]) +
      p[, 3] * dnorm(y, mean = mu[, 3], sd = sigma[, 3])
  )) / length(y))
}

EM_gaussian <-
  function(X, Y, Xval, Yval, M0, model, structure, trace) {
    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_val <- dat_Norm(K = 3, n = nrow(Xval))

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

      par_mat_val$mu1 <- predict(mu1_glm, newdata = Xval)
      par_mat_val$mu2 <- predict(mu2_glm, newdata = Xval)
      par_mat_val$mu3 <- predict(mu3_glm, newdata = Xval)

      par_mat_val[, c("p1", "p2", "p3")] <-
        predict(p_glm, newdata = Xval, type = "probs")

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
  function(X, Y, Xval, Yval, M0, structure, trace) {
    par_mat <- dat_Norm(K = 3, n = nrow(X))
    par_mat_val <- dat_Norm(K = 3, n = nrow(Xval))

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
                       bag.fraction = 1, train.fraction = 1, trace=T)
        mu2_bst <- gbm(Y ~ ., data = X, weights = par_mat$weight2,
                       distribution = "gaussian", n.trees = 50,
                       interaction.depth = 4, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 1, trace=T)
        mu3_bst <- gbm(Y ~ ., data = X, weights = par_mat$weight3,
                       distribution = "gaussian", n.trees = 50,
                       interaction.depth = 4, shrinkage = 0.1,
                       bag.fraction = 1, train.fraction = 1, trace=T)
      }

      if (structure == "p" |structure == "both") {
        p_bst <-
          BST(
            X = X,
            Y = par_mat[,c("z1","z2","z3")],
            Xval = Xval,
            Yval = par_mat_val[, c("z1", "z2", "z3")],
            M = 20,
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



