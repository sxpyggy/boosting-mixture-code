library(caret)
library(nnet)
library(rpart)
library(partykit)
library(mboost)

negLL<-function(y,p){
  # y and p are of dimension of K, i.e., a data point
  if (is.vector(y)==T){
    y<-t(y);p<-t(p)
  }
  # y and p are of dimension of n*K;
  if (is.vector(y)==F){
    y<-as.matrix(y)
    p<-as.matrix(p)
  }
  -mean(diag(y%*%t(log(p))))
}

node_pre<-function(y,K){
  # optimal updates of node in boosting
  (K-1)/K * sum(y)/sum(abs(y)*(1-abs(y)))
}

FtoP<-function(FF){
  # F has dimension of n*K
  EF<-exp(FF)
  SEF<-apply(EF,1,sum)
  P<-EF/SEF
  P
}

PtoF<-function(P){
  logP<-log(P)
  logPa<-apply(logP,1,mean)
  FF<-logP-logPa
  FF
}

dat_BST<-function(K,n){
  # put into a suitable format
  dat_bst<-data.frame(BST=matrix(NA,ncol=3*K,nrow=n))
  bst_names<-NULL
  for (k in 1:K){
    vnames<-c(paste("BST","_y",k,sep=""),paste("BST","_f",k,sep=""),paste("BST","_p",k,sep=""))
    bst_names<-c(bst_names,vnames)
  }
  names(dat_bst)<-bst_names
  dat_bst
}

dat_Norm<-function(K,n){
  # put into a suitable format
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

sim_data<-function(n,size,seed){
  set.seed(seed)
  X1 = rnorm(n, 2, 1)
  X2 = rexp(n, 2)
  X3 = rbinom(n, 1, 0.5)
  X4 = rgamma(n,0.5,0.5)
  F1 = X1+log(X2)
  F2 = 1-X1+2*X3-X1*X2
  F3 = log(X2)+X3*X1+2*sin(X1)

  P<-FtoP(FF=cbind(F1,F2,F3))
  dat<-data.frame(X1,X2,X3,X4,F1,F2,F3,P1=P[,1],P2=P[,2],P3=P[,3],Y1=NA,Y2=NA,Y3=NA)
  for (i in 1:nrow(dat)){
    dat[i,c("Y1","Y2","Y3")]<-t(rmultinom(1,size=size,dat[i,c("P1","P2","P3")]))/size
  }
  dat[,c("F1","F2","F3")]<-PtoF(dat[,c("P1","P2","P3")])
  dat
}

sim_gaussian <- function(n, seed) {
  set.seed(seed)
  X1 = rnorm(n, 2, 1)
  X2 = rexp(n, 2)
  X3 = rbinom(n, 1, 0.5)
  X4 = rgamma(n, 0.5, 0.5)
  F1 = X1 + log(X2)
  F2 = 1 - X1 + 2 * X3 - X1 * X2
  F3 = log(X2) + X3 * X1 + 2 * sin(X1)
  P <- FtoP(FF = cbind(F1, F2, F3))
  MU1 <- 7 + 2 * X1 + exp(0.3 * X2)
  MU2 <- 1 - 0.1 * X1 ^ 2 + 5 * X3 - X1 * X2
  MU3 <- -5 - 2 * X1 + 4 * sin(X2) + 0.2 * log(X4)
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
      SG1 = 3,
      SG2 = 1,
      SG3 = 3,
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

BST <- function(X, Y, Xval, Yval, M, cp, maxdepth, lr, trace) {
  K <- ncol(Y)
  n <- nrow(Y)
  n_val <- nrow(Yval)
  val_bst <- NULL
  if (is.null(Yval)!=TRUE){
  val_bst <- dat_BST(K, n_val)
  }
  dat_bst <- dat_BST(K, n)
  Ind_y <- which(grepl("y", names(dat_bst)) == T)
  Ind_f <- which(grepl("f", names(dat_bst)) == T)
  Ind_p <- which(grepl("p", names(dat_bst)) == T)

  # initialization
  dat_bst[, Ind_p] <- 1 / K
  dat_bst[, Ind_f] <- PtoF(dat_bst[, Ind_p])
  dat_bst[, Ind_y] <- Y - dat_bst[, Ind_p]
  if (is.null(Yval)!=TRUE){
  val_bst[, Ind_p] <- 1 / K
  val_bst[, Ind_f] <- PtoF(val_bst[, Ind_p])
  }
  Train_loss <- NULL
  Valid_loss <- NULL
  Tree_save <- list()

  # boosting
  for (m in 1:M) {
    dat_bst[, Ind_y] <- Y - dat_bst[, Ind_p]
    tree_save<-list()
    for (k in 1:K) {
      fit_tree <-
        rpart(
          dat_bst[, Ind_y[k]] ~ .,
          control = rpart.control(
            #minsplit = 20,
            #minbucket = round(20 / 3),
            cp = cp,
            maxcompete = 4,
            maxsurrogate = 0,
            usesurrogate = 2,
            xval = 0,
            surrogatestyle = 0,
            maxdepth = maxdepth
          ),
          method = "anova",
          data = X
        )
      fit_party <- as.party(fit_tree)
      tree_save[[k]] <- fit_party
      fitted_node <-
        fitted(fit_party) # node and responses in the node
      names(fitted_node) = c("node", "response")
      if (is.null(Yval)!=TRUE){
      valid_node <-
        predict(fit_party, newdata = Xval, type = "node")
      }
      f_pred <-
        tapply(fitted_node$response, fitted_node$node, node_pre, K = K)
      f_pred <-
        data.frame(f = as.vector(f_pred), node = as.numeric(names(f_pred))) # best updates to F, i.e., f
      for (i in 1:nrow(dat_bst)) {
        dat_bst[i, Ind_f[k]] <-
          dat_bst[i, Ind_f[k]] + lr * f_pred$f[f_pred$node == fitted_node$node[i]]
        if (is.null(Yval)!=TRUE) {
          if (i <= nrow(val_bst)){
          val_bst[i, Ind_f[k]] <-
            val_bst[i, Ind_f[k]] + lr * f_pred$f[f_pred$node == valid_node[i]]
          }
        }
      }
    }
    dat_bst[, Ind_p] <- FtoP(dat_bst[, Ind_f])
    if(is.null(Yval)!=TRUE){
    val_bst[, Ind_p] <- FtoP(val_bst[, Ind_f])
    }
    Train_loss[m] <- negLL(Y, dat_bst[, Ind_p])
    if(is.null(Yval)!=TRUE){
    Valid_loss[m] <- negLL(Yval, val_bst[, Ind_p])
    }
    Tree_save[[m]]<- tree_save
    if(is.null(Yval)!=TRUE&trace==T){
    print(paste("iteration:",round(m,0),";","train loss:",
                round(Train_loss[m],4),";","validation loss:",
                round(Valid_loss[m],4)))
    }
    if(is.null(Yval)==TRUE&trace==T){
      print(paste("iteration:",round(m,0),";","train loss:",
                  round(Train_loss[m],4),";","validation loss: NULL"))

    }
  }
  list(Train_loss=Train_loss, Valid_loss=Valid_loss,
       fitted=dat_bst,
       valid=val_bst,
       Tree_save=Tree_save, lr=lr)
}

predict_BST <- function(X, BST_model, M_best, lr) {
  # number of iterations
  M <- length(BST_model)
  # number of levels
  K <- unique(sapply(BST_model, length))
  # number of samples
  n <- nrow(X)
  test_bst <- dat_BST(K, n)
  Ind_y <- which(grepl("y", names(test_bst)) == T)
  Ind_f <- which(grepl("f", names(test_bst)) == T)
  Ind_p <- which(grepl("p", names(test_bst)) == T)
  test_bst[,Ind_p]<-1/K
  test_bst[,Ind_f]<-PtoF(test_bst[,Ind_p])
  Test_loss<-NULL

  for (m in 1:M_best) {
    for (k in 1:K) {
      fit_party <- BST_model[[m]][[k]]
      fitted_node <-
        fitted(fit_party) # node and responses in the node
      names(fitted_node) = c("node", "response")
      test_node <-
        predict(fit_party, newdata = X, type = "node")
      f_pred <-
        tapply(fitted_node$response, fitted_node$node, node_pre, K = K)
      f_pred <-
        data.frame(f = as.vector(f_pred), node = as.numeric(names(f_pred))) # best updates to F, i.e., f
      for (i in 1:nrow(test_bst)) {
        test_bst[i, Ind_f[k]] <-
          test_bst[i, Ind_f[k]] + lr * f_pred$f[f_pred$node == test_node[i]]
      }
    }
    test_bst[, Ind_p] <- FtoP(test_bst[, Ind_f])
  }
  data.frame(test_bst)
}
