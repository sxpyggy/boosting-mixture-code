library(caret)
library(nnet)
library(rpart)
library(partykit)
library(mboost)

negLL<-function(y,p){
  # the negative log-likelihood for the multinomial distribution
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
  # optimal updates of node in multinomial boosting
  # see equation 32 in Friedman (2001)
  (K-1)/K * sum(y)/sum(abs(y)*(1-abs(y)))
}

FtoP<-function(FF){
  # F has dimension of n*K
  # From linear predictor to probs
  EF<-exp(FF)
  SEF<-apply(EF,1,sum)
  P<-EF/SEF
  P
}

PtoF<-function(P){
  # From probs to linear predictors
  logP<-log(P)
  logPa<-apply(logP,1,mean)
  FF<-logP-logPa
  FF
}

dat_BST<-function(K,n){
  # put into a suitable format
  # y, gradients (observed - prediction); f, linear predictor; p, prediction
  dat_bst<-data.frame(BST=matrix(NA,ncol=3*K,nrow=n))
  bst_names<-NULL
  for (k in 1:K){
    vnames<-c(paste("BST","_y",k,sep=""),paste("BST","_f",k,sep=""),paste("BST","_p",k,sep=""))
    bst_names<-c(bst_names,vnames)
  }
  names(dat_bst)<-bst_names
  dat_bst
}

sim_data<-function(n,size,seed){

  # simulation of the multinomial distributed variables
  set.seed(seed)
  # X1 = rnorm(n, 2, 1)
  # X2 = rexp(n, 2)
  # X3 = rbinom(n, 1, 0.5)
  # X4 = rgamma(n,0.5,0.5)
  # F1 = X1+log(X2)
  # F2 = 1-0.5*X1+2*X3-X1*X2
  # F3 = log(X2)+X3*X1+2*sin(X1)

  X1 = rnorm(n)
  X2 = rnorm(n)
  X3 = rnorm(n)
  X4 = rbinom(n,1,0.5)
  F1 = tanh(X1^2)
  F2 = tanh(X2^2)
  F3 = tanh(X3^2*X4)

  #pairs(cbind(F1,F2,F3,X1,X2,X3,X4))

  P<-FtoP(FF=cbind(F1,F2,F3))

  #pairs(cbind(P,X1,X2,X3,X4))
  #boxplot(P)

  dat<-
    data.frame(X1,X2,X3,X4,F1,F2,F3,P1=P[,1],P2=P[,2],P3=P[,3],Y1=NA,Y2=NA,Y3=NA)
  for (i in 1:nrow(dat)){
    dat[i,c("Y1","Y2","Y3")]<-
      t(rmultinom(1,size=size,dat[i,c("P1","P2","P3")]))/size
  }
  dat[,c("F1","F2","F3")]<-PtoF(dat[,c("P1","P2","P3")])
  dat

}

fit_tree <- function(X,dat_bst,Xval,val_bst,maxdepth,cp,
                      lr,Ind_y,Ind_f,K,k){
  # this function is used to perform one gradient boosting step, i.e., training
  # one weaker learner of tree for one class.
set.seed(1)
  tree <-rpart(dat_bst[, Ind_y[k]] ~ .,control = rpart.control(
    #minsplit = 20,
    minbucket = 50,
    cp = cp, maxcompete = 4, maxsurrogate = 0, usesurrogate = 2,xval = 0,
    surrogatestyle = 0, maxdepth = maxdepth), method = "anova", data = X)
  # dat_bst contains gradients(y), linear predictor(f), estimated probs (p)
  # Ind_y[k] is the gradient for the k-th class
  # cp is the complexity parameter, maxdepth is the depth of the weaker learners
  # X contains all the covariates

  fit_party <- as.party(tree) # the node index in as.party is different from the tree
  fitted_node <-fitted(fit_party) # node and responses in the node
  names(fitted_node) = c("node", "response")
  f_pred <-tapply(fitted_node$response, fitted_node$node, node_pre, K = K)
  # equation 32 in Firedman（2001）
  f_pred <-data.frame(f = as.vector(f_pred), node = as.numeric(names(f_pred)))

  if (is.null(Xval)!=TRUE){
    valid_node <-predict(fit_party, newdata = Xval, type = "node")
  }
  for (i in 1:nrow(dat_bst)) {
    dat_bst[i, Ind_f[k]] <-
      dat_bst[i, Ind_f[k]] + lr * f_pred$f[f_pred$node == fitted_node$node[i]]
    if (is.null(Xval)!=TRUE) {
      if (i <= nrow(val_bst)){
        # normally the size of validation is smaller than the training data
        val_bst[i, Ind_f[k]] <-
          val_bst[i, Ind_f[k]] + lr * f_pred$f[f_pred$node == valid_node[i]]
      }
    }
  }
  # update the linear predictor as f_{m-1}+ lr * f_{m}

  return(list(dat_bst=dat_bst,val_bst=val_bst,tree_save=fit_party,tree_0=tree))
  # tree_save is the party version of the tree
  # tree_0 is the rpart tree with different node index from fit_party
}

BST <- function(X, Y, Pinit, Xval, Yval, Pvalinit,
                         M, cp, maxdepth, lr, trace, patience,parallel) {
  # X, Xval: covariates. Y, Yval: multinomial response.
  # Pinit, Pvalinit: initial values for probs.
  # M: number of iterations.
  # cp, maxdepth: complex parameter, depth of weak learners of tree.
  # lr: learning rate in the boosting.
  # trace: show loss or not.
  # patience: early stop if validation loss is not decreased.
  set.seed(1)
  K <- ncol(Y)
  n <- nrow(Y)
  n_val <- nrow(Yval)
  val_bst <- NULL
  if (is.null(Yval)!=TRUE){
    val_bst <- dat_BST(K, n_val)
  }
  dat_bst <- dat_BST(K, n)

  # the index of gradients, linear predictor, estimated probabilities
  Ind_y <- which(grepl("y", names(dat_bst)) == T)
  Ind_f <- which(grepl("f", names(dat_bst)) == T)
  Ind_p <- which(grepl("p", names(dat_bst)) == T)

  # initialization
  dat_bst[, Ind_p]<- ifelse(is.null(Pinit),1/K ,Pinit)
  dat_bst[, Ind_f] <- PtoF(dat_bst[, Ind_p])
  dat_bst[, Ind_y] <- Y - dat_bst[, Ind_p]
  if (is.null(Yval)!=TRUE){
    val_bst[, Ind_p] <- ifelse(is.null(Pvalinit),1/K ,Pvalinit)
    val_bst[, Ind_f] <- PtoF(val_bst[, Ind_p])
  }
  Train_loss <- NULL
  Valid_loss <- NULL
  Valid_imp<- NULL # incremental improvement of validation loss
  Valid_imp[1] <- 0
  Valid_impT<-NULL
  Valid_impT[1]<-T
  Tree_save <- list()
  Tree_0<-list()
  # boosting
  for (m in 1:M) {
    dat_bst[, Ind_y] <- Y - dat_bst[, Ind_p]
    tree_save<-list()
    tree_0<-list()

    if (parallel==T){
      clusterExport(c1,varlist = list("node_pre"))
      res<- clusterApplyLB(c1,1:K,fun=fit_tree,X=X,dat_bst=dat_bst,Xval=Xval,
                          val_bst=val_bst,maxdepth=maxdepth,cp=cp,
                            lr=lr,Ind_y=Ind_y,Ind_f=Ind_f,K=K)}

    if (parallel==F){
      res<-list()
      for (k in 1:K){
      res[[k]]<-fit_tree(X,dat_bst,Xval,val_bst,maxdepth,cp,
                          lr,Ind_y,Ind_f,K,k)
      }
    }

    for (k in 1:K){
      tree_save[[k]]=res[[k]][["tree_save"]] # one tree
      tree_0[[k]]=res[[k]][["tree_0"]]
      dat_bst[,Ind_f[k]]<- res[[k]][["dat_bst"]][,Ind_f[k]]
      val_bst[,Ind_f[k]]<- res[[k]][["val_bst"]][,Ind_f[k]]
    }

    dat_bst[, Ind_p] <- FtoP(dat_bst[, Ind_f])
    if(is.null(Yval)!=TRUE){
      val_bst[, Ind_p] <- FtoP(val_bst[, Ind_f])
    }
    Train_loss[m] <- negLL(Y, dat_bst[, Ind_p])
    if(is.null(Yval)!=TRUE){
      Valid_loss[m] <- negLL(Yval, val_bst[, Ind_p])
    }
    if(is.null(Yval)!=TRUE&m>1){
      Valid_imp[m] <- Valid_loss[m-1]-Valid_loss[m]
      Valid_impT[m]<- Valid_imp[m]>10^-4
    }
    Tree_save[[m]]<- tree_save # K trees
    Tree_0[[m]]<- tree_0

    if(is.null(Yval)!=TRUE&trace==T){
      print(paste("iteration:",round(m,0),"; ","train loss:",
                  round(Train_loss[m],4),"; ","validation loss:",
                  round(Valid_loss[m],4),"; ","validation improve (10^-4):",
                  round(Valid_imp[m]*10^4,4),"; ", Valid_impT[m],sep=""
      ))
    }
    if(is.null(Yval)==TRUE&trace==T){
      print(paste("iteration:",round(m,0),"; ","train loss:",
                  round(Train_loss[m],4),"; ","validation loss: NULL",sep=""))

    }
    if(is.null(Yval)!=TRUE&m>patience){
      if (sum(Valid_impT[(m-patience+1):m])==0){
        break
      }
    }
  }
  list(Train_loss=Train_loss, Valid_loss=Valid_loss,
       fitted=dat_bst,
       valid=val_bst,
       Tree_save=Tree_save,
       Tree_0=Tree_0,
       lr=lr)
}

predict_BST <- function(X, bst, Pinit, M_best, type) {
  # X is the test covariates; bst is the output from function BST_parallel
  # Pinit is the initial values for the test data
  # M_best is the number of weak learners used.
  # type can be response ("response") or the linear predictor ("link")
  BST_model<-bst$Tree_save
  lr <- bst$lr
  M <- length(BST_model)
  M_best<-ifelse(is.null(M_best),M,M_best)
  # number of levels
  K <- unique(sapply(BST_model, length))
  # number of samples
  n <- nrow(X)
  test_bst <- dat_BST(K, n)
  Ind_y <- which(grepl("y", names(test_bst)) == T)
  Ind_f <- which(grepl("f", names(test_bst)) == T)
  Ind_p <- which(grepl("p", names(test_bst)) == T)
  test_bst[,Ind_p]<-ifelse(is.null(Pinit), 1/K ,Pinit)
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
  if (type=="response"){
    return(data.frame(test_bst[,Ind_p]))
  }
  if (type=="link"){
    return(data.frame(test_bst[,Ind_f]))
  }
}
