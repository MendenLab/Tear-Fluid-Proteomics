library(glmnet)
library(beeswarm)
library(pROC)
# generic cross-validation function to genrate splits
crossvalidation <- function(entity, nFold=5, random=T) {
  #-----------------------------------------------------------------------------------------------------------
  #DESC:    This function performs a "nFold" crossvalidation on a given input called "entity". 
  #         It is testing the created list "cv" and throws ERROR?s, if those are not passed(see
  #         ERROR-section below).
  #IN:      entity        => vector, which should be crossvalidated
  #         nFold         => numeric value in charge of which tpye of "nFold" crossvalidation will 
  #                          be applied to "entity"
  #         random        => this parameter decides if a seed will be set or not 
  #OUT:     Returns the list "cv", which contains as many lists as the size of "nFold". 
  #         Each element in "cv" is called "fold_XXX".
  #         Each of these "fold_XXX" contains three sets: "test", "xTrain" and "train"   
  #-----------------------------------------------------------------------------------------------------------
  
  # set seed to fixed value if not randomly selected
  if(!random){
    set.seed(754)
  }
  
  #creating a shuffled vector from the values of "entity"
  vec <- entity[sample(length(entity))]
  
  #defining the size of each set
  size <- floor(length(vec)/nFold)
  binSize <- rep(size, nFold)
  modulo <- length(vec) %% size
  if (modulo > 0) {
    binSize[1:modulo] <- binSize[1:modulo] +1
  }
  
  #setting the start indicies and the end indicies for the sets
  binEndIdx <- cumsum(binSize)
  binStardIdx <- c(1, binEndIdx[1:(nFold-1)]+1)
  
  #iterating over all the values of "vec" and putting them into a set
  #creating the list "cv" with the sublists and sets in them
  cv <- list()
  for (i in 1:nFold){
    
    testIdx <- binStardIdx[i]:binEndIdx[i]
    test <- vec[testIdx] #test
    xTrainIdx <- binStardIdx[(i %% nFold) + 1]:binEndIdx[(i %% nFold) + 1]
    xTrain <- vec[xTrainIdx] #validation
    train <- vec[setdiff(1:length(vec), union(testIdx, xTrainIdx))]
    
    cv[[paste("fold_", i, sep="")]] <- list(test=test, xTrain=xTrain, train=train)
  }

  return(cv) 
}


# cross-validation function for linear regression lambda parameter 
runCV <- function(feat, drug, method="elasticNet", plotPath="", mixPara=NA) {
  #-----------------------------------------------------------------------------------------------------------
  #DESC:    This function performs a 5 fold crossvalidation of a linear regression.
  #IN:      feat          => matrix containing molecular features (cell line times features)
  #         drug          => vector with drug response (names are cell line IDs)
  #         method        => either ridge, elasticNet or lasso
  #         plotPath      => optional path for cross-validation plot
  #OUT:     Returns the cross-validated lambda value
  #-----------------------------------------------------------------------------------------------------------
  if (method == "ridge")
    mixPara = 0
  if (method == "lasso")
    mixPara = 1
  if (method == "elasticNet" & is.na(mixPara))
    mixPara = 0.5
  
  cv <- crossvalidation(names(drug), 10)
  Rp <- list()
  lambda <- list()
  pred <- c()
  true <- c()
  for (i in 1:length(cv)) {
    # 1) train model
    fit <- glmnet(feat[cv[[i]]$train,], 
                  drug[cv[[i]]$train], alpha=mixPara, family = "binomial")
    lambda[[i]] <- fit$lambda
    # 2) predict on validation set
    xTrain_pred <- predict(fit,
                           type="response",
                           newx=feat[cv[[i]]$xTrain,])

    # 3) calculate performance on validaiton set
    # Check if both classes are present in the fold
    if (length(unique(drug[cv[[i]]$xTrain])) == 2) {
      Rp[[i]] <- apply(xTrain_pred, 2, function(x) if (sd(x)>0) {cor(x, drug[cv[[i]]$xTrain])} else{NA})
    } else {
      # Assign NA if only one class is present
      Rp[[i]] <- rep(0, ncol(xTrain_pred))
    }
    # 4) prediction on independent test set
    pred <- c(pred, predict(fit, feat[cv[[i]]$test,],
                            type="response")[,which.max(Rp[[i]])])
    true <- c(true, drug[cv[[i]]$test])
  }
  
  nMax <- max(unlist(lapply(Rp, length)))
  Rp <- lapply(Rp, function(x) c(x, rep(NA, nMax-length(x))))
  Rp <- matrix(unlist(Rp), nrow=length(Rp), byrow=T, dimnames = list(NULL,names(Rp[[1]])))
  
  # 5) chose representative L1-Norm and/or L2-Norm
  Rp_mean <- apply(Rp, 2, mean)
  sIDX_overall <- which.max(Rp_mean)  #This is the best lambda of all the cross validation?
  
  # # optional: plot of cross-validation
  # if (plotPath != "")
  #   pdf(plotPath, height=10)
  # layout(matrix(1:2, 2, 1))
  # l <- lambda[[which.max(unlist(lapply(lambda, length)))]]
  # Rp_error <- apply(Rp, 2, var)
  # plot(log(l), Rp_mean,
  #       main=paste(method, "10-fold cross-validation"),
  #       xlab="ln(lambda)",
  #       ylab="mean Pearson correlation", frame=F,
  #       ylim=range(c(Rp_mean-Rp_error, Rp_mean+Rp_error), na.rm=T))
  #  arrows(log(l), Rp_mean-Rp_error,
  #         log(l), Rp_mean+Rp_error,
  #         code=3,angle=90,length=0.1,col="lightgrey");
  #  abline(v=log(fit$l)[sIDX_overall], col="red")
  #  axis(3, log(fit$l)[sIDX_overall], names(sIDX_overall), col="red")
  # 
  # plot(drug[names(pred)], pred, frame=F, col="#0000FF55", pch=16,
  #      main="compare prediction on independent test set",
  #      xlab="observed drug response", ylab="predicted drug response",
  #      sub=paste("Pearson correlation =", formatC(cor(drug[names(pred)], pred, use="complete.obs"), digits=2)))

  # 6) # Generate and plot the ROC curve
  roc_obj <- roc(true, pred)
  auc_value <- auc(roc_obj)
  auc_ci <- ci.auc(roc_obj)
  
  # ROCAUC and CI
  # Combine AUC and CI into a single string
  auc_summary <- paste0("AUC = ", round(auc_value, 3), " ± (", 
                        round(auc_ci[3] - auc_ci[2], 3), ")")
  
  return(list(s=names(sIDX_overall),
              pred=pred[names(drug)],
              method=method,
              auc=c(auc_value, auc_ci[3] - auc_ci[2])))
}


# retrain and visualise the used features
reTrainModel <- function(feat, drug, cvBuild, plotPath="") {
  if (cvBuild$method == "ridge")
    mixPara = 0
  if (cvBuild$method == "lasso")
    mixPara = 1
  if (cvBuild$method == "elasticNet")
    mixPara = 0.5
  
  cell <- intersect(names(drug), rownames(feat))
  fit <- glmnet(feat[cell,], 
                drug[cell], alpha=mixPara, family = "binomial")
  
  beta <- fit$beta[,cvBuild$s]
  beta <- sort(beta[beta!=0])
  
  # optional: plot of cross-validation
  if (plotPath != "")
    pdf(plotPath, height=10)
  layout(matrix(1:2, 2, 1))
  plot(fit, xvar='lambda', col="grey",
       xlim=c(min(log(fit$lambda))-1, max(log(fit$lambda))))
  labels <- fit$beta[,ncol(fit$beta)]
  text(min(log(fit$lambda))-0.5, labels, names(labels), cex=0.5)
  abline(v=log(fit$lambda[as.numeric(gsub("s", "", cvBuild$s)) + 1]), lwd=2, col="red")
  
  barplot(beta, las=T, col=c("red", "blue")[as.factor(beta>0)], 
          border=F, ylab="weight", las=2, cex.names=0.7)
  
  if (plotPath != "")
    dev.off()
  return(beta = beta)
}


# check if models can be better than random
checkRandom <- function(feat, drug, method, plotPath="", n=10) {
  
  Rm <- c()
  Rr <- c()
  for (i in 1:n) {
    print(paste(" <- ", i))
    Rm <- c(Rm, cor(drug, runCV(feat, drug, method)$pred))
    Rr <- c(Rr, cor(drug, runCV(feat, sample(drug), method)$pred))
  }

  # calculate bayes factor
  probLead <- sum(Rm > Rr)
  probLead <- probLead / length(Rm)
  bayesFact <- probLead / (1-probLead)
  
  # plot histogram
  # dev.off()

  # Set layout and margins
  layout(1)
  par(mar = c(5, 5, 4, 2) + 0.1) 
  
  if (plotPath != "")
    pdf(plotPath)

  hist(Rm, xlim=range(Rm, Rr), col="#0000FF55",
       xlab="Pearson correlation", main="Compare random with predicted model",
       sub=paste("t-test p-val =", formatC(t.test(Rm, Rr)$p.value, digits=2),
                 "/ bayes factor =", formatC(bayesFact, digits=2)))
  hist(Rr, add=T, col="#FF000055")
  legend("top", c("random", "model"), pch = 15, col=c("#FF000055", "#0000FF55"), cex=0.7)

  if (plotPath != "")
       dev.off()
    
    return(bayesFact)
}


# zoom in on features with largest weight
plotWeights <- function(feat, drug, weights, plotPath="", n=5) {
  if (length(weights) < n)
    n <- length(weights)
  
  # optional: plot of cross-validation
  if (plotPath != "")
    pdf(plotPath, height=8)
  
  # mutual exclusivity plot
  if (n > 1) {
    layout(matrix(c(1,1,2,3,2,3), 2, 3, byrow=F))
    
    f <- feat[,names(sort(abs(weights), decreasing=T)[1:n])]
    f <- f[names(drug),]
    f <- t(f[,sort(colSums(f), index.return=T, decreasing=F)$ix])
    for (i in rev(1:nrow(f))) {
      f <- f[,sort(f[i,], index.return=T, decreasing=F)$ix]
    }
    image(t(f[,colSums(f)>0]), col=c("lightgrey", "black"), 
          xaxt='n', yaxt='n', frame=F)
    axis(2, seq(0,1,1/(n-1)), rownames(f), las=2, cex.axis=0.6)
    
    # weight plot 
    barplot(sort(abs(weights), 
                 decreasing=T)[1:n], las=2,
            col=c("red", "blue")[factor(sign(weights)>0)][sort(abs(weights), 
                                                               index.return=T,
                                                               decreasing=T)$ix][1:5],
            border=F, ylab="|weight|", cex.names=0.7)
    legend("topright", c("neg. weight", "pos. weight"), pch = 15, col=c("red", "blue"))
    
   if (plotPath != "")
     dev.off()
  }
}
