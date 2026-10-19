# SSV summary
# Author: YingXiang
# Date: 28/05/2024

# =========================
# Load packages
# =========================
library(readr)
library(dplyr)
library(BiocManager)
library(pamr)
library(preprocessCore)
library(glmnet)
library(randomForest)
library(e1071)
library(HiDimDA)
library(ipred)
library(vsn)
library(MGSDA)
library(xgboost)
library(umap)
library(data.table)
library(Matrix)
library(lightgbm)
library(limma)
library(reshape)
library(ggplot2)
library(grid)
library(caret)
library(readxl)
library(sva)
library(EDASeq)
library(edgeR)
library(RUVSeq)
library(Biobase)
library(BiocGenerics)
library(tidyverse)
library(readxlsb)
library(ranger)
library(zoo)

# =========================
# Necessary parameters
# =========================
# raw_data: numeric dataframe (rows = features, columns = samples)
# log2_transform: whether the dataset needs log2 transform: TRUE or FALSE
# class_id: factor vector of class labels for each sample
#
# batch_id: batch id vector for each sample
# if there isn't batch information, use umap to classify different batches
#
# normalization function
# non-RNASeq: QN, MED, VSN
# RNASeq/miRNASeq: TMM, TC, UQ, MED, DESeq, QN, RUV*3

# =========================
# Frozen quantile normalization
# =========================
quant.norm <- function(data = NULL, train = NULL, test = NULL, ref.dis = NULL) {
  train <- as.matrix(data[, train, drop = FALSE])
  test  <- as.matrix(data[, test, drop = FALSE])
  
  if (!is.null(train) && !is.null(test)) {
    stopifnot(nrow(train) == nrow(test))
  }
  if (is.null(train)) {
    stopifnot(!is.null(ref.dis))
  }
  
  if (!is.null(train) && ncol(train) > 0) {
    train.qn <- preprocessCore::normalize.quantiles(train, copy = TRUE)
    dimnames(train.qn) <- dimnames(train)
    ref.dis <- as.numeric(sort(train.qn[, 1]))
  } else {
    train.qn <- NULL
  }
  
  if (is.null(test) || ncol(test) == 0) {
    test.fqn <- NULL
  } else {
    test.fqn <- apply(test, 2, function(x) {
      ord <- rank(x, ties.method = "average")
      ref.dis[ord]
    })
    test.fqn <- as.matrix(test.fqn)
    if (ncol(test.fqn) != ncol(test)) {
      test.fqn <- matrix(test.fqn, nrow = nrow(test), ncol = ncol(test))
    }
    rownames(test.fqn) <- rownames(test)
    colnames(test.fqn) <- colnames(test)
  }
  
  return(list(
    train.qn = train.qn,
    test.fqn = test.fqn,
    ref.dis = ref.dis
  ))
}

# =========================
# Median normalization
# =========================
med.norm <- function(data = NULL, train = NULL, test = NULL, ref.dis = NULL) {
  train <- as.matrix(data[, train, drop = FALSE])
  test  <- as.matrix(data[, test, drop = FALSE])
  
  if (!is.null(train) && !is.null(test)) {
    stopifnot(nrow(train) == nrow(test))
  }
  if (is.null(train)) {
    stopifnot(!is.null(ref.dis))
  }
  
  if (!is.null(train) && ncol(train) > 0) {
    ref.dis <- median(train, na.rm = TRUE)
    train_median <- apply(train, 2, median, na.rm = TRUE)
    train_shift <- train_median - ref.dis
    shifts.train <- matrix(
      rep(train_shift, each = nrow(train)),
      nrow = nrow(train),
      ncol = ncol(train)
    )
    train.mn <- train - shifts.train
    rownames(train.mn) <- rownames(train)
    colnames(train.mn) <- colnames(train)
  } else {
    train.mn <- NULL
  }
  
  if (is.null(test) || ncol(test) == 0) {
    test.fmn <- NULL
  } else {
    test_median <- apply(test, 2, median, na.rm = TRUE)
    test_shift <- test_median - ref.dis
    shifts.test <- matrix(
      rep(test_shift, each = nrow(test)),
      nrow = nrow(test),
      ncol = ncol(test)
    )
    test.fmn <- test - shifts.test
    rownames(test.fmn) <- rownames(test)
    colnames(test.fmn) <- colnames(test)
  }
  
  return(list(
    train.mn = train.mn,
    test.fmn = test.fmn,
    ref.dis = ref.dis
  ))
}

# =========================
# VSN normalization
# =========================
vs.norm <- function(data = NULL, train = NULL, test = NULL, ref.dis = NULL) {
  train <- as.matrix(data[, train, drop = FALSE])
  test  <- as.matrix(data[, test, drop = FALSE])
  
  if (!is.null(train) && !is.null(test)) {
    stopifnot(nrow(train) == nrow(test))
  }
  if (is.null(train)) {
    stopifnot(!is.null(ref.dis))
  }
  
  if (!is.null(train) && ncol(train) > 0) {
    train0 <- 2^train
    ref.dis <- vsn::vsn2(train0)
    train.vsn <- log2(exp(as.matrix(ref.dis)))
    rownames(train.vsn) <- rownames(train)
    colnames(train.vsn) <- colnames(train)
  } else {
    train.vsn <- NULL
  }
  
  if (is.null(test) || ncol(test) == 0) {
    test.fvsn <- NULL
  } else {
    test0 <- 2^test
    test.fvsn0 <- vsn::vsn2(test0, ref.dis)
    test.fvsn <- log2(exp(as.matrix(test.fvsn0)))
    rownames(test.fvsn) <- rownames(test)
    colnames(test.fvsn) <- colnames(test)
  }
  
  return(list(
    train.vsn = train.vsn,
    test.fvsn = test.fvsn,
    ref.dis = ref.dis
  ))
}

# =========================
# RNA-seq normalization class and methods
# =========================
precision <- methods::setClass(
  "precision",
  slots = c(
    raw.train.data = "ANY",
    raw.test.data = "ANY",
    harmon.train.data = "list",
    harmon.test.data = "list",
    classification.result = "list",
    cluster.result = "list"
  )
)

create_precision.classification <- function(traindata, testdata, trainlabel, testlabel) {
  object <- methods::new(
    Class = "precision",
    raw.train.data = list(data = traindata, label = trainlabel),
    raw.test.data = list(data = testdata, label = testlabel)
  )
  return(object)
}

harmon.method.TMM <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw, group = factor(groups), genes = rownames(raw))
  d <- calcNormFactors(dat.DGE, method = "TMM")
  scaling.factor <- d$samples$norm.factors * d$samples$lib.size / 1e6
  dat.harmonized <- cpm(d)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.method.TMM.frozen <- function(rawtrain, train_groups, rawtest, test_groups) {
  dat.DGE <- DGEList(
    counts = cbind(rawtrain[, 1, drop = FALSE], rawtest),
    group = factor(c(train_groups[1], test_groups)),
    genes = rownames(rawtest)
  )
  d <- calcNormFactors(dat.DGE, method = "TMM", refColumn = 1)
  scaling.factor <- (d$samples$norm.factors * d$samples$lib.size / 1e6)[-1]
  dat.harmonized <- cpm(d)[, -1, drop = FALSE]
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.TMM <- function(object) {
  object@harmon.train.data$TMM <- harmon.method.TMM(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$TMM <- harmon.method.TMM(
      raw = object@raw.test.data$data,
      groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.TMM.frozen <- function(object) {
  object@harmon.train.data$TMM.frozen <- harmon.method.TMM(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$TMM.frozen <- harmon.method.TMM.frozen(
      rawtrain = object@raw.train.data$data,
      train_groups = object@raw.train.data$label,
      rawtest = object@raw.test.data$data,
      test_groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.method.TC <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw, group = factor(groups), genes = rownames(raw))
  scaling.factor <- dat.DGE$samples$lib.size / 1e6
  dat.harmonized <- cpm(dat.DGE, normalized.lib.sizes = FALSE)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.TC <- function(object) {
  object@harmon.train.data$TC <- harmon.method.TC(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$TC <- harmon.method.TC(
      raw = object@raw.test.data$data,
      groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.method.UQ <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw, group = factor(groups), genes = rownames(raw))
  q.factor <- apply(dat.DGE$counts, 2, function(x) quantile(x[x != 0], probs = 0.75))
  scaling.factor <- q.factor / 1e6
  dat.harmonized <- t(t(raw) / scaling.factor)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.UQ <- function(object) {
  object@harmon.train.data$UQ <- harmon.method.UQ(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$UQ <- harmon.method.UQ(
      raw = object@raw.test.data$data,
      groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.method.med <- function(raw, groups) {
  dat.DGE <- DGEList(counts = raw, group = factor(groups), genes = rownames(raw))
  m.factor <- apply(dat.DGE$counts, 2, function(x) median(x[x != 0]))
  scaling.factor <- m.factor / 1e6
  dat.harmonized <- t(t(raw) / scaling.factor)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.med <- function(object) {
  object@harmon.train.data$med <- harmon.method.med(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$med <- harmon.method.med(
      raw = object@raw.test.data$data,
      groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.method.DESeq <- function(raw, groups) {
  colnames(raw) <- paste0("V", 1:ncol(raw))
  condition <- data.frame(SampleName = colnames(raw), Condition = factor(groups))
  rownames(condition) <- colnames(raw)
  dat.DGE <- DESeq2::DESeqDataSetFromMatrix(countData = raw, colData = condition, design = ~ Condition)
  dat.DGE <- DESeq2::estimateSizeFactors(dat.DGE)
  scaling.factor <- DESeq2::sizeFactors(dat.DGE)
  dat.harmonized <- DESeq2::counts(dat.DGE, normalized = TRUE)
  res <- list(dat.harmonized = dat.harmonized, scaling.factor = scaling.factor)
  return(res)
}

harmon.DESeq <- function(object) {
  object@harmon.train.data$DESeq <- harmon.method.DESeq(
    raw = object@raw.train.data$data,
    groups = object@raw.train.data$label
  )
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$DESeq <- harmon.method.DESeq(
      raw = object@raw.test.data$data,
      groups = object@raw.test.data$label
    )
  }
  return(object)
}

harmon.method.QN <- function(raw) {
  dat.harmonized <- preprocessCore::normalize.quantiles(as.matrix(raw))
  colnames(dat.harmonized) <- colnames(raw)
  rownames(dat.harmonized) <- rownames(raw)
  res <- list(dat.harmonized = dat.harmonized)
  return(res)
}

harmon.method.QN.frozen <- function(rawtrain, rawtest) {
  rawtrain <- as.matrix(rawtrain)
  rawtest  <- as.matrix(rawtest)
  
  dat.harmonized <- preprocessCore::normalize.quantiles(rawtrain)
  rownames(dat.harmonized) <- rownames(rawtrain)
  colnames(dat.harmonized) <- colnames(rawtrain)
  
  ref.dis <- as.numeric(sort(dat.harmonized[, 1]))
  
  dat.harmonized.test <- apply(rawtest, 2, function(x) {
    ord <- rank(x, ties.method = "average")
    ref.dis[ord]
  })
  dat.harmonized.test <- as.matrix(dat.harmonized.test)
  if (ncol(dat.harmonized.test) != ncol(rawtest)) {
    dat.harmonized.test <- matrix(dat.harmonized.test, nrow = nrow(rawtest), ncol = ncol(rawtest))
  }
  rownames(dat.harmonized.test) <- rownames(rawtest)
  colnames(dat.harmonized.test) <- colnames(rawtest)
  
  res <- list(dat.harmonized = dat.harmonized.test)
  return(res)
}

harmon.QN <- function(object) {
  object@harmon.train.data$QN <- harmon.method.QN(raw = object@raw.train.data$data)
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$QN <- harmon.method.QN(raw = object@raw.test.data$data)
  }
  return(object)
}

harmon.QN.frozen <- function(object) {
  object@harmon.train.data$QN.frozen <- harmon.method.QN(raw = object@raw.train.data$data)
  if (!is.null(object@raw.test.data)) {
    object@harmon.test.data$QN.frozen <- harmon.method.QN.frozen(
      rawtrain = object@raw.train.data$data,
      rawtest = object@raw.test.data$data
    )
  }
  return(object)
}

# =========================
# ML functions
# pam, svm, lasso, xgboost, knn, random forest
# =========================

# PAM
pam.intcv <- function(X, y, vt.k = NULL, n.k = 30, kfold = 5, folds = NULL, seed) {
  ptm <- proc.time()
  set.seed(seed)
  
  X <- as.matrix(X)
  y <- factor(as.vector(y))
  
  if (ncol(X) != length(y)) {
    stop(paste0("PAM input mismatch: ncol(X)=", ncol(X), ", length(y)=", length(y)))
  }
  
  if (is.null(rownames(X))) {
    rownames(X) <- paste0("feature_", seq_len(nrow(X)))
  }
  
  data.pam <- list(
    x = X,
    y = y,
    geneids = rownames(X),
    genenames = rownames(X)
  )
  
  fit.pam <- pamr::pamr.train(
    data.pam,
    threshold = vt.k,
    n.threshold = n.k
  )
  
  fit.cv <- new.pamr.cv(
    fit = fit.pam,
    data = data.pam,
    nfold = kfold
  )
  
  best.threshold <- fit.cv$threshold[max(which(fit.cv$error == min(fit.cv$error, na.rm = TRUE)))]
  mc <- fit.cv$error[which.min(fit.cv$error)]
  
  model <- pamr::pamr.train(
    data.pam,
    threshold = best.threshold,
    n.threshold = n.k
  )
  
  # Do not call pamr.listgenes() inside the ML loop.
  # In small/imbalanced splits, pamr.listgenes() can throw a dimnames error
  # even when pamr.train() and pamr.predict() are usable for error estimation.
  coefs <- NULL
  
  time <- proc.time() - ptm
  return(list(mc = mc, time = time, model = model, cfs = coefs))
}

new.pamr.cv <- function(fit, data, nfold = 5, ...) {
  x <- data$x[fit$gene.subset, fit$sample.subset, drop = FALSE]
  
  if (is.null(fit$newy)) {
    y <- factor(data$y[fit$sample.subset])
  } else {
    y <- factor(data$newy[fit$sample.subset])
  }
  
  this.call <- match.call()
  nsccv2 <- get("nsccv", envir = asNamespace("pamr"))
  balanced.folds <- get("balanced.folds", envir = asNamespace("pamr"))
  folds <- balanced.folds(y, nfolds = nfold)
  
  junk <- nsccv2(
    x, y,
    object = fit,
    folds = folds,
    survival.time = data$survival.time,
    censoring.status = data$censoring.status,
    ngroup.survival = fit$ngroup.survival,
    problem.type = fit$problem.type,
    ...
  )
  
  junk$call <- this.call
  junk$newy <- fit$newy
  junk$sample.subset <- fit$sample.subset
  return(junk)
}

trycatch.func <- function(expr, msg = "") {
  out <- tryCatch(
    {
      expr
    },
    warning = function(cond) {
      message("warning: ")
      message(cond)
      return(NULL)
    },
    error = function(cond) {
      message("error: ")
      message(cond)
      return(NA)
    },
    finally = {
      message(msg)
    }
  )
  return(out)
}

pam.predict <- function(pam.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- pamr::pamr.predict(
    pam.intcv.model$model,
    newx = pred.obj,
    threshold = pam.intcv.model$model$threshold,
    type = "class"
  )
  
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  
  prob <- pamr::pamr.predict(
    pam.intcv.model$model,
    newx = pred.obj,
    threshold = pam.intcv.model$model$threshold,
    type = "posterior"
  )
  
  return(list(pred = pred, mc = mc, prob = prob))
}

# SVM
svm.intcv <- function(kfold = 5, X, y, seed) {
  ptm <- proc.time()
  set.seed(seed)
  
  svm_tune <- tune.svm(
    x = data.matrix(X),
    y = factor(y),
    tunecontrol = tune.control(cross = kfold)
  )
  
  time <- proc.time() - ptm
  return(list(mc = svm_tune$best.performance, time = time, model = svm_tune$best.model))
}

svm.predict <- function(svm.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(svm.intcv.model$model, newdata = t(pred.obj))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  return(list(pred = pred, mc = mc))
}

# LASSO
lasso.intcv <- function(kfold = 5, X, y, seed, alp = 1) {
  ptm <- proc.time()
  set.seed(seed)
  
  cv.fit <- glmnet::cv.glmnet(
    x = data.matrix(t(X)),
    y = factor(y),
    family = "binomial",
    type.measure = "class",
    alpha = alp,
    nfold = kfold
  )
  
  #mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.1se)]
  #coefs <- trycatch.func(coef(cv.fit, s = "lambda.1se"))
  
  #consistent with the PRECISION package’s CV threshold method
  mc <- cv.fit$cvm[which(cv.fit$lambda == cv.fit$lambda.min)]
  coefs <- trycatch.func(coef(cv.fit, s = "lambda.min"))
  
  time <- proc.time() - ptm
  
  return(list(mc = mc, time = time, model = cv.fit, cfs = coefs))
}

lasso.predict <- function(lasso.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(
    lasso.intcv.model$model,
    newx = t(pred.obj),
    s = lasso.intcv.model$model$lambda.1se,
    type = "class"
  )
  
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  
  prob <- predict(
    lasso.intcv.model$model,
    newx = t(pred.obj),
    s = lasso.intcv.model$model$lambda.1se
  )
  
  return(list(pred = pred, mc = mc, prob = prob))
}

# XGBoost
cvxgb <- function(datacv, fold) {
  minerror <- Inf
  best_param <- NULL
  index <- 1
  
  for (min_child_weightvalue in seq(1, 9, 4)) {
    param <- list(
      min_child_weight = min_child_weightvalue,
      objective = "binary:logistic",
      eval_metric = "error"
    )
    
    cvout <- xgb.cv(
      data = datacv,
      nrounds = 100,
      nfold = fold,
      params = param,
      maximize = FALSE,
      verbose = FALSE,
      early_stopping_rounds = 20
    )
    
    err <- cvout$evaluation_log$test_error_mean
    
    if (min(err, na.rm = TRUE) < minerror) {
      minerror <- min(err, na.rm = TRUE)
      best_param <- param
      index <- which.min(err)
    }
  }
  
  return(list(param = best_param, nround = index, cverror = minerror))
}

xgboost.intcv <- function(kfold = 5, X, y, seed) {
  ptm <- proc.time()
  set.seed(seed)
  
  id <- as.numeric(as.character(y))
  dtrain <- xgb.DMatrix(data = t(X), label = id)
  cvresult <- cvxgb(datacv = dtrain, fold = kfold)
  bestparam <- cvresult$param
  print("best param is")
  print(bestparam)
  
  nround <- cvresult$nround
  trainresult <- xgb.train(
    params = bestparam,
    data = dtrain,
    nrounds = nround,
    verbose = 0
  )
  
  time <- proc.time() - ptm
  return(list(mc = cvresult$cverror, time = time, model = trainresult))
}

xgboost.predict <- function(xgboost.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(xgboost.intcv.model$model, t(pred.obj))
  prediction <- as.numeric(pred > 0.5)
  pred_class <- factor(prediction, levels = c(0, 1))
  mc <- tabulate.ext.err.func(pred_class, pred.obj.group.id)
  return(list(pred = pred_class, mc = mc, prob = pred))
}

# KNN
knn.intcv <- function(kfold = 5, X, y, seed) {
  ptm <- proc.time()
  set.seed(seed)
  
  ctrl <- trainControl(
    method = "repeatedcv",
    repeats = 3,
    number = kfold
  )
  
  knn <- train(
    x = data.matrix(t(X)),
    y = factor(y),
    method = "knn",
    tuneLength = 9,
    trControl = ctrl
  )
  
  time <- proc.time() - ptm
  return(list(mc = 1 - max(knn$results$Accuracy), time = time, model = knn, cfs = NULL))
}

knn.predict <- function(knn.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(knn.intcv.model$model, newdata = data.matrix(t(pred.obj)))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  return(list(pred = pred, mc = mc))
}

# Random Forest
ranfor.intcv <- function(kfold = 5, X, y, seed) {
  ptm <- proc.time()
  set.seed(seed)
  
  control <- trainControl(
    method = "cv",
    number = 5,
    search = "random"
  )
  
  rf <- train(
    x = data.matrix(t(X)),
    y = factor(y),
    method = "ranger",
    metric = "Accuracy",
    tuneLength = 5,
    preProcess = c("center", "scale"),
    trControl = control
  )
  
  time <- proc.time() - ptm
  return(list(mc = 1 - max(rf$results$Accuracy), time = time, model = rf, cfs = NULL))
}

ranfor.predict <- function(ranfor.intcv.model, pred.obj, pred.obj.group.id) {
  pred <- predict(ranfor.intcv.model$model, newdata = t(pred.obj))
  mc <- tabulate.ext.err.func(pred, pred.obj.group.id)
  return(list(pred = pred, mc = mc))
}

tabulate.ext.err.func <- function(pred.obj, obs.grp) {
  return(1 - sum(diag(table(pred.obj, obs.grp))) / length(obs.grp))
}

# =========================
# ML loop
# =========================
# n: times of ml
# data: log2 transformed dataframe (row:features, col:samples)
# class: a vector containing class id (factor 0 and 1)
# batch: a vector containing batch id
# batch_time: whether to divide batches by time
# early_batch_id: if batch_time is not FALSE, plug in the early batch id
# batch_balance: should the two batches be trimmed to the same length
# ml: ml methods
# split_combination: a vector containing splitting strategy, e.g. c(50,70,80,90,100)
# data_type_rnaseq: whether the data is RNA-seq
# file: prefix of the result file name

# make sure that the colnames of data is the same as the names of batch
ml_sum <- function(
    n = times,
    data,
    class,
    batch,
    batch_time = FALSE,
    early_batch_id,
    batch_balance = FALSE,
    ml,
    split_combination,
    data_type_rnaseq = FALSE,
    file,
    norm_custom = NULL
) {
  dim1 <- as.character(split_combination)
  
  if (data_type_rnaseq == TRUE) {
    dim2 <- c("non", "TC", "UQ", "med", "TMM", "TMM.frozen", "DESeq", "QN", "QN.frozen")
  } else {
    dim2 <- c("non", "qn", "mn", "vsn")
  }
  if (!is.null(norm_custom)) {
    dim2 <- norm_custom
  }
  
  dim3 <- as.character(seq_len(n))
  
  for (ml_method in ml) {
    eval(parse(
      text = paste0(
        ml_method,
        "_result = array(NA, dim = c(length(dim1), length(dim2), length(dim3)), dimnames = list(dim1, dim2, dim3))"
      )
    ))
  }
  
  class <- unlist(class)
  
  for (k in seq_len(n)) {
    set.seed(k + 2024)
    unique_batch <- unique(batch)
    
    if (batch_time) {
      batch1 <- which(batch %in% early_batch_id)
      batch2 <- which(!batch %in% early_batch_id)
    } else {
      unique_chars <- unique(batch)
      shuffled_unique_chars <- sample(unique_chars)
      half_length <- ceiling(length(shuffled_unique_chars) / 2)
      group1_chars <- shuffled_unique_chars[1:half_length]
      group2_chars <- shuffled_unique_chars[(half_length + 1):length(shuffled_unique_chars)]
      batch1 <- which(batch %in% group1_chars)
      batch2 <- which(batch %in% group2_chars)
    }
    
    for (split in split_combination) {
      i <- which(split == split_combination)
      
      if (batch_balance == TRUE) {
        length_diff <- abs(length(batch1) - length(batch2))
        if (length(batch1) > length(batch2)) {
          batch1_id <- batch1[-sample(seq_along(batch1), length_diff)]
          batch2_id <- batch2
        } else {
          batch1_id <- batch1
          batch2_id <- batch2[-sample(seq_along(batch2), length_diff)]
        }
      } else {
        batch1_id <- batch1
        batch2_id <- batch2
      }
      
      train_2 <- sample(batch2_id, length(batch2_id) * (100 - split) / 100)
      train_1 <- sample(batch1_id, length(batch1_id) * split / 100)
      train_id <- sort(c(train_2, train_1))
      
      test_2 <- batch2_id[!(batch2_id %in% train_2)]
      test_1 <- batch1_id[!(batch1_id %in% train_1)]
      test_id <- sort(c(test_2, test_1))
      
      train_label <- class[train_id]
      test_label <- class[test_id]
      
      if (data_type_rnaseq == TRUE) {
        object <- create_precision.classification(
          traindata = data[, train_id],
          testdata = data[, test_id],
          train_label,
          test_label
        )
        object <- harmon.TC(object)
        object <- harmon.UQ(object)
        object <- harmon.med(object)
        object <- harmon.TMM(object)
        object <- harmon.TMM.frozen(object)
        object <- harmon.DESeq(object)
        object <- harmon.QN(object)
        object <- harmon.QN.frozen(object)
      }
      
      for (norm in dim2) {
        # Use the actual position of the current normalization in dim2.
        # This prevents out-of-bound indexing when norm_custom is used.
        j <- match(norm, dim2)
        if (is.na(j)) {
          stop(paste("Normalization", norm, "not found in dim2"))
        }
        
        if (data_type_rnaseq == TRUE) {
          if (norm == "non") {
            train <- object@raw.train.data$data
            test <- object@raw.test.data$data
          } else {
            eval(parse(text = paste0("train <- object@harmon.train.data$", norm, "$dat.harmonized")))
            eval(parse(text = paste0("test <- object@harmon.test.data$", norm, "$dat.harmonized")))
          }
          
          train <- as.matrix(train)
          test <- as.matrix(test)
          train_label <- object@raw.train.data$label
          test_label <- object@raw.test.data$label
        } else {
          if (norm == "non") {
            test <- data[, test_id, drop = FALSE]
            train <- data[, train_id, drop = FALSE]
          }
          
          if (norm == "qn") {
            qn_result <- quant.norm(data, train_id, test_id)
            test <- qn_result$test.fqn
            train <- qn_result$train.qn
          }
          
          if (norm == "mn") {
            mn_result <- med.norm(data, train_id, test_id)
            test <- mn_result$test.fmn
            train <- mn_result$train.mn
          }
          
          if (norm == "vsn") {
            vsn_result <- vs.norm(data, train_id, test_id)
            test <- vsn_result$test.fvsn
            train <- vsn_result$train.vsn
          }
          
          train <- as.matrix(train)
          test <- as.matrix(test)
        }
        
        for (ml_method in ml) {
          if (ml_method == "knn") {
            knn.intcv.model <- knn.intcv(kfold = 5, X = train, y = train_label, seed = k)
            pred.obj <- test
            pred <- predict(knn.intcv.model$model, newdata = data.matrix(t(pred.obj)))
            error <- mean(pred != test_label)
          }
          
          if (ml_method == "lasso") {
            lasso.intcv.model <- lasso.intcv(X = train, y = train_label, seed = k)
            error <- mean(
              predict(
                lasso.intcv.model$model,
                newx = t(test),
                s = lasso.intcv.model$model$lambda.1se,
                type = "class"
              ) != test_label
            )
          }
          
          if (ml_method == "pam") {
            pam.intcv.model <- pam.intcv(kfold = 5, X = as.matrix(train), y = train_label, seed = k)
            error <- mean(
              pamr.predict(
                pam.intcv.model$model,
                newx = test,
                threshold = pam.intcv.model$model$threshold,
                type = "class"
              ) != test_label
            )
          }
          
          if (ml_method == "rf") {
            ranfor.intcv.model <- ranfor.intcv(kfold = 5, X = train, y = train_label, seed = k)
            pred.obj <- test
            pred <- predict(ranfor.intcv.model$model, newdata = t(pred.obj))
            error <- mean(pred != test_label)
          }
          
          if (ml_method == "svm") {
            svm.intcv.model <- svm.intcv(kfold = 5, X = t(train), y = train_label, seed = k)
            pred.obj <- test
            error <- mean(predict(svm.intcv.model$model, newdata = t(pred.obj)) != test_label)
          }
          
          if (ml_method == "xgb") {
            xgboost.intcv.model <- xgboost.intcv(
              X = train,
              y = as.numeric(train_label) - 1,
              seed = k
            )
            pred <- predict(xgboost.intcv.model$model, t(test))
            pred_class <- factor(as.numeric(pred > 0.5), levels = c(0, 1))
            error <- mean(pred_class != factor(test_label, levels = c(0, 1)))
          }
          
          eval(parse(text = paste0(ml_method, "_result[i, j, k] = error")))
        }
      }
    }
  }
  
  for (ml_method in ml) {
    eval(parse(text = paste0("result <<- ", ml_method, "_result")))
    filename <- paste(
      file,
      n,
      ml_method,
      "batch_balance",
      batch_balance,
      "result.csv",
      sep = "_"
    )
    write.table(result, file = filename, sep = "\t", col.names = NA, quote = FALSE)
  }
}

# =========================
# Result plot
# deprecate this function.
# =========================
result_plot <- function(
    n,
    ml,
    file,
    batch_balance = FALSE,
    data_type_rnaseq = FALSE,
    split_combination
) {
  error <- vector("numeric")
  se <- vector("numeric")
  
  for (i in ml) {
    if (data_type_rnaseq == FALSE) {
      norm <- c("None", "QN", "MN", "VSN")
      filename <- paste(file, n, i, "batch_balance", batch_balance, "result.csv", sep = "_")
      eval(parse(text = paste0(i, "_result = read.csv('", filename, "', sep = '\t')")))
      eval(parse(text = paste0(i, "_result_non = rowMeans((", i, "_result[, -1])[, c(T, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_qn = rowMeans((", i, "_result[, -1])[, c(F, T, F, F)])")))
      eval(parse(text = paste0(i, "_result_mn = rowMeans((", i, "_result[, -1])[, c(F, F, T, F)])")))
      eval(parse(text = paste0(i, "_result_vsn = rowMeans((", i, "_result[, -1])[, c(F, F, F, T)])")))
      
      eval(parse(text = paste0(i, "_result_non_se = apply((", i, "_result[, -1])[, c(T, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_qn_se = apply((", i, "_result[, -1])[, c(F, T, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_mn_se = apply((", i, "_result[, -1])[, c(F, F, T, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_vsn_se = apply((", i, "_result[, -1])[, c(F, F, F, T)], 1, sd) / 10")))
      
      eval(parse(text = paste0("error = c(error, ", i, "_result_non, ", i, "_result_qn, ", i, "_result_mn, ", i, "_result_vsn)")))
      eval(parse(text = paste0("se = c(se, ", i, "_result_non_se, ", i, "_result_qn_se, ", i, "_result_mn_se, ", i, "_result_vsn_se)")))
    } else {
      norm <- c("non", "TC", "UQ", "med", "TMM", "TMM.frozen", "DESeq", "QN", "QN.frozen")
      filename <- paste(file, n, i, "batch_balance", batch_balance, "result.csv", sep = "_")
      eval(parse(text = paste0(i, "_result = read.csv('", filename, "', sep = '\t')")))
      eval(parse(text = paste0(i, "_result_non = rowMeans((", i, "_result[, -1])[, c(T, F, F, F, F, F, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_TC = rowMeans((", i, "_result[, -1])[, c(F, T, F, F, F, F, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_UQ = rowMeans((", i, "_result[, -1])[, c(F, F, T, F, F, F, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_med = rowMeans((", i, "_result[, -1])[, c(F, F, F, T, F, F, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_TMM = rowMeans((", i, "_result[, -1])[, c(F, F, F, F, T, F, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_TMM.frozen = rowMeans((", i, "_result[, -1])[, c(F, F, F, F, F, T, F, F, F)])")))
      eval(parse(text = paste0(i, "_result_DESeq = rowMeans((", i, "_result[, -1])[, c(F, F, F, F, F, F, T, F, F)])")))
      eval(parse(text = paste0(i, "_result_QN = rowMeans((", i, "_result[, -1])[, c(F, F, F, F, F, F, F, T, F)])")))
      eval(parse(text = paste0(i, "_result_QN.frozen = rowMeans((", i, "_result[, -1])[, c(F, F, F, F, F, F, F, F, T)])")))
      
      eval(parse(text = paste0(i, "_result_non_se = apply((", i, "_result[, -1])[, c(T, F, F, F, F, F, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_TC_se = apply((", i, "_result[, -1])[, c(F, T, F, F, F, F, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_UQ_se = apply((", i, "_result[, -1])[, c(F, F, T, F, F, F, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_med_se = apply((", i, "_result[, -1])[, c(F, F, F, T, F, F, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_TMM_se = apply((", i, "_result[, -1])[, c(F, F, F, F, T, F, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_TMM.frozen_se = apply((", i, "_result[, -1])[, c(F, F, F, F, F, T, F, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_DESeq_se = apply((", i, "_result[, -1])[, c(F, F, F, F, F, F, T, F, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_QN_se = apply((", i, "_result[, -1])[, c(F, F, F, F, F, F, F, T, F)], 1, sd) / 10")))
      eval(parse(text = paste0(i, "_result_QN.frozen_se = apply((", i, "_result[, -1])[, c(F, F, F, F, F, F, F, F, T)], 1, sd) / 10")))
      
      eval(parse(text = paste0("error = c(error, ", i, "_result_non, ", i, "_result_TC, ", i, "_result_UQ, ", i, "_result_med, ", i, "_result_TMM, ", i, "_result_TMM.frozen, ", i, "_result_DESeq, ", i, "_result_QN, ", i, "_result_QN.frozen)")))
      eval(parse(text = paste0("se = c(se, ", i, "_result_non_se, ", i, "_result_TC_se, ", i, "_result_UQ_se, ", i, "_result_med_se, ", i, "_result_TMM_se, ", i, "_result_TMM.frozen_se, ", i, "_result_DESeq_se, ", i, "_result_QN_se, ", i, "_result_QN.frozen_se)")))
    }
  }
  
  strategy <- factor(split_combination)
  
  result_result <- data.frame(
    error = error,
    se = se,
    ml = rep(ml, each = length(norm) * length(split_combination)),
    norm = rep(norm, each = length(split_combination), length(ml)),
    x = rep(seq_along(split_combination), length(norm) * length(ml))
  )
  
  p <- ggplot(result_result, aes(x = x, y = error, color = ml, linetype = norm)) +
    geom_line() +
    scale_x_continuous(breaks = seq_along(split_combination), labels = split_combination) +
    geom_errorbar(aes(ymin = error - se * 1.96, ymax = error + se * 1.96), width = 0.1) +
    labs(x = "", title = paste0(toupper(substring(file, 1, 1)), substring(file, 2)))
  
  print(p)
}

# =========================
# Result plot: triple splitting
# deprecate this function.
# =========================
result_plot_all <- function(
    n,
    ml,
    file_contrast1,
    file_contrast2,
    batch_balance,
    data_type_rnaseq,
    split_comination
) {
  error <- vector("numeric")
  se <- vector("numeric")
  
  for (i in ml) {
    filename <- paste(file_contrast2, n, i, "batch_balance", batch_balance, "result.csv", sep = "_")
    eval(parse(text = paste0(i, "_result = read.csv('", filename, "', sep = '\t')")))
    eval(parse(text = paste0(i, "_result_non = rowMeans((", i, "_result[, -1])[, c(T, F, F, F)])")))
    eval(parse(text = paste0(i, "_result_qn = rowMeans((", i, "_result[, -1])[, c(F, T, F, F)])")))
    eval(parse(text = paste0(i, "_result_mn = rowMeans((", i, "_result[, -1])[, c(F, F, T, F)])")))
    eval(parse(text = paste0(i, "_result_vsn = rowMeans((", i, "_result[, -1])[, c(F, F, F, T)])")))
    
    eval(parse(text = paste0(i, "_result_non_se = apply((", i, "_result[, -1])[, c(T, F, F, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_qn_se = apply((", i, "_result[, -1])[, c(F, T, F, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_mn_se = apply((", i, "_result[, -1])[, c(F, F, T, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_vsn_se = apply((", i, "_result[, -1])[, c(F, F, F, T)], 1, sd) / 10")))
    
    eval(parse(text = paste0("error = c(error, ", i, "_result_non, ", i, "_result_qn, ", i, "_result_mn, ", i, "_result_vsn)")))
    eval(parse(text = paste0("se = c(se, ", i, "_result_non_se, ", i, "_result_qn_se, ", i, "_result_mn_se, ", i, "_result_vsn_se)")))
  }
  
  result_1 <- data.frame(
    error = error,
    se = se,
    ml = rep(ml, each = 4 * length(1)),
    norm = rep(c("None", "QN", "MN", "VSN"), each = length(1), length(ml)),
    x = rep(length(split_comination) + 1, 4 * length(ml))
  )
  
  error <- vector("numeric")
  se <- vector("numeric")
  
  for (i in ml) {
    filename <- paste(file_contrast1, n, i, "batch_balance", batch_balance, "result.csv", sep = "_")
    eval(parse(text = paste0(i, "_result = read.csv('", filename, "', sep = '\t')")))
    eval(parse(text = paste0(i, "_result_non = rowMeans((", i, "_result[, -1])[, c(T, F, F, F)])")))
    eval(parse(text = paste0(i, "_result_qn = rowMeans((", i, "_result[, -1])[, c(F, T, F, F)])")))
    eval(parse(text = paste0(i, "_result_mn = rowMeans((", i, "_result[, -1])[, c(F, F, T, F)])")))
    eval(parse(text = paste0(i, "_result_vsn = rowMeans((", i, "_result[, -1])[, c(F, F, F, T)])")))
    
    eval(parse(text = paste0(i, "_result_non_se = apply((", i, "_result[, -1])[, c(T, F, F, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_qn_se = apply((", i, "_result[, -1])[, c(F, T, F, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_mn_se = apply((", i, "_result[, -1])[, c(F, F, T, F)], 1, sd) / 10")))
    eval(parse(text = paste0(i, "_result_vsn_se = apply((", i, "_result[, -1])[, c(F, F, F, T)], 1, sd) / 10")))
    
    eval(parse(text = paste0("error = c(error, ", i, "_result_non, ", i, "_result_qn, ", i, "_result_mn, ", i, "_result_vsn)")))
    eval(parse(text = paste0("se = c(se, ", i, "_result_non_se, ", i, "_result_qn_se, ", i, "_result_mn_se, ", i, "_result_vsn_se)")))
  }
  
  result_2 <- data.frame(
    error = error,
    se = se,
    ml = rep(ml, each = 4 * length(split_comination)),
    norm = rep(c("None", "QN", "MN", "VSN"), each = length(split_comination), length(ml)),
    x = rep(seq_along(split_comination), 4 * length(ml))
  )
  
  result <- rbind(result_1, result_2)
  
  p <- ggplot(result, aes(x = x, y = error, color = ml, linetype = norm)) +
    geom_line() +
    scale_x_continuous(
      breaks = 1:(length(split_comination) + 1),
      labels = c(
        paste("contrast 1 (", split_comination, ")", sep = ""),
        "contrast 2 (100)"
      )
    ) +
    geom_errorbar(aes(ymin = error - se * 1.96, ymax = error + se * 1.96), width = 0.1) +
    labs(x = "", title = paste(file_contrast1, "VS", file_contrast2))
  
  print(p)
}

result_plot2 <- function(
    n,
    ml,
    file,
    batch_balance = FALSE,
    data_type_rnaseq = FALSE,
    split_combination
) {
  
  # 预设 normalization 顺序
  norm_order <- if (!data_type_rnaseq) {
    c("non", "qn", "mn", "vsn")
  } else {
    c("non", "TC", "UQ", "med", "TMM", "TMM.frozen", "DESeq", "QN", "QN.frozen")
  }
  
  # 读取并汇总单个文件
  summarize_result_file <- function(filename, norm_order) {
    dat <- read.csv(filename, sep = "\t", check.names = FALSE)
    
    # 第一列通常是 split 比例
    split_vals <- dat[, 1]
    mat <- dat[, -1, drop = FALSE]
    
    res_list <- list()
    
    for (norm_name in norm_order) {
      # 按列名前缀匹配，比如 "^non\\."、"^qn\\."
      norm_pattern <- gsub("([.|()\\^{}+$*?]|\\[|\\]|\\\\)", "\\\\\\1", norm_name)
      cols <- grep(paste0("^", norm_pattern, "\\.?"), colnames(mat), ignore.case = TRUE)
      
      if (length(cols) == 0) next
      
      submat <- mat[, cols, drop = FALSE]
      
      # 去掉全 NA 列
      keep_col <- colSums(!is.na(submat)) > 0
      submat <- submat[, keep_col, drop = FALSE]
      
      if (ncol(submat) == 0) next
      
      res_list[[norm_name]] <- data.frame(
        split = split_vals,
        error = rowMeans(submat, na.rm = TRUE),
        se = apply(submat, 1, sd, na.rm = TRUE) / sqrt(ncol(submat)),
        norm = norm_name,
        stringsAsFactors = FALSE
      )
    }
    
    if (length(res_list) == 0) return(NULL)
    
    out <- do.call(rbind, res_list)
    rownames(out) <- NULL
    out
  }
  
  all_res <- list()
  
  for (method in ml) {
    filename <- paste(file, n, method, "batch_balance", batch_balance, "result.csv", sep = "_")
    
    if (!file.exists(filename)) {
      warning("File not found: ", filename)
      next
    }
    
    tmp <- summarize_result_file(filename, norm_order)
    
    if (is.null(tmp) || nrow(tmp) == 0) next
    
    tmp$ml <- method
    all_res[[method]] <- tmp
  }
  
  if (length(all_res) == 0) {
    stop("No valid result data found.")
  }
  
  result_result <- do.call(rbind, all_res)
  rownames(result_result) <- NULL
  
  # 如果 split_combination 提供了，就用它控制横轴顺序
  if (!missing(split_combination) && !is.null(split_combination)) {
    result_result$split <- factor(result_result$split, levels = split_combination)
  } else {
    result_result$split <- factor(result_result$split, levels = unique(result_result$split))
  }
  
  result_result$norm <- factor(result_result$norm, levels = norm_order)
  result_result$ml <- factor(result_result$ml, levels = ml)
  
  p <- ggplot(
    result_result,
    aes(
      x = split,
      y = error,
      color = ml,
      linetype = norm,
      group = interaction(ml, norm)
    )
  ) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_errorbar(
      aes(ymin = error - 1.96 * se, ymax = error + 1.96 * se),
      width = 0.1
    ) +
    labs(
      x = "",
      y = "Error",
      title = paste0(toupper(substr(file, 1, 1)), substr(file, 2, nchar(file)))
    ) +
    theme_bw()
  
  print(p)
}
