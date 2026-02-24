suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(CoxBoost)
  library(randomForestSRC)
  library(gbm)
  library(plsRcox)
  library(MASS)
  library(survivalsvm)
  library(superpc)
  library(dplyr)
  library(yaml)
})

# ------------------------------------------------------------
# 101-model framework (10 algorithms × 10 feature selectors = 100) + 1 baseline
# Algorithms:
# CoxBoost, Lasso, Ridge, RSF, Enet, GBM, plsRcox, stepwise Cox, survival-SVM, SuperPC
#
# Inputs (must exist):
#  - results/survival/features_matrix.csv  (samples x genes/features)
#  - results/survival/survival_table.csv   (sample,time,event) with no missing values
#
# Optional external validation cohorts:
#  - data/external/<COHORT_NAME>_features.csv  (samples x features; same feature names)
#  - data/external/<COHORT_NAME>_survival.csv  (sample,time,event)
# The script will evaluate C-index per cohort and rank models by mean C-index.
# ------------------------------------------------------------

cfg <- yaml::read_yaml("config/config.yaml")
set.seed(cfg$project$seed)

out_dir <- file.path(cfg$project$out_dir, "survival_101")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

X <- read.csv("results/survival/features_matrix.csv", row.names = 1, check.names = FALSE)
y <- read.csv("results/survival/survival_table.csv", check.names = FALSE)

stopifnot(all(c("sample","time","event") %in% colnames(y)))
y <- y %>% mutate(sample = as.character(sample))
rownames(y) <- y$sample
y <- y[rownames(X), , drop = FALSE]

if (any(is.na(y$time)) || any(is.na(y$event))) stop("survival_table.csv has NA in time/event. Fill it first.")
y$time <- as.numeric(y$time)
y$event <- as.integer(y$event)

# train/test split (internal)
test_size <- ifelse(is.null(cfg$survival$test_size), 0.3, cfg$survival$test_size)
n <- nrow(X)
test_idx <- sample(seq_len(n), size = floor(test_size * n))
train_idx <- setdiff(seq_len(n), test_idx)

X_train <- as.matrix(X[train_idx, , drop = FALSE])
X_test  <- as.matrix(X[test_idx,  , drop = FALSE])
y_train <- y[train_idx, , drop = FALSE]
y_test  <- y[test_idx,  , drop = FALSE]

surv_train <- Surv(y_train$time, y_train$event)
surv_test  <- Surv(y_test$time,  y_test$event)

# helper: C-index
cindex <- function(time, event, risk) {
  survConcordance(Surv(time, event) ~ risk)$concordance
}

# Univariate Cox ranking (fallback selector)
uni_rank <- function(Xm, ydf) {
  pvals <- apply(Xm, 2, function(g) {
    fit <- coxph(Surv(ydf$time, ydf$event) ~ g)
    summary(fit)$coefficients[,"Pr(>|z|)"][1]
  })
  names(sort(pvals, decreasing = FALSE))
}

# Feature selection wrappers (return ordered feature names)
select_features <- function(method, Xm, ydf, top_k = 30) {
  feats <- colnames(Xm)
  ranked <- NULL

  safe_top <- function(vec) {
    vec <- vec[vec %in% feats]
    if (length(vec) == 0) vec <- uni_rank(Xm, ydf)
    head(unique(vec), top_k)
  }

  if (method == "CoxBoost") {
    cb <- CoxBoost(time = ydf$time, status = ydf$event, x = Xm, stepno = 100)
    co <- coef(cb, at.step = cb$stepno)
    ranked <- names(sort(abs(co), decreasing = TRUE))
    ranked <- ranked[abs(co[ranked]) > 0]
    return(safe_top(ranked))
  }

  if (method %in% c("Lasso","Ridge","Enet")) {
    alpha <- ifelse(method=="Lasso", 1, ifelse(method=="Ridge", 0, 0.5))
    fit <- cv.glmnet(x = Xm, y = surv_train, family = "cox", alpha = alpha, nfolds = 10)
    co <- as.matrix(coef(fit, s = "lambda.min"))
    co <- co[,1]
    ranked <- names(sort(abs(co), decreasing = TRUE))
    ranked <- ranked[abs(co[ranked]) > 0]
    return(safe_top(ranked))
  }

  if (method == "RSF") {
    rf <- rfsrc(Surv(time, event) ~ ., data = data.frame(time=ydf$time, event=ydf$event, Xm),
                ntree = 500, importance = TRUE)
    imp <- rf$importance
    ranked <- names(sort(imp, decreasing = TRUE))
    return(safe_top(ranked))
  }

  if (method == "GBM") {
    df <- data.frame(time=ydf$time, event=ydf$event, Xm)
    g <- gbm(Surv(time, event) ~ ., data=df, distribution="coxph", n.trees=2000,
             interaction.depth=3, shrinkage=0.01, n.minobsinnode=10, bag.fraction=0.7, verbose=FALSE)
    s <- summary(g, plotit=FALSE)
    ranked <- as.character(s$var)
    return(safe_top(ranked))
  }

  if (method == "plsRcox") {
    df <- data.frame(time=ydf$time, event=ydf$event, Xm)
    pls <- plsRcox::plsRcox(Surv(time, event) ~ ., data=df, nt=5, verbose=FALSE)
    # rank by absolute coefficients
    co <- pls$coefficients
    if (is.null(co)) return(safe_top(NULL))
    co <- as.numeric(co)
    names(co) <- colnames(Xm)
    ranked <- names(sort(abs(co), decreasing = TRUE))
    return(safe_top(ranked))
  }

  if (method == "stepCox") {
    df <- data.frame(time=ydf$time, event=ydf$event, Xm)
    # start with null model, then stepwise forward; to keep it stable, cap max vars
    null <- coxph(Surv(time, event) ~ 1, data=df)
    full <- coxph(Surv(time, event) ~ ., data=df)
    st <- suppressWarnings(stepAIC(null, scope = list(lower=null, upper=full), direction="both", trace=FALSE))
    vars <- names(st$coefficients)
    return(safe_top(vars))
  }

  if (method == "survivalSVM") {
    # no stable variable importance; fallback to univariate ranking
    return(safe_top(uni_rank(Xm, ydf)))
  }

  if (method == "SuperPC") {
    dat <- list(
      x = t(Xm),
      y = ydf$time,
      censoring.status = 1-ydf$event, # superpc uses 1=censored
      featurenames = colnames(Xm)
    )
    spc <- superpc.train(dat, type = "survival")
    # choose a reasonable threshold by cross-val and take top features from scores
    cv <- superpc.cv(spc, dat)
    # select threshold giving max score
    best_i <- which.max(cv$scor)
    thr <- cv$thresholds[best_i]
    idx <- which(abs(spc$feature.scores) >= thr)
    ranked <- dat$featurenames[order(abs(spc$feature.scores), decreasing=TRUE)]
    ranked <- ranked[ranked %in% dat$featurenames[idx]]
    return(safe_top(ranked))
  }

  # default fallback
  safe_top(uni_rank(Xm, ydf))
}

# Model training/prediction wrappers: return risk for new data
fit_predict <- function(alg, Xm_train, ydf_train, Xm_new) {
  surv_obj <- Surv(ydf_train$time, ydf_train$event)

  if (alg == "CoxBoost") {
    cb <- CoxBoost(time=ydf_train$time, status=ydf_train$event, x=Xm_train, stepno=200)
    lp <- predict(cb, newdata=Xm_new, type="lp", at.step=cb$stepno)
    return(as.numeric(lp))
  }

  if (alg %in% c("Lasso","Ridge","Enet")) {
    alpha <- ifelse(alg=="Lasso", 1, ifelse(alg=="Ridge", 0, 0.5))
    fit <- cv.glmnet(x = Xm_train, y = surv_obj, family="cox", alpha=alpha, nfolds=10)
    lp <- predict(fit, newx = Xm_new, s="lambda.min", type="link")
    return(as.numeric(lp))
  }

  if (alg == "RSF") {
    df <- data.frame(time=ydf_train$time, event=ydf_train$event, Xm_train)
    rf <- rfsrc(Surv(time, event) ~ ., data=df, ntree=1000)
    pred <- predict(rf, newdata=data.frame(Xm_new))
    # Use -mortality as risk proxy (higher risk => worse)
    risk <- pred$predicted
    return(as.numeric(risk))
  }

  if (alg == "GBM") {
    df <- data.frame(time=ydf_train$time, event=ydf_train$event, Xm_train)
    g <- gbm(Surv(time, event) ~ ., data=df, distribution="coxph", n.trees=3000,
             interaction.depth=3, shrinkage=0.01, n.minobsinnode=10, bag.fraction=0.7, verbose=FALSE)
    best <- gbm.perf(g, method="OOB", plot.it=FALSE)
    lp <- predict(g, newdata=data.frame(Xm_new), n.trees=best, type="link")
    return(as.numeric(lp))
  }

  if (alg == "plsRcox") {
    df <- data.frame(time=ydf_train$time, event=ydf_train$event, Xm_train)
    pls <- plsRcox::plsRcox(Surv(time, event) ~ ., data=df, nt=5, verbose=FALSE)
    lp <- predict(pls, newdata=data.frame(Xm_new), type="lp")
    return(as.numeric(lp))
  }

  if (alg == "stepCox") {
    df <- data.frame(time=ydf_train$time, event=ydf_train$event, Xm_train)
    null <- coxph(Surv(time, event) ~ 1, data=df)
    full <- coxph(Surv(time, event) ~ ., data=df)
    st <- suppressWarnings(stepAIC(null, scope=list(lower=null, upper=full), direction="both", trace=FALSE))
    lp <- predict(st, newdata=data.frame(Xm_new), type="lp")
    return(as.numeric(lp))
  }

  if (alg == "survivalSVM") {
    df <- data.frame(time=ydf_train$time, event=ydf_train$event, Xm_train)
    svm <- survivalsvm(Surv(time, event) ~ ., data=df, type="regression",
                      kernel="rbf", gamma.mu=0.1, opt.meth="quadprog")
    # survivalsvm returns prognostic index via predict()
    pred <- predict(svm, newdata=data.frame(Xm_new))
    return(as.numeric(pred$predicted))
  }

  if (alg == "SuperPC") {
    dat <- list(
      x = t(Xm_train),
      y = ydf_train$time,
      censoring.status = 1-ydf_train$event,
      featurenames = colnames(Xm_train)
    )
    spc <- superpc.train(dat, type="survival")
    cv <- superpc.cv(spc, dat)
    best_i <- which.max(cv$scor)
    thr <- cv$thresholds[best_i]
    # SuperPC prediction needs test dat list
    dat_new <- list(
      x = t(Xm_new),
      y = rep(1, nrow(Xm_new)),
      censoring.status = rep(1, nrow(Xm_new)),
      featurenames = colnames(Xm_new)
    )
    pr <- superpc.predict(spc, dat, dat_new, threshold=thr, n.components=1)
    # use v.pred as risk-like score
    return(as.numeric(pr$v.pred))
  }

  stop(paste("Unknown algorithm:", alg))
}

selectors <- c("CoxBoost","Lasso","Ridge","RSF","Enet","GBM","plsRcox","stepCox","survivalSVM","SuperPC")
models    <- c("CoxBoost","Lasso","Ridge","RSF","Enet","GBM","plsRcox","stepCox","survivalSVM","SuperPC")

# Load optional external cohorts
load_external <- function() {
  ext_dir <- "data/external"
  if (!dir.exists(ext_dir)) return(list())
  feat_files <- list.files(ext_dir, pattern="_features\.csv$", full.names=TRUE)
  cohorts <- list()
  for (ff in feat_files) {
    name <- sub("_features\.csv$", "", basename(ff))
    sf <- file.path(ext_dir, paste0(name, "_survival.csv"))
    if (!file.exists(sf)) next
    Xc <- read.csv(ff, row.names=1, check.names=FALSE)
    yc <- read.csv(sf, check.names=FALSE)
    if (!all(c("sample","time","event") %in% colnames(yc))) next
    rownames(yc) <- as.character(yc$sample)
    common <- intersect(rownames(Xc), rownames(yc))
    Xc <- as.matrix(Xc[common, , drop=FALSE])
    yc <- yc[common, , drop=FALSE]
    cohorts[[name]] <- list(X=Xc, y=yc)
  }
  cohorts
}

external <- load_external()

results <- list()

# Baseline model (no selection) as the +1 to make 101
baseline_alg <- "Ridge"
baseline_feats <- colnames(X_train)
risk_test <- fit_predict(baseline_alg, X_train, y_train, X_test)
c_test <- cindex(y_test$time, y_test$event, risk_test)

ext_scores <- c()
if (length(external) > 0) {
  for (nm in names(external)) {
    Xc <- external[[nm]]$X
    yc <- external[[nm]]$y
    # align features
    commonf <- intersect(colnames(X_train), colnames(Xc))
    if (length(commonf) < 5) next
    rt <- fit_predict(baseline_alg, X_train[,commonf,drop=FALSE], y_train, Xc[,commonf,drop=FALSE])
    ext_scores[nm] <- cindex(yc$time, yc$event, rt)
  }
}
results[[1]] <- data.frame(
  selector="None",
  model=baseline_alg,
  n_features=length(baseline_feats),
  cindex_internal=c_test,
  cindex_external_mean=ifelse(length(ext_scores)>0, mean(ext_scores), NA_real_),
  cindex_mean=ifelse(length(ext_scores)>0, mean(c(c_test, ext_scores)), c_test)
)

# 100 combos: selector × model
counter <- 2
for (sel in selectors) {
  feats <- select_features(sel, X_train, y_train, top_k = 30)
  Xm_tr <- X_train[, feats, drop=FALSE]
  Xm_te <- X_test[,  feats, drop=FALSE]

  for (mdl in models) {
    # predict internal
    risk <- fit_predict(mdl, Xm_tr, y_train, Xm_te)
    ci_in <- cindex(y_test$time, y_test$event, risk)

    # external
    ext <- c()
    if (length(external) > 0) {
      for (nm in names(external)) {
        Xc <- external[[nm]]$X
        yc <- external[[nm]]$y
        commonf <- intersect(feats, colnames(Xc))
        if (length(commonf) < 5) next
        rxc <- fit_predict(mdl, Xm_tr[,commonf,drop=FALSE], y_train, Xc[,commonf,drop=FALSE])
        ext[nm] <- cindex(yc$time, yc$event, rxc)
      }
    }

    results[[counter]] <- data.frame(
      selector=sel,
      model=mdl,
      n_features=length(feats),
      cindex_internal=ci_in,
      cindex_external_mean=ifelse(length(ext)>0, mean(ext), NA_real_),
      cindex_mean=ifelse(length(ext)>0, mean(c(ci_in, ext)), ci_in)
    )
    counter <- counter + 1
  }
}

res_df <- bind_rows(results) %>% arrange(desc(cindex_mean))

write.csv(res_df, file.path(out_dir, "model_ranking_101.csv"), row.names = FALSE)

best <- res_df[1,]
writeLines(c(
  "Best model (by mean C-index):",
  paste0("Selector: ", best$selector),
  paste0("Model: ", best$model),
  paste0("n_features: ", best$n_features),
  paste0("C-index internal: ", round(best$cindex_internal, 4)),
  paste0("C-index external mean: ", round(best$cindex_external_mean, 4)),
  paste0("C-index mean: ", round(best$cindex_mean, 4))
), con = file.path(out_dir, "best_model.txt"))

message("Saved ranking to results/survival_101/model_ranking_101.csv")
message("Best model summary saved to results/survival_101/best_model.txt")
