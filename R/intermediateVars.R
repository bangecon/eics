##' Checking for collider bias using random forests.
##'
##' This function checks whether control variables that a process (e.g.
##' researcher "intuition", a model-selection algorithm) selected might create
##' collider bias with respect to the current treatment.
##'
##' @param data is the universe of all variables considered, which may include
##' both controls and excluded variables (e.g. variables eliminated by a model
##' selection algorithm).
##' @param treatment is the name of the treatment variable contained in
##' \code{data}.
##' @param controls is a list() of names of selected controls contained in
##' \code{data}.
##' @param ntree is the numerical number of trees to use in the random forest
##' models.
##' @param zthresh is the z-score used as a threshold for flagging a variable
##' for suspicion as an alternative treatment, feedback, or intermediate bias.
##' @param out is a logical indicating whether to create a file named
##' "alternateTreatments.csv" in the current working directory that contains
##' the variables flagged as alternate treatments.
##' @param parallel is a logical indicating whether to use parallel processing
##' @return \code{intermediateVars} returns a list that contains the following
##' components:
##' \item{treatment}{name of the treatment variable.}
##' \item{controls}{the list of control variables.}
##' \item{intermediateVariables}{the list of possible intermediate variables.}
##' \item{zThresh}{the z-score threshold for suspicion.}
##' \item{zMatrix}{the scaled variable importance matrix between the treatment,
##' controls, and outcome.}
##' @export
intermediateVars <- function(data,
                             treatment,
                             controls,
                             ntree = 100,
                             zthresh = 1,
                             out = FALSE,
                             parallel = FALSE) {
  X <- data[, which(names(data) %in% c(treatment, controls))]
  X <- X[, -which(names(X) %in% treatment)]
  TX <- cbind(data[, which(names(data) %in% treatment)], X)
  names(TX)[1] <- treatment
  # Build predictor importance z-matrix -------------------
  if (parallel == TRUE) {
    mdl <-
      foreach(
        nodetrees = rep(ceiling(ntree / 5), 5),
        .combine = randomForest::combine,
        .multicombine = TRUE,
        .packages = 'randomForest'
      ) %dorng% {
        .libPaths(libs)
        randomForest(TX[,-1], TX[, 1], ntree = nodetrees, importance = TRUE)
      }
    zmat <- data.frame(c(0, importance(
      mdl, type = 1, scale = TRUE
    )))
    rownames(zmat) <- colnames(TX)
    for (i in 2:(ncol(TX) - 1)) {
      mdl <-
        foreach(
          nodetrees = rep(ceiling(ntree / 5), 5),
          .combine = randomForest::combine,
          .multicombine = TRUE,
          .packages = 'randomForest'
        ) %dorng% {
          .libPaths(libs)
          randomForest(TX[,-i],
                       TX[, i],
                       ntree = nodetrees,
                       importance = TRUE)
        }
      predImpZ <- importance(mdl, type = 1, scale = TRUE)
      zmat <- cbind(zmat,
                    c(predImpZ[1:(i - 1)], 0, predImpZ[i:length(predImpZ)]))
    }
    mdl <-
      foreach(
        nodetrees = rep(ceiling(ntree / 5), 5),
        .combine = randomForest::combine,
        .multicombine = TRUE,
        .packages = 'randomForest'
      ) %dorng% {
        .libPaths(libs)
        randomForest(TX[, -ncol(TX)], TX[, ncol(TX)], ntree = nodetrees, importance = TRUE)
      }
    zmat <- cbind(zmat, c(importance(
      mdl, type = 1, scale = TRUE
    ), 0))
    unregister()
  } else {
    mdl <-
      randomForest(TX[,-1], TX[, 1], ntree = ntree, importance = TRUE)
    zmat <- data.frame(c(0, importance(
      mdl, type = 1, scale = TRUE
    )))
    rownames(zmat) <- colnames(TX)
    for (i in 2:(ncol(TX) - 1)) {
      mdl <-
        randomForest(TX[,-i], TX[, i], ntree = ntree, importance = TRUE)
      predImpZ <- importance(mdl, type = 1, scale = TRUE)
      zmat <- cbind(zmat,
                    c(predImpZ[1:(i - 1)], 0, predImpZ[i:length(predImpZ)]))
    }
    mdl <-
      randomForest(TX[, -ncol(TX)], TX[, ncol(TX)], ntree = ntree, importance = TRUE)
    zmat <- cbind(zmat, c(importance(
      mdl, type = 1, scale = TRUE
    ), 0))
  }
  zmat <- (zmat - mean(unlist(zmat))) / sd(unlist(zmat))
  colnames(zmat) <- colnames(TX)
  # Check for downstream intermediate variables -
  TpredictX <-
    data.frame(t(as.matrix(zmat[1, 2:ncol(zmat)])))
  colnames(TpredictX) <- 'zscore'
  rownames(TpredictX) <- colnames(zmat)[2:nrow(zmat)]
  intVars <- rownames(subset(TpredictX, zscore > zthresh))
  if (out == TRUE) {
    write.csv(intVars,
              paste0(treatment, 'intermediateVariables.csv'))
  }
  out <- list(
    treatment = treatment,
    controls = controls,
    intermediateVariables = intVars,
    zMatrix = zmat,
    zThresh = ztrhesh
  )
  out
}
