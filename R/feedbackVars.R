##' Checking for feedback loops using random forests.
##'
##' This function checks whether control variables that a process (e.g.
##' researcher "intuition", a model-selection algorithm) selected might create
##' feedback loops with the current treatment.
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
##' @return \code{feedbackVars} returns a list that contains the following
##' components:
##' \item{treatment}{name of the treatment variable.}
##' \item{controls}{the list of control variables.}
##' \item{feedbackVariables}{the list of possible feedback variables.}
##' \item{zThresh}{the z-score threshold for suspicion.}
##' \item{zMatrix}{the scaled variable importance matrix between the treatment,
##' controls, and outcome.}
##' @export
feedbackVars <- function(data,
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
  # Check for "feedback" between T and X's --------------------
  fbCheck <- logical(length = nrow(zmat))
  names(fbCheck) <- rownames(zmat)
  for (i in 2:ncol(TX)) {
    fbCheck[i] <- (zmat[1, i] > zthresh & zmat[i, 1] > zthresh)
  }
  fbCheck <- data.frame(fbCheck)
  names(fbCheck) <- 'zVzthresh'
  fbVars <- rownames(subset(fbCheck, zVzthresh == TRUE))
  if (out == TRUE) {
    write.csv(fbVars, paste0(treatment, 'feedbackVariables.csv'))
  }
  out <- list(
    treatment = treatment,
    controls = controls,
    feedbackVariables = fbVars,
    zThresh = zthresh,
    impZ = impZ,
    zMatrix = zmat
  )
  out
}
