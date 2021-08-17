##' Checking for alternate treatments using random forests.
##'
##' This function checks whether other variables that a previous process (e.g.
##' researcher "intuition", a model-selection algorithm) excluded might in
##' fact be a more useful treatment. In essence, it checks whether the current
##' treatment acted as a "collider" between eliminated variables and the
##' outcome (and hence may have influenced its elimination).
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
##' @return \code{altTreat} returns a list that contains the following components:
##' \item{treatment}{name of the treatment variable.}
##' \item{controls}{the list of control variables.}
##' \item{excludedVariables}{the list of variables in data not included as
##' controls.}
##' \item{alternateTreatments}{the list of possible alternate treatments.}
##' \item{impZ}{the scaled variable importance matrix between the excluded
##' variables and the treatment,}
##' \item{zThresh}{the z-score threshold for suspicion.}
##' @export
altTreat <- function(data,
                     treatment,
                     controls,
                     ntree = 100,
                     zthresh = 1,
                     out = FALSE,
                     parallel = TRUE) {
  X <- data[, which(names(data) %in% c(treatment, controls))]
  X <- X[, -which(names(X) %in% treatment)]
  otherX <-
    data[, -which(names(data) %in% c(treatment, controls))]
  TX <- cbind(data[, which(names(data) %in% treatment)], X)
  names(TX)[1] <- treatment
  # Check  any "other predictors" might be a better treatment ------
  if (parallel == TRUE) {
    unregister <- function() {
      env <- foreach:::.foreachGlobals
      rm(list = ls(name = env), pos = env)
    }
    no_cores <- ceiling(detectCores() - 1)
    libs <- .libPaths()
    cl <- makePSOCKcluster(no_cores, outfile = '')
    registerDoParallel(cl, no_cores)
    clusterEvalQ(cl)
    clusterExport(cl, c('data', 'TX', 'otherX', 'ntree'))
    mdlone <-
      foreach(
        nodetrees = rep(ceiling(ntree / 5), 5),
        .combine = randomForest::combine,
        .multicombine = TRUE,
        .packages = 'randomForest'
      ) %dorng% {
        .libPaths(libs)
        randomForest(otherX,
                     TX[, 1],
                     ntree = nodetrees,
                     importance = T)
      }
    unregister()
  } else {
    mdlone <-
      randomForest(otherX, TX[, 1], ntree = ntree, importance = T)
  }
  impZ <-
    data.frame(scale(importance(
      mdlone, type = 1, scale = T
    )))
  names(impZ) <- 'impZ'
  altTreat <- rownames(subset(impZ, impZ > zthresh))
  if (out == T) {
    write.csv(altTreat,
              paste0(treatment, 'alternateTreatments.csv'))
  }
  out <- list(
    treatment = treatment,
    controls = controls,
    excludedVariables = names(otherX),
    alternateTreatments = altTreat,
    impZ = impZ,
    zThresh = zthresh
  )
  out
}
