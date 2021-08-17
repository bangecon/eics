##' Empirically-Informed Covariate Selection (EICS) using random forests.
##'
##' This function checks all three possible violations of the backdoor
##' criterion by:
##' 1. Checking whether an eliminated variable might be superior to the
##' current treatment (separate function: altTreat);
##' 2. Checking whether one or more control variables selected for the model
##' create a feedback loop with the current treatment;
##' 3. Checking whether one or more control variables selected create intermediate
##' bias with the current treatment.
##'
##' @param data is the universe of all variables considered, which may include
##' both controls and excluded variables (e.g. variables eliminated by a model
##' selection algorithm).
##' @param outcome is the name of the outcome variable contained in \code{data}.
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
##' @return \code{eics} returns a list that contains the following components:
##' \item{treatment}{name of the treatment variable.}
##' \item{controls}{the list of control variables.}
##' \item{excludedVariables}{the list of variables in data not included as
##' controls.}
##' \item{alternateTreatments}{the list of possible alternate treatments.}
##' \item{feedbackVariables}{the list of possible feedback variables.}
##' \item{intermediateVariables}{the list of possible intermediate variables.}
##' \item{alternateOutcomes}{the list of possible alternate outcomes.}
##' \item{zThresh}{the z-score threshold for suspicion.}
##' \item{impZ}{the scaled variable importance matrix between the excluded
##' variables and the treatment,}
##' \item{zMatrix}{the scaled variable importance matrix between the treatment,
##' controls, and outcome.}
##' @export
eics <- function(data,
                 outcome,
                 treatment,
                 controls,
                 ntree = 100,
                 zthresh = 1,
                 out = FALSE,
                 parallel = FALSE) {
    Y <- data[, which(names(data) %in% outcome)]
    X <- data[, which(names(data) %in% c(treatment, controls))]
    X <- X[,-which(names(X) %in% treatment)]
    TX <- cbind(data[, which(names(data) %in% treatment)], X)
    otherX <-
        data[,-which(names(data) %in% c(treatment, controls))]
    names(TX)[1] <- treatment
    # Check for Alternate Treatments ------------------------------------------
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
        clusterExport(cl, c('data', 'Y', 'TX', 'otherX', 'ntree'))
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
    if (out == TRUE) {
        write.csv(altTreat,
                  paste0(treatment, 'alternateTreatments.csv'))
    }
    # Build Importance Vector from Y to X -------------------------------------
    impY <- NULL
    if (parallel == TRUE) {
        for (i in 2:ncol(TX)) {
            mdltwo <-
                foreach(
                    nodetrees = rep(ceiling(ntree / 5), 5),
                    .combine = randomForest::combine,
                    .multicombine = TRUE,
                    .packages = 'randomForest'
                ) %dorng% {
                    .libPaths(libs)
                    randomForest(cbind(Y, TX[, -i], otherX),
                                 TX[, i],
                                 ntree = nodetrees,
                                 importance = TRUE)
                }
            impY[i - 1] <-
                importance(mdltwo, type = 1, scale = TRUE)[1]
        }
    } else {
        mdltwo <-
            randomForest(cbind(Y, TX[, -i], otherX),
                         TX[, i],
                         ntree = ntree,
                         importance = TRUE)
        impY[i - 1] <- importance(mdltwo, type = 1, scale = TRUE)[1]
    }
    names(impY) <- colnames(TX)[-1]
    # Check for Alternate Outcomes --------------------------------------------
    altOutcomes <- rownames(subset(impY, impY > zthresh))
    if (out == TRUE) {
        write.csv(altOutcomes,
                  paste0(treatment, 'alternateOutcomes.csv'))
    }
    # Build 2-Way Predictor Importance Matrix ---------------------------------
    if (parallel == TRUE) {
        mdl <-
            foreach(
                nodetrees = rep(ceiling(ntree / 5), 5),
                .combine = randomForest::combine,
                .multicombine = TRUE,
                .packages = 'randomForest'
            ) %dorng% {
                .libPaths(libs)
                randomForest(TX[, -1], TX[, 1], ntree = nodetrees, importance = TRUE)
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
                    randomForest(TX[, -i],
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
                randomForest(TX[,-ncol(TX)], TX[, ncol(TX)], ntree = nodetrees, importance = TRUE)
            }
        zmat <- cbind(zmat, c(importance(
            mdl, type = 1, scale = TRUE
        ), 0))
    } else {
        mdl <-
            randomForest(TX[, -1], TX[, 1], ntree = ntree, importance = TRUE)
        zmat <- data.frame(c(0, importance(
            mdl, type = 1, scale = TRUE
        )))
        rownames(zmat) <- colnames(TX)
        for (i in 2:(ncol(TX) - 1)) {
            mdl <-
                randomForest(TX[, -i], TX[, i], ntree = ntree, importance = TRUE)
            predImpZ <- importance(mdl, type = 1, scale = TRUE)
            zmat <- cbind(zmat,
                          c(predImpZ[1:(i - 1)], 0, predImpZ[i:length(predImpZ)]))
        }
        mdl <-
            randomForest(TX[,-ncol(TX)], TX[, ncol(TX)], ntree = ntree, importance = TRUE)
        zmat <- cbind(zmat, c(importance(
            mdl, type = 1, scale = TRUE
        ), 0))
    }
    zmat <- (zmat - mean(unlist(zmat))) / sd(unlist(zmat))
    colnames(zmat) <- colnames(TX)
    # Check for Feedback ------------------------------------------------------
    fbCheck <- logical(length = nrow(zmat))
    names(fbCheck) <- rownames(zmat)
    for (i in 2:ncol(TX)) {
        fbCheck[i] <- (zmat[1, i] > zthresh & zmat[i, 1] > zthresh)
    }
    fbCheck <- data.frame(fbCheck)
    fbVars <- rownames(subset(fbCheck, fbCheck == TRUE))
    if (out == TRUE) {
        write.csv(fbVars, paste0(treatment, 'feedbackVariables.csv'))
    }
    # Check for Intermediate Variables ----------------------------------------
    TpredictX <-
        data.frame(t(as.matrix(zmat[1, 2:ncol(zmat)])))
    colnames(TpredictX) <- 'zscore'
    rownames(TpredictX) <- colnames(zmat)[2:nrow(zmat)]
    intVars <- rownames(subset(TpredictX, zscore > zthresh))
    if (out == TRUE) {
        write.csv(intVars,
                  paste0(treatment, 'intermediateVariables.csv'))
    }
    # Check for Alternate Outcomes --------------------------------------------
    names(impY) <- colnames(TX)[-1]
    altOutcomes <- rownames(subset(impY, impY > zthresh))
    if (out == TRUE) {
        write.csv(altOutcomes,
                  paste0(treatment, 'alternateOutcomes.csv'))
    }
    unregister()
    out <- list(
        treatment = treatment,
        controls = controls,
        excludedVariables = names(otherX),
        alternateTreatments = altTreat,
        intermediateVariables = intVars,
        feedbackVariables = fbVars,
        alternateOutcomes = altOutcomes,
        zThresh = ztrhesh,
        impZ = impZ,
        zMatrix = zmat
    )
    class(out) <- c('eics', class(out))
    out
}
