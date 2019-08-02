#ROBUST Adaptive Eigenvalues Decomposition Discriminant Analysis :
# with this version we can jointly handle:
#noise labels (between existing groups in both test and training set)
#noise labels (with unobserved groups wrongly labelled as belonging to existing group
#in both test and training set)
#outliers in both test and training

#TRANDSUCTIVE APPROACH: robust semisupervised EM algorithm from robust parameters estimation
# and unobserved classes detection

#' Robust Apadtive Eigenvalues Decomposition Discriminant Analysis (transductive approach)
#'
#' @param X_train A numeric matrix of observations where rows correspond to observations and columns correspond to variables. The group membership of each observation is known - labeled data.
#' @param class_train A vector (if numeric it will be coerced to factor) with distinct entries representing a classification of the corresponding observations in X_train.
#' @param X_test A numeric matrix of observations where rows correspond to observations and columns correspond to variables. The group membership of each observation may not be known - unlabeled data.
#' @param model_names A character string indicating the desired models to be tested. With default NULL, all available models are tested. The models available for univariate and multivariate data are described in \code{\link{mclust::mclustModelNames}}
#' @param G A numeric vector indicating the number of expected classes in X_test. With default NULL, models with G = length(unique(cltrain)): (length(unique(cltrain))+2) are tested.
#' @param alpha_train The proportion of observations to be trimmed in X_test.
#' @param alpha_test The proportion of observations to be trimmed in X_train.
#' @param restr_factor The constant >= 1 that constrains the allowed differences among group scatters. Larger values imply larger differences of group scatters, a value of 1 specifies the strongest restriction (equivalent to set model_names=EII)
#' @param ctrl_EM A list of control parameters for the EM algorithm for estimation of model parameters; see also code{\link{control_EM}}
#' @param ctrl_restr A list of control parameters for the algorithm performing constrains on the differences among group scatters; see also code{\link{control_restr}}
#' @param ctrl_init A list of control parameters for the algorithm initialization; see also code{\link{control_init}}
#' @param verbose A logical argument controlling if a text progress bar is displayed during the fitting procedure. By default is TRUE if the session is interactive, and FALSE otherwise.
#' @param ... Additional internal arguments not to be provided by the user.
#' @return An object of class "raeddat" providing a list of output components for each model in model_names, with the Best model (according to BIC) first
#' @export
#'
#' @examples

raedda_t <-
  function (X_train,
            class_train,
            X_test,
            model_names = NULL,
            G = NULL,
            # the number of possible groups considered in the estimation
            alpha_train = 0,
            alpha_test = 0.05,
            restr_factor = NULL,
            ctrl_EM = control_EM(),
            ctrl_restr = control_restr(),
            ctrl_init = control_init(),
            verbose = interactive(),
            ...) {

    X_train <- data.matrix(X_train)
    X_test <- data.matrix(X_test)
    
    if(is.factor(class_train)) { # I do this for avoiding problems in naming objects afterwards
      class_train <- droplevels(class_train)
    }

    if (is.null(model_names)) {
      if (ncol(X_train) == 1) {
        model_names <- c("E", "V")
      } else {
        model_names <- mclust::mclust.options("emModelNames")
      }
    }
    if (is.null(G)) {
      G <-
        length(unique(class_train)):(length(unique(class_train)) + 2) #if I do not specify a specific G, I will consider the G current labels of the training, G+1 and G+2
    }
    if (verbose) {
      cat("fitting ...\n")
      flush.console()
      pbar <- txtProgressBar(min = 0,
                             max = length(model_names),
                             style = 3)
      on.exit(close(pbar))
      ipbar <- 0
    }
    RES <- list()
    bestBIC <- -Inf
    RES[["Best"]] <- list()
    for (model_name in model_names) {
      for (g in G) {
        RES[[model_name]][[as.character(g)]] <- list()
        RES[[model_name]][[as.character(g)]] <-
          raedda_t_model(
            X_train = X_train,
            class_train = class_train,
            X_test = X_test,
            model_name = model_name,
            G = g,
            alpha_train = alpha_train,
            alpha_test = alpha_test,
            restr_factor = restr_factor,
            ctrl_EM = ctrl_EM,
            ctrl_restr = ctrl_restr,
            ctrl_init = ctrl_init,
            verbose = verbose,
            ...
          )
        if (!is.na(RES[[model_name]][[as.character(g)]]$bic)) {
          if (RES[[model_name]][[as.character(g)]]$bic > bestBIC) {
            RES[["Best"]] <- RES[[model_name]][[as.character(g)]]
            bestBIC <- RES[[model_name]][[as.character(g)]]$bic
          }
        }
      }
      if (verbose) {
        ipbar <- ipbar + 1
        setTxtProgressBar(pbar, ipbar)
      }
    }
    RES$Best$train$X_train <-
      X_train #I also return input training and test data
    RES$Best$test$X_test <- X_test
    if (length(RES$Best$G) == 0) {
      warning("none of the selected models could be fitted")
    } else if (verbose) {
      cat(
        "\nA",
        RES$Best$model_name,
        "patterned model with",
        RES$Best$G,
        "groups was selected"
      )
    }
    class(RES) <- "raedda_t"
    RES
  }
