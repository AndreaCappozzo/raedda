raedda_l <-
  function(X_train,
           class_train,
           alpha_train = 0,
           model_names = NULL,
           ctrl_init = control_init(),
           verbose = interactive(),
           ...) {
    X_train <- data.matrix(X_train)
    if (is.null(model_names)) {
      if(ncol(X_train)==1){
        model_names <- c("E", "V")
      } else {
        model_names <- mclust::mclust.options("emModelNames")
      }
    }
    if (verbose) {
      cat("fitting Learning Phase...\n")
      flush.console()
      pbarL <- txtProgressBar(min = 0,
                              max = length(model_names),
                              style = 3)
      on.exit(close(pbarL))
      ipbarL <- 0
    }
    RES <- list()
    bestBIC <- -Inf
    RES[["Best"]] <- list()
    for (model_name in model_names) {
      RES[[model_name]] <- list()
      RES[[model_name]] <- raedda_l_model(X_train,
                                           class_train,
                                           alpha_train,
                                           model_name,
                                           ctrl_init = ctrl_init,
                                           ...)
      if (!is.na(RES[[model_name]]$bic)) {
        if (RES[[model_name]]$bic > bestBIC) {
          RES[["Best"]] <- RES[[model_name]]
          bestBIC <- RES[[model_name]]$bic
        }
      }
      if (verbose) {
        ipbarL <- ipbarL + 1
        setTxtProgressBar(pbarL, ipbarL)
      }
    }
    RES$Best$train$X_train <-
      X_train #I also return input training data
    if (verbose) {
      cat("\nA",
          RES$Best$model_name,
          "patterned model was selected in the Learning Phase")
    }
    class(RES) <- "redda"
    RES
  }

#ROBUST Adaptive Eigenvalues Decomposition Discriminant Analysis :
# with this version we can jointly handle:
#noise labels (between existing groups in both test and training set)
#noise labels (with unobserved groups wrongly labelled as belonging to existing group
#in both test and training set)
#outliers in both test and training

#INDUCTIVE APPROACH
#1st step- Learning phase: robust parameter estimation employing robust discriminant analysis
#2nd step- Discovery phase: I perform an EM algorithm on the augmented test set (test + trimmed obs)

raedda_d <- function(fit_learning,
                     #the result from the robust learning phase
                     X_test,
                     model_names = NULL,
                     G = NULL,
                     #the number of possible groups considered in the discovery phase
                     alpha_discovery = .05,
                     #alpha in the discovery phase
                     restr_factor_d = NULL,
                     ctrl_EM = control_EM(),
                     ctrl_restr = control_restr(),
                     ctrl_init = control_init(),
                     verbose = interactive(),
                     ...) {
  X_test <- data.matrix(X_test)
  #check if the selected model can be estimated according to the partial ordering defined in the igraph
  #NOTE: for the **I models some structures are redundant and henceforth I do not consider them.
  #For example if we learnt EEI the models EEI and EEE are equivalent and need not be estimated
  modelscope_check <- switch(
    fit_learning$Best$model_name,
    "E" = c("E", "V"),
    "V" = "V",
    "EII" = c(
      "EII",
      "VII",
      "VEI",
      "EVI",
      "VVI",
      "EVE",
      "EEV",
      "VEV",
      "EVV",
      "VVV"
    ),
    "VII" = c("VII", "VVI",
              "VEV", "VVV"),
    "EEI" = c(
      "EEI",
      "VEI",
      "EVI",
      "VVI",
      "EVE",
      "VEE",
      "EEV",
      "VEV",
      "EVV",
      "VVV"
    ),
    "VEI" = c("VEI", "VVI", "VEV", "VVV"),
    "EVI" = c("EVI", "VVI", "EVE", "VVE", "EVV", "VVV"),
    "VVI" =  c("VVI", "VVE", "VVV"),
    "EEE" = c("EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV"),
    "EVE" = c("EVE", "EVV", "VVE", "VVV"),
    "VEE" = c("VEE", "VEV", "VVE", "VVV"),
    "VVE" = c("VVE", "VVV"),
    "EEV" = c("EEV", "VEV", "EVV", "VVV"),
    "VEV" = c("VEV", "VVV"),
    "EVV" = c("EVV", "VVV"),
    "VVV" = c("VVV")
  )

  if (is.null(model_names)) {
    model_names <- modelscope_check
  } else if (any(!(model_names %in% modelscope_check))) {
    stop(
      'One or more of the selected models cannot be estimated due to order incompatibility with the parsimonious structure selected in the learning phase'
    )
  }
  G <- sort(G)
  if (is.null(G)) {
    G <-
      fit_learning$Best$G:(fit_learning$Best$G + 2) #if I do not specify a specific G, I will consider the G current labels of the learning phase, G+1 and G+2
  }

  if (G[1] < fit_learning$Best$G) {
    stop(
      "The expected number of groups must be equal or greater than the number of groups in the learning phase"
    )
  }

  if(is.null(restr_factor_d)){

    # the learning phase gives an idea of the magnitude of the eigenvalue-ratio
    # in the known groups, this information can be exploited
    # for avoiding setting in advance the value of restr_factor_d

    eigenvalues_known_groups <-
      tryCatch(
        apply(fit_learning$Best$parameters$variance$sigma, 3, function(s2)
          eigen(s2, only.values = T)$val),
        error = function(e)
          1
      )

    restr_factor_train <-
      abs(max(eigenvalues_known_groups) / min(eigenvalues_known_groups)) #in order to avoid spurious solution, I can set restr.factor to be no larger than the eigenvalue-ratio in the known groups. Abs is added for avoiding numerical problems

    restr_factor_d <- restr_factor_train
  }

  if (verbose) {
    cat("fitting Discovery Phase...\n")
    flush.console()
    pbarL <- txtProgressBar(min = 0,
                            max = length(model_names),
                            style = 3)
    on.exit(close(pbarL))
    ipbarL <- 0
  }

  out <- list()
  bestBIC <- -Inf
  out[["Best"]] <- list()
  for (model_name in model_names) {
    for (g in setdiff(G, fit_learning$Best$G)) {
      # setdiff(G, object$Best$G) means I do not loop for G==object$Best$G since the parameters won't be updated
      out[[model_name]][[as.character(g)]] <- list()
      out[[model_name]][[as.character(g)]] <-
        raedda_d_model(
          fit_learning = fit_learning,
          X_test = X_test,
          model_name = model_name,
          G = g,
          alpha_discovery = alpha_discovery,
          restr_factor_d = restr_factor_d,
          ctrl_EM = ctrl_EM,
          ctrl_restr = ctrl_restr,
          ctrl_init = ctrl_init
        )

      if (!is.na(out[[model_name]][[as.character(g)]]$bic)) {
        if (out[[model_name]][[as.character(g)]]$bic > bestBIC) {
          out[["Best"]] <- out[[model_name]][[as.character(g)]]
          bestBIC <- out[[model_name]][[as.character(g)]]$bic
        }
      }
    }
    if (verbose) {
      ipbarL <- ipbarL + 1
      setTxtProgressBar(pbarL, ipbarL)
    }
  }
  #I check whether the best model is the one with the same number of groups as in the learning phase
  if (G[1] == fit_learning$Best$G) {
    out[[fit_learning$Best$model_name]][[as.character(fit_learning$Best$G)]] <-
      raedda_d_model(
        fit_learning = fit_learning,
        X_test = X_test,
        model_name = fit_learning$Best$model_name,
        G =
          fit_learning$Best$G,
        alpha_discovery = alpha_discovery,
        restr_factor_d = restr_factor_d,
        restr_factor = restr_factor,
        ctrl_EM = ctrl_EM,
        ctrl_restr = ctrl_restr,
        ctrl_init = ctrl_init,
        verbose = verbose
      )
    if (out[[fit_learning$Best$model_name]][[as.character(fit_learning$Best$G)]]$bic > bestBIC) {
      out[["Best"]] <-
        out[[fit_learning$Best$model_name]][[as.character(fit_learning$Best$G)]]
    }
  }

  if (length(out$Best$G) == 0) {
    warning("none of the selected models could be fitted")
  } else if (verbose) {
    cat(
      "\nA",
      out$Best$model_name,
      "patterned model with",
      out$Best$G,
      "groups was selected in the Discovery Phase"
    )
  }
  class(out) <- "raedda_d"
  out
}
# RAEDDA_i is just a wrapper around REDDA and RAEDDA_d: it simultaneously fits both learning and discovery phase
# via an inductive approach, giving both Xtrain and Xtest as inputs
raedda_i <- function (X_train,
                      class_train,
                      X_test,
                      model_names_l = NULL,
                      model_names_d = NULL,
                      G = NULL,
                      # the number of possible groups considered in the discovery phase
                      alpha_train = 0,
                      alpha_discovery = 0.05,
                      restr_factor_d = NULL,
                      ctrl_EM = control_EM(),
                      ctrl_restr = control_restr(),
                      ctrl_init = control_init(),
                      verbose = interactive(),
                      ...)
{

  learning_phase <-
    raedda_l(
      X_train = X_train,
      class_train = class_train,
      alpha_train = alpha_train,
      model_names = model_names_l,
      ctrl_init = ctrl_init,
      verbose = verbose,
      ...
    )
  if (verbose) cat("\n")
  discovery_phase <-
    raedda_d(
      fit_learning = learning_phase,
      X_test = X_test,
      model_names = model_names_d,
      G = G,
      alpha_discovery = alpha_discovery,
      restr_factor_d = restr_factor_d,
      ctrl_EM = ctrl_EM,
      ctrl_restr = ctrl_restr,
      ctrl_init = ctrl_init,
      verbose = verbose,
      ...
    )

  RES <- list()
  RES$learning_phase <- learning_phase
  RES$discovery_phase <- discovery_phase
  class(RES) <- "raedda_i"
  RES
}
