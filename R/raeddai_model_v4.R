raedda_l_model <- function(X_train,
                           class_train,
                           alpha_train,
                           # The proportion of obs in the Xtrain to be trimmed
                           model_name,
                           ctrl_init,
                           ...) {

  N_train <- nrow(X_train)
  ltrain <- mclust::unmap(class_train)
  d <- ncol(X_train)
  G <- ncol(ltrain)
  class_train <- as.factor(class_train)
  classLabel <- levels(class_train)

  # Core algorithm ----------------------------------------------------------


  if (alpha_train != 0) {

    nsamp <- ctrl_init$n_samp
    max_iter_init <- ctrl_init$max_iter

    N_train_trim <- N_train - ceiling(N_train * alpha_train)
    robust_result <- patterned_MCD(nsamp = nsamp, # I perform the estimation starting from nsamp J_g subsets of sample size (d+1), inspired to what done in \cite{Hubert2018}
                                   X_train,
                                   class_train,
                                   ltrain,
                                   alpha_train,
                                   model_name, max_iter_init)
    ll <- robust_result$ll
    fitm <- robust_result$fitm

  } else {
    N_train_trim <- NULL
    fitm <- mclust::mstep(modelName = model_name,
                          data = X_train,
                          z = ltrain)
    mTau.train <-
      matrix(log(fitm$parameters$pro),
             nrow(X_train),
             fitm$G,
             byrow = TRUE)
    lDensity.train <- do.call(mclust::cdens, c(list(
      data = X_train,
      logarithm = TRUE
    ), fitm))

    sum.train <- mTau.train + lDensity.train
    mat.train <- ltrain * sum.train
    ll <- sum(mat.train)
  }


  # Checking if errors in the procedure -------------------------------------

  fitetrain <-
    tryCatch(
      do.call(mclust::estep, c(list(data = X_train), fitm)),
      error = function(e) {
        list(z = NA)
      }
    )
  emptyz <- ifelse(all(!is.na(fitetrain$z)), yes = FALSE, no = TRUE)

  # Results Collection ------------------------------------------------------

  if (!emptyz) {
    res <- list()
    res$N_train <- N_train
    res$N_train_after_trimming <- N_train_trim
    res$alpha_train <- alpha_train
    res$d <- d
    res$G <- G
    res$model_name <- model_name
    res$parameters <- fitm$parameters
    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if (d == 1) {
      names(res$parameters$mean) <- classLabel
    } else {
      colnames(res$parameters$mean) <- classLabel
    }

    ztrain <- fitetrain$z
    cltrain <-
      factor(sapply(map(ztrain), function(i)
        classLabel[i]), levels = classLabel) # I classify a posteriori also the trimmed units
    pos_trimmed_train <- NULL
    cltrain_after_trimming <- NULL
    if (alpha_train != 0) {
      D_Xtrain_cond <- do.call(mclust::cdens, c(list(
        data = X_train, # computing the component density, this is done because I am interested in trimming also obs that might have
        logarithm = T
      ), fitm)) # been WRONGLY assinged to a class
      ind_D_Xtrain_cdens <-
        cbind(1:N_train, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
      D_Xtrain <-
        D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
      pos_trimmed_train <-
        # I trim the conditional density \phi(x_n; \mu_g, \Sigma_g) when x_n comes from group g
        which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                           [[ceiling(N_train * alpha_train)]]))
      cltrain_after_trimming <- cltrain
      cltrain_after_trimming <-
        factor(cltrain, levels = c(classLabel, "0"))
      cltrain_after_trimming[pos_trimmed_train] <- "0"
    }

    res$train <- list()
    res$train$z <- ztrain
    res$train$cl <- cltrain
    res$train$cl_after_trimming <- cltrain_after_trimming
    res$train$obs_trimmed <- pos_trimmed_train
    res$train$alpha_train <- alpha_train
    bic_all <-
      2 * ll - mclust::nMclustParams(
        modelName = model_name,
        d = d,
        G = G,
        noise = FALSE,
        equalPro = FALSE
      ) * log(fitm$n) # if alpha_train !=0 this is trimmed BIC
    res$ll <- ll
    res$bic <- bic_all
  }
  else {
    res <- list()
    res$error <-
      "Either training groups too small or trimming level too large for this model type"
    res$ll <- NA
    res$bic <- NA
  }
  res
}


raedda_d_model <- function(fit_learning,
                           X_test,
                           model_name,
                           G,
                           alpha_discovery,
                           #alpha in the discovery phase
                           restr_factor_d,
                           # the eigenvalue restriction to be applied to the (C-G) extra groups only
                           ctrl_EM,
                           ctrl_restr,
                           ctrl_init,
                           ...) {

  N_test_original <- nrow(X_test)
  if (!is.null(fit_learning$Best$N_train_after_trimming)) {
    # if trimming was applied in the learning phase, I create an augmented test set including the trimmed obs from the train
    X_test <-
      rbind(X_test, fit_learning$Best$train$X_train[fit_learning$Best$train$obs_trimmed, , drop =
                                              F])
  }
  N_test <- nrow(X_test)
  d <- ncol(X_test)
  N_train_trimmed <- N_test - N_test_original
  fitm_learning <- fit_learning$Best
  fitm_learning$n <- N_test
  fitm_learning$modelName <- fitm_learning$model_name # for compatibility with Mclust functions
  if (alpha_discovery != 0) {
    N_test_trim <- N_test - ceiling(N_test * alpha_discovery)
  } else {
    N_test_trim <- NULL
  }
  X_test_fit <- X_test

  iter <- 0

  # Same group as Learning Phase, nothing to estimate -----------------------

  if (G == fit_learning$Best$G) {
    # if G is set to be the same as number of grps in the learning phase I do not have to estimate anything,
    # so I simply save the parameters obtained in the learning phase and I compute the loglik
    fitm <- fitm_learning
    fite <- do.call(mclust::estep, c(list(data = X_test), fitm))
    z <- fite$z
    z_fit <- fite$z
    emptyz <- TRUE
    if (all(!is.na(fite$z)))
    {
      emptyz <- FALSE
    }
    if (alpha_discovery != 0) {
      D <-
        do.call(mclust::dens, c(list(
          data = X_test, #compute the Density for Parameterized MVN Mixtures for each obs
          logarithm = F
        ), fitm)) #I temporarily discard those alpha_test% of obs whose density is lowest
      pos_trimmed_test <-
        which(D <= (sort(D, decreasing = F)[[ceiling(N_test * alpha_discovery)]]))
      z_fit <- z[-pos_trimmed_test, , drop = F]
      X_test_fit <- X_test[-pos_trimmed_test, , drop = F]
    }
    ll <-
      sum(do.call(mclust::dens, c(
        list(data = X_test_fit, logarithm = TRUE), # value of the trimmed log-likelihood
        fitm
      )))
  } else {

    # G_discovery > G_learning ------------------------------------------------

    index <-
      1:fit_learning$Best$G

    # Initialization for the H hidden classes and EM algorithm ---------------------------------

    EM_tol <- ctrl_EM$tol
    EM_max_iter <- ctrl_EM$max_iter
    aitken <- ctrl_EM$aitken
    nstart_EM <- ctrl_EM$nstart_EM
    extra_groups <- setdiff(1:G, index)
    rep_EM_results_collection <- vector("list", length = nstart_EM)

    for (rep_EM in 1:nstart_EM) {
    # The estimation needs to be carried out only for C-G groups,
    # since I am deploying an inductive approach

      fitm <-
        init_H_hidden_classes(
          fitm_TRAIN=fitm_learning,
          X_test = X_test,
          model_name = model_name,
          G = G,
          N_test = N_test,
          d = d
        )

    # For such groups I check if the eigenvalue-ratio is satisfied, if not I enforce
    # the constraints

      if (!any(is.na(fitm$parameters$variance$sigma))) {
      # if there are NA or NULL something went wrong and I will not compute the constrained maximization
      if((fitm$modelName=="VVE"|fitm$modelName=="EVE") & (G != fit_learning$Best$G)){
        fitm$X <- X_test # I add the data on which the M-step is computed since I need them for the MM
      }
      suppressWarnings(fitm <-
                         constr_Sigma(
                           fitm = fitm,
                           restr_factor = restr_factor_d,
                           extra_groups = extra_groups, # constraints enforced only on the extra groups (Inductive approach)
                           ctrl_restr = ctrl_restr
                         )) # this performs constrained estimation of Sigma according to the selected model, that is, the initial values in the EM algorithm satisfy the eigenvalues-ratio
    }

    # EM algorithm ------------------------------------------------------------

    llold <- -Inf
    ll <- -Inf
    llstore <- rep(-Inf, 3) # for Aitken
    criterion <- TRUE

    while (criterion) {
      iter <- iter + 1
      fite <-
        tryCatch(
          do.call(mclust::estep, c(list(data = X_test), fitm)),
          error = function(e)
            list(z = NA)
        ) #expectation step
      emptyz <- TRUE
      if (all(!is.na(fite$z)))
      {
        emptyz <- FALSE
        z <- fite$z
        z_fit <- fite$z
        if (alpha_discovery != 0) {
          D <-
            do.call(mclust::dens, c(list(
              data = X_test, #compute the Density for Parameterized MVN Mixtures for each obs
              logarithm = F
            ), fitm)) #I temporarily discard those alpha_test% of obs whose density is lowest
          pos_trimmed_test <-
            which(D <= (sort(D, decreasing = F)[[ceiling(N_test * alpha_discovery)]]))
          z_fit <- z[-pos_trimmed_test, , drop = F]
          X_test_fit <- X_test[-pos_trimmed_test, , drop = F]
        }

        #For extra_groups I manually update pro, mean and variance
        pro_extra <- colMeans(z_fit[, extra_groups, drop = F])
        fitm$parameters$pro <-
          c(fit_learning$Best$parameters$pro * (1 - sum(pro_extra)), pro_extra) # proportions are re-estimated for each class
        if (fitm$d == 1) {
          fitm$parameters$mean[extra_groups] <-
            mclust::covw(X = X_test_fit,
                         Z = z_fit,
                         normalize = F)$mean[, extra_groups] # I update the mu vector just for the extra groups
        } else {
          fitm$parameters$mean[, extra_groups] <-
            mclust::covw(X = X_test_fit,
                         Z = z_fit,
                         normalize = F)$mean[, extra_groups] # I update the mu vector just for the extra groups
        }

        fitm$parameters$variance <-
          UPDATE_SIGMA(
            fitm = fitm,
            extra_groups = extra_groups,
            z = z_fit,
            X = X_test_fit
          )

        if (!(any(is.na(fitm$parameters$variance$sigma))|any(is.na(fitm$parameters$variance$cholsigma)))) {
          # if there are NA something went wrong and I will not compute the constrained maximization
          if(fitm$modelName=="VVE"|fitm$modelName=="EVE"){
            fitm$X <- X_test_fit # I add the data on which the M-step is computed since I need them for the MM
          }

          suppressWarnings(fitm <-
                             constr_Sigma(
                               fitm = fitm,
                               restr_factor = restr_factor_d,
                               extra_groups = extra_groups, # constraints enforced only on the extra groups (Inductive approach)
                               ctrl_restr = ctrl_restr
                             ) )# this performs constrained estimation of Sigma according to the selected model
        }

        ll <-
          suppressWarnings(tryCatch(
            sum(do.call(mclust::dens, c(
              list(data = X_test_fit, logarithm = TRUE), #value of the trimmed log-likelihood
              fitm
            ))),
            error = function(e)
              - Inf
          ))
        ll <-
          ifelse(is.nan(ll) |
                   is.na(ll), -Inf, ll) #If llis NA or Nan it proceeds till the next estep and then the loop breaks
        llstore <- c(llstore[-1], ll)
        if (aitken) {
          criterion <- (Aitken(llstore)$linf - ll) > EM_tol
        } else {
          criterion <- (ll - llold) > EM_tol
        }
        criterion <- (criterion) & (iter < EM_max_iter)
        llold <- ll
      }
      else {
        criterion <- FALSE
      }
    }
    rep_EM_results_collection[[rep_EM]] <-
      list(ll = ll,
           fitm = fitm,
           iter = iter)
    }

    # Best solution out of the nstart_EM initializations

    best_EM_iteration <- which.max(sapply(1:nstart_EM, function(EM_rep) rep_EM_results_collection[[EM_rep]]$ll))
    ll <- rep_EM_results_collection[[best_EM_iteration]]$ll
    fitm <- rep_EM_results_collection[[best_EM_iteration]]$fitm
    iter <- rep_EM_results_collection[[best_EM_iteration]]$iter

}
  # Results Collection ------------------------------------------------------


  classLabel <- names(fit_learning$Best$parameters$pro)
  if (G > fit_learning$Best$G) {
    classLabel <-
      c(classLabel, paste0("HIDDEN_GROUP_", 1:(G - fit_learning$Best$G)))
  }
  res <- list() # is this an empty list
  if (!emptyz) {
    res$N_test_augmented <- N_test
    res$N_test <- N_test_original
    res$N_train_trimmed <-
      N_train_trimmed #these obs come from the tr but were not used in the robust learning phase since discarded by the trimming approach
    res$N_test_augmented_after_trimming <- N_test_trim
    res$d <- d
    res$G <- G
    if(G != fit_learning$Best$G) {
    res$iter <- iter
    res$converged <- (iter < EM_max_iter)
    }
    res$model_name <- model_name
    res$parameters <- fitm$parameters
    #Name the different groups
    names(res$parameters$pro) <- classLabel
    if (d == 1) {
      names(res$parameters$mean) <- classLabel
    } else {
      colnames(res$parameters$mean) <- classLabel
    }
    
    fite_test <- do.call(mclust::estep, c(list(data = X_test),
                                          fitm))
    z_test <- fite_test$z
    
    cl <-
      factor(sapply(mclust::map(z_test), function(i)
        classLabel[i]), levels = classLabel)
    pos_trimmed_test = NULL
    cltest_after_trimming = NULL
    if (alpha_discovery != 0) {
      D <-
        do.call(mclust::dens, c(list(
          data = X_test, #compute the Density for Parameterized MVN Mixtures for each obs
          logarithm = F
        ), fitm))
      pos_trimmed_test <-
        which(D <= (sort(D, decreasing = F)[[ceiling(N_test * alpha_discovery)]]))
      cltest_after_trimming <-
        factor(cl, levels = c(classLabel, "0"))
      cltest_after_trimming[pos_trimmed_test] <- "0"
      fitm$n <- N_test_trim # used later in the computation of the TBIC
    }
    res$test_augmented <- list()
    res$test_augmented$z <- z_test
    res$test_augmented$cl <- cl
    res$test_augmented$cl_after_trimming <- cltest_after_trimming
    res$test_augmented$obs_trimmed <- pos_trimmed_test
    res$test <- list()
    res$test$z <- z_test[1:N_test_original, ]
    res$test$cl <- cl[1:N_test_original]
    res$test$cl_after_trimming <-
      cltest_after_trimming[1:N_test_original]
    res$test$obs_trimmed <-
      pos_trimmed_test[pos_trimmed_test <= N_test_original]
    res$test$alpha_discovery <- alpha_discovery
    res$training_trimmed <- list()
    res$training_trimmed$z <- z_test[-(1:N_test_original), ]
    res$training_trimmed$cl <- cl[-(1:N_test_original)]
    res$training_trimmed$cl_after_trimming <-
      cltest_after_trimming[-(1:N_test_original)]
    res$training_trimmed$obs_trimmed <-
      pos_trimmed_test[pos_trimmed_test > N_test_original]
    bic_all <-
      robust_bic_raedda_d(
        modelName = model_name,
        loglik = ll,
        n = fitm$n,
        d = fitm$d,
        H = fitm$G - fit_learning$Best$G,
        restr_factor_d = restr_factor_d
      )
    res$ll <- ll
    res$bic <- bic_all
  } else {
    res$error <- "Groups too small for this model type"
    res$ll <- NA
    res$bic <- NA
  }
  res
}
