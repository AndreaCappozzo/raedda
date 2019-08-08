

# Robust Initialization for the G known groups ----------------------------------------------------------

patterned_MCD <-
  # main function that will be repeated nsamp times, run only if alpha_Xtrain!=0
  function(nsamp,
           Xtrain,
           cltrain,
           ltrain,
           alpha_Xtrain,
           modelName,
           max_iter_init) {
    robust_start <-
      function(Xtrain,
               cltrain,
               ltrain,
               alpha_Xtrain,
               modelName,
               max_iter_init) {
        Ntrain <- nrow(Xtrain)
        J_ind <-
          tryCatch(
            c(sapply(
              split(seq(cltrain), cltrain),
              sample,
              size = ncol(Xtrain) + 1,
              replace = FALSE
            )),
            error = function(e)
              seq(cltrain)
          )

        #Init with the random (p + 1)-subset J_g
        fitm <- tryCatch(
          mclust::mstep(
            modelName = modelName,
            data = Xtrain[J_ind, , drop = F],
            z = ltrain[J_ind,]
          ),
          error = function(e) {
            list(parameters = NA)
          }
        )
        Xtrain_fit <-
          NULL # I initialize its value for preventing erros when the EM algorithm fails
        emptyz <- TRUE
        pos_old <-
          1:ceiling(Ntrain * alpha_Xtrain)  # used in the while loop
        # initial random value just using the first ceiling(Ntrain * alpha_Xtrain) obs
        # the estimation will continue running until for two consecutive iterations exactly
        # the same obs are trimmed
        criterion <- TRUE
        iter_init <- 0
        while (criterion) {
          iter_init <- iter_init + 1
          # in performing robust estimation I discard those obs with smallest value of
          # component density given their label
          fitetrain <-
            tryCatch(
              do.call(mclust::estep, c(list(data = Xtrain), fitm)),
              error = function(e) {
                list(z = NA)
              }
            )
          emptyz <- TRUE
          if (all(!is.na(fitetrain$z))) {
            emptyz <- FALSE
            D_Xtrain_cond <- do.call(mclust::cdens, c(list(
              data = Xtrain, # computing the component density, this is done because I am interested in trimming also obs that might have
              logarithm = T
            ), fitm)) # been WRONGLY assinged to a class
            ind_D_Xtrain_cdens <-
              cbind(1:Ntrain, mclust::map(ltrain)) # matrix with obs and respective group from ltrain
            D_Xtrain <-
              D_Xtrain_cond[ind_D_Xtrain_cdens] # I compute D_g conditioning ON the fact that I know the supposed true class
            pos_trimmed_train <-
              which(D_Xtrain <= (sort(D_Xtrain, decreasing = F)
                                 [[ceiling(Ntrain * alpha_Xtrain)]]))
            if (length(pos_trimmed_train) != ceiling(Ntrain * alpha_Xtrain)) {
              #condition due to the fact that for some models (the simplest ones usually) it might happen that 2 obs have exactly the same comp prob, leading to errors
              pos_trimmed_train <-
                pos_trimmed_train[1:ceiling(Ntrain * alpha_Xtrain)]
            }
            ltrain_fit <- ltrain[-pos_trimmed_train, , drop = F]
            Xtrain_fit <- Xtrain[-pos_trimmed_train, , drop = F]
            fitm <-
              tryCatch(
                mclust::mstep(
                  modelName = modelName,
                  data = Xtrain_fit,
                  z = ltrain_fit
                ),
                error = function(e) {
                  list(parameters = NA)
                }
              )
            if (all(pos_old == pos_trimmed_train)) {
              criterion <- FALSE
            } else {
              pos_old <- pos_trimmed_train
              criterion <- (criterion) & (iter_init < max_iter_init)
            }
          }
          else {
            criterion <- FALSE
          }
        }
        if (!emptyz) {
          mTau.train <-
            matrix(log(fitm$parameters$pro),
                   nrow(Xtrain_fit),
                   fitm$G,
                   byrow = TRUE)
          lDensity.train <- do.call(mclust::cdens, c(list(
            data = Xtrain_fit,
            logarithm = TRUE
          ), fitm))

          sum.train <- mTau.train + lDensity.train
          mat.train <- ltrain_fit * sum.train
          ll <- sum(mat.train)
        } else {
          ll <- NA
        }
        fitm$X <-
          Xtrain_fit # this is used in the initialization of raedda_transductive_model: for VVE and EVE models for which the constraint is not satisfied the data from which the M-step is computed are needed for the MM algorithm
        return(list(ll = ll, fitm = fitm))
      }
    n_init <-
      replicate(
        nsamp,
        expr = robust_start(
          Xtrain,
          cltrain,
          ltrain,
          alpha_Xtrain,
          modelName,
          max_iter_init
        ),
        simplify = T
      )

    if (!all(is.na(n_init[1,]))) {
      ind <- which.max(n_init[1, ])
    } else {
      ind <-
        1 #if no initial subsample worked, get the first one and then the alg will break down further in the code
    }
    return(n_init[, ind])

  }


# Initialization for the H hidden classes ---------------------------------


init_H_hidden_classes <-
  function(fitm_TRAIN, X_test,
           model_name,
           G, N_test, d) {
    G_train <- fitm_TRAIN$G
    H <- G - G_train
    index <-
      1:G_train #for the known groups I will use the estimated parameters obtained in the learning phase

    extra_groups <- setdiff(1:G, index)

    J_h_ind <-
      # draw a random (p + 1) subset for each hidden class
      tryCatch(
        sample(
          x = 1:N_test,
          size = H * (d + 1),
          replace = FALSE
        ),
        error = function(e){
          sample(
            x = 1:N_test,
            size = (d + 1), # if N_test is small, draw at least a sample of d+1 units
            replace = FALSE
          )
        }
      )

    z_fake <-
      matrix(0, nrow = N_test, ncol = G) # construct z_fake cos I need to get the same structure of model_name in the m-step
    z_fake[matrix(c(J_h_ind, rep(extra_groups, each = length(J_h_ind) /
                                   H)), ncol = 2)] <-
      1 # out of the H(p+1) units I have randomly chosen from Y, H goes to the first hidden group, another H goes for the second hidden group etc...
    # this is a workaround in order to make it work regardless of the number of obs J_hidden extracted and number of hidden groups H
    z_fake[-J_h_ind, index] <-
      1 / G_train # I then initialized the other obs to have all(rowSums(z_fake)==1)==TRUE: the estimates for the known group will not be retained

    fitm <-
      #Note that this step is done only for obtaining the correct structure for the variance covariance matrix according to modelName, parameters estimated here will not be retained
      mclust::mstep(modelName = model_name,
                    data = X_test,
                    z = z_fake)


    #For extra_groups I manually update pro, mean and variance
    unif_sample <- runif(n = H)
    pro_extra <- unif_sample / sum(unif_sample) * H / G
    fitm$parameters$pro <-
      c(fitm_TRAIN$parameters$pro * G_train / G, pro_extra)
    if (fitm$d == 1) {
      fitm$parameters$mean <-
        c(
          fitm_TRAIN$parameters$mean,
          mclust::covw(
            X = X_test,
            Z = z_fake,
            normalize = F
          )$mean[, extra_groups]
        )
    } else {
      fitm$parameters$mean <-
        cbind(
          fitm_TRAIN$parameters$mean,
          mclust::covw(
            X = X_test,
            Z = z_fake,
            normalize = F
          )$mean[, extra_groups]
        ) # I update the mu vector just for the extra groups
    }
    fitm$parameters$variance <-
      INITIALIZE_SIGMA(object = fitm_TRAIN, fitm = fitm) #this step is done in order to have the correct parsimonious structure for the variance component, keeping the estimations in the learning phase fixed


    fitm$parameters$variance <-
      UPDATE_SIGMA(
        fitm = fitm,
        extra_groups = extra_groups,
        z = z_fake,
        X = X_test
      )

    fitm
  }

#   n_init <-
#     replicate(
#       nstart,
#       expr = RANDOM_START(z_known_groups, fitm_TRAIN, G, Xtest, modelName),
#       simplify = T
#     )
#   ind <- which.max(n_init[1,])
#   return(n_init[2, ind]$fitm)
# }
