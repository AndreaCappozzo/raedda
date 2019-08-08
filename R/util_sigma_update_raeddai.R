

INITIALIZE_SIGMA <-
  function(object, fitm) {
    #Function for initializing the variance component according to the desired parsimonious structure, keeping the estimations in the learning phase fixed
    modelName_NEW <- fitm$modelName
    modelName_OLD <- object$modelName
    index <-
      1:object$G
    NEW <- unlist(strsplit(modelName_NEW, split = ""))
    OLD <- unlist(strsplit(modelName_OLD, split = ""))
    if (fitm$d == 1) {
      if (NEW == "E") {
        fitm$parameters$variance$sigmasq <-
          object$parameters$variance$sigmasq
      } else {
        fitm$parameters$variance$sigmasq[index] <-
          object$parameters$variance$sigmasq
        fitm$parameters$variance$scale <-
          fitm$parameters$variance$sigmasq
      }
      return(fitm$parameters$variance)
    }
    DECOMP_OLD <-
      SIGMA_COMP(object$parameters$variance)
    fitm$parameters$variance <-
      SIGMA_COMP(fitm$parameters$variance) # done for models that do not return volume shape orientation automatically
    #Scale
    if (NEW[1] == "E") {
      fitm$parameters$variance$scale <- DECOMP_OLD$scale
    } else {
      fitm$parameters$variance$scale[index] <- DECOMP_OLD$scale
    }
    #Shape
    if (NEW[2] == "E" | NEW[2] == "I") {
      fitm$parameters$variance$shape <- DECOMP_OLD$shape
    } else {
      fitm$parameters$variance$shape[, index] <- DECOMP_OLD$shape
    }
    #Orientation
    if (NEW[3] == "E" | NEW[3] == "I") {
      fitm$parameters$variance$orientation <- DECOMP_OLD$orientation
    } else {
      fitm$parameters$variance$orientation[, , index] <-
        DECOMP_OLD$orientation
    }

    #sigmasq for EII, VII & VVI
    if (modelName_NEW == "EII" |
        modelName_NEW == "VII" | modelName_NEW == "VVI") {
      fitm$parameters$variance$sigmasq <- fitm$parameters$variance$scale
    }

    #Sigma for EII EEI and EEE
    if (modelName_NEW == "EII" |
        modelName_NEW == "EEI" | modelName_NEW == "EEE") {
      fitm$parameters$variance$Sigma <-
        object$parameters$variance$Sigma
    }
    #cholSigma for EEE
    if (modelName_NEW == "EEE") {
      fitm$parameters$variance$cholSigma <-
        chol(fitm$parameters$variance$Sigma)
    }
    #sigma for all models
    fitm$parameters$variance$sigma <- mclust::decomp2sigma(
      fitm$d,
      fitm$G,
      scale =  fitm$parameters$variance$scale,
      shape =  fitm$parameters$variance$shape,
      orientation =  fitm$parameters$variance$orientation
    )

    #cholsigma for VVV
    if (modelName_NEW == "VVV") {
      fitm$parameters$variance$cholsigma <-
        tryCatch(array(as.vector(apply(
          fitm$parameters$variance$sigma,
          3, chol
        )),
        dim = c(fitm$d, fitm$d, fitm$G)), error= function(e) NA)
    }
    return(fitm$parameters$variance)
  }

UPDATE_SIGMA <-
  function(fitm, extra_groups, z, X = Xtest) {
    #extra_groups identifies for which groups sigma needs to be updated
    modelName <- fitm$modelName
    #I compute the necessary quantities to be used for the estimation under different parsimonious structures
    W_extra <-
      mclust::covw(X = X,
           Z = z,
           normalize = F)$W[, , extra_groups, drop = F]
    n_extra <- colSums(z[, extra_groups, drop = F])
    if (fitm$d == 1) {
      if (modelName == "V") {
        sigmasq_extra <- exp(log(as.vector(W_extra)) - log(n_extra))
        #I update the parameters
        fitm$parameters$variance$sigmasq[extra_groups] <-
          sigmasq_extra
        fitm$parameters$variance$scale[extra_groups] <-
          sigmasq_extra
      }
      return(fitm$parameters$variance)
    }
    if (modelName == "VII") {
      scale_extra <-
        apply(W_extra, 3, function(x)
          sum(diag(x))) / (fitm$d * n_extra) # Model \lambda_kI Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
      fitm$parameters$variance$sigmasq[extra_groups] <-
        scale_extra
    }

    if (modelName == "VEI") {
      invB <- diag(1 / fitm$parameters$variance$shape)
      scale_extra <-
        apply(W_extra, 3, function(x)
          sum(diag(x %*% invB))) / (fitm$d * n_extra) # Model \lambda_kB Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
    }

    if (modelName == "EVI") {
      B_extra <-
        apply(W_extra, 3, function(x)
          diag(x) / exp(1 / fitm$d * sum(log(diag(
            x
          ))))) # Model \lambdaB_k Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$shape[, extra_groups] <-
        B_extra
    }

    if (modelName == "VVI") {
      B_extra <-
        apply(W_extra, 3, function(x)
          diag(x) / exp(1 / fitm$d * sum(log(diag(
            x
          ))))) # Model \lambda_kB_k Celeux Govaert
      scale_extra <-
        apply(W_extra, 3, function(x)
          exp(1 / fitm$d * sum(log(diag(
            x
          ))))) / n_extra
      #I update the parameters
      fitm$parameters$variance$shape[, extra_groups] <-
        B_extra
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
      fitm$parameters$variance$sigmasq[extra_groups] <-
        scale_extra
    }

    if (modelName == "EVE") {
      D <- fitm$parameters$variance$orientation
      numA <- apply(W_extra, 3, function(x)
        diag(t(D) %*% x %*% D))
      denA <- apply(numA, 2, function(x)
        exp(1 / fitm$d * sum(log(x))))
      shape_extra <- exp(sweep(
        x = log(numA),
        MARGIN = 2,
        STATS = log(denA),
        FUN = "-"
      )) # Model \lambda_DA_kD' Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$shape[, extra_groups] <-
        shape_extra
    }

    if (modelName == "VEE") {
      C <-
        fitm$parameters$variance$orientation %*% diag(fitm$parameters$variance$shape) %*%
        t(fitm$parameters$variance$orientation)
      scale_extra <-
        apply(W_extra, 3, function(x)
          sum(diag(x %*% solve(C)))) / (fitm$d * n_extra) # Model \lambda_kDAD' Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
    }

    if (modelName == "EEV") {
      orientation_extra <-
        array(apply(W_extra, 3, function(x)
          eigen(x)$vec), dim = dim(W_extra)) # Model \lambdaD_kAD_k' Celeux Govaert
      #I update the parameters
      fitm$parameters$variance$orientation[, , extra_groups] <-
        orientation_extra
    }

    if (modelName == "VEV") {
      invA <- diag(1 / fitm$parameters$variance$shape)
      orientation_extra <-
        array(as.vector(apply(W_extra, 3, function(x)
          eigen(x)$vectors)), dim = dim(W_extra)) # Model \lambda_kD_kAD_k' Celeux Govaert
      scale_extra <-
        sapply(1:length(extra_groups), function(k)
          exp(log(sum(
            diag(
              W_extra[, , k] %*% orientation_extra[, , k] %*% invA %*% t(orientation_extra[, , k])
            )
          )) - log(fitm$d) - log(n_extra[k])))
      #I update the parameters
      fitm$parameters$variance$orientation[, , extra_groups] <-
        orientation_extra
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
    }

    if (modelName == "EVV") {
      C <-
        array(as.vector(apply(W_extra, 3, function(x)
          x / (
            det(x) ^ (1 / fitm$d)
          ))), dim = dim(W_extra)) # Model \lambdaD_kA_kD_k' Celeux Govaert
      orientation_extra <-
        array(as.vector(apply(C, 3, function(x)
          eigen(x)$vectors)), dim = dim(W_extra))
      shape_extra <-
        apply(C, 3, function(x)
          eigen(x, only.values = T)$val)
      #I update the parameters
      fitm$parameters$variance$orientation[, , extra_groups] <-
        orientation_extra
      fitm$parameters$variance$shape[, extra_groups] <-
        shape_extra
    }

    if (modelName == "VVE") {
      D <- fitm$parameters$variance$orientation
      C_k <- #it is named A_k in the paper
        sapply(1:length(extra_groups), function(g)
          diag(t(D) %*% W_extra[, , g] %*% D) / n_extra[g]) # Model \lambda_kDA_kD' Celeux Govaert
      scale_extra <-
        apply(C_k, 2, function(x)
          exp(1 / fitm$d * sum(log(x))))
      shape_extra <-
        exp(sweep(
          log(C_k),
          2,
          STATS = log(scale_extra),
          FUN = "-"
        ))
      #I update the parameters
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
      fitm$parameters$variance$shape[, extra_groups] <-
        shape_extra
    }

    if (modelName == "VVV") {
      eigenval <-
        abs(apply(sweep(W_extra, 3, n_extra, "/"), 3, function(x)
          eigen(x, only.values = TRUE)$val)) # Model \lambda_kD_kA_kD_k' Celeux Govaert
      # NOTE: I add abs() for avoiding numerical problems with almost singular Sample cov matrices for the extra groups
      scale_extra <-
        tryCatch(
          expr = exp(1 / fitm$d * colSums(log(eigenval))),
          warning = function(w)
            NA
        )
      shape_extra <-
        exp(sweep(log(eigenval), 2, log(scale_extra), "-"))
      orientation_extra <- array(apply(W_extra, 3, function(x)
        eigen(x)$vec),
        dim = dim(W_extra))
      #I update the parameters
      fitm$parameters$variance$scale[extra_groups] <- scale_extra
      fitm$parameters$variance$shape[, extra_groups] <-
        shape_extra
      fitm$parameters$variance$orientation[, , extra_groups] <-
        orientation_extra
    }

    #I update the sigma array
    fitm$parameters$variance$sigma <- mclust::decomp2sigma(
      fitm$d,
      fitm$G,
      scale =  fitm$parameters$variance$scale,
      shape =  fitm$parameters$variance$shape,
      orientation =  fitm$parameters$variance$orientation
    )
    if (modelName == "VVV") {
      fitm$parameters$variance$cholsigma <-
        tryCatch(array(as.vector(apply(
          fitm$parameters$variance$sigma,
          3, chol
        )),
        dim = c(fitm$d, fitm$d, fitm$G)), error = function(e)
          NA)
    }
    return(fitm$parameters$variance)
  }
