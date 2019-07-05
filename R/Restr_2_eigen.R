#General function that applies restr maximization according to modelName,
# to be used in the transductive and inductive approach via the formal extra_groups
# If extra groups == C (all groups: known + hidden) transductive approach
# If extra groups == C-G == H (only hidden) inductive approach

constr_Sigma <-
  function(fitm, restr_factor, extra_groups, ctrl_restr) {
    if (fitm$modelName == "EVI" |
        fitm$modelName == "EVE" |
        fitm$modelName == "EVV") {
      fitm <- restr_2_EV_(fitm, restr_factor, extra_groups, ctrl_restr)
    } else if (fitm$modelName == "VII" | fitm$modelName == "VEI" | fitm$modelName == "VEE" |
               fitm$modelName == "VEV") {
      fitm <- restr_2_VE_(fitm, restr_factor, extra_groups, ctrl_restr)
    } else if (fitm$modelName == "V" | fitm$modelName == "VVI" | fitm$modelName == "VVE" |
               fitm$modelName == "VVV") {
      fitm <- restr_2_VV_(fitm, restr_factor, extra_groups, ctrl_restr)
    }
    return(fitm)
  }

# Functions from TCLUST package -------------------------------------------


restr2_eigenv <-
  function (autovalues, ni.ini, restr.fact, zero.tol) {
    # original function for performing eigenvalues-ratio restrictions from TCLUST package
    ev <- autovalues

    ###### function parameters:
    ###### ev: matrix containin eigenvalues									#n proper naming - changed autovalue to ev (eigenvalue)
    ###### ni.ini: current sample size of the clusters						#n proper naming - changed ni.ini to cSize (cluster size)
    ###### factor: the factor parameter in tclust program
    ###### init: 1 if we are applying restrictions during the smart inicialization
    ######       0 if we are applying restrictions during the C steps execution

    ###### some inicializations

    if (!is.matrix (ev))
      #n	checking for malformed ev - argument (e.g. a vector of determinants instead of a matrix)
      if (is.atomic (ev))
        #n
        ev <- t (as.numeric (ev))									#n
      else
        #n
        ev <- as.matrix (ev)										#n

      stopifnot (ncol (ev) == length (ni.ini))							#n	check wether the matrix autovalues and ni.ini have the right dimension.

      d <- t (ev)

      p <- nrow (ev)
      K <- ncol (ev)

      n <- sum(ni.ini)

      nis <- matrix(data = ni.ini,
                    nrow = K,
                    ncol = p)

      #m	MOVED: this block has been moved up a bit. see "old position of "block A" for it's old occurrence
      idx.nis.gr.0 <-
        nis > zero.tol										#n	as nis is in R we have to be carefull when checking against 0
      used.ev <- ni.ini > zero.tol										#n
      ev.nz <- ev[, used.ev]												#n	non-zero eigenvalues
      #m
      #o	if ((max (d[nis > 0]) <= zero.tol))									#m
      #i	if ((max (d[idx.nis.gr.0]) <= zero.tol))							#n	"idx.nis.gr.0" instead of (nis > 0)
      if ((max (ev.nz) <= zero.tol))
        #n	simplify syntax
        return (matrix (0, nrow = p, ncol = K))							#m
      #m
      ###### we check if the  eigenvalues verify the restrictions			#m
      #m
      #o	if (max (d[nis > 0]) / min (d[nis > 0]) <= restr.fact)				#m
      #i	if (max (d[idx.nis.gr.0]) / min (d[idx.nis.gr.0]) <= restr.fact)	#n	"idx.nis.gr.0" instead of (nis > 0)
      if (max (ev.nz) / min (ev.nz) <= restr.fact)
        #n	simplify syntax

      {
        #m
        #o		d[!idx.nis.gr.0] <- mean (d[idx.nis.gr.0])						#n	"idx.nis.gr.0" instead of (nis > 0)
        ev[,!used.ev] <- mean (ev.nz)									#n	simplify syntax
        return (ev)														#m
      }																	#m

      ###### d_ is the ordered set of values in which the restriction objective function change the definition
      ###### points in d_ correspond to  the frontiers for the intervals in which this objective function has the same definition
      ###### ed is a set with the middle points of these intervals

      #o	d_ <- sort (c (d, d / restr.fact))
      d_ <-
        sort (c (ev, ev / restr.fact))								#n	using ev instead of d
      dim <-
        length (d_)													##2do: continue here cleaning up restr2_eigenv
      d_1 <- d_
      d_1[dim + 1] <- exp(log(d_[dim]) + log(2))
      d_2 <- c (0, d_)
      ed <- exp(log(d_1 + d_2) - log(2))
      dim <- dim + 1


      ##o
      ##o	old position of "block A"
      ##o

      ###### the only relevant eigenvalues are those belong to a clusters with sample size greater than 0.
      ###### eigenvalues corresponding to a clusters with 0 individuals has no influence in the objective function
      ###### if all the eigenvalues are 0 during the smart initialization we assign to all the eigenvalues the value 1

      ###### we build the sol array
      ###### sol[1],sol[2],.... this array contains the critical values of the interval functions which defines the m objective function
      ###### we use the centers of the interval to get a definition for the function in each interval
      ###### this set with the critical values (in the array sol) contains the optimum m value

      t <- s <- r <- array(0, c(K, dim))
      sol <- sal <- array(0, c(dim))

      for (mp_ in 1:dim)
      {
        for (i in 1:K)
        {
          r[i, mp_] <-
            sum ((d[i,] < ed[mp_])) + sum((d[i,] > ed[mp_] * restr.fact))
          s[i, mp_] <- sum (d[i,] * (d[i,] < ed[mp_]))
          t[i, mp_] <- sum (d[i,] * (d[i,] > ed[mp_] * restr.fact))
        }

        sol[mp_] <-
          exp(
            log(sum (exp(log(ni.ini) - log(n) + log((s[, mp_] + t[, mp_] / restr.fact))))) -
              log(sum(exp(log(ni.ini) - log(n) + log (r[, mp_])))))

        e <-	sol[mp_] * (d < sol[mp_]) +
          d * (d >= sol[mp_]) * (d <= restr.fact * sol[mp_]) +
          (restr.fact * sol[mp_]) * (d > restr.fact * sol[mp_])
        #o <- -1 / 2 *exp(log(nis) - log(n) + log(log(e) + exp(log(d)-log(e))))
        o <- -1 / 2 *exp(log(nis) - log(n)) * (log(e) + exp(log(d)-log(e)))

        sal[mp_] <- sum(o)
      }

      ###### m is the optimum value for the eigenvalues procedure
      #o	eo <- which.max (c (sal))						## remove c ()
      #	m <- sol[eo]
      m <- sol[which.max (sal)]						#n
      ###### based on the m value we get the restricted eigenvalues

      t (m * (d < m) + d * (d >= m) * (d <= restr.fact * m) + (restr.fact * m) * (d > restr.fact * m))	##	the return value
  }
restr2_deter_ <-
  function (autovalues, ni.ini, restr.fact, zero.tol = 1e-16){
    ###### function parameters:
    ###### autovalues: matrix containing eigenvalues
    ###### ni.ini: current sample size of the clusters
    ###### factor: the factor parameter in tclust program
    ###### some initializations

    p = nrow (autovalues)

    if (p == 1)
      return  (restr2_eigenv (autovalues, ni.ini, restr.fact, zero.tol))

    K <-  ncol (autovalues)

    es <-  apply (autovalues, 2, function(x) exp(sum(log(x)))) # determinants for each cluster
    #log_es <-  apply (autovalues, 2, function(x) sum(log(x))) # log of the determinants
    idx.ni.ini.gr.0 <-
      ni.ini > zero.tol		#groups with size larger than zero.tol (checking for empty clusters) #n	as ni.ini is in R we have to be carefull when checking against 0

    ######	we check if all the determinants in no empty populations are 0
    #o	if (max(es[ni.ini > 0]) <= zero.tol)	##	all eigenvalues are somehow zero.
    if (max(es[idx.ni.ini.gr.0]) <= zero.tol)
      #n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)
      return (matrix (0, p, K))												##		-> return zero mat

    ######  we put in d the determinants of the populations (converting es into a matrix of dim 1 x K)
    d = t(es)	#### --> dim (d) = 1 x K (has once been "d <- matrix (es, nrow = 1)")
    #d = t(log_es)	#### --> dim (d) = 1 x K (has once been "d <- matrix (es, nrow = 1)")
    #n	put this block into a function (for improved readability)
    autovalues_det <-
      HandleSmallEv (autovalues, zero.tol)						##	handling close to zero eigenvalues here
    #autovalues_det are the normalized eigenvalues
    #cat ("\n1d^(1/p):\t", d^(1/p), "\n")

    ######	we check if all the determinants verify the restrictions
    #o	if (max (d[ni.ini > 0]) / min (d[ni.ini > 0]) <= restr.fact)
    if (max (d[idx.ni.ini.gr.0]) / min (d[idx.ni.ini.gr.0]) <= restr.fact)
      #n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)
    {
      #o		d [ni.ini == 0] <-  mean (d[ni.ini > 0])								## and get the mean - determinants for all clusters without observations.
      d [!idx.ni.ini.gr.0] <-
        mean (d[idx.ni.ini.gr.0])						#n	"idx.ni.ini.gr.0" instead of (ni.ini > 0)

      dfin <- exp(1/p*log(d))
    }
    else
      dfin <-
      restr2_eigenv (autovalues = exp(1/p*log(d)), ni.ini, restr.fact = restr.fact ^ (1 / p), zero.tol) #restricted determinats
    ######  we apply the restriction to the determinants by using the restr2_eigenv function
    ######  In order to apply this function is neccessary to transform d and factor with the power (1/p)
    #cat ("\nfin:\t", dfin, "\n")
    #.multbyrow (autovalues_det, dfin)											## autovalues_det %*% diag (dfin)
    exp(sweep(
      x = log(autovalues_det),
      MARGIN = 2,
      STATS = log(dfin),
      FUN = "+"
    ))
  }

HandleSmallEv <-
  function (autovalues, zero.tol){
    #n	a part of restr2_deter_, which handles almost zero eigenvalues
    ##	handling close to zero eigenvalues here
    ######  populations with one eigenvalue close to 0 are very close to be contained in a hyperplane
    ######  autovalues2 is equal to autovalues except for the columns corresponding to populations close to singular
    ######  for these populations we put only one eigenvalue close to 0 and the rest far from 0
    #K <- nrow (autovalues)													#n

    #o		autovalues[autovalues < zero.tol] = zero.tol
    autovalues[autovalues <= zero.tol] <-
      zero.tol							#n	"<= zero.tol" for checking for zero

    mi <-
      apply(autovalues, 2, min)											##	the minimum eigenvalue of each cluster
    ma <-
      apply(autovalues, 2, max)											##	the maximum eigenvalue of each cluster

    #o		idx.iter <- (1:K) [ma/mi>1 / zero.tol]									#o	all clusters which have almost zero eigenvalues
    idx.iter <-
      which (mi / ma <= zero.tol)									#n	making more obvious for what to check!

    for (i in idx.iter)
      ##	for each of these clusters. set all "normal" eigenvalues to a high value
      autovalues[autovalues[, i] > mi[i] / zero.tol, i] <-
      mi[i] / zero.tol

    #o		es2 = apply(autovalues, 2, prod)
    det <-  apply (autovalues, 2, function(x) exp(sum(log(x))))											#n	the new determinants
    #logdet <-  apply (autovalues, 2, function(x) sum(log(x)))											#n	the new determinants

    ######	autovalues_det contains the autovalues corrected by the determinant
    ######	the product of the eigenvalues of each column in autovalues_det is equal to 1

    #o		autovalues_det <- .multbyrow (autovalues, es2^(-1/p))

    p <- nrow (autovalues)													#n
    #autovalues_det <- .multbyrow (autovalues, det^(-1/p))					#n
    autovalues_det <-
      exp(sweep(
        x = log(autovalues),
        MARGIN = 2,
        STATS = (1 / p)*log(det),
        #STATS = (1 / p)*logdet,
        FUN = "-"
      )) #no using additional function and more efficient
    #autovalues_det are the normalized eigenvalues
    return (autovalues_det)
  }

# Majorization Minimization Algorithm Browne2014

MM1 <-
  function(D_0 = diag(dim(W)[1]),
           A,
           W,
           w,
           ctrl_restr) {
    D <- D_0
    F_obj <- Inf
    F_obj_old <- Inf
    iter <- 0
    criterion <- TRUE
    while (criterion) {
      iter <- iter + 1
      F_t <-
        matrix(apply(
          sapply(1:dim(W)[3], function(g)
            diag(1 / A[, g]) %*% t(D) %*% W[, , g] - w[g] * diag(1 / A[, g]) %*% t(D)),
          1,
          sum
        ), dim(W)[1], dim(W)[1], byrow = F)
      SVD_F <- svd(F_t)
      P <- SVD_F$u
      R <- SVD_F$v
      D <- R %*% t(P)
      F_obj  <-
        sum(sapply(1:dim(W)[3], function(g)
          sum(diag(
            W[, , g] %*% D %*% diag(1 / A[, g]) %*% t(D)
          ))))
      criterion <-
        ((F_obj_old- F_obj) > ctrl_restr$MM_tol) &
        (iter < ctrl_restr$MM_max_iter)
      F_obj_old <- F_obj
    }
    D
  }

# Restrictions for EV_, VE_ and VV_ models -------------------------------------

restr_2_EV_ <- function(fitm, restr_factor, extra_groups,
                           ctrl_restr) {
  eigenvalues <-
    exp(log(fitm$parameters$variance$shape[,extra_groups, drop=FALSE]) + log(fitm$parameters$variance$scale))
  if (ifelse(
    is.na(max(eigenvalues) / min(eigenvalues) <= restr_factor) |
    is.infinite(max(eigenvalues) / min(eigenvalues)),
    TRUE,
    max(eigenvalues) / min(eigenvalues) <= restr_factor
  )) {
    return(fitm)
  }
  iter <- 0
  while ((max(eigenvalues) / min(eigenvalues) - restr_factor > ctrl_restr$tol) &
         iter <= ctrl_restr$max_iter) {
    iter <- iter + 1
    eigenvalues <- restr2_eigenv(
      autovalues = eigenvalues,
      ni.ini = colSums(fitm$z[,extra_groups, drop=FALSE]),
      restr.fact = restr_factor,
      zero.tol = .Machine$double.eps
    )

    #2nd step: .restr2_deter
    if(length(extra_groups) == fitm$G) { # Transductive approach
      eigenvalues <- restr2_deter_(
        autovalues = eigenvalues,
        ni.ini = colSums(fitm$z[,extra_groups, drop=FALSE]),
        #ni.ini = rep(1, fitm$G),
        restr.fact = 1,
        #I force the determinants to be equal
        zero.tol = .Machine$double.eps
      )
    } else { # Inductive approach

      scale_extra_group <- apply(eigenvalues, 2, function(x)
        exp(1 / fitm$d*sum(log(x))))

      shape_extra_group <- exp(sweep(
        x = log(eigenvalues),
        2,
        log(scale_extra_group),
        FUN = "-"
      ))

      eigenvalues <- exp(sweep(
        x = log(shape_extra_group),
        2,
        log(fitm$parameters$variance$scale), # same volume of the known groups
        FUN = "+"
      ))
    }
  }

  scale_RESTR <- apply(eigenvalues, 2, function(x)
    exp(1 / fitm$d*sum(log(x))))[1]
  shape_RESTR <- exp(log(eigenvalues) - log(scale_RESTR))

  if(fitm$modelName == "EVE" & length(extra_groups) == fitm$G) { # Transductive approach
    W <-
      mclust::covw(fitm$X, Z = fitm$z, normalize = F)$W #sample scatter matrices for the un-trimmed obs
    E <-
      apply(W, 3, function(x)
        eigen(x, only.values = T)$val) #eigenvalues of W
    w <- apply(E, 2, max) #biggest eigenv of W_g, g=1, ..., G
    orientation_RESTR <- MM1(A = shape_RESTR, W = W, w = w, ctrl_restr = ctrl_restr)
    fitm$parameters$variance$orientation <- orientation_RESTR
  }
  #I update the variance components in the output of the M-step
  fitm$parameters$variance$shape[, extra_groups] <- shape_RESTR
  fitm$parameters$variance$scale <- scale_RESTR
  fitm$parameters$variance$sigma <- mclust::decomp2sigma(
    d = fitm$d,
    G = fitm$G,
    scale = fitm$parameters$variance$scale,
    shape = fitm$parameters$variance$shape,
    orientation = fitm$parameters$variance$orientation
  )
  return(fitm)
}

restr_2_VE_ <- function(fitm, restr_factor, extra_groups, ctrl_restr) {

  if(fitm$modelName == "VII"){ #done in order to obtain scale shape orientation components for VII model
    fitm$parameters$variance <- SIGMA_COMP(fitm$parameters$variance)
  }
  eigenvalues <-
    fitm$parameters$variance$shape %o% fitm$parameters$variance$scale[extra_groups]

  if (ifelse(is.na(max(eigenvalues)/min(eigenvalues) <= restr_factor)|is.infinite(max(eigenvalues)/min(eigenvalues)), TRUE, max(eigenvalues)/min(eigenvalues) <= restr_factor)) {
    return(fitm)
  }

  scale_0 <- scale_RESTR <- fitm$parameters$variance$scale[extra_groups] # original scale
  shape_0 <- shape_RESTR <- fitm$parameters$variance$shape # original shape

  iter <- 0
  while ((max(eigenvalues) / min(eigenvalues) - restr_factor > ctrl_restr$tol) &
         iter <= ctrl_restr$max_iter) {
    iter <- iter + 1
    eigenvalues <-
      restr2_eigenv(
        autovalues = eigenvalues,
        ni.ini = colSums(fitm$z[, extra_groups, drop=FALSE]), # the groups size
        restr.fact = restr_factor,
        zero.tol = .Machine$double.eps
      )

    scale_0 <- scale_RESTR
    shape_0 <- shape_RESTR

    if(length(extra_groups) == fitm$G) { # Transductive approach
      shape_RESTR <-
        apply(sapply(1:fitm$G, function(g)
          exp(log(eigenvalues[, g]) - log(scale_RESTR[g]))), 1, sum) /
        (exp(1 / fitm$d*sum(log((apply(sapply(1:fitm$G, function(g)
          eigenvalues[, g] / scale_RESTR[g]), 1, sum))))))
    } else { # Inductive approach
      shape_RESTR <- shape_0
    }
    scale_RESTR <-
      sapply(1:length(extra_groups), function(g)
        sum(exp(log(eigenvalues[, g])-log(shape_RESTR))) / fitm$d)
    eigenvalues <-  shape_RESTR %o% scale_RESTR
  }
  #I update the variance components in the output of the M-step
  fitm$parameters$variance$shape <- shape_RESTR
  fitm$parameters$variance$scale[extra_groups] <- scale_RESTR
  if (fitm$modelName == "VII") { #I also update sigmasq for VII model
    fitm$parameters$variance$sigmasq[extra_groups] <- fitm$parameters$variance$scale[extra_groups]
  }
  fitm$parameters$variance$sigma <-
    mclust::decomp2sigma(
      d = fitm$d,
      G = fitm$G,
      scale = fitm$parameters$variance$scale,
      shape = fitm$parameters$variance$shape,
      orientation = fitm$parameters$variance$orientation
    )
  return(fitm)
}

restr_2_VV_ <-
  function(fitm,
           restr_factor,
           extra_groups,
           ctrl_restr) {
    if (fitm$modelName == "V") {
      max_sigmasq <- max(fitm$parameters$variance$sigmasq)
      min_sigmasq <- min(fitm$parameters$variance$sigmasq)
      #univariate (d=1) case
      if (max_sigmasq / min_sigmasq <= restr_factor) {
        return(fitm)
      } else if (length(extra_groups) == 1) {
        # when there is just one extra group
        # and it has variance smaller than every other known group,
        # I heuristically protect against spurious solutions
        # truncating the extra sigmasq to the smallest known sigmasq
        if (all(fitm$parameters$variance$sigmasq[extra_groups] <
                fitm$parameters$variance$sigmasq[-extra_groups])) {
          sigmasq_RESTR <-
            min(fitm$parameters$variance$sigmasq[-extra_groups])
        } else { # if it is not the smallest I keep it as it is
          # since for sure it will not cause spurious solutions
          sigmasq_RESTR <- fitm$parameters$variance$sigmasq[extra_groups]
        }
      } else {
        sigmasq_RESTR <- as.vector(
          restr2_eigenv(
            fitm$parameters$variance$sigmasq[extra_groups],
            ni.ini = colSums(fitm$z[, extra_groups, drop = FALSE]),
            restr.fact = restr_factor,
            zero.tol = .Machine$double.eps
          )
        )
      }
      fitm$parameters$variance$sigmasq[extra_groups] <-
        sigmasq_RESTR
      fitm$parameters$variance$scale[extra_groups] <- sigmasq_RESTR
      return(fitm)
    }

    if (fitm$modelName == "VVV") {
      #done in order to obtain scale shape orientation components for VVV model
      fitm$parameters$variance <-
        SIGMA_COMP(fitm$parameters$variance)
    }

    eigenvalues <-
      exp(sweep(
        x = log(fitm$parameters$variance$shape[, extra_groups, drop = FALSE]),
        2,
        log(fitm$parameters$variance$scale[extra_groups]),
        FUN = "+"
      ))

    if (ifelse(
      is.na(max(eigenvalues) / min(eigenvalues) <= restr_factor) |
      is.infinite(max(eigenvalues) / min(eigenvalues)),
      TRUE,
      max(eigenvalues) / min(eigenvalues) <= restr_factor
    )) {
      # the ifelse is needed to prevent errors due to almost singular covariance matrices
      return(fitm)
    }

    eigenvalues_RESTR <-
      restr2_eigenv(
        eigenvalues,
        ni.ini = colSums(fitm$z[, extra_groups, drop = FALSE]),
        restr.fact = restr_factor,
        zero.tol = .Machine$double.eps
      )

    scale_RESTR <- apply(eigenvalues_RESTR, 2, function(x)
      exp(1 / fitm$d * sum(log(x))))
    shape_RESTR <- exp(sweep(
      x = log(eigenvalues_RESTR),
      2,
      log(scale_RESTR),
      FUN = "-"
    ))

    if (fitm$modelName == "VVE" &
        length(extra_groups) == fitm$G) { # Transductive approach
      W <-
        mclust::covw(fitm$X, Z = fitm$z, normalize = F)$W #sample scatter matrices for the un-trimmed obs
      E <-
        apply(W, 3, function(x)
          eigen(x, only.values = T)$val) #eigenvalues of W
      w <- apply(E, 2, max) #biggest eigenv of W_g, g=1, ..., G
      orientation_RESTR <-
        MM1(
          A = shape_RESTR,
          W = W,
          w = w,
          ctrl_restr = ctrl_restr
        )
      fitm$parameters$variance$orientation <- orientation_RESTR
    }
    #I update the variance components in the output of the M-step
    fitm$parameters$variance$shape[, extra_groups] <- shape_RESTR
    fitm$parameters$variance$scale[extra_groups] <- scale_RESTR
    fitm$parameters$variance$sigma <-
      mclust::decomp2sigma(
        d = fitm$d,
        G = fitm$G,
        scale = fitm$parameters$variance$scale,
        shape = fitm$parameters$variance$shape,
        orientation = fitm$parameters$variance$orientation
      )
    #cholsigma for VVV
    if (fitm$modelName == "VVV") {
      fitm$parameters$variance$cholsigma <-
        array(as.vector(apply(
          fitm$parameters$variance$sigma,
          3, chol
        )),
        dim = c(fitm$d, fitm$d, fitm$G))
    }
    fitm
  }

