

# Utility control functions -----------------------------------------------

control_init <- function(n_samp = 50,
                         n_start_extra_classes = 30,
                         max_iter = 1e02) {
  # Initialization control paramaters
  list(
    n_samp = n_samp,
    n_start_extra_classes = n_start_extra_classes,
    max_iter = max_iter
  )
}

control_EM <- function(nstart_EM = 5,
                       tol = 1e-05,
                       max_iter = 1e02,
                       aitken = TRUE)
  # EM control parameters
{
  list(
    nstart_EM = nstart_EM,
    tol = tol,
    max_iter = max_iter,
    aitken = aitken
  )
}

control_restr <-
  function(tol = 1e-10,
           max_iter = 1e02,
           MM_tol = 1e-10,
           MM_max_iter = 1e04)
    # Constrained maximization control parameters
    #MM is related to Majorization - Minimization algorithm (used for VVE, EVE models)
  {
    list(
      tol = tol,
      max_iter = max_iter,
      MM_tol = MM_tol,
      MM_max_iter = MM_max_iter
    )
  }
