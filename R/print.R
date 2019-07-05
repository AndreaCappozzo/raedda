# Printing Methods for the different estimation procedures ------------------------------

print.raedda_t <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust and Adaptive Gaussian finite mixture model for classification \n")
  txt_1 <- paste("Transductive Approach")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat(txt_1, "\n")
  cat( paste0("  ", "Model = ", x$Best$model_name, "\n") )
  cat( paste0("  ", "G = ", x$Best$G, "\n") )
  cat( paste0("  ", "Training trimming level = ", x$Best$train$alpha_train, "\n") )
  cat( paste0("  ", "Test trimming level = ", x$Best$test$alpha_test, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$Best$bic, 3), "\n") )
}

print.raedda_i <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust and Adaptive Gaussian finite mixture model for classification \n")
  txt_1 <- paste("Inductive Approach")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat(txt_1, "\n")
  cat( paste0("  ", "Model = ", x$discovery_phase$Best$model_name, "\n") )
  cat( paste0("  ", "G = ", x$discovery_phase$Best$G, "\n") )
  cat( paste0("  ", "Training trimming level = ", x$learning_phase$Best$train$alpha_train, "\n") )
  cat( paste0("  ", "Discovery trimming level = ", x$discovery_phase$Best$test$alpha_discovery, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$discovery_phase$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$discovery_phase$Best$bic, 3), "\n") )
}

print.redda <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust and Adaptive Gaussian finite mixture model for classification \n")
  txt_1 <- paste("Inductive Approach: Learning phase")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat(txt_1, "\n")
  cat( paste0("  ", "Model = ", x$Best$model_name, "\n") )
  cat( paste0("  ", "G = ", x$Best$G, "\n") )
  cat( paste0("  ", "Training trimming level = ", x$Best$train$alpha_train, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$Best$bic, 3), "\n") )
}

print.raedda_d <- function(x, ...)
{
  cat("\n")
  txt <- paste(" ", "Robust and Adaptive Gaussian finite mixture model for classification \n")
  txt_1 <- paste("Inductive Approach: Discovery phase")
  sep <- paste0(rep("-", max(nchar(txt)) + 1),
                collapse = "")
  cat(sep, "\n")
  cat(txt)
  cat(sep, "\n")
  cat("\n")
  cat(txt_1, "\n")
  cat( paste0("  ", "Model = ", x$Best$model_name, "\n") )
  cat( paste0("  ", "G = ", x$Best$G, "\n") )
  cat( paste0("  ", "Discovery trimming level = ", x$Best$test$alpha_discovery, "\n") )
  cat( paste0("  ", "Log-likelihood= ", round(x$Best$ll, 3), "\n") )
  cat( paste0("  Robust BIC= ", round(x$Best$bic, 3), "\n") )
}
