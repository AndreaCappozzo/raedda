
# Predict function to be used with raedda_l output ------------------------

# It's simply a wrapper around raedda_d where no extra classes are sought
# Consistent with Conventions for R Modeling Packages
# https://tidymodels.github.io/model-implementation-principles/index.html

predict.redda <- function(object, new_data, ...){
    if (!inherits(object, "redda"))
      stop("object not of class \"redda\"")
    classLabel <- names(object$Best$parameters$pro)
    fite <- estep(
      modelName = object$Best$model_name,
      data = new_data,
      parameters = object$Best$parameters
    )
    z <- fite$z
    cl <-
      factor(sapply(mclust::map(z), function(i)
        classLabel[i]), levels = classLabel)
    out <- list(class = cl, prob = z)
    out
  }
