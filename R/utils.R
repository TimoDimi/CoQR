

loss_model_VaR <- function(theta, df, model, prob_level){
  # The parameter theta only contains the VaR parameters here for the M-estimator!
  beta <- prob_level$beta

  m <- model_fun(theta=theta, df=df, prob_level=prob_level, model=model, model_type="first")
  v <- m$m1
  loss <- (v - df$x) * ((df$x <= v) - beta)
  return(loss)
}



loss_model_CoVaR <- function(theta, df, m1, model, prob_level){
  # The parameter theta only contains the CoVaR parameters here for the second step M-estimator!
  # One has to pass predictions m1 here!
  alpha <- prob_level$alpha

  m <- model_fun(theta=theta, df=df, prob_level=prob_level, model=model, SRM="CoVaR", model_type="second", m1=m1)
  v <- m1
  c <- m$m2

  loss_c <- (df$x > v) * ((c - df$y) * ((df$y <= c) - alpha))

  return(loss_c)
}



collect_data <- function(data, x, y, z=NULL){
  # This function collects the data from either data, or x,y,z and saves it as a data frame

  # Checks if everything is entered correctly!
  if (is.null(data) & (is.null(x) | is.null(y) )){stop("data and (x or y) are not specified")}

  # if (!is.numeric(y) | !is.numeric(x) | !is.numeric(z))
  #   stop("There is a value in x,y, or z that is not numeric!")

  if (any(is.na(y)) | any(is.na(x)) | any(is.na(z)))
    stop("Data contains NAs!")

  if (!(all(is.finite(y)) & all(is.finite(x)) & all(is.finite(z))))
    stop("Not all values are finite!")


  # check dimensions of x,y and z
  # if ( (length(x)!=length(y)) || (length(x)!=dim(z)[1]) ) {stop("Dimensions of x,y, or z do not match!")}

  # Stop if columns of z are multicolinear
  # if (det(t(z) %*% z) < 10^(-8)) {stop("The matrix z is (close to) perfect multicolinearity!")}

  # Collect data
  if (is.null(data)){
    # Case that x,y (and z) are specified
    if (is.null(z)){
      if ( (length(x)!=length(y))) {stop("Dimensions of x and y do not match!")}
      data <- data.frame(x=x, y=y)
    } else {
      # If z is a vector, convert to a nx1 matrix
      if (!is.matrix(z)) { z <- matrix(z, ncol=1, dimnames=list(NULL, c("z1"))) }

      # Delete any constant predictors (i.e., intercepts) of z; add later on
      if (any(apply(z, 2, var)==0)) {
        intercept_index <- which(apply(z, 2, var)==0)
        z <- z[,-intercept_index]
      }

      # Empty colnames of z?
      if(is.null(colnames(z))){
        colnames(z) <- paste0("z", 1:(dim(z)[2]))
      }

      # Add intercept
      z <- cbind("Intercept"=1, z)

      data <- data.frame(x=x, y=y, z) %>% as_tibble()
    }
  } else {
  data <- data %>%
    # dplyr::select(x,y, contains("z"), contains("Date"), contains("Date_index")) %>%
    tibble::as_tibble()
  }

  # Treat empty Date_index and Date columns
  # if(!("Date_index" %in% colnames(data))) {data <- data %>% tibble::add_column(Date_index=1:nrow(data))}
  if(!("Date" %in% colnames(data))) {data <- data %>% dplyr::mutate(Date=lubridate::as_date(1:nrow(data)))}

  data %>% dplyr::select(Date, x, y, everything())
}





#
# collect_data2 <- function(x,y,z=NULL){
#   # This function collects the data and saves it as a data frame
#
#   # if (any(is.na(y)) | any(is.na(xq)) | any(is.na(xe)))
#   #   stop("Data contains NAs!")
#   # if (!(all(is.finite(y)) & all(is.finite(xq)) & all(is.finite(xe))))
#   #   stop("Not all values are finite!")
#
#
#
#   # Add an option for z=NULL
#   if (is.null(z)){
#     if ( (length(x)!=length(y))) {stop("Dimensions of x and y do not match!")}
#     df <- data.frame(x=x, y=y)
#   } else {
#     # If z is a vector, convert to a matrix
#     # ToDo: Add intercept here or later?
#     if (is.numeric(z)) {
#       z <- cbind("(Intercept)"=1,
#                  matrix(z, ncol=1, dimnames=list(NULL, c("z1"))))
#     }
#
#     if (is.matrix(z)) {
#     # check dimensions of x,y and z
#     if ( (length(x)!=length(y)) || (length(x)!=dim(z)[1]) ) {stop("Dimensions of x,y, or z do not match!")}
#
#     # Stop if columns of z are multicolinear
#     if (det(t(z) %*% z) < 10^(-8)) {stop("The matrix z is (close to) perfect multicolinearity!")}
#
#     # Find intercept, and move to the first column
#     if (any(apply(z, 2, var)==0)) {
#       intercept_index <- which(apply(z, 2, var)==0)
#       z <- z[-intercept_index]
#       }
#
#     if(is.null(colnames(z))){ colnames(z) <- paste0("z", 1:(dim(z)[2]))
#
#     df <- data.frame(x=x, y=y, z=z)
#   }
#
#   return(df)
# }
#
#
#
#
#
#
#         if ( !is.null(x) & !is.null(y) ){
#           data <- collect_data(x,y,z) %>%
#             mutate(Date_index=1:n(), Date=lubridate::as_date(Date_index))
#         } else error("Specify either a data frame or numerical vectors x and y.")
#       } else {
#         # Case that x,y (and z) are specified
#
#         data <- data %>%
#           dplyr::select(x,y, contains("z"), contains("Date"), contains("Date_index")) %>%
#           tibble::as_tibble()
#
#
#
#       # if (any(is.na(y)) | any(is.na(xq)) | any(is.na(xe)))
#       #   stop("Data contains NAs!")
#       # if (!(all(is.finite(y)) & all(is.finite(xq)) & all(is.finite(xe))))
#       #   stop("Not all values are finite!")
#
#
#
#       # Add an option for z=NULL
#       if (is.null(z)){
#         if ( (length(x)!=length(y))) {stop("Dimensions of x and y do not match!")}
#         df <- data.frame(x=x, y=y)
#       } else {
#         # If z is a vector, convert to a matrix
#         # ToDo: Add intercept here or later?
#         if (is.numeric(z)) {
#           z <- cbind("(Intercept)"=1,
#                      matrix(z, ncol=1, dimnames=list(NULL, c("z1"))))
#         }
#
#         if (is.matrix(z)) {
#           # check dimensions of x,y and z
#           if ( (length(x)!=length(y)) || (length(x)!=dim(z)[1]) ) {stop("Dimensions of x,y, or z do not match!")}
#
#           # Stop if columns of z are multicolinear
#           if (det(t(z) %*% z) < 10^(-8)) {stop("The matrix z is (close to) perfect multicolinearity!")}
#
#           # Find intercept, and move to the first column
#           if (any(apply(z, 2, var)==0)) {
#             intercept_index <- which(apply(z, 2, var)==0)
#             z <- z[-intercept_index]
#           }
#
#           if(is.null(colnames(z))){ colnames(z) <- paste0("z", 1:(dim(z)[2]))
#
#           df <- data.frame(x=x, y=y, z=z)
#           }
#
#           return(df)
#         }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # Old function with "wo_intercept" option. But I think we do not need it anymore
# collect_data_old <- function(x,y,z, wo_intercept=FALSE){
#   # This function collects the data and saves it as a data frame
#   if (!(is.matrix(z))) {z <- matrix(z, ncol=1)} # If z is a vector, convert to matrix
#   if (det(t(z) %*% z) < 10^(-8)) {stop("The matrix z is (close to) perfect multicolinearity!")}
#
#   if ( (length(x)!=length(y)) || (length(x)!=dim(z)[1]) ) {stop("Dimensions of x,y, or z do not match!")}
#
#
#   # ToDo: Somehow identify the intercept of z! For now, it must be the first component
#   if (wo_intercept){
#     if (var(z[,1])==0) { z <- as.matrix(z[,-1], ncol=dim(z)[2]-1)} # Only delete the intercept if there is one
#     colnames_z <- paste0("z", 1:dim(z)[2])
#     colnames(z) <- colnames_z
#   } else {
#     if (var(z[,1])!=0) { z <- cbind(1,z)} # Only add an intercept if there is none
#     colnames_z <- paste0("z", 1:(dim(z)[2]-1))
#     colnames_z <- c("intercept", colnames_z)
#     colnames(z) <- colnames_z
#   }
#
#   df <- data.frame(x=x, y=y, z)
#
#   return(df)
# }
#
