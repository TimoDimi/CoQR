


#' Joint Dynamic Models for the VaR and CoVaR
#'
#' Estimates a joint dynamic semiparametric model for the pair (VaR, CoVaR):
#' \deqn{VaR_\beta(X_t | Z_t) = v_t(\theta^v)
#' \deqn{CoVaR_\alpha(Y_t | Z_t) = c_t(\theta^c)
#'
#' @param data A tsibble that holds the variables x, y, possibly covariates and an index columns with the name Date
#' @param model Specify the model type; see \link{model_fun} for details
#' @param SRM The systemic risk measure under consideration; currently only the "CoVaR" is implemented
#' @param beta Probability level for the VaR
#' @param alpha Probability level for the CoVaR
#' @param theta0 Starting values for the model parameters. If NULL, then standard values are used
#' @param optim_replications A vector with two integer entries indicating how often the M-estimator of the (VaR, CoVaR) model will be restarted. The default is c(1,3)
#' @param ... Further arguments (does not apply here)
#'
#' @return A 'CoQR' object
#'
#' @seealso
#'
#' @examples
# (1) Simulate bivariate data (x,y) with covariate z
#' eps <- mvtnorm::rmvt(n=1000, sigma=matrix(c(1,0.5,0.5,1),2), df = 8)
#' z <- rnorm(1000)
#' xy <- c(1,1) + cbind(2*z, 2.5*z) + eps
#'
#' # Collect data as tsibble
#' data <- tsibble(Date=1:length(z),
#'                 x=xy[,1],
#'                 y=xy[,2],
#'                 Intercept=1,
#'                 z=z,
#'                 index=Date)
#'
#' # Estimate the 'joint_linear' CoQR regression model
#' obj <- CoQR(data=data,
#'             model = "joint_linear",
#'             beta=0.95, alpha=0.95)
#'
#' # Estimate the standard errors of the parameters
#' summary(obj)
#'
#' # Plot the times series
#' plot(obj)
#'
#'
#'
#' # (2) Simulate bivariate GARCH data
#' library(rmgarch)
#' library(dplyr)
#' data(dji30retw)
#'
#' # Estimate the "CoCAViaR_SAV_diag" model on the negative percentage log-returns
#'
#' data <- dji30retw %>%
#' tibble::rownames_to_column("Date") %>%
#'   mutate(Date=lubridate::as_date(Date)) %>%
#'   select(Date, JPM, BAC) %>%
#'   as_tsibble(index=Date) %>%
#'   rename(x=JPM, y=BAC)


#' obj <- CoQR(data=data,
#'            model="CoCAViaR_SAV_diag")
#'
#' # Covariance estimation and display the parameter estimates with standard errors
#' summary(obj)
#'
#'
#' @references \href{https://arxiv.org/abs/.....}{A Dynamic Co-Quantile Regression}
#' @rdname CoQR
#' @export
CoQR <- function(...) {
  UseMethod('CoQR')
}



#' @rdname CoQR
#' @export
CoQR.default <- function(data=NULL,
                         model="CoCAViaR_SAV_diag", SRM="CoVaR", beta=0.95, alpha=0.95,
                         theta0=NULL, optim_replications=c(1,3)){

  fit <- CoQR.fit(data=data,
                  model=model, SRM=SRM, beta=beta, alpha=alpha,
                  theta0=theta0, optim_replications=optim_replications)

  fit$call <- match.call()
  fit
}




#' @rdname CoQR
#' @export
CoQR.fit <- function(data,
                     model, SRM, beta, alpha,
                     theta0, optim_replications){

  # Collect input data as a tibble
  # data <- collect_data(data=data, x=x, y=y, z=z)

  # Check inputs:
  if (!is_tsibble(data)) stop("Error: Please enter a 'tsibble' object for the argument 'data'.")


  models_implemented <- c("joint_linear",
                          "CoCAViaR_SAV_diag", "CoCAViaR_SAV_fullA", "CoCAViaR_SAV_full",
                          "CoCAViaR_AS_pos", "CoCAViaR_AS_signs", "CoCAViaR_AS_mixed",
                          "CoCAViaR_8pCrossCoVaR", "CoCAViaR_10p", "CoCAViaR_7pVaRViolation")

  if ( !( is.list(model) | (is.character(model) & (model %in% models_implemented)) ) ){
    stop("Please provide for 'model' either a list of functions or a string matching one of the pre-implemented models!")
  }

  TT <- dim(data)[1]
  SRM <- match.arg(SRM, c("MES","CoVaR"))
  prob_level <- switch(SRM,
                       MES = {list(beta=beta, alpha=NA)},
                       CoVaR = {list(beta=beta, alpha=alpha)})

  # Assign default implemented starting values if theta0==NULL
  if (is.null(theta0)){
    theta0_list <- theta_fun(model=model, theta=theta0, df=data)
    theta0 <- theta0_list$theta_start_default
  }
  # Split the starting value!
  theta0_list <- theta_fun(model=model, theta=theta0, df=data)
  theta01 <- theta0_list$theta1
  theta02 <- theta0_list$theta2


  ### VaR Optimization
  thetav_est_rep <- matrix(NA, nrow=optim_replications[1], ncol=length(theta01))
  Mest_obj1 <- rep(NA, optim_replications[1])

  for (rep in 1:optim_replications[1]){
    opt_v <- optim(par=theta01,
                   fn=function(theta,...){
                     mean(loss_model_VaR(theta,...),na.rm=TRUE)},
                   df=data, model=model, prob_level=prob_level)

    thetav_est_rep[rep, ] <- opt_v$par
    Mest_obj1[rep] <- opt_v$value

    theta01 <- thetav_est_rep[rep, ] +
      MASS::mvrnorm(n=1, mu=rep(0,length(theta01)), Sigma=0.1*diag(length(theta01)))
    while(!is.finite(mean(loss_model_VaR(theta01,df=data, model=model, prob_level=prob_level), na.rm=TRUE))){
      theta01 <- thetav_est_rep[rep, ] +
        MASS::mvrnorm(n=1, mu=rep(0,length(theta01)), Sigma=0.1*diag(length(theta01)))
    }
  }
  thetav_est <- thetav_est_rep[which.min(Mest_obj1),]
  m1_est <- model_fun(thetav_est, df=data, prob_level=prob_level, model=model, model_type="first")$m1


  ### Second step: CoVaR/MES optimization
  theta2_est_rep <- matrix(NA, nrow=optim_replications[2], ncol=length(theta02))
  Mest_obj2 <- rep(NA, optim_replications[2])

  for (rep in 1:optim_replications[2]){
    opt_m2 <- optim(par=theta02,
                    fn=function(theta,...){
                      switch(SRM,
                             MES = {mean(loss_model_MES(theta,...),na.rm=TRUE)},
                             CoVaR = {mean(loss_model_CoVaR(theta,...),na.rm=TRUE)})},
                    df=data, m1=m1_est, model=model, prob_level=prob_level)

    theta2_est_rep[rep, ] <- opt_m2$par
    Mest_obj2[rep] <- opt_m2$value


    theta02 <- theta2_est_rep[rep, ] +
      MASS::mvrnorm(n=1, mu=rep(0,length(theta02)), Sigma=0.1*diag(length(theta02)))
    while(!is.finite(switch(SRM,
                            MES = {mean(loss_model_MES(theta02, df=data, m1=m1_est, model=model, prob_level=prob_level), na.rm=TRUE)},
                            CoVaR = {mean(loss_model_CoVaR(theta02, df=data, m1=m1_est, model=model, prob_level=prob_level), na.rm=TRUE)}))){
      theta02 <- theta2_est_rep[rep, ] +
        MASS::mvrnorm(n=1, mu=rep(0,length(theta02)), Sigma=0.1*diag(length(theta02)))
    }
  }
  theta2_est <- theta2_est_rep[which.min(Mest_obj2),]

  theta_est <- c(thetav_est, theta2_est)

  # Create a data.frame with in sample predictions
  m_est  <- model_fun(theta_est, data, prob_level, model, SRM)
  data <- data %>%
    dplyr::mutate(VaR=as.numeric(m_est$m1), !!SRM:=as.numeric(m_est$m2)) %>%
    dplyr::select(Date, x, y, tidyselect::everything(), VaR, CoVaR)

  # create a variable with colnames
  if (model == "joint_linear"){
    colnames_help <- data %>%
      as_tibble() %>%
      dplyr::select(-c(Date, x, y, VaR, CoVaR)) %>%
      colnames()
    colnames <- list(VaR=colnames_help, CoVaR=colnames_help)
  } else {
    colnames <- theta0_list$theta_names
  }

  # Return an object of class "CoQR"
  fit <- list(
    theta=theta_est,
    data=data,
    colnames=colnames,
    SRM=SRM,
    model=model,
    prob_level=prob_level,
    call=match.call()
  )

  class(fit) <- "CoQR"
  fit
}


#' CoQR summary method
#'
#' @param CoQR_object A CoQR object
#' @param method The method to compute the standard errors. Either "asymptotic" or "boot"
#' @param ... Further input parameters
#'
#' @return A summary.CoQR object
#' @export
summary.CoQR <- function(CoQR_object, method='asymptotic',...){
  cov <- cov.CoQR(CoQR_object, method,...)
  se <- sqrt(diag(cov))
  tval <- CoQR_object$theta/se
  coef_mat <- cbind(
    Estimate     = CoQR_object$theta,
    `Std. Error` = se,
    `t value`    = tval,
    `Pr(>|t|)`   = 2 * stats::pt(abs(tval), dim(CoQR_object$data)[1] - length(CoQR_object$theta) , lower.tail = FALSE)
  )

  # VaR and CoVaR coef mats
  coef_mat_VaR <- coef_mat[1:length(CoQR_object$colnames$VaR),]
  rownames(coef_mat_VaR) <- CoQR_object$colnames$VaR

  coef_mat_CoVaR <- coef_mat[-(1:length(CoQR_object$colnames$VaR)),]
  rownames(coef_mat_CoVaR) <- CoQR_object$colnames$CoVaR

  # Return an object of class "summary.CoQR"
  object <- list(cov=cov,
                 se=se,
                 coef_mat_VaR=coef_mat_VaR,
                 coef_mat_CoVaR=coef_mat_CoVaR)
  class(object) <- "summary.CoQR"
  object
}


#' CoQR p arameter covariance estimation
#'
#' @param CoQR_object A CoQR object
#' @param method The method to compute the standard errors. Either "asymptotic" or "boot"
#' @param ... Further input parameters
#'
#' @return covariance matrix
#' @export
cov.CoQR <- function(CoQR_object, method, ...) {
  if (method == 'asymptotic') {
    cov <- vcovA(CoQR_object, ...)
  } else if(method == 'boot') {
    cov <- vcovB(CoQR_object, ...)
  } else if(method == 'boot_fast') {
    cov <- vcovBfast(CoQR_object, ...)
  } else {
    stop('method can be asymptotic, boot, or boot_fast')
  }
  cov
}



#' print method for the sum.CoQR object
#'
#' @param sum.CoQR_object A sum.CoQR object
#'
#' @return Nothing
#' @export
print.summary.CoQR <- function(sum.CoQR_object){
  # Print the VaR and CoVaR coefficients in separate calls:
  cat("\nVaR Coefficients:\n")
  stats::printCoefmat(sum.CoQR_object$coef_mat_VaR,
                      signif.legend=FALSE)
  cat("\nCoVaR Coefficients:\n")
  stats::printCoefmat(sum.CoQR_object$coef_mat_CoVaR)
}



#' print method for the CoQR class
#'
#' @param obj CoQR object
#' @param digits printed digits after the comma
#'
#' @export
print.CoQR <- function(obj, digits=4){
  theta_info <- theta_fun(model=obj$model, theta=obj$theta, df=obj$data %>% dplyr::select(-c("VaR", "CoVaR")))
  q1 <- theta_info$length_theta1
  q2 <- theta_info$length_theta2

  cat("Call:\n")
  cat(deparse(obj$call), "\n")
  cat("\nVaR Parameter Estimates:\n")
  print(format(obj$theta[1:q1], digits = digits), quote = FALSE)
  cat("\nCoVaR Parameter Estimates:\n")
  print(format(obj$theta[(q1+1):(q1+q2)], digits = digits), quote = FALSE)
}


#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot



#' Autoplot a CoQR object
#'
#' @param obj CoQR object
#' @param facet_names Titles of the facets
#'
#' @return A ggplot object
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr `%>%`
autoplot.CoQR <- function(obj, facet_names=c("X / VaR","Y / CoVaR")){
  df_tmp <- obj$data %>%
    dplyr::mutate(VaR_violation=(x > VaR))

  # ToDo: Use dplyr functions instead of reshape2!
  df_long <- dplyr::left_join(
    df_tmp %>%
      dplyr::select(Date, VaR_violation, x, y) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="NegativeReturns"),
    df_tmp %>%
      dplyr::select(Date, VaR_violation, VaR, obj$SRM) %>%
      dplyr::rename(x=VaR, y=obj$SRM) %>%
      reshape2::melt(id.vars=c("Date","VaR_violation"), variable.name="Symbol", value.name="SRMForecasts"),
    by=c("Date", "VaR_violation", "Symbol")) %>%
    tibble::as_tibble()

  # Rename for the facet names
  levels(df_long$Symbol) <- facet_names

  # ggplot2
  p <- ggplot(df_long %>% arrange(VaR_violation)) +
    ggplot2::geom_point(aes(x=Date, y=NegativeReturns, color=VaR_violation)) +
    ggplot2::scale_colour_manual(values = c("grey", "black")) +
    ggnewscale::new_scale_color() + # For two color scales!!!
    ggplot2::geom_line(aes(x=Date, y=SRMForecasts, color=Symbol)) +
    ggplot2::scale_colour_manual(values = c("red", "blue")) +
    ggplot2::facet_wrap(~Symbol, ncol=1) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position="bottom") +
    ggplot2::ylab("Negative Returns")

  p
}



#' CoQR plot method
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.CoQR <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}


#' CoQR fitted method
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
fitted.CoQR <- function(obj){
  m <- model_fun(obj$theta, obj$data, obj$prob_level, obj$model, obj$SRM)
  cbind(m$m1, m$m2)
}



#' CoQR residuals
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
residuals.CoQR <- function(obj){
  cbind(obj$data$x, obj$data$y) - fitted(obj)
}

# Not sure which one to use. Maybe both and use predict for cross-sectional predictions, and forecast for time-series forecasting!
predict.CoQR <- function(obj){
}



#' Generic forecast method
#'
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
forecast <- function(...) {
  UseMethod('forecast')
}


#' CoQR forecast method
#'
#' @param obj CoQR object
#' @param newdata new data set to use for forecasting (if not re-estimated in this time step)
#'
#' @return Vector of forecasts
#' @export
forecast.CoQR <- function(obj, newdata=NULL){
  if (is.null(newdata)) {
    df <- obj$data %>% dplyr::select(x,y)
  } else {
    # ToDo: Check if newdata is reasonable!
    df <- newdata
  }

  m <- model_fun(obj$theta, df=df, prob_level=obj$prob_level, model=obj$model, SRM=obj$SRM, forecast=TRUE)
  forecast <- c(m1_FC=tail(m$m1,1), m2_FC=tail(m$m2,1))
  forecast
}
