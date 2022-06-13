
#' Title
#'
#' @param data
#' @param model
#' @param length_IS
#' @param refit_freq
#' @param SRM
#' @param beta
#' @param alpha
#' @param theta0
#' @param optim_replications
#'
#' @return
#' @export
#'
#' @examples
CoQRroll <- function(data, model="CoCAViaR_SAV_diag",
                     length_IS=1000, refit_freq=100,
                     SRM="CoVaR", beta=0.95, alpha=0.95,
                     theta0=NULL, optim_replications=c(1,3)){

  # data must be a data.frame with columns Date, Date_index, x, y
  TT <- dim(data)[1]
  length_OOS <- TT - length_IS
  refit_points <- seq(length_IS+1,TT,by=refit_freq)

  FC_df <- data.frame()
  CoQR_objects <- list()

  # Loop over all OOS days
  for (tt in (length_IS+1):TT){
    data_tt <- data[(tt-length_IS):(tt-1),]

    # Only refit at certain points!
    if (tt %in% refit_points){

      # Iterative starting values from the previous fit
      if (tt==refit_points[1]) theta_start <- theta0 else theta_start <- CoQR_obj$theta

      # ToDo: Error handling if fit fails?
      CoQR_obj <- CoQR(data=data_tt, model=model,
                       SRM=SRM, beta=beta, alpha=alpha,
                       theta0=theta_start, optim_replications=optim_replications)
      CoQR_objects <- append(CoQR_objects, list(CoQR_obj))
    }
    FCs <- forecast(CoQR_obj,
                    newdata=data_tt %>% dplyr::select(x,y))

    FC_df <- dplyr::bind_rows(FC_df,
                       data.frame(Date_index=tt, m1_FC=FCs[1], m2_FC=FCs[2],row.names = NULL))
  }

  rownames(FC_df) <- NULL
  FC_df <- FC_df %>%
    dplyr::rename(VaR=m1_FC, !!SRM:=m2_FC) %>%
    dplyr::left_join(data, by="Date_index") %>%
    tibble::as_tibble() %>%
    dplyr::select(Date, Date_index, x, y, VaR, !!SRM)


  # Return an object of class "CoQRroll"
  obj <- list(
    FC_df=FC_df,
    data=df,
    SRM=SRM,
    alpha=alpha,
    beta=beta,
    CoQR_objects=CoQR_objects)

  class(obj) <- "CoQRroll"

  obj
}


#' Title
#'
#' @param obj
#'
#' @return
#' @export
#'
#' @examples
print.CoQRroll <- function(obj){
  print(obj$FC_df)
}





#' Title
#'
#' @param obj
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.CoQRroll <- function(obj, ...){
  p <- autoplot(obj, ...)
  print(p)
}





#' Title
#'
#' @param obj
#' @param facet_names
#'
#' @return
#' @export
#'
#' @examples
#' @import ggplot2
#' @importFrom magrittr `%>%`
autoplot.CoQRroll <- function(obj, facet_names=c("X / VaR","Y / CoVaR")){
  # Add options etc! this is a little lame so far...
  # Shall we include the estimates into the df?
  # ToDo: Improve facet and legend labels

  df_tmp <- obj$FC_df %>%
    dplyr::mutate(VaR_violation=(x > VaR))

  df_long <- dplyr::left_join(
    df_tmp %>%
      dplyr::select(Date, Date_index, VaR_violation, x, y) %>%
      reshape2::melt(id.vars=c("Date","Date_index","VaR_violation"), variable.name="Symbol", value.name="NegativeReturns"),
    df_tmp %>%
      dplyr::select(Date, Date_index, VaR_violation, VaR, obj$SRM) %>%
      dplyr::rename(x=VaR, y=obj$SRM) %>%
      reshape2::melt(id.vars=c("Date","Date_index","VaR_violation"), variable.name="Symbol", value.name="SRMForecasts"),
    by=c("Date", "Date_index", "VaR_violation", "Symbol")) %>%
    tibble::as_tibble()

  levels(df_long$Symbol) <- facet_names

  p <- ggplot2::ggplot(df_long %>% arrange(VaR_violation)) +
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


