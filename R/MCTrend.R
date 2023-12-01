#' Monte Carlo Trend Analysis
#'
#' This function performs Monte Carlo trend analysis on input data and generates plots.
#'
#' @export
#'
#' @param x A data frame containing the input data. The first raw expected to contain model names or time series names.
#' @param n_rep Number of replications for the Monte Carlo simulation.
#' @param plot_title Title for the plot.
#' @param int A number indicating lower threshold value of the interval within which no trend is defined, the upper value is calculated based on this value, by default a lower value of 0.25 is considered.
#' @param opt A number indicating type of results, for opt = 1 returns test result, opt = 2 returns plot
#'
#'
#' @return A data frame and a plot containing results of the trend analysis.
#'
#'
#' @importFrom dplyr select relocate mutate
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot geom_rect geom_line geom_point geom_text facet_wrap theme_light labs aes theme
#' @importFrom lmomco parwei lmoms cdfwei pp quawei
#' @importFrom magrittr %<>% %>%
#' @importFrom trend sens.slope
#'
#' @examples
#'
#'   # file for example
#'   file <- MCTrend::example
#'
#'   # Apply the test
#'   MCTrend::MCTrend(x = file, n_rep = 100, plot_title = 'Precipitaciones', int = 0.1, opt = 1)
#'
#'   # plot of the result of the test
#'   MCTrend::MCTrend(x = file, n_rep = 100, plot_title = 'Precipitaciones', int = 0.1, opt = 2)
#'

MCTrend <- function(x, n_rep, plot_title, int = 0.25, opt) {

  X_obs <- Y_obs <- eq_x <- eq_y <- NULL # Evitar problemas con "Undefined global functions or variables"  al chequear paquete

  Model <- x %>% select(-1) %>% colnames() %>% as.data.frame() # Extrae nombre de los modelos o de las columnas
  MC_trend1 <- data.frame(Model, Model) # Crea un DF auxiliar de 2 columnas, 1a para los nombres y 2a para los valores de la tendencia
  colnames(MC_trend1) <- c('Model','MC_trend') # Cambia nombre de los encabezados del DF


  T_S_MC_slope <- list() # Lista auxiliar vacia
  WEI_T_S <- list() # Lista auxiliar vacia
  Obs_sen_slope <- list() # Lista auxiliar vacia
  Quawei_T_S <- list() # Lista auxiliar vacia

  for (i in 2:dim(x)[2]){
    v <- x[,i]
    set.seed(i)
    result <- as.data.frame(replicate(n_rep, sample(v, dim(x)[1], replace = FALSE)))  # n_rep es el numero de replicas o repeticiones desordenadas del vector original
    T_S_MC <- result %>% apply(2,sens.slope) # Se obtienen las pendientes de sen de cada replica
    T_S_MC_slope[i-1] <- data.frame(Sen_slope = sapply(T_S_MC, function(x){as.numeric(x[1])}))
    WEI_T_S[[i-1]] <- lmomco::parwei(lmomco::lmoms(unlist(T_S_MC_slope[i-1]))) # Se ajusta una Weibull para las pendientes de las replicas
    Obs_sen_slope[[i-1]] <- sens.slope(v)
    MC_trend1[i-1,2]  <- cdfwei(unlist(Obs_sen_slope[[i-1]][1]),lmomco::parwei(lmomco::lmoms(unlist(T_S_MC_slope[i-1]))))
    PP_T_S <- lmomco::pp(unlist(T_S_MC_slope[i-1])) # tienen todos la misma probb
    Quawei_T_S[[i-1]] <- lmomco::quawei(PP_T_S, WEI_T_S[[i-1]])
  }

  ## To data frame and sort
  Probabilidad = data.frame(lmomco::pp(unlist(T_S_MC_slope[1]))) # tienen todos la misma probb
  colnames(Probabilidad) <- 'Probabilidad'

  T_S_MC_slope <- as.data.frame(do.call(cbind, T_S_MC_slope))
  colnames(T_S_MC_slope) <- colnames(x)[-1]
  CDF_Empirica <- T_S_MC_slope %>% apply(2, sort, decreasing=F) %>% as.data.frame()

  CDF_Ajustada <- as.data.frame(do.call(cbind, Quawei_T_S))
  colnames(CDF_Ajustada) <- colnames(x)[-1]


  # Plots
  long_format_CDF_emp <- reshape2::melt(CDF_Empirica, variable.name = "Model", value.name = 'CDF_Empirica')
  long_format_CDF_ajus <- reshape2::melt(CDF_Ajustada, variable.name = "Model", value.name = 'CDF_Ajustada')

  long_format_MC_trend <- matrix(replicate(dim(x)[2]-1,t(Probabilidad)),ncol=1)
  long_format_MC_trend <- cbind(Probabilidad = long_format_MC_trend, long_format_CDF_emp,long_format_CDF_ajus)
  long_format_MC_trend %<>% select(-4) %>% relocate(Probabilidad, .after = Model)

  Obs_sen_slope_Aux = sapply(Obs_sen_slope, function(x){as.numeric(x[1])})
  long_format_Obs <- cbind(MC_trend1,X_obs = Obs_sen_slope_Aux)
  colnames(long_format_Obs)[2] <- 'Y_obs'
  long_format_Obs$Y_obs <- as.double(long_format_Obs$Y_obs)

  long_format_Obs %<>% mutate(eq_y = paste('Yobs == ',round(Y_obs,2)))
  long_format_Obs %<>% mutate(eq_x = paste('Xobs == ',round(X_obs,2)))

  labels <- c("Ajustada", "Empirica")
  names(labels) <- c("Model", "Real")

  plot1 <- ggplot(long_format_MC_trend, aes(x=CDF_Empirica, y=Probabilidad, group = 1)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = (1 - int), ymax = Inf),
              fill = "pink", alpha = 0.01) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = int),
              fill = "pink", alpha = 0.01) +

    geom_line(aes(color = Model),linetype = 'solid', linewidth = 2)+
    geom_line(aes(x=CDF_Ajustada, y=Probabilidad),color='black',linetype = 'dashed') +
    geom_point(long_format_Obs, mapping=aes(x=X_obs, y=Y_obs),
               color='blue', shape=0, size=2.5, stroke = 1) +

    geom_text(data=long_format_Obs,
              aes(x = -1.0,y = 0.8,
                  label = eq_x),
              size = 3, parse = TRUE, inherit.aes=FALSE, color = 'blue') +
    geom_text(data=long_format_Obs,
              aes(x = -1.0,y = 0.65,
                  label = eq_y),
              size = 3, parse = TRUE, inherit.aes=FALSE, color = 'blue') +

    facet_wrap(~Model) +
    theme_light() +
    labs(title = plot_title,
         subtitle='Monte-Carlo trend test',
         caption = 'Dashed line corresponds to real data, solid color line corresponds to the model',
         y = 'F(x)', x = 'Slope [mm/yr]') +
    theme(legend.position="none")

  # Se entregan resultados
  if(opt == 1) {
    result_DF <- data.frame(MC_trend1, Obs_sen_slope_Aux)
    result_DF %<>% mutate(Test = ifelse((result_DF$MC_trend < int) | (result_DF$MC_trend > (1 - int)), 'Trend', 'No trend'))
    colnames(result_DF)[2:3] <- c("Y_obs", "X_obs")
    return(result_DF)
  } else {
    return(plot1)
  }

}
