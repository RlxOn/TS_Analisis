#################### Series de tiempo ARIMA ####################################
pacman::p_load(tseries, forecast, timeSeries, astsa, lmtest, nortest, car, ggplot2,
               tidyverse, lubridate, gridExtra)

set.seed(1989)  #Taylor's Version

#Funcion para graficar serie, descomposicion y correlaciones ===================
ts_analyze <- function(ts, type = "additive", lag_max = NULL){
  #Inicializamos los parametros
  ts <- ts
  type <- type
  lag_max <- lag_max
  
  #Descomponemos la serie en tendencia, ciclos y comp aleatoria
  descomposicion <- decompose(ts, type = type)
    
  descomp_df <- data.frame(
    tiempo = getTime(descomposicion$x),
    serie = descomposicion$x,
    tendencia = descomposicion$trend,
    ciclos = descomposicion$seasonal,
    aleatoria = descomposicion$random)
  
  #Guardamos los plots
  gg_serie <- ggplot(data = descomp_df, aes(x = tiempo, y = serie)) + 
    geom_line(col = "#e30b5d", linewidth = 1.5) + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + ggtitle("Serie original")
  
  gg_tendencia <- ggplot(data = descomp_df, aes(x = tiempo, y = tendencia)) + 
    geom_line(col = "#0b5de3", linewidth = 1.5) + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + ggtitle("Tendencia")
  
  gg_ciclos <- ggplot(data = descomp_df, aes(x = tiempo, y = ciclos)) + 
    geom_line(col = "#0b5de3", linewidth = 1.5) + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + ggtitle("Ciclos") 
  
  gg_aleatoria <- ggplot(data = descomp_df, aes(x = tiempo, y = aleatoria)) + 
    geom_line(col = "#e30b5d", linewidth = 1.5) + xlab(NULL) + ylab(NULL) + 
    theme_minimal() + ggtitle("Componente aleatoria") + 
    geom_hline(yintercept = 0, col = "red", linewidth = 1, linetype = "dashed")
  
  print(
  grid.arrange(gg_serie, gg_tendencia, gg_ciclos, gg_aleatoria, ncol = 2))
  
  #Hacemos el ACF y el PACF
  confianza <- 1.96 / sqrt(length(ts))
  
  #ACF
  acf_ts <- acf(ts, plot = FALSE, lag.max = lag_max)
  acf_df <- data.frame(lag = acf_ts$lag, acf = acf_ts$acf)
  
  gg_acf <- ggplot(acf_df, aes(x = lag, y = acf)) + geom_bar(stat = "identity", fill = "#0b5de3") +
    geom_hline(yintercept = c(-confianza, confianza), linetype = "dashed", color = "red",
               linewidth = 1) +
    labs(title = "ACF", x = "Lag", y = NULL) + theme_minimal() + 
    geom_hline(yintercept = 0, col = "black")
  
  #PACF
  pacf_ts <- pacf(ts, plot = FALSE, lag.max = lag_max)
  pacf_df <- data.frame(lag = pacf_ts$lag, pacf = pacf_ts$acf)
  
  gg_pacf <- ggplot(pacf_df, aes(x = lag, y = pacf)) + geom_bar(stat = "identity", fill = "#0be391") +
    geom_hline(yintercept = c(-confianza, confianza), linetype = "dashed", color = "red",
               linewidth = 1) +
    labs(title = "PACF", x = "Lag", y = NULL) + theme_minimal() + 
    geom_hline(yintercept = 0, col = "black")
  
  print(
    grid.arrange(gg_acf, gg_pacf, ncol = 1)
  )
  
  #Pruebas de estacionariedad
  dfuller <- adf.test(ts, alternative = "stationary")$p.value
  kpss <- kpss.test(ts)$p.value
  pp <- pp.test(ts, alternative = "stationary")$p.value
  
  Pruebas_df <- data.frame(
  Pruebas = c("Dickey-Fuller", "KPSS", "Phillips-Perron"),
  P_values = c(dfuller, kpss, pp),
  H_Alt = c("Estacionaria", "No estacionaria", "Estacionaria")
  )
  
  print(Pruebas_df)
  
  #Lista de retorno
  output_list <- list(descomp_df, acf_ts, pacf_ts, Pruebas_df) 
  names(output_list) <- c("Descomposicion", "ACF", "PACF", "Pruebas")
  
  return(output_list)
}




#Funcion para ajustar un modelo ARIMA ==========================================
arima_fit <- function(ts, seasonal = TRUE, stationary = FALSE, auto = TRUE,
                      order_coef = NULL, seasonal_coef = NULL, drift = TRUE, mean = TRUE, h = 12){

  #Inicializar parametros
  ts <- ts
  seasonal <- seasonal
  starionary <- stationary
  auto <- auto
  order_coef <- order_coef
  seasonal_coef <- seasonal_coef
  drift <- drift
  mean <- mean
  h <- h
  
  if(auto == TRUE){
    fit_arima <- auto.arima(ts, stationary = stationary, seasonal = seasonal,
                            allowdrift = drift, allowmean = mean)
    }else{
    fit_arima <- Arima(ts, order = order_coef, seasonal = seasonal_coef, 
                       include.mean = mean, include.drift = drift)
    }
    
    #Coeficientes
    coef_tab <- broom::tidy(fit_arima)
    conf <- confint(fit_arima)
    
    coef_tab <- coef_tab %>%
      mutate(Lower = conf[,1],
             Upper = conf[,2],
             Is_Sign = ifelse(sign(Lower) == sign(Upper), TRUE, FALSE))
    
    #Performance y bondad de ajuste
    metricas <- broom::glance(fit_arima)
    
    #Supuestos de los residuales
    resd <- fit_arima$residuals
    
    ##Normalidad
    jbt <- jarque.bera.test(resd)$p.value
    ad <- ad.test(resd)$p.value
    
    ##Media cero
    ttest <- t.test(resd, mu = 0)$p.value
    wtest <- wilcox.test(resd, mu = 0)$p.value
    
    ##Homocedasticidad
    aux <- 1:length(resd)
    bp <- bptest(resd ~ aux)$p.value
    
    pruebas <- data.frame(
      Pruebas = c("Jarque.Bera", "Anderson-Darling", "Prueta t",
                  "Wilcoxon", "Breusch-Pagan"),
      p_values = c(jbt, ad, ttest, wtest, bp)
    )
    
    print(pruebas)
    
    ##Correlacion nula
    LB <- c()
    BP <- c()
    
    for(i in 1:12){
      x <- Box.test(resd, type = "Ljung-Box", lag = i)$p.value
      y <- Box.test(resd, type = "Box-Pierce", lag = i)$p.value
      
      LB <- append(LB, x)
      BP <- append(BP, y)
    }
    
    
    #Graficos necesarios
    ##Observados vs Ajustados
    aj <- fit_arima$fitted
    
    obs_plot <- ggplot() + 
      geom_line(data = NULL, aes(x = getTime(ts), y = ts, color = "Observados"), linewidth = 1.2) + 
      geom_line(data = NULL, aes(x = getTime(aj), y = aj, color = "Ajustados"), linewidth = 1.2) + 
      labs(x = NULL, y = NULL, title = "Observados vs Ajustados", color = NULL) + 
      scale_color_manual(values = c("Observados" = "#e30b5d", "Ajustados" = "navy")) + 
      theme_minimal() + theme(legend.position = "bottom")
    
    ##LjungBox y BoxPierce p-values
    corr_plot <- ggplot() + geom_point(data = NULL, aes(x = 1:12, y = LB, color = "Ljung-Box"), size = 3) + 
      geom_point(data = NULL, aes(x = 1:12, y = BP, color = "Box-Pierce"), size = 3) +
      geom_hline(yintercept = 0.05, col = "red", linetype = "dashed") + 
      labs(x = NULL, y = NULL, title = "Correlacion nula", subtitle = "p-values", color = NULL) + 
      scale_color_manual(values = c("Ljung-Box" = "#e30b5d", "Box-Pierce" = "#0b5de3")) + 
      theme_minimal() + theme(legend.position = "bottom") + 
      scale_x_continuous(n.breaks = 12)
    
    ##Serie residuales
    res_plot <- ggplot() + geom_line(data = NULL, aes(x = getTime(resd), y = resd), col = "navy", linewidth = 1.2) + 
      geom_hline(yintercept = mean(resd), col = "red", linetype = "dashed", linewidth = 1) + 
      theme_minimal() + labs(x = NULL, y = NULL, title = "Residuales")
    
    ##Normalidad residuales
    dens_plot <- ggplot(data = NULL, aes(x = resd)) + geom_density(fill = "#000059", col = "#000059") + 
      theme_minimal() + labs(title = "Densidad Residuales", x = NULL, y = NULL)
    
    print(
      grid.arrange(obs_plot, corr_plot, res_plot, dens_plot, ncol = 2)
    )
    
    #Prediccion
    pred <- forecast::forecast(fit_arima, h = h, level = 95)
    
    print(
    ggplot() + geom_line(data = NULL, aes(x = getTime(ts), y = ts, color = "Serie"), linewidth = 1.2) +
      geom_line(data = NULL, aes(x = getTime(pred$mean), y = pred$mean, color = "Prediccion"), linewidth = 1.2) + 
      geom_ribbon(data = NULL, aes(x = getTime(pred$mean), ymin = pred$lower,
                                   ymax = pred$upper), alpha = 0.6, fill = "light blue") + 
      labs(title = "Prediccion", x = NULL, y = NULL, color = NULL) + 
      scale_color_manual(values = c("Serie" = "#000059", "Prediccion" = "#0a0aff")) + 
      theme_minimal() + theme(legend.position = "bottom")
    )
    
    #Lista de salida
    output_list <- list(fit_arima, coef_tab, metricas, pruebas, pred)
    
    names(output_list) <- c("Modelo", "Coeficientes", "Metricas", "Pruebas",
                            "Prediccion")
    
    return(output_list)
}
