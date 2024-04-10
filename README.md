# TS_Analisis
Script que contiene funciones para el análisis de series de tiempo y ajuste de un modelo arima.

La función ts_analyze() recibe tres parámetros: ts (serie de tiempo), type (descomposición de la serie "additive" o "multiplicative") y lag_max (cantidad de lags para el acf y pacf).

La función arima_fit() recibe a la serie de tiempo, y parámetros relacionados con las funciones auto.arima(), Arima() y forecast().
