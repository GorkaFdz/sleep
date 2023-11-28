# Índices de valores atípicos con el paquete anomalyze
anomalia <- function(tiempo, valor){
  ts <- tibble(tiempo, valor)
  ts_anomalia <- ts%>%
    time_decompose(valor)%>%
    anomalize(remainder)
  index <- which(ts_anomalia$anomaly == "Yes")
  return (unname(index))
}

# Índices de valores atípicos fuera de un mínimo y máximo
anomalia_valor <- function(datos, min, max){
  index <- which(datos >= max | datos <= min)
  return (index)
}

# Índices de valores atípicos con percentiles
anomalia_percentil <- function(datos, min, max){
  inf <- quantile(datos, min)
  sup <- quantile(datos, max)
  index <- which(datos <= inf | datos >= sup)
  return (index)
}

# Cálculo de la mejor combinación de nudos para el modelo de splines (devuelve el modelo)
cv <- function(x, y){
  # crear la lista de los tipos de nodos
  nodos <- list()
  for (i in 3:5) {
    knots <- seq(min(x), max(x), len = i + 2)[-c(1, i + 2)]
    nodos[i-2] <- list(knots)
  }
  knots3 <- unname(quantile(x, c(0.1, 0.5, 0.9)))
  knots4 <- unname(quantile(x, c(0.05, 0.35, 0.65, 0.95)))
  knots5 <- unname(quantile(x, c(0.05, 0.275, 0.5, 0.725, 0.95)))
  nodos[4] <- list(knots3)
  nodos[5] <- list(knots4)
  nodos[6] <- list(knots5)
  
  # particiones de los datos para entrenar el modelo
  cv <- createDataPartition(y, times = 10, p = 0.8)
  
  mejor <- Inf
  
  datos <- data.frame(x = x, y = y)
  # entrenar modelos con los distintos nodos
  for (nodo in nodos){
    media <- c()
    # usar todas las muestras
    for (sample in cv){
      test <- datos[-unique(sample),]
      train <- datos[unique(sample),]
      independiente <- train$x
      model <- lm(train$y ~ bs(independiente, knots = nodo, degree = 3))
      suppressWarnings({
        predictions <- predict(model, data.frame(independiente = test$x))
      })
      RMSE <- RMSE(predictions, test$y)
      media[length(media)+1] <- RMSE 
    }
    # calcular la media de los errores
    if (mean(media) < mejor){
      mejor <- mean(media)
      modelo <- model
    }
  }
  return (modelo)
}

# Cálcula valor con el que se alcanza el máximo del modelo devuelto por cv()
optimo <- function(modelo){
  piece <- RegSplineAsPiecePoly(modelo, "bs(independiente, knots = nodo, degree = 3)")
  
  raiz <- solve(piece, deriv = 1)
  extremos <- predict(piece, raiz)
  
  max <- which.max(extremos)
  return (raiz[max])
}

# Valor maximo a partir de un dataframe
maximo <- function(df){
  bootstrap <- bootstraps(df, times = 50)
  media <- c()
  for (sample in bootstrap$splits){
    data <- as.data.frame(sample)
    # cambiar nombre de las columnas
    modelo <- train(accel ~ times, data = data, method = "gam", trControl = trainControl(method = "cv", number = 10))
    newx <- seq(min(data$times), max(data$times), len = 200)
    predicciones <- predict(modelo, newdata = data.frame(times = newx))
    max <- newx[which.max(predicciones)]
    media[length(media)+1] <- max
  }
  return (mean(media))
}