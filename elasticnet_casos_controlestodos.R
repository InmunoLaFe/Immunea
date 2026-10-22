# ============================================================
# ELASTIC NET PARA SELECCIÓN DE VARIABLES
# Outcome: Factor_caso_1_EA
#   1 = caso EA
#   0 = no caso
# ============================================================

# ------------------------------------------------------------
# 1. Instalar y cargar paquetes
# ------------------------------------------------------------

paquetes <- c(
  "readxl",
  "dplyr",
  "glmnet",
  "pROC",
  "ggplot2"
)

paquetes_nuevos <- paquetes[
  !paquetes %in% rownames(installed.packages())
]

if (length(paquetes_nuevos) > 0) {
  install.packages(paquetes_nuevos)
}

library(readxl)
library(dplyr)
library(glmnet)
library(pROC)
library(ggplot2)


# ------------------------------------------------------------
# 2. Importar la base de datos
# ------------------------------------------------------------

ruta_archivo <- paste0(
  "C:/Users/victo/OneDrive/Desktop/Research/INMUNEA/",
  "Estadistica Victor/Clean_INMUNEA_ALZHEIMER_DB.xlsx"
)

datos_originales <- read_excel(
  path = ruta_archivo,
  sheet = 1,
  .name_repair = "minimal"
)

dim(datos_originales)
names(datos_originales)


# ------------------------------------------------------------
# 3. Definir outcome y variables predictoras
# ------------------------------------------------------------

variable_resultado <- "Factor_caso_1_EA"

variables_predictoras <- c(
  "postCG_unsw",
  "postCG_sw",
  "CD4+_Naive",
  "CD4+_Mem_central",
  "CD4+_Mem_efectora",
  "CD4+_efectora_terminal",
  "CD8_totales",
  "CD8+_naive",
  "CD8+_Mem_central",
  "CD8+_Mem_efectora_ CD27+dim",
  "CD8+_Efectora",
  "CD8+_Efectora_terminal",
  "doble_positivas",
  "doble_negativas_no_γ/δ",
  "doble_negativas_γ/δ",
  "Tregs",
  "Treg_memoria_efectora",
  "Treg_Th1-like",
  "Treg_Th2-like",
  "Treg_Th17-like",
  "Treg_Th22-like",
  "Th_reguladoras_foliculares",
  "Th_convencionales",
  "Otros_Th_convencionales",
  "Th9",
  "Th1",
  "Th2",
  "Th17",
  "Th22",
  "Th_foliculares",
  "Th_foliculares_otros",
  "Tfh1",
  "Tfh2",
  "Tfh17",
  "Celulas_NK",
  "Células_NKT",
  "age",
  "sex"
)


# ------------------------------------------------------------
# 4. Comprobar que todas las variables existen
# ------------------------------------------------------------

variables_necesarias <- c(
  variable_resultado,
  variables_predictoras
)

variables_no_encontradas <- setdiff(
  variables_necesarias,
  names(datos_originales)
)

if (length(variables_no_encontradas) > 0) {
  stop(
    paste0(
      "Estas variables no se han encontrado en la base de datos:\n",
      paste(variables_no_encontradas, collapse = "\n")
    )
  )
}

cat("Todas las variables se han localizado correctamente.\n")


# ------------------------------------------------------------
# 5. Seleccionar las variables del análisis
# ------------------------------------------------------------

datos_elastic <- datos_originales[, variables_necesarias]


# ------------------------------------------------------------
# 6. Convertir variables a formato numérico
# ------------------------------------------------------------

# Outcome binario 0/1
datos_elastic[[variable_resultado]] <- as.numeric(
  as.character(datos_elastic[[variable_resultado]])
)

# Sex:
# 0 = hombre
# 1 = mujer
datos_elastic[["sex"]] <- as.numeric(
  as.character(datos_elastic[["sex"]])
)

# Convertir el resto de predictores cuantitativos a numérico
variables_continuas <- setdiff(
  variables_predictoras,
  "sex"
)

datos_elastic[variables_continuas] <- lapply(
  datos_elastic[variables_continuas],
  function(x) {
    as.numeric(as.character(x))
  }
)


# ------------------------------------------------------------
# 7. Comprobaciones básicas
# ------------------------------------------------------------

# Distribución del outcome
print(
  table(
    datos_elastic[[variable_resultado]],
    useNA = "ifany"
  )
)

# Distribución por sexo
print(
  table(
    datos_elastic[["sex"]],
    useNA = "ifany"
  )
)

# Comprobar valores válidos del outcome
valores_outcome <- unique(
  na.omit(datos_elastic[[variable_resultado]])
)

if (!all(valores_outcome %in% c(0, 1))) {
  stop(
    "Factor_caso_1_EA contiene valores diferentes de 0 y 1."
  )
}

# Comprobar valores válidos de sex
valores_sex <- unique(
  na.omit(datos_elastic[["sex"]])
)

if (!all(valores_sex %in% c(0, 1))) {
  stop(
    "La variable sex contiene valores diferentes de 0 y 1."
  )
}


# ------------------------------------------------------------
# 8. Examinar datos faltantes
# ------------------------------------------------------------

missing_por_variable <- data.frame(
  variable = names(datos_elastic),
  n_missing = sapply(
    datos_elastic,
    function(x) sum(is.na(x))
  ),
  porcentaje_missing = round(
    100 * sapply(
      datos_elastic,
      function(x) mean(is.na(x))
    ),
    2
  )
)

missing_por_variable <- missing_por_variable %>%
  arrange(desc(n_missing))

print(missing_por_variable)


# ------------------------------------------------------------
# 9. Análisis por casos completos
# ------------------------------------------------------------

n_inicial <- nrow(datos_elastic)

datos_completos <- datos_elastic %>%
  filter(complete.cases(.))

n_final <- nrow(datos_completos)

cat("\nNúmero inicial de participantes:", n_inicial, "\n")
cat("Número con datos completos:", n_final, "\n")
cat("Participantes excluidos por missing:", n_inicial - n_final, "\n")

# Distribución del outcome después de excluir missing
print(
  table(datos_completos[[variable_resultado]])
)


# ------------------------------------------------------------
# 10. Eliminar predictores sin variabilidad
# ------------------------------------------------------------

numero_valores_unicos <- sapply(
  datos_completos[variables_predictoras],
  function(x) length(unique(x))
)

variables_sin_variabilidad <- names(
  numero_valores_unicos[numero_valores_unicos <= 1]
)

if (length(variables_sin_variabilidad) > 0) {
  cat(
    "\nVariables eliminadas por no tener variabilidad:\n",
    paste(variables_sin_variabilidad, collapse = "\n"),
    "\n"
  )
  
  variables_predictoras_finales <- setdiff(
    variables_predictoras,
    variables_sin_variabilidad
  )
} else {
  variables_predictoras_finales <- variables_predictoras
}


# ------------------------------------------------------------
# 11. Crear matriz X y vector y
# ------------------------------------------------------------

# Extraer directamente los predictores como matriz numérica.
# Esto evita problemas con nombres que contienen +, espacios,
# tildes, guiones o caracteres griegos.

X <- datos_completos[, variables_predictoras_finales, drop = FALSE]

# Asegurar que todas las columnas sean numéricas
X[] <- lapply(
  X,
  function(x) as.numeric(as.character(x))
)

# Convertir a matriz, formato requerido por glmnet
X <- as.matrix(X)

# Vector de resultado binario
y <- as.numeric(
  datos_completos[[variable_resultado]]
)

# Comprobaciones
dim(X)
length(y)
colnames(X)

# Confirmar que la matriz sea numérica
str(X)
storage.mode(X)

# ------------------------------------------------------------
# 12. Configuración de validación cruzada
# ------------------------------------------------------------

set.seed(20260722)

# Con muestras pequeñas no conviene utilizar más folds
# que sujetos haya en la categoría minoritaria.

n_clase_minoritaria <- min(table(y))

numero_folds <- min(
  10,
  n_clase_minoritaria
)

if (numero_folds < 3) {
  stop(
    paste0(
      "Hay muy pocos sujetos en una de las categorías para ",
      "realizar validación cruzada: n = ",
      n_clase_minoritaria
    )
  )
}

cat(
  "\nNúmero de folds empleado:",
  numero_folds,
  "\n"
)

# Se utiliza la misma partición para todos los valores de alpha,
# de modo que puedan compararse de manera homogénea.

set.seed(20260722)

numero_folds <- 5

foldid <- integer(length(y))

foldid[y == 0] <- sample(
  rep(
    seq_len(numero_folds),
    length.out = sum(y == 0)
  )
)

foldid[y == 1] <- sample(
  rep(
    seq_len(numero_folds),
    length.out = sum(y == 1)
  )
)

# Comprobar distribución de clases por fold
table(foldid, y)


# ------------------------------------------------------------
# 13. Optimización del parámetro alpha
# ------------------------------------------------------------

# alpha = 0: Ridge
# alpha = 1: LASSO
# 0 < alpha < 1: Elastic Net

grid_alpha <- seq(
  from = 0,
  to = 1,
  by = 0.1
)

lista_modelos <- vector(
  mode = "list",
  length = length(grid_alpha)
)

resultados_alpha <- data.frame(
  alpha = grid_alpha,
  lambda_min = NA_real_,
  lambda_1se = NA_real_,
  deviance_cv_minima = NA_real_,
  n_variables_lambda_min = NA_integer_,
  n_variables_lambda_1se = NA_integer_
)

for (i in seq_along(grid_alpha)) {
  
  alpha_actual <- grid_alpha[i]
  
  modelo_cv <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = alpha_actual,
    foldid = foldid,
    type.measure = "deviance",
    standardize = TRUE,
    intercept = TRUE,
    nlambda = 200
  )
  
  lista_modelos[[i]] <- modelo_cv
  
  indice_lambda_min <- which.min(modelo_cv$cvm)
  
  coef_min <- coef(
    modelo_cv,
    s = "lambda.min"
  )
  
  coef_1se <- coef(
    modelo_cv,
    s = "lambda.1se"
  )
  
  resultados_alpha$lambda_min[i] <- modelo_cv$lambda.min
  resultados_alpha$lambda_1se[i] <- modelo_cv$lambda.1se
  
  resultados_alpha$deviance_cv_minima[i] <-
    modelo_cv$cvm[indice_lambda_min]
  
  resultados_alpha$n_variables_lambda_min[i] <-
    sum(as.vector(coef_min)[-1] != 0)
  
  resultados_alpha$n_variables_lambda_1se[i] <-
    sum(as.vector(coef_1se)[-1] != 0)
}

print(resultados_alpha)


# ------------------------------------------------------------
# 14. Seleccionar el mejor alpha
# ------------------------------------------------------------

indice_mejor_alpha <- which.min(
  resultados_alpha$deviance_cv_minima
)

mejor_alpha <- resultados_alpha$alpha[
  indice_mejor_alpha
]

modelo_final_cv <- lista_modelos[[indice_mejor_alpha]]

cat("\nMejor alpha:", mejor_alpha, "\n")
cat(
  "Lambda mínimo:",
  modelo_final_cv$lambda.min,
  "\n"
)
cat(
  "Lambda 1-SE:",
  modelo_final_cv$lambda.1se,
  "\n"
)


# ------------------------------------------------------------
# 15. Gráfico de rendimiento según alpha
# ------------------------------------------------------------

grafico_alpha <- ggplot(
  resultados_alpha,
  aes(
    x = alpha,
    y = deviance_cv_minima
  )
) +
  geom_line() +
  geom_point(size = 2) +
  geom_vline(
    xintercept = mejor_alpha,
    linetype = "dashed"
  ) +
  labs(
    title = "Selección del parámetro alpha",
    subtitle = paste0(
      "Alpha óptimo = ",
      mejor_alpha
    ),
    x = "Alpha",
    y = "Deviance binomial media en validación cruzada"
  ) +
  theme_classic()

print(grafico_alpha)


# ------------------------------------------------------------
# 16. Gráfico de validación cruzada del modelo seleccionado
# ------------------------------------------------------------

plot(modelo_final_cv)

title(
  paste0(
    "Elastic Net: alpha = ",
    mejor_alpha
  )
)


# ------------------------------------------------------------
# 17. Extraer coeficientes con lambda.min
# ------------------------------------------------------------

coeficientes_lambda_min <- as.matrix(
  coef(
    modelo_final_cv,
    s = "lambda.min"
  )
)

tabla_coeficientes_min <- data.frame(
  variable = rownames(coeficientes_lambda_min),
  coeficiente = as.numeric(coeficientes_lambda_min[, 1]),
  row.names = NULL
) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    seleccionada = coeficiente != 0,
    OR_penalizada = exp(coeficiente)
  ) %>%
  arrange(desc(abs(coeficiente)))

cat("\nCoeficientes con lambda.min:\n")
print(tabla_coeficientes_min)


# ------------------------------------------------------------
# 18. Variables seleccionadas con lambda.min
# ------------------------------------------------------------

variables_seleccionadas_min <- tabla_coeficientes_min %>%
  filter(seleccionada) %>%
  select(
    variable,
    coeficiente,
    OR_penalizada
  )

cat("\nVariables seleccionadas con lambda.min:\n")
print(variables_seleccionadas_min)


# ------------------------------------------------------------
# 19. Extraer coeficientes con lambda.1se
# ------------------------------------------------------------

coeficientes_lambda_1se <- as.matrix(
  coef(
    modelo_final_cv,
    s = "lambda.1se"
  )
)

tabla_coeficientes_1se <- data.frame(
  variable = rownames(coeficientes_lambda_1se),
  coeficiente = as.numeric(coeficientes_lambda_1se[, 1]),
  row.names = NULL
) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    seleccionada = coeficiente != 0,
    OR_penalizada = exp(coeficiente)
  ) %>%
  arrange(desc(abs(coeficiente)))

cat("\nCoeficientes con lambda.1se:\n")
print(tabla_coeficientes_1se)


# ------------------------------------------------------------
# 20. Variables seleccionadas con lambda.1se
# ------------------------------------------------------------

variables_seleccionadas_1se <- tabla_coeficientes_1se %>%
  filter(seleccionada) %>%
  select(
    variable,
    coeficiente,
    OR_penalizada
  )

cat("\nVariables seleccionadas con lambda.1se:\n")
print(variables_seleccionadas_1se)


# ------------------------------------------------------------
# 21. Predicciones obtenidas mediante validación cruzada
# ------------------------------------------------------------

# Se repite el modelo óptimo conservando las predicciones
# out-of-fold de la validación cruzada.

set.seed(20260722)

modelo_final_predicciones <- cv.glmnet(
  x = X,
  y = y,
  family = "binomial",
  alpha = mejor_alpha,
  foldid = foldid,
  type.measure = "deviance",
  standardize = TRUE,
  intercept = TRUE,
  nlambda = 200,
  keep = TRUE
)

# Identificar la columna correspondiente a lambda.min
indice_lambda_pred <- which.min(
  abs(
    modelo_final_predicciones$lambda -
      modelo_final_predicciones$lambda.min
  )
)

probabilidades_cv <- modelo_final_predicciones$fit.preval[
  ,
  indice_lambda_pred
]

# Eliminar posibles valores no disponibles
indices_validos <- complete.cases(
  probabilidades_cv,
  y
)

roc_cv <- roc(
  response = y[indices_validos],
  predictor = probabilidades_cv[indices_validos],
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc_cv <- as.numeric(
  auc(roc_cv)
)

cat(
  "\nAUC basada en predicciones out-of-fold:",
  round(auc_cv, 3),
  "\n"
)

plot(
  roc_cv,
  main = paste0(
    "ROC Elastic Net mediante validación cruzada\n",
    "AUC = ",
    round(auc_cv, 3)
  )
)


# ------------------------------------------------------------
# 22. Guardar los resultados en archivos CSV
# ------------------------------------------------------------

directorio_resultados <- paste0(
  "C:/Users/victo/OneDrive/Desktop/Research/INMUNEA/",
  "Estadistica Victor/Resultados_Elastic_Net"
)

if (!dir.exists(directorio_resultados)) {
  dir.create(
    directorio_resultados,
    recursive = TRUE
  )
}

write.csv(
  resultados_alpha,
  file = file.path(
    directorio_resultados,
    "comparacion_alpha.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tabla_coeficientes_min,
  file = file.path(
    directorio_resultados,
    "coeficientes_lambda_min.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  variables_seleccionadas_min,
  file = file.path(
    directorio_resultados,
    "variables_seleccionadas_lambda_min.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tabla_coeficientes_1se,
  file = file.path(
    directorio_resultados,
    "coeficientes_lambda_1se.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  variables_seleccionadas_1se,
  file = file.path(
    directorio_resultados,
    "variables_seleccionadas_lambda_1se.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nResultados guardados en:\n",
  directorio_resultados,
  "\n"
)