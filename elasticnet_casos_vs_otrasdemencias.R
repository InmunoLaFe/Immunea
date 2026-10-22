# ============================================================
# ELASTIC NET PARA SELECCIÓN DE VARIABLES
#
# Comparación:
#   grupo original = 0  -> outcome_modelo = 0
#   grupo original = 2  -> outcome_modelo = 1
#
# Se excluye:
#   grupo original = 1
#
# Interpretación:
#   coeficiente positivo -> mayor probabilidad de grupo 2
#   coeficiente negativo -> mayor probabilidad de grupo 0
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

cat(
  "Dimensiones de la base original:",
  nrow(datos_originales),
  "filas y",
  ncol(datos_originales),
  "columnas\n"
)


# ------------------------------------------------------------
# 3. Definir outcome y variables predictoras
# ------------------------------------------------------------

variable_grupo_original <- "grupo"
variable_resultado <- "grupo_2_vs_0"

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

variables_necesarias_originales <- c(
  variable_grupo_original,
  variables_predictoras
)

variables_no_encontradas <- setdiff(
  variables_necesarias_originales,
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
# 5. Seleccionar grupos 0 y 2 y recodificar el outcome
# ------------------------------------------------------------

datos_elastic <- datos_originales[
  ,
  variables_necesarias_originales
]

datos_elastic[[variable_grupo_original]] <- as.numeric(
  as.character(datos_elastic[[variable_grupo_original]])
)

cat("\nDistribución original de la variable grupo:\n")

print(
  table(
    datos_elastic[[variable_grupo_original]],
    useNA = "ifany"
  )
)

# Mantener únicamente grupos 0 y 2
datos_elastic <- datos_elastic %>%
  filter(
    !is.na(.data[[variable_grupo_original]]),
    .data[[variable_grupo_original]] %in% c(0, 2)
  )

cat("\nDistribución tras conservar únicamente grupos 0 y 2:\n")

print(
  table(
    datos_elastic[[variable_grupo_original]],
    useNA = "ifany"
  )
)

# Recodificación:
# grupo original 0 -> 0
# grupo original 2 -> 1

datos_elastic[[variable_resultado]] <- ifelse(
  datos_elastic[[variable_grupo_original]] == 2,
  1,
  0
)

cat("\nDistribución del outcome recodificado:\n")
cat("0 = grupo original 0\n")
cat("1 = grupo original 2\n")

print(
  table(
    datos_elastic[[variable_resultado]],
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 6. Convertir predictores a formato numérico
# ------------------------------------------------------------

# Sexo:
# 0 = hombre
# 1 = mujer

datos_elastic[["sex"]] <- as.numeric(
  as.character(datos_elastic[["sex"]])
)

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

valores_outcome <- unique(
  na.omit(datos_elastic[[variable_resultado]])
)

if (!all(valores_outcome %in% c(0, 1))) {
  stop(
    "El outcome recodificado contiene valores diferentes de 0 y 1."
  )
}

if (length(valores_outcome) != 2) {
  stop(
    "Después del filtrado no están presentes ambos grupos."
  )
}

valores_sex <- unique(
  na.omit(datos_elastic[["sex"]])
)

if (!all(valores_sex %in% c(0, 1))) {
  stop(
    "La variable sex contiene valores diferentes de 0 y 1."
  )
}

cat("\nDistribución por sexo:\n")

print(
  table(
    datos_elastic[["sex"]],
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# 8. Examinar datos faltantes
# ------------------------------------------------------------

variables_analisis <- c(
  variable_resultado,
  variables_predictoras
)

datos_analisis <- datos_elastic[
  ,
  variables_analisis
]

missing_por_variable <- data.frame(
  variable = names(datos_analisis),
  
  n_missing = sapply(
    datos_analisis,
    function(x) sum(is.na(x))
  ),
  
  porcentaje_missing = round(
    100 * sapply(
      datos_analisis,
      function(x) mean(is.na(x))
    ),
    2
  ),
  
  row.names = NULL
)

missing_por_variable <- missing_por_variable %>%
  arrange(desc(n_missing))

cat("\nDatos ausentes por variable:\n")

print(missing_por_variable)


# ------------------------------------------------------------
# 9. Análisis por casos completos
# ------------------------------------------------------------

n_inicial <- nrow(datos_analisis)

datos_completos <- datos_analisis %>%
  filter(complete.cases(.))

n_final <- nrow(datos_completos)

cat(
  "\nParticipantes pertenecientes a grupos 0 o 2:",
  n_inicial,
  "\n"
)

cat(
  "Participantes con datos completos:",
  n_final,
  "\n"
)

cat(
  "Excluidos por datos ausentes:",
  n_inicial - n_final,
  "\n"
)

cat("\nDistribución final del outcome:\n")
cat("0 = grupo original 0\n")
cat("1 = grupo original 2\n")

print(
  table(datos_completos[[variable_resultado]])
)

if (length(unique(datos_completos[[variable_resultado]])) != 2) {
  stop(
    paste0(
      "Tras excluir los datos ausentes no permanecen ",
      "ambas categorías del outcome."
    )
  )
}


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
    "\nVariables eliminadas por falta de variabilidad:\n",
    paste(variables_sin_variabilidad, collapse = "\n"),
    "\n"
  )
  
  variables_predictoras_finales <- setdiff(
    variables_predictoras,
    variables_sin_variabilidad
  )
  
} else {
  
  variables_predictoras_finales <- variables_predictoras
  
  cat(
    "\nNo se encontraron predictores sin variabilidad.\n"
  )
}

cat(
  "Número final de predictores candidatos:",
  length(variables_predictoras_finales),
  "\n"
)


# ------------------------------------------------------------
# 11. Crear matriz X y vector y
# ------------------------------------------------------------

X <- datos_completos[
  ,
  variables_predictoras_finales,
  drop = FALSE
]

X[] <- lapply(
  X,
  function(x) as.numeric(as.character(x))
)

X <- as.matrix(X)

# Outcome:
# 0 = grupo original 0
# 1 = grupo original 2

y <- as.numeric(
  datos_completos[[variable_resultado]]
)

if (anyNA(X)) {
  stop("La matriz X contiene valores ausentes.")
}

if (anyNA(y)) {
  stop("El vector y contiene valores ausentes.")
}

if (!all(y %in% c(0, 1))) {
  stop(
    "El vector y contiene valores diferentes de 0 y 1."
  )
}

cat("\nDimensiones de X:", dim(X), "\n")
cat("Longitud de y:", length(y), "\n")
cat("Formato de X:", storage.mode(X), "\n")


# ------------------------------------------------------------
# 12. Validación cruzada estratificada de 5 folds
# ------------------------------------------------------------

set.seed(20260722)

n_clase_minoritaria <- min(table(y))

if (n_clase_minoritaria < 5) {
  stop(
    paste0(
      "La categoría minoritaria contiene únicamente ",
      n_clase_minoritaria,
      " participantes. No pueden generarse ",
      "5 folds estratificados."
    )
  )
}

numero_folds <- 5

foldid <- integer(length(y))

# Distribuir participantes del grupo original 0
foldid[y == 0] <- sample(
  rep(
    seq_len(numero_folds),
    length.out = sum(y == 0)
  )
)

# Distribuir participantes del grupo original 2
foldid[y == 1] <- sample(
  rep(
    seq_len(numero_folds),
    length.out = sum(y == 1)
  )
)

cat("\nDistribución de los grupos en cada fold:\n")
cat("0 = grupo original 0; 1 = grupo original 2\n")

print(
  table(foldid, y)
)


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
  error_estandar_lambda_min = NA_real_,
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
    nlambda = 100,
    maxit = 1000000
  )
  
  lista_modelos[[i]] <- modelo_cv
  
  indice_lambda_min <- which.min(
    modelo_cv$cvm
  )
  
  coef_min <- coef(
    modelo_cv,
    s = "lambda.min"
  )
  
  coef_1se <- coef(
    modelo_cv,
    s = "lambda.1se"
  )
  
  resultados_alpha$lambda_min[i] <-
    modelo_cv$lambda.min
  
  resultados_alpha$lambda_1se[i] <-
    modelo_cv$lambda.1se
  
  resultados_alpha$deviance_cv_minima[i] <-
    modelo_cv$cvm[indice_lambda_min]
  
  resultados_alpha$error_estandar_lambda_min[i] <-
    modelo_cv$cvsd[indice_lambda_min]
  
  # Tolerancia numérica:
  # evita considerar seleccionados coeficientes
  # infinitesimales en modelos Ridge.
  
  resultados_alpha$n_variables_lambda_min[i] <-
    sum(
      abs(as.vector(coef_min)[-1]) > 1e-6
    )
  
  resultados_alpha$n_variables_lambda_1se[i] <-
    sum(
      abs(as.vector(coef_1se)[-1]) > 1e-6
    )
}

cat("\nResultados para cada valor de alpha:\n")

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
cat("Lambda mínimo:", modelo_final_cv$lambda.min, "\n")
cat("Lambda 1-SE:", modelo_final_cv$lambda.1se, "\n")


# ------------------------------------------------------------
# 15. Gráfico del rendimiento según alpha
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
      "Comparación grupo 0 frente a grupo 2; alpha óptimo = ",
      mejor_alpha
    ),
    x = "Alpha",
    y = "Deviance binomial media en validación cruzada"
  ) +
  theme_classic()

print(grafico_alpha)


# ------------------------------------------------------------
# 16. Gráfico de validación cruzada
# ------------------------------------------------------------

plot(modelo_final_cv)

title(
  paste0(
    "Elastic Net: grupo 0 frente a grupo 2; alpha = ",
    mejor_alpha
  )
)


# ------------------------------------------------------------
# 17. Coeficientes con lambda.min
# ------------------------------------------------------------

coeficientes_lambda_min <- as.matrix(
  coef(
    modelo_final_cv,
    s = "lambda.min"
  )
)

tabla_coeficientes_min <- data.frame(
  variable = rownames(coeficientes_lambda_min),
  
  coeficiente = as.numeric(
    coeficientes_lambda_min[, 1]
  ),
  
  row.names = NULL
) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    seleccionada = abs(coeficiente) > 1e-6,
    
    OR_penalizada = exp(coeficiente),
    
    direccion = case_when(
      coeficiente > 1e-6 ~
        "Mayor probabilidad de pertenecer al grupo 2",
      
      coeficiente < -1e-6 ~
        "Mayor probabilidad de pertenecer al grupo 0",
      
      TRUE ~
        "No seleccionada"
    )
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
    OR_penalizada,
    direccion
  )

cat("\nVariables seleccionadas con lambda.min:\n")

print(variables_seleccionadas_min)


# ------------------------------------------------------------
# 19. Coeficientes con lambda.1se
# ------------------------------------------------------------

coeficientes_lambda_1se <- as.matrix(
  coef(
    modelo_final_cv,
    s = "lambda.1se"
  )
)

tabla_coeficientes_1se <- data.frame(
  variable = rownames(coeficientes_lambda_1se),
  
  coeficiente = as.numeric(
    coeficientes_lambda_1se[, 1]
  ),
  
  row.names = NULL
) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    seleccionada = abs(coeficiente) > 1e-6,
    
    OR_penalizada = exp(coeficiente),
    
    direccion = case_when(
      coeficiente > 1e-6 ~
        "Mayor probabilidad de pertenecer al grupo 2",
      
      coeficiente < -1e-6 ~
        "Mayor probabilidad de pertenecer al grupo 0",
      
      TRUE ~
        "No seleccionada"
    )
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
    OR_penalizada,
    direccion
  )

cat("\nVariables seleccionadas con lambda.1se:\n")

print(variables_seleccionadas_1se)


# ------------------------------------------------------------
# 21. Predicciones out-of-fold y curva ROC
# ------------------------------------------------------------

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
  nlambda = 100,
  maxit = 1000000,
  keep = TRUE
)

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

indices_validos <- complete.cases(
  probabilidades_cv,
  y
)

# ROC:
# 0 = grupo original 0
# 1 = grupo original 2

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

ci_auc_cv <- ci.auc(
  roc_cv,
  method = "delong"
)

cat(
  "\nAUC basada en predicciones out-of-fold:",
  round(auc_cv, 3),
  "\n"
)

cat(
  "IC 95% del AUC:",
  round(ci_auc_cv[1], 3),
  "-",
  round(ci_auc_cv[3], 3),
  "\n"
)

plot(
  roc_cv,
  main = paste0(
    "ROC Elastic Net: grupo 0 frente a grupo 2\n",
    "AUC = ",
    round(auc_cv, 3),
    " (IC 95%: ",
    round(ci_auc_cv[1], 3),
    "-",
    round(ci_auc_cv[3], 3),
    ")"
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)


# ------------------------------------------------------------
# 22. Guardar predicciones individuales
# ------------------------------------------------------------

predicciones_individuales <- data.frame(
  grupo_modelo = y,
  
  grupo_original = ifelse(
    y == 0,
    0,
    2
  ),
  
  probabilidad_cv_grupo_2 = probabilidades_cv
)


# ------------------------------------------------------------
# 23. Guardar resultados
# ------------------------------------------------------------

directorio_resultados <- paste0(
  "C:/Users/victo/OneDrive/Desktop/Research/INMUNEA/",
  "Estadistica Victor/Resultados_Elastic_Net_Grupo_0_vs_2"
)

if (!dir.exists(directorio_resultados)) {
  dir.create(
    directorio_resultados,
    recursive = TRUE
  )
}

write.csv(
  missing_por_variable,
  file = file.path(
    directorio_resultados,
    "datos_ausentes_por_variable.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  resultados_alpha,
  file = file.path(
    directorio_resultados,
    "comparacion_alpha_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tabla_coeficientes_min,
  file = file.path(
    directorio_resultados,
    "coeficientes_lambda_min_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  variables_seleccionadas_min,
  file = file.path(
    directorio_resultados,
    "variables_seleccionadas_lambda_min_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  tabla_coeficientes_1se,
  file = file.path(
    directorio_resultados,
    "coeficientes_lambda_1se_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  variables_seleccionadas_1se,
  file = file.path(
    directorio_resultados,
    "variables_seleccionadas_lambda_1se_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  predicciones_individuales,
  file = file.path(
    directorio_resultados,
    "predicciones_out_of_fold_grupo_0_vs_2.csv"
  ),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat(
  "\nResultados guardados en:\n",
  directorio_resultados,
  "\n"
)