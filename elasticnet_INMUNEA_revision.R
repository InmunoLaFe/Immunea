# ============================================================
# INMUNEA — ELASTIC NET Y MODELO LOGÍSTICO EXPLORATORIO
# Versión revisada para la respuesta a los revisores
# (Reviewer 1, Comment 3: IC de sensibilidad/especificidad y
#  magnitud/contracción de los coeficientes del Elastic Net)
#
# Contrastes (outcome = 1 / referencia = 0):
#   A. AD vs NDC      grupo 1 vs grupo 0
#   B. AD vs OD       grupo 2 (OD) = 1, grupo 1 (AD) = 0
#   C. OD vs NDC      grupo 2 (OD) = 1, grupo 0 (NDC) = 0
#   D. AD vs no-AD    Factor_caso_1_EA (1 = AD; 0 = OD + NDC)
#
# Cambios respecto a la versión original (elasticnet_casos_vs_*.R):
#   1. Un único script para los cuatro contrastes (misma semilla,
#      mismos folds estratificados, misma rejilla de alpha y mismos
#      argumentos de cv.glmnet que los scripts originales).
#   2. Coeficientes estandarizados (log-odds por 1 DE), para comparar
#      su magnitud entre variables (Supplementary Table 5).
#   3. Recuento de variables retenidas con tolerancia numérica: con
#      Ridge (alpha = 0) los coeficientes no son exactamente cero
#      (~1e-37) y antes se contaban como "seleccionados".
#   4. Se marca cuando lambda.min coincide con el lambda máximo de la
#      secuencia (modelo nulo, solo intercepto).
#   5. fit.preval de cv.glmnet (glmnet >= 4.0) está en escala del
#      predictor lineal (logit): se transforma con plogis() antes de
#      guardarlo como probabilidad. El AUC no cambia.
#   6. Se guardan las trayectorias de los coeficientes (paths) y la
#      curva de validación cruzada del alpha seleccionado.
#   7. Elastic Net restringido a los 4 predictores del modelo
#      logístico exploratorio, comparando coeficientes penalizados y
#      no penalizados (contracción; Supplementary Table 5C).
#   8. Modelo logístico exploratorio AD vs NDC: AUC con IC de DeLong
#      y IC 95% bootstrap de sensibilidad y especificidad en el
#      umbral de Youden, con semilla fija (reproducible).
#   9. Supplementary Figure 3.
# ============================================================


# ------------------------------------------------------------
# 1. Instalar y cargar paquetes
# ------------------------------------------------------------

paquetes <- c("readxl", "dplyr", "glmnet", "pROC", "ggplot2", "patchwork", "ggrepel")

paquetes_nuevos <- paquetes[!paquetes %in% rownames(installed.packages())]
if (length(paquetes_nuevos) > 0) {
  install.packages(paquetes_nuevos)
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(glmnet)
  library(pROC)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})


# ------------------------------------------------------------
# 2. Rutas e importación de la base de datos
# ------------------------------------------------------------

# Se pueden pasar como argumentos: Rscript script.R <base.xlsx> <carpeta_resultados>
argumentos <- commandArgs(trailingOnly = TRUE)

ruta_archivo <- if (length(argumentos) >= 1) argumentos[1] else paste0(
  "C:/Users/victo/OneDrive/Desktop/Research/INMUNEA/",
  "Estadistica Victor/Clean_INMUNEA_ALZHEIMER_DB.xlsx"
)

directorio_resultados <- if (length(argumentos) >= 2) argumentos[2] else paste0(
  "C:/Users/victo/OneDrive/Desktop/Research/INMUNEA/",
  "Estadistica Victor/Resultados_revision_R1C3"
)

if (!dir.exists(directorio_resultados)) {
  dir.create(directorio_resultados, recursive = TRUE)
}

datos_originales <- read_excel(
  path = ruta_archivo,
  sheet = 1,
  .name_repair = "minimal"
)

datos_originales$grupo <- as.numeric(as.character(datos_originales$grupo))

cat(
  "Dimensiones de la base original:",
  nrow(datos_originales), "filas y",
  ncol(datos_originales), "columnas\n"
)


# ------------------------------------------------------------
# 3. Variables predictoras (36 inmunológicas + edad + sexo)
# ------------------------------------------------------------

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
  "sex"            # 0 = hombre, 1 = mujer
)

variables_no_encontradas <- setdiff(
  c("grupo", "Factor_caso_1_EA", variables_predictoras),
  names(datos_originales)
)

if (length(variables_no_encontradas) > 0) {
  stop(paste0(
    "Estas variables no se han encontrado en la base de datos:\n",
    paste(variables_no_encontradas, collapse = "\n")
  ))
}


# ------------------------------------------------------------
# 4. Funciones auxiliares
# ------------------------------------------------------------

# Tolerancia para considerar un coeficiente distinto de cero
tolerancia <- 1e-10

# Desviación estándar con denominador n (la que usa glmnet al estandarizar)
de_glmnet <- function(x) sqrt(mean((x - mean(x))^2))

# Selección de casos completos y construcción de X e y
preparar_datos <- function(datos, variable_resultado) {

  datos_modelo <- datos[, c(variable_resultado, variables_predictoras)]

  datos_modelo[] <- lapply(datos_modelo, function(x) as.numeric(as.character(x)))

  datos_completos <- datos_modelo %>% filter(complete.cases(.))

  y <- datos_completos[[variable_resultado]]

  if (!all(y %in% c(0, 1)) || length(unique(y)) != 2) {
    stop("El outcome debe contener exclusivamente 0 y 1, con ambas categorías.")
  }

  # Eliminar predictores sin variabilidad
  n_valores <- sapply(datos_completos[variables_predictoras], function(x) length(unique(x)))
  predictores_finales <- variables_predictoras[n_valores > 1]

  X <- as.matrix(datos_completos[, predictores_finales, drop = FALSE])

  list(
    X = X,
    y = y,
    n_inicial = nrow(datos_modelo),
    n_final = nrow(datos_completos)
  )
}

# Folds estratificados de 5 particiones (idéntico a los scripts originales)
crear_folds <- function(y, numero_folds = 5, semilla = 20260722) {

  set.seed(semilla)

  foldid <- integer(length(y))

  foldid[y == 0] <- sample(rep(seq_len(numero_folds), length.out = sum(y == 0)))
  foldid[y == 1] <- sample(rep(seq_len(numero_folds), length.out = sum(y == 1)))

  foldid
}

# Elastic Net con optimización de alpha y lambda
ajustar_elastic_net <- function(X, y, etiqueta,
                                grid_alpha = seq(0, 1, by = 0.1),
                                nlambda = 100, maxit = 1000000) {

  foldid <- crear_folds(y)

  lista_modelos <- lapply(grid_alpha, function(alpha_actual) {
    cv.glmnet(
      x = X,
      y = y,
      family = "binomial",
      alpha = alpha_actual,
      foldid = foldid,
      type.measure = "deviance",
      standardize = TRUE,
      intercept = TRUE,
      nlambda = nlambda,
      maxit = maxit,
      keep = TRUE
    )
  })

  resultados_alpha <- data.frame(
    contraste = etiqueta,
    alpha = grid_alpha,
    lambda_min = sapply(lista_modelos, function(m) m$lambda.min),
    lambda_1se = sapply(lista_modelos, function(m) m$lambda.1se),
    lambda_max = sapply(lista_modelos, function(m) max(m$lambda)),
    deviance_cv_minima = sapply(lista_modelos, function(m) min(m$cvm)),
    error_estandar_lambda_min = sapply(lista_modelos, function(m) m$cvsd[which.min(m$cvm)]),
    n_variables_lambda_min = sapply(lista_modelos, function(m)
      sum(abs(as.vector(coef(m, s = "lambda.min"))[-1]) > tolerancia)),
    n_variables_lambda_1se = sapply(lista_modelos, function(m)
      sum(abs(as.vector(coef(m, s = "lambda.1se"))[-1]) > tolerancia))
  )

  # lambda.min en el extremo de la secuencia = modelo nulo (solo intercepto)
  resultados_alpha$modelo_nulo <-
    abs(resultados_alpha$lambda_min - resultados_alpha$lambda_max) < 1e-12

  # Mejor alpha: menor deviance en validación cruzada
  indice_mejor_alpha <- which.min(resultados_alpha$deviance_cv_minima)
  mejor_alpha <- grid_alpha[indice_mejor_alpha]
  modelo_final_cv <- lista_modelos[[indice_mejor_alpha]]

  # Coeficientes en escala original y estandarizados (por 1 DE)
  de_x <- apply(X, 2, de_glmnet)

  coef_min <- as.vector(coef(modelo_final_cv, s = "lambda.min"))[-1]
  coef_1se <- as.vector(coef(modelo_final_cv, s = "lambda.1se"))[-1]

  tabla_coeficientes <- data.frame(
    contraste = etiqueta,
    variable = colnames(X),
    de = de_x,
    coeficiente_lambda_min = coef_min,
    coeficiente_estandarizado_lambda_min = coef_min * de_x,
    OR_por_DE_lambda_min = exp(coef_min * de_x),
    seleccionada_lambda_min = abs(coef_min) > tolerancia,
    coeficiente_lambda_1se = coef_1se,
    coeficiente_estandarizado_lambda_1se = coef_1se * de_x,
    OR_por_DE_lambda_1se = exp(coef_1se * de_x),
    seleccionada_lambda_1se = abs(coef_1se) > tolerancia,
    row.names = NULL
  ) %>%
    arrange(desc(abs(coeficiente_estandarizado_lambda_min)))

  # Trayectorias de los coeficientes estandarizados (paths)
  trayectorias <- as.matrix(modelo_final_cv$glmnet.fit$beta) * de_x

  # Predicciones out-of-fold en lambda.min
  indice_lambda_pred <- which.min(abs(modelo_final_cv$lambda - modelo_final_cv$lambda.min))
  predictor_lineal_cv <- modelo_final_cv$fit.preval[, indice_lambda_pred]

  roc_cv <- roc(
    response = y,
    predictor = predictor_lineal_cv,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  ci_auc_cv <- ci.auc(roc_cv, method = "delong")

  predicciones_individuales <- data.frame(
    resultado_real = y,
    fold = foldid,
    predictor_lineal_cv = predictor_lineal_cv,
    probabilidad_cv_grupo_1 = plogis(predictor_lineal_cv)
  )

  list(
    etiqueta = etiqueta,
    n_outcome = sum(y == 1),
    n_referencia = sum(y == 0),
    resultados_alpha = resultados_alpha,
    mejor_alpha = mejor_alpha,
    modelo = modelo_final_cv,
    coeficientes = tabla_coeficientes,
    trayectorias = trayectorias,
    predicciones = predicciones_individuales,
    auc_cv = as.numeric(auc(roc_cv)),
    ci_auc_cv = as.numeric(ci_auc_cv)[c(1, 3)],
    modelo_nulo = resultados_alpha$modelo_nulo[indice_mejor_alpha]
  )
}

guardar_csv <- function(tabla, nombre) {
  write.csv(tabla, file.path(directorio_resultados, nombre),
            row.names = FALSE, fileEncoding = "UTF-8")
}


# ------------------------------------------------------------
# 5. Elastic Net en los cuatro contrastes (38 predictores)
# ------------------------------------------------------------

# A. AD (1) vs NDC (0)
datos_A <- datos_originales %>% filter(!is.na(grupo), grupo %in% c(0, 1))
prep_A <- preparar_datos(datos_A, "grupo")

# B. AD vs OD: OD (grupo 2) = 1, AD (grupo 1) = 0
datos_B <- datos_originales %>%
  filter(!is.na(grupo), grupo %in% c(1, 2)) %>%
  mutate(grupo_2_vs_1 = as.numeric(grupo == 2))
prep_B <- preparar_datos(datos_B, "grupo_2_vs_1")

# C. OD vs NDC: OD (grupo 2) = 1, NDC (grupo 0) = 0
datos_C <- datos_originales %>%
  filter(!is.na(grupo), grupo %in% c(0, 2)) %>%
  mutate(grupo_2_vs_0 = as.numeric(grupo == 2))
prep_C <- preparar_datos(datos_C, "grupo_2_vs_0")

# D. AD vs no-AD (OD + NDC); el script original usaba nlambda = 200
#    y maxit por defecto (100000)
prep_D <- preparar_datos(datos_originales, "Factor_caso_1_EA")

resultados_en <- list(
  ajustar_elastic_net(prep_A$X, prep_A$y, "AD vs NDC (outcome: AD)"),
  ajustar_elastic_net(prep_B$X, prep_B$y, "AD vs OD (outcome: OD)"),
  ajustar_elastic_net(prep_C$X, prep_C$y, "OD vs NDC (outcome: OD)"),
  ajustar_elastic_net(prep_D$X, prep_D$y, "AD vs no-AD (outcome: AD)",
                      nlambda = 200, maxit = 100000)
)

for (res in resultados_en) {
  cat("\n------------------------------------------------------------\n")
  cat(res$etiqueta, "| n =", res$n_outcome + res$n_referencia,
      "(", res$n_outcome, "/", res$n_referencia, ")",
      "| mejor alpha =", res$mejor_alpha,
      "| modelo nulo:", res$modelo_nulo, "\n")
  print(res$resultados_alpha)
  cat("Variables retenidas (lambda.min / lambda.1se):\n")
  print(res$coeficientes %>%
          filter(seleccionada_lambda_min | seleccionada_lambda_1se) %>%
          select(variable, coeficiente_estandarizado_lambda_min,
                 coeficiente_estandarizado_lambda_1se))
}

# Supplementary Table 5A: resumen por contraste
tabla_suplementaria_5A <- do.call(rbind, lapply(resultados_en, function(res) {
  fila <- res$resultados_alpha[res$resultados_alpha$alpha == res$mejor_alpha, ]
  data.frame(
    contraste = res$etiqueta,
    n_outcome = res$n_outcome,
    n_referencia = res$n_referencia,
    alpha = res$mejor_alpha,
    lambda_min = fila$lambda_min,
    lambda_1se = fila$lambda_1se,
    modelo_nulo = fila$modelo_nulo,
    variables_lambda_min = fila$n_variables_lambda_min,
    variables_lambda_1se = fila$n_variables_lambda_1se,
    # en un modelo nulo el AUC out-of-fold no es interpretable
    auc_out_of_fold = ifelse(fila$modelo_nulo, NA, res$auc_cv),
    auc_ic95_inf = ifelse(fila$modelo_nulo, NA, res$ci_auc_cv[1]),
    auc_ic95_sup = ifelse(fila$modelo_nulo, NA, res$ci_auc_cv[2])
  )
}))

cat("\nSupplementary Table 5A:\n")
print(tabla_suplementaria_5A)

# Supplementary Table 5B: coeficientes estandarizados retenidos
tabla_suplementaria_5B <- do.call(rbind, lapply(resultados_en, function(res) {
  res$coeficientes %>% filter(seleccionada_lambda_min | seleccionada_lambda_1se)
}))

cat("\nSupplementary Table 5B:\n")
print(tabla_suplementaria_5B)

guardar_csv(tabla_suplementaria_5A, "SupplTable5A_resumen_contrastes.csv")
guardar_csv(tabla_suplementaria_5B, "SupplTable5B_coeficientes_retenidos.csv")
guardar_csv(do.call(rbind, lapply(resultados_en, `[[`, "resultados_alpha")),
            "comparacion_alpha_todos_contrastes.csv")
guardar_csv(do.call(rbind, lapply(resultados_en, `[[`, "coeficientes")),
            "coeficientes_estandarizados_todos_contrastes.csv")

for (i in seq_along(resultados_en)) {
  res <- resultados_en[[i]]
  sufijo <- c("AD_vs_NDC", "AD_vs_OD", "OD_vs_NDC", "AD_vs_noAD")[i]
  guardar_csv(res$predicciones, paste0("predicciones_out_of_fold_", sufijo, ".csv"))
  guardar_csv(
    data.frame(log_lambda = log(res$modelo$glmnet.fit$lambda), t(res$trayectorias),
               check.names = FALSE),
    paste0("trayectorias_coeficientes_", sufijo, ".csv")
  )
}


# ------------------------------------------------------------
# 6. Modelo logístico exploratorio AD vs NDC
#    (age, sex, DN γδ y DN no-γδ; misma especificación que el
#     modelo publicado)
# ------------------------------------------------------------

datos_logistico <- datos_originales %>%
  filter(grupo %in% c(0, 1)) %>%
  transmute(
    grupo = grupo,
    age = as.numeric(age),
    sex = as.numeric(sex),
    dn_gd = as.numeric(`doble_negativas_γ/δ`),
    dn_no_gd = as.numeric(`doble_negativas_no_γ/δ`)
  ) %>%
  na.omit()

cat("\nModelo logístico: n =", nrow(datos_logistico),
    "| AD =", sum(datos_logistico$grupo == 1),
    "| NDC =", sum(datos_logistico$grupo == 0), "\n")

modelo_logistico <- glm(
  grupo ~ age + dn_gd + dn_no_gd + sex,
  data = datos_logistico,
  family = binomial
)

tabla_OR <- data.frame(
  variable = names(coef(modelo_logistico)),
  OR = exp(coef(modelo_logistico)),
  exp(suppressMessages(confint(modelo_logistico))),
  row.names = NULL,
  check.names = FALSE
)
names(tabla_OR)[3:4] <- c("IC95_inf", "IC95_sup")
print(tabla_OR, digits = 3)

probabilidades <- predict(modelo_logistico, type = "response")

roc_logistico <- roc(
  response = datos_logistico$grupo,
  predictor = probabilidades,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

ic_auc_delong <- ci.auc(roc_logistico, method = "delong")

set.seed(20261022)
ic_auc_bootstrap <- ci.auc(roc_logistico, method = "bootstrap",
                           boot.n = 2000, boot.stratified = TRUE)

umbral_youden <- coords(
  roc_logistico, "best", best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity", "tp", "fn", "tn", "fp"),
  transpose = FALSE
)

# IC 95% bootstrap (2.000 réplicas estratificadas, percentil),
# re-estimando el umbral de Youden en cada réplica
set.seed(20261022)
ic_youden <- ci.coords(
  roc_logistico,
  x = "best",
  best.method = "youden",
  ret = c("sensitivity", "specificity"),
  boot.n = 2000,
  boot.stratified = TRUE
)

tabla_rendimiento <- data.frame(
  AUC = as.numeric(auc(roc_logistico)),
  AUC_IC95_DeLong_inf = ic_auc_delong[1],
  AUC_IC95_DeLong_sup = ic_auc_delong[3],
  AUC_IC95_bootstrap_inf = ic_auc_bootstrap[1],
  AUC_IC95_bootstrap_sup = ic_auc_bootstrap[3],
  umbral_Youden = umbral_youden$threshold,
  VP = umbral_youden$tp, FN = umbral_youden$fn,
  VN = umbral_youden$tn, FP = umbral_youden$fp,
  sensibilidad = umbral_youden$sensitivity,
  sensibilidad_IC95_inf = ic_youden$sensitivity[1],
  sensibilidad_IC95_sup = ic_youden$sensitivity[3],
  especificidad = umbral_youden$specificity,
  especificidad_IC95_inf = ic_youden$specificity[1],
  especificidad_IC95_sup = ic_youden$specificity[3]
)

cat("\nRendimiento del modelo logístico exploratorio:\n")
print(t(round(tabla_rendimiento, 3)))

guardar_csv(tabla_OR, "modelo_logistico_OR.csv")
guardar_csv(tabla_rendimiento, "modelo_logistico_AUC_sensibilidad_especificidad.csv")


# ------------------------------------------------------------
# 7. Elastic Net con los 4 predictores del modelo exploratorio
#    (contracción frente a los coeficientes no penalizados;
#     Supplementary Table 5C)
# ------------------------------------------------------------

X_4 <- as.matrix(datos_logistico[, c("age", "sex", "dn_gd", "dn_no_gd")])
colnames(X_4) <- c("age", "sex", "doble_negativas_γ/δ", "doble_negativas_no_γ/δ")

resultado_4 <- ajustar_elastic_net(X_4, datos_logistico$grupo,
                                   "AD vs NDC, 4 predictores (outcome: AD)")

print(resultado_4$resultados_alpha)

de_4 <- apply(X_4, 2, de_glmnet)
coef_no_penalizados <- coef(modelo_logistico)[c("age", "sex", "dn_gd", "dn_no_gd")]
ic_no_penalizados <- suppressMessages(confint(modelo_logistico))[c("age", "sex", "dn_gd", "dn_no_gd"), ]

coef_en_4 <- resultado_4$coeficientes[match(colnames(X_4), resultado_4$coeficientes$variable), ]

tabla_suplementaria_5C <- data.frame(
  variable = colnames(X_4),
  de = de_4,
  beta_estandarizado_no_penalizado = coef_no_penalizados * de_4,
  IC95_inf = ic_no_penalizados[, 1] * de_4,
  IC95_sup = ic_no_penalizados[, 2] * de_4,
  beta_estandarizado_EN_lambda_min = coef_en_4$coeficiente_estandarizado_lambda_min,
  beta_estandarizado_EN_lambda_1se = coef_en_4$coeficiente_estandarizado_lambda_1se,
  row.names = NULL
)

tabla_suplementaria_5C$contraccion_lambda_min_pct <- 100 * (
  1 - tabla_suplementaria_5C$beta_estandarizado_EN_lambda_min /
    tabla_suplementaria_5C$beta_estandarizado_no_penalizado
)

cat("\nSupplementary Table 5C (alpha =", resultado_4$mejor_alpha,
    "; lambda.min =", round(resultado_4$modelo$lambda.min, 3),
    "; lambda.1se =", round(resultado_4$modelo$lambda.1se, 3), "):\n")
print(tabla_suplementaria_5C, digits = 3)
cat("AUC out-of-fold (lambda.min):", round(resultado_4$auc_cv, 3),
    "(IC 95%", round(resultado_4$ci_auc_cv, 3), ")\n")

guardar_csv(tabla_suplementaria_5C, "SupplTable5C_contraccion_4_predictores.csv")
guardar_csv(resultado_4$resultados_alpha, "comparacion_alpha_4_predictores.csv")


# ------------------------------------------------------------
# 8. Supplementary Figure 3: curvas de validación cruzada y
#    trayectorias de los coeficientes estandarizados
# ------------------------------------------------------------

tinta <- "#0b0b0b"; tinta2 <- "#52514e"; rejilla <- "#e4e3df"; gris <- "#c9c8c2"

colores <- c(
  "Age" = "#2a78d6", "Sex (female)" = "#eb6834", "DN γδ T cells" = "#1baf7a",
  "DN non-γδ T cells" = "#eda100", "Th17-like Tregs" = "#e87ba4", "Th9-like cells" = "#008300"
)

etiquetas <- c(
  "age" = "Age", "sex" = "Sex (female)", "doble_negativas_γ/δ" = "DN γδ T cells",
  "doble_negativas_no_γ/δ" = "DN non-γδ T cells", "Treg_Th17-like" = "Th17-like Tregs",
  "Th9" = "Th9-like cells"
)

tema <- theme_classic(base_size = 8.5) +
  theme(
    text = element_text(colour = tinta), axis.text = element_text(colour = tinta2),
    axis.line = element_line(colour = tinta2, linewidth = 0.3),
    axis.ticks = element_line(colour = tinta2, linewidth = 0.3),
    panel.grid.major.y = element_line(colour = rejilla, linewidth = 0.25),
    plot.title = element_text(size = 8.5, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 7.5, colour = tinta2, hjust = 0),
    legend.position = "none", plot.margin = margin(4, 8, 4, 4)
  )

panel_cv <- function(m, titulo, subtitulo, xmin = -Inf) {
  d <- data.frame(ll = log(m$lambda), cvm = m$cvm, lo = m$cvm - m$cvsd, hi = m$cvm + m$cvsd)
  d <- d[d$ll >= xmin, ]
  iguales <- isTRUE(all.equal(m$lambda.min, m$lambda.1se))
  etq <- if (iguales) {
    data.frame(x = log(m$lambda.min), y = max(d$hi), t = "λmin = λ1se", h = 1.05)
  } else {
    data.frame(x = log(c(m$lambda.min, m$lambda.1se)), y = max(d$hi),
               t = c("λmin", "λ1se"), h = c(1.1, -0.1))
  }
  ggplot(d, aes(ll, cvm)) +
    geom_text(data = etq, aes(x, y, label = t, hjust = h), size = 2.4, colour = tinta, vjust = 1) +
    geom_linerange(aes(ymin = lo, ymax = hi), colour = gris, linewidth = 0.35) +
    geom_point(colour = tinta2, size = 0.7) +
    geom_vline(xintercept = log(m$lambda.min), linetype = "dashed", linewidth = 0.35, colour = tinta) +
    geom_vline(xintercept = log(m$lambda.1se), linetype = "dotted", linewidth = 0.45, colour = tinta) +
    labs(x = expression(log(lambda)), y = "CV binomial deviance (±1 SE)",
         title = titulo, subtitle = subtitulo) +
    tema
}

panel_trayectorias <- function(m, trayectorias, destacar, no_penalizados = NULL, xmin = -Inf) {
  mantener <- log(m$glmnet.fit$lambda) >= xmin
  ll <- log(m$glmnet.fit$lambda)[mantener]
  trayectorias <- trayectorias[, mantener, drop = FALSE]
  d <- do.call(rbind, lapply(rownames(trayectorias), function(v)
    data.frame(var = v, ll = ll, b = trayectorias[v, ])))
  d$etiqueta <- ifelse(d$var %in% destacar, etiquetas[d$var], "Other")
  dd <- d[d$etiqueta != "Other", ]
  de <- do.call(rbind, lapply(split(dd, dd$var), function(x) x[which.min(x$ll), ]))
  rango <- range(ll)
  p <- ggplot() +
    geom_hline(yintercept = 0, colour = tinta2, linewidth = 0.3) +
    geom_line(data = d[d$etiqueta == "Other", ], aes(ll, b, group = var), colour = gris, linewidth = 0.35) +
    geom_line(data = dd, aes(ll, b, group = var, colour = etiqueta), linewidth = 0.7) +
    geom_vline(xintercept = log(m$lambda.min), linetype = "dashed", linewidth = 0.35, colour = tinta) +
    geom_vline(xintercept = log(m$lambda.1se), linetype = "dotted", linewidth = 0.45, colour = tinta) +
    scale_colour_manual(values = colores) +
    labs(x = expression(log(lambda)), y = "Standardised coefficient (per 1 SD)") +
    tema
  if (!is.null(no_penalizados)) {
    x_ml <- rango[1] - 0.08 * diff(rango)
    d_ml <- data.frame(ll = x_ml, b = no_penalizados$b, etiqueta = etiquetas[no_penalizados$var])
    p <- p +
      geom_point(data = d_ml, aes(ll, b, colour = etiqueta), shape = 17, size = 1.8) +
      annotate("text", x = x_ml, y = max(no_penalizados$b) + 0.12, label = "ML",
               size = 2.4, colour = tinta2)
    de$ll <- x_ml
  }
  p +
    geom_text_repel(data = de, aes(ll, b, label = etiqueta), colour = tinta, size = 2.4,
                    direction = "y", hjust = 1, nudge_x = -0.04 * diff(rango),
                    segment.colour = tinta2, segment.size = 0.2, min.segment.length = 0,
                    box.padding = 0.15, seed = 1) +
    scale_x_continuous(expand = expansion(mult = c(0.32, 0.03)))
}

res_A <- resultados_en[[1]]
res_B <- resultados_en[[2]]
destacar_A <- c("age", "sex", "doble_negativas_γ/δ", "doble_negativas_no_γ/δ")
destacar_B <- c("Treg_Th17-like", "Th9", "age")

figura <-
  (panel_cv(res_A$modelo, "A  AD vs NDC: 38 candidate predictors",
            sprintf("n = %d (%d AD, %d NDC); α = %g (ridge)",
                    res_A$n_outcome + res_A$n_referencia, res_A$n_outcome,
                    res_A$n_referencia, res_A$mejor_alpha)) |
     panel_trayectorias(res_A$modelo, res_A$trayectorias, destacar_A)) /
  (panel_cv(res_B$modelo, "B  AD vs OD: 38 candidate predictors",
            sprintf("n = %d (%d AD, %d OD); α = %g (LASSO)",
                    res_B$n_outcome + res_B$n_referencia, res_B$n_referencia,
                    res_B$n_outcome, res_B$mejor_alpha), xmin = -4.5) |
     panel_trayectorias(res_B$modelo, res_B$trayectorias, destacar_B, xmin = -4.5)) /
  (panel_cv(resultado_4$modelo, "C  AD vs NDC: exploratory-model predictors",
            sprintf("n = %d (%d AD, %d NDC); α = %g (ridge)",
                    resultado_4$n_outcome + resultado_4$n_referencia, resultado_4$n_outcome,
                    resultado_4$n_referencia, resultado_4$mejor_alpha)) |
     panel_trayectorias(resultado_4$modelo, resultado_4$trayectorias, destacar_A,
                        no_penalizados = data.frame(
                          var = tabla_suplementaria_5C$variable,
                          b = tabla_suplementaria_5C$beta_estandarizado_no_penalizado)))

ggsave(file.path(directorio_resultados, "Supplementary_Figure_3_ElasticNet.png"), figura,
       width = 180, height = 210, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(directorio_resultados, "Supplementary_Figure_3_ElasticNet.tiff"), figura,
       width = 180, height = 210, units = "mm", dpi = 300, compression = "lzw", bg = "white")
ggsave(file.path(directorio_resultados, "Supplementary_Figure_3_ElasticNet.pdf"), figura,
       width = 180, height = 210, units = "mm", device = cairo_pdf, bg = "white")


# ------------------------------------------------------------
# 9. Información de la sesión
# ------------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(directorio_resultados, "sessionInfo.txt"))

cat("\nResultados guardados en:\n", directorio_resultados, "\n")
