suppressPackageStartupMessages({
  library(tercen)
  library(dplyr)
  library(dtplyr)
  library(drc)
  library(data.table)
})

ctx <- tercenCtx()

if(!ctx$hasNumericXAxis) stop("x axis is missing.")
if(length(ctx$yAxis) < 1) stop("y axis is missing.")

function.type <- ctx$op.value('function.type', as.character, "Four-parameter log-logistic")
n.predictions <- ctx$op.value('n.predictions', as.double, 100)
response.output <- ctx$op.value('response.output', as.character, "50, 90, 99")
response.output <- as.numeric(trimws(strsplit(response.output, ",")[[1]]))
relative.response <- ctx$op.value('relative.response', as.logical, TRUE)
maximum.x <- ctx$op.value('maximum.x', as.double, 1e6)

dose.transformation <- ctx$op.value('dose.transformation', as.character, "None") # log10, none
dt <- switch(dose.transformation,
             Log = exp(1),
             Log10 = 10,
             None = NULL)

model.function <- switch(
  function.type,
  "Three-parameter log-logistic" = "LL.3",
  "Four-parameter log-logistic" = "LL.4",
  "Michaelis-Menten" = "MM.2"
)

par_names <- switch(
  function.type,
  "Three-parameter log-logistic" = c("b", "d", "e"),
  "Four-parameter log-logistic" = c("b", "c", "d", "e"),
  "Michaelis-Menten" = c("d", "e")
)

dt_in <- ctx %>% 
  dplyr::select(.x, .y, .ri, .ci) %>%
  data.table::as.data.table()

get_pseudo_r2 <- function(mod) {
  predicted <- mod$predres[, "Predicted values"]
  actual <- mod$predres[, "Residuals"] + predicted
  rss <- sum((predicted - actual) ^ 2)
  tss <- sum((actual - mean(actual)) ^ 2)
  1 - rss/tss
}

# NOTE: Global limits calculation removed from here. 
# It is safer to calculate limits relative to the specific subset of data in the loop.

df_result <- dt_in[, 
  {
      # 1. Fit the model
      mod <- try(drm(
        .y ~ .x, fct = match.fun(model.function)(), logDose = dt
      ), silent = TRUE)
      
      if(!inherits(mod, 'try-error')) {
        coef <- mod$coefficients
        names(coef) <- gsub(pattern = ":\\(Intercept\\)", "", names(coef))
        out <- as.data.frame(t(coef))
        
        # Generate Prediction Curve
        x.pred <- seq(min(.x), max(.x), length.out = n.predictions)
        y.pred <- predict(mod, newdata = data.frame(x.pred))
        out <- cbind(out, x.pred, y.pred)
        
        out["pseudo_R2"] <- get_pseudo_r2(mod)
        
        # 2. Calculate Inverse Predictions (ED values)
        if(model.function %in% c("LL.3", "LL.4", "MM.2")) {
          
          # Optimization: Define search limits locally based on this group's data range
          # Searching 0 to 1,000,000 for a curve spanning 0.1-10 is computationally dangerous
          
          rx <- range(.x, na.rm = TRUE)
          # Ensure lower bound is > 0 for Log models to avoid log(0)=-Inf
          lower_limit <- if(rx[1] <= 0) 1e-9 else rx[1] / 100
          upper_limit <- min(rx[2] * 100, maximum.x) # Cap at Maximum X input
          
          local_limits <- c(lower_limit, upper_limit)
          
          f <- function(x, y) y - predict(mod, data.frame(.x = x))[1]
          
          for(i in response.output) {
            y_ed <- ifelse(
              relative.response & model.function == "LL.4",
              ((out$d[1] - out$c[1]) * i / 100) + out$c[1],
              out$d[1] * i / 100
            )
            
            # Use local_limits instead of global limits
            x <- try(uniroot(f, local_limits, y = y_ed)$root, silent = TRUE)
            
            if(inherits(x, 'try-error')) x <- NA_real_
            vn <- paste0("X", i)
            out[[paste0("X", i)]] <- x
            out[[paste0("Y", i)]] <- as.double(y_ed)
          }
        }
      } else {
        # Handle Failures
        if(length(unique(.y)) == 1) {
           x.pred <- seq(min(.x), max(.x), length.out = n.predictions)
           out <- data.frame(x.pred = x.pred, y.pred = .y[1])
           out[paste0("X", response.output)] <- NA_real_
           out[paste0("Y", response.output)] <- .y[1]
        } else {
           out <- data.frame(x.pred = NA_real_, y.pred = NA_real_)
           out[paste0("X", response.output)] <- NA_real_
           out[paste0("Y", response.output)] <- NA_real_
        }
        out[par_names] <- NA_real_
        out["pseudo_R2"] <- NA_real_
        # Ensure columns exist even on failure
        out[paste0("X", response.output)] <- NA_real_
        out[paste0("Y", response.output)] <- NA_real_
        
        out <- out[c(par_names, "x.pred", "y.pred", "pseudo_R2", paste0("X", response.output), paste0("Y", response.output))]
      }
  out
  }, by = c(".ri", ".ci")
]

if(model.function == "LL.4") {
  if("d" %in% names(df_result) && "c" %in% names(df_result)){
    df_result <- df_result %>% mutate(Span = d - c)
  }
} 

sum.table <- df_result %>%
  dplyr::select(-x.pred, -y.pred) %>%
  unique() %>% 
  arrange(.ri, .ci) %>%
  as_tibble() %>%
  ctx$addNamespace() 

pred.table <- df_result %>%
  dplyr::select(.ri, .ci, x.pred, y.pred) %>%
  dplyr::rename(x_pred = x.pred, y_pred = y.pred) %>%
  arrange(.ri, .ci) %>%
  as_tibble() %>%
  ctx$addNamespace()

ctx$save(list(sum.table, pred.table))
