#' @title Rice Detection via Polynomial Regression
#' @description Checks if a time series segment fits a rice-like growth curve using quadratic regression.
#' @param x Vector of time indices.
#' @param y Vector of vegetation index values.
#' @param alpha Significance level for the regression coefficients.
#' @return Logical; TRUE if the segment matches rice growth criteria.
#' @export
regression.rice <- function(x, y, alpha = 0.05) {
  dat <- data.frame(x, y)
  dat <- dat[!is.na(dat$x), ]
  lm.poly2 <- lm(y ~ poly(x, 2), data = dat)
  lm.coeffs <- summary(lm.poly2)$coefficients
  lm.coeffs[2, 1] > 0 & lm.coeffs[2, 4] <= alpha & lm.coeffs[3, 1] < 0 & lm.coeffs[3, 4] <= alpha
}

#' @title Get Model P-value
#' @description Extracts the overall p-value from a linear model object.
#' @param model A linear model object (class `lm`).
#' @param ... Additional arguments passed to `pf`.
#' @return The overall p-value of the model.
#' @export
lm.pvalue <- function(model, ...) {
  lm.summary <- summary(model)
  fstat <- lm.summary$fstatistic
  pf(q = fstat[1], df1 = fstat[2], df2 = fstat[3], ...)
}

#' @title Get Model R-squared
#' @description Extracts the R-squared value from a linear model object.
#' @param model A linear model object.
#' @return The R-squared value.
#' @export
lm.rsquare <- function(model) {
  lm.summary <- summary(model)
  lm.summary$r.squared
}

#' @title Get Model Coefficient or P-value
#' @description Extracts a specific coefficient or its p-value from a linear model.
#' @param model A linear model object.
#' @param coef.id Index of the coefficient (0 for intercept, 1 for first term, etc.).
#' @param pvalue Logical; if TRUE, returns the p-value instead of the estimate.
#' @return The coefficient estimate or p-value.
#' @export
lm.coeff <- function(model, coef.id = 1, pvalue = FALSE) {
  lm.summary <- summary(model)
  lm.summary$coefficients[coef.id + 1, ifelse(pvalue, 4, 1)]
}

#' @title Get Regression Intercept P-value (Experimental)
#' @description Extracts the p-value of the intercept from a simple linear regression.
#' @param x Vector of X values.
#' @param y Vector of Y values.
#' @param alpha Significance level.
#' @return The p-value of the intercept.
#' @export
reg.intercept <- function(x, y, alpha = 0.05) {
  lm.model <- lm(y ~ x, data = data.frame(x = 1:length(x)))
  summary(lm.model)$coefficients[1, 4]
}

#' @title Simple Regression Model (Experimental)
#' @description Returns a linear model object for a simple regression.
#' @param coeff Coefficient index.
#' @param alpha Significance level.
#' @return An `lm` object.
#' @export
reg.fstat <- function(coeff = 1, alpha = 0.05) {
  lm.model <- lm(y ~ x, data = data.frame(x = 1:length(y)))
  lm.model
}

#' @title T-test P-value
#' @description Performs a t-test and returns the p-value.
#' @param pixel.data Numeric vector of data.
#' @param ... Additional arguments passed to `t.test`.
#' @return The p-value, or 1 if data points are fewer than 25.
#' @export
ttest.p <- function(pixel.data, ...) {
  pixel.data <- pixel.data[!is.na(pixel.data)]
  if (length(pixel.data) >= 25) {
    tres <- t.test(x = pixel.data, ...)
    result <- round(tres$p.value, 3)
  } else {
    result <- 1
  }
  result
}

# z.test <- function(z){
#   pnorm
# }

#' @title Safe Savitzky-Golay Filter
#' @description Applies a Savitzky-Golay filter to a vector, ensuring enough data points are present.
#' @param data Numeric vector.
#' @param filt.len Filter length.
#' @param ... Additional arguments passed to `signal::sgolayfilt`.
#' @return The filtered numeric vector.
#' @export
safe.sg <- function(data, filt.len, ...) {
  nna <- which(!is.na(data))
  if (length(nna) > filt.len) {
    data[nna] <- signal::sgolayfilt(data[nna], n = filt.len, ...)
  }
  data
}

#' @title Linear Interpolation (Integer)
#' @description Performs linear interpolation and rounds results to the nearest integer.
#' @param ... Arguments passed to `approx`.
#' @return Rounded interpolated values.
#' @export
approx_int <- function(...) {
  result <- approx(...)$y
  result <- round(result)
  result
}

#' @title Label Result
#' @description Prepends a label to a data frame or matrix result.
#' @param x Data frame or matrix.
#' @param label Label string to prepend.
#' @param ... Additional arguments.
#' @return The modified data structure with the label prepended.
#' @export
labelResult <- function(x, label, ...) {
  if (length(x) > 0) {
    if (class(x) == "data.frame") x <- data.frame(label, x, ...) else x <- cbind(label, x)
  }
  x
}
