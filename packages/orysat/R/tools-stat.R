regression.rice <- function(x,y, alpha=0.05){
  dat <- data.frame(x,y)
  dat <- dat[!is.na(dat$x),]  
  lm.poly2 <- lm(y~poly(x,2), data = dat)
  lm.coeffs <- summary(lm.poly2)$coefficients
  return(lm.coeffs[2,1]>0 & lm.coeffs[2,4] <= alpha & lm.coeffs[3,1]<0 & lm.coeffs[3,4] <= alpha)
}

lm.pvalue <- function(model,...){
  lm.summary <- summary(model)
  fstat <- lm.summary$fstatistic
  return(pf(q=fstat[1],df1=fstat[2], df2=fstat[3], ...))
} 

lm.rsquare <- function(model){
  lm.summary <- summary(model)
  return(lm.summary$r.squared)
} 

lm.coeff <- function(model, coef.id=1, pvalue=FALSE){
  lm.summary <- summary(model)
  return(lm.summary$coefficients[coef.id+1, ifelse(pvalue,4,1)])
}

reg.intercept <- function(x, y, alpha=0.05){
  lm.model <- lm(y~x, data=data.frame(x=1:length))
  summary(lm.model)$coefficients[1,4]
}

reg.fstat <- function(coeff=1, alpha=0.05){
  lm.model <- lm(y~x, data=data.frame(x=1:length))
   
}

ttest.p <- function(pixel.data, ...){
   pixel.data<- pixel.data[!is.na(pixel.data)]
  if(length(pixel.data)>=25) {
    tres <- t.test(x=pixel.data, ...)
    result <- round(tres$p.value,3)
  } else result  <- 1
  return(result)
}

# z.test <- function(z){
#   pnorm
# }

safe.sg <- function(data, filt.len, ...){
  nna <- which(!is.na(data))
  if(length(nna)>filt.len){
    data[nna] <- signal::sgolayfilt(data[nna], n=filt.len, ...)  
  }
  return(data)
}

approx_int <- function(...){
  result <- approx(...)$y
  result <- round(result)
  return(result)
  
}

labelResult <- function(x, label, ...){
  
  if(length(x)>0){
    if (class(x)=="data.frame") x <- data.frame(label,x, ...) else x <- cbind(label,x)  
    
  } 
  return(x)
}
