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

safe.sg <- function(data,...){
  nna <- which(!is.na(data))
  data[nna] <- signal::sgolayfilt(data[nna],...)
  return(data)
}

approx_int <- function(...){
  result <- approx(...)$y
  result <- round(result)
  return(result)
  
}

polregIT.rice <- function(pix.evi, evi.date, dates.covered, width = 16, increments = 5, eos.evires=0.5, alpha=0.5){
  idx.starts <- seq(1,length(pix.evi)-width, by=increments)
  idx.ends <- idx.starts+width-1
  dat.rice <- vector()
  dat.evi <- data.frame(x=evi.date, y=pix.evi)
  dat.evi <- na.omit(dat.evi)
  last.rice <- 0 #as.numeric(dates.covered[length(dates.covered)])
  for(i in 1:length(idx.starts)){
    newx <- data.frame(x=as.numeric(seq(dates.covered[idx.starts[i]]-7, dates.covered[idx.ends[i]], by="day")))
    if(last.rice>newx$x[1]) next
    this.inc <- dat.evi[which(dat.evi$x %in% newx$x),]
    if(nrow(this.inc)<5) next
    this.inc$x <- as.numeric(this.inc$x)
    inc.lm <- lm(y~poly(x,2), data=this.inc)
    inc.smry <- summary(inc.lm)
    inc.fstat <- inc.smry$fstatistic
    pval <- pf(inc.fstat[1],inc.fstat[2], inc.fstat[3], lower.tail = FALSE)
    if(pval<=alpha &
    inc.smry$coefficients[2,1]>0 & inc.smry$coefficients[2,4]<=alpha  & # Increasing at the start
    inc.smry$coefficients[3,1]<0 & inc.smry$coefficients[3,4]<=alpha) # Decreasing after the peak
    {
      new.evi <- predict(inc.lm, newdata=newx)
      max.evi <- max(new.evi)
      eos.evi <- max.evi - ((max.evi-new.evi[1])*eos.evires)
      if(new.evi[length(new.evi)]<=eos.evi){
        cnt.rec <- nrow(this.inc)
        rsq <- inc.smry$r.squared
        rice.sos <- newx$x[1]
        pkevi <- which(new.evi==max.evi)
        rice.eos <- newx$x[pkevi+min(which(new.evi[pkevi:length(new.evi)]<=eos.evi))]
        dat.rice <- rbind(dat.rice, c(rice.sos, rice.eos, cnt.rec, rsq))
        last.rice <- rice.eos
      }
      # else {
      #   new.evi <- data.frame(x=newx$x, y=new.evi)
      # }
      #
    } else next
  }
  return(dat.rice)
}


