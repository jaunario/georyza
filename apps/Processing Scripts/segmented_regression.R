x <- as.Date(colnames(bgd_gt_evi_2009305.2011081)[-(1:3)],"A%Y%j")
y <- t(bgd_gt_evi_2009305.2011081[,-(1:3)])
y <- as.data.frame(y)
ipr <- mapply(polregIT.rice, pix.evi = y, evi.date = data.frame(x), dates.covered = data.frame(x), alpha = 0.1)
bgd_gt_evi_2009305.2011081$newrice <- sapply(ipr, length)
aa <- bgd_gt_evi_2009305.2011081[,c(1:3,69)]


dat.evi <- data.frame(x,y)
dat.evi <- data.frame(x=mat.jday[1,],y=mat.evi[1,])

polregIT.rice(pix.evi=y[,67], evi.date = x, dates.covered = x, alpha=0.05)
ipr <- mapply(polregIT.rice, pix.evi = mat.evi, evi.date = mat.jday, dates.covered = data.frame(date=date.acqdates), alpha = 0.1)

fft(dat.evi$y)

ipr <- mapply(polregIT.rice, pix.evi = mat.evi, evi.date = mat.jday, dates.covered = data.frame(date=date.acqdates), alpha = 0.1)
stk.evi[3310]


seg1.lm <- lm(y~poly(x,2), data = dat.evi[1:30,])
dat.evi <- data.frame(x,y)
dat.evi <- na.omit(dat.evi)
plot(dat.evi)
lines(y=signal::sgolayfilt(x=dat.evi$y,ts = dat.evi$x),x=dat.evi$x)
x <- mat.jday[1,]
y <- mat.evi[1,]

rm.noise <- function(x,y)
  
  fillx <- approx(x,y, xout=x)
  dy <- diff(fillx$y) 
  decr <- c(FALSE, dy<0)
  cleany <- vector()
  for(i in 1:nrow(fillx)){
    if(i==1) {
      cleany <- fillx$y[i]
    } else {
      if(cleany[i-1]<fillx$y[i]) decr <- TRUE
    }
    
  }
  repeat{
  dify <- diff(y)
  grp.neg <- consecutive.groups(which(dify<0))
  singles <- which(sapply(grp.neg, length)==1)
  if(length(singles)>0) {
  y[unlist(grp.neg[singles])] <- NA
  y <- round(approx(y, x=x, xout=x)$y)
  
  dify <- diff(y)
  grp.pos <- consecutive.groups(which(dify>0))
  singles <- which(sapply(grp.pos, length)==1)
  singles <- singles[singles!=1]
  y[unlist(grp.pos[singles])] <- NA
  y <- round(approx(y, x=x, xout=x)$y)
  } else break
}

repeat{
  dify <- diff(y)
  grp.neg <- consecutive.groups(which(dify>0))
  singles <- which(sapply(grp.neg, length)==1)
  if(length(singles)>0) {
    y[unlist(grp.neg[singles])] <- NA
    y <- round(approx(y, x=x, xout=x)$y)
  } else break
}


y_clean <- y[1]
y <- y[-1]
while (length(y)>0){
  if(is.na(y[1])){
    y_clean <- c(y_clean, y[i])
  } else {
    dify1 <- y[-1]-y[1]
    dify[1]>0
    
  }  
  
}

diff(mat.evi[,1])

ts(y, frequency = 46, start = c(YEAR-1,36))
                  segment.size=8
segment.increment = 2
dat <- data.frame(x,y=y/10000)
dat <- na.omit(dat)

y_hat <- mean(dat$y, na.rm=TRUE)
consecutive.groups(which(dat$y-y_hat<0))
alpha = 0.05
i = 1 # Iterator
seg.start=1
sections <- vector()
section.models <- list()
repeat {
  inc.seq <-seg.start:((seg.start+segment.size-1)+(i-1)*segment.increment)
  thisseg <- data.frame(x=x[inc.seq], y=y[inc.seq])
  segment.lm <- lm(y~x, thisseg)
  results.seglm <- summary(segment.lm)
  fstat <- results.seglm$fstatistic
  pval <- pf(fstat[1],fstat[2], fstat[3], lower.tail = FALSE)
  if(pval<=alpha) {
    i <- i+1 
    prev.result <- segment.lm
  } else {
    last.inc <- (max(inc.seq)-(segment.increment-1)):max(inc.seq)
    sections <- rbind(sections, c(seg.start, max(inc.seq[!inc.seq %in% last.inc])))
    seg.start <- last.inc
    section.models[[nrow(sections)]] <- prev.result
    i <- 1
  }
  if(max(inc.seq)==nrow(dat)) {
    sections <- rbind(sections, c(seg.start, max(inc.seq)))
    section.models[[nrow(sections)]] <- segment.lm
    break
  }
}

# Fit an overall linear regression
lm1 <- lm(y~x,dat)
# If p-value 
library(segmented)
o <- segmented(lm1, seg.Z = ~x, psi = list(x=c(dat$x[28],dat$x[37])), control = seg.control(display = FALSE))
summary(o)
segment.lm$effects


