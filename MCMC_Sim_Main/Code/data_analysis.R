args = commandArgs(trailingOnly=TRUE)

filename <- ''

if (length(args)==0) {
  stop('One argument must be supplied - (input filename).csv', call.=FALSE)
} else if (length(args)==1) {
  filename <- args[1]
} else{
  stop('Too many arguments supplied!', call.=FALSE)
}

data <- file(filename, open ='r')
line1 <- readLines(data)[[1]]
msp <- as.list(strsplit(line1, ',')[[1]])
invisible(unlist(msp))

mspvector <- c()
for (i in 1:length(msp)) {
  mspvector <- c(mspvector, as.numeric(msp[[i]]))
}

acorr <- acf(mspvector, length(mspvector), plot=F)
k <- 1
repeat{
    if(acorr$acf[k] < 0){
        break
    }
    k = k+1
}
acorr <- acf(mspvector, k, plot=F)

plot(acorr, main = 'Autocorrelation of mean squared position values', xlab  = 'Monte Carlo Time')

repeat{

  acorr <- acf(mspvector, k, plot=F)
  acfdata <- data.frame(lag=acorr$lag, acf=acorr$acf)
  acorr.log <- log(acfdata[,2])
  x <- seq(length(acorr.log))
  #loess_fit <- loess(acorr.log ~ x)
  fit <- lm(acorr.log ~ x)
  if((summary(fit)$r.squared) > 0.99){
        break
  }
  k = k - 1
}

acorr <- acf(mspvector, k, plot=F)
acfdata <- data.frame(lag=acorr$lag, acf=acorr$acf)
acorr.log <- log(acfdata[,2])
x <- seq(length(acorr.log))
fit <- lm(acorr.log ~ x)
plot(acorr.log, type='o' , pch=1, cex=.2 , col = 'steelblue2', main = 'Log autocorrelation of mean squared position values', ylab = 'Log(ACF)', xlab  = 'Monte Carlo Time')
abline(coef(fit)[1:2])
cf <- round(coef(fit), 2)
eq <- paste0("log(acorr) = ", cf[1], ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " x ")
mtext(eq, 3, line=-2)

n <- length(mspvector)
B <- 100
boot.samples = matrix(sample(mspvector,size=n*B,replace=TRUE), B, n)
boot.statistics = apply(boot.samples,1,mean)
se = sd(boot.statistics)
print(se)

#lines(predict(loess_fit ), col='red', lwd=2)
# write.table(acfdata, 'acfdata.txt', sep='\t')
