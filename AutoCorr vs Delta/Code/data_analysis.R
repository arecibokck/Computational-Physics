args = commandArgs(trailingOnly=TRUE)

filename <- args[1]

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

    if(acorr$acf[k] < 0 ){
        break
    }
    k = k+1
}

repeat{

  acorr <- acf(mspvector, k, plot=F)
  acfdata <- data.frame(lag=acorr$lag, acf=acorr$acf)
  acorr.log <- log(acfdata[,2])
  x <- seq(length(acorr.log))
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
if((summary(fit)$r.squared) > 0.90){
	cf <- round(coef(fit), 3)
	print(cf[2])
}
