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
unlist(msp)

mspvector <- c()
for (i in 1:length(msp)) {
  mspvector <- c(mspvector, as.numeric(msp[[i]]))
}

acorr <- acf(mspvector, 200, plot=F)

plot(acorr, main = 'Autocorrelation of mean squared position values', xlab  = 'Monte Carlo Time')

acfdata <- data.frame(lag=acorr$lag, acf=acorr$acf)

acorr.log <- log(acfdata[,2])

plot(acorr.log, main = 'Log autocorrelation of mean squared position values', ylab = 'Log ACF', xlab  = 'Monte Carlo Time')

write.table(acfdata, 'acfdata.txt', sep='\t')
