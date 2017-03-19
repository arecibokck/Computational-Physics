args = commandArgs(trailingOnly=TRUE)

filename = ''

if (length(args)==0) {
  stop('One argument must be supplied - (input filename).csv', call.=FALSE)
} else if (length(args)==1) {
  filename <- args[1]
} else{
  stop('Too many arguments supplied!', call.=FALSE)
}

data <- file(filename, open ='r')
line1 <- readLines(data)[[1]]
energy <- as.list(strsplit(line1, ',')[[1]])
invisible(unlist(energy))

energyvector <- c()
for (i in 1:length(energy)) {
  num <- as.numeric(energy[[i]])
  if(!is.nan(num)){
      energyvector <- c(energyvector,num)
  }
}

n <- length(energyvector)
B <- 100
boot.samples = matrix(sample(energyvector,size=n*B,replace=TRUE), B, n)
boot.statistics = apply(boot.samples,1,mean)
se = sd(boot.statistics)
print(se)
