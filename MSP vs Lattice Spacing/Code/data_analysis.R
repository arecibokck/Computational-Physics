args = commandArgs(trailingOnly=TRUE)

filename <- ''

if (length(args)==0) {
  stop('One argument must be supplied - (input filename).csv', call.=FALSE)
} else if (length(args)==1) {
  msp <- args[1]
} else{
  stop('Too many arguments supplied!', call.=FALSE)
}

msp <- unlist(strsplit(msp, split=" "))
mspvector <- c()
for (i in 1:length(msp)) {
  num <- as.numeric(msp[[i]])
  if(!is.nan(num)){
      mspvector <- c(mspvector,num)
  }
}

n <- length(mspvector)
B <- 100
boot.samples = matrix(sample(mspvector,size=n*B,replace=TRUE), B, n)
boot.statistics = apply(boot.samples,1,mean)
se = sd(boot.statistics)
print(se)
