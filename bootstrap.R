## sample size
N <- 1000
## sample data
x <- rnorm(n=N, mean=0, sd=0.1)
## the empirical distribution functions visualised
hist(x)
## estimator t is mean
mean(x)
## number of bootstrap samples R
R <- 200
## first for a single sample
## random indices sampled from 1 to N with replacement
ii <- sample.int(N, replace=TRUE)
## the single bootstrap sample
xstar <- x[ii]

## no we generate R samples
## first allocate some memory in a container
xstar <- array(NA, dim=c(R, N))
## now loop over the samples
for(r in c(1:R)) {
  ii <- sample.int(N, replace=TRUE)
  xstar[r,] <- x[ii]
}
## any of the container not filled yet?
any(is.na(xstar))
## bootstrap replicates of the mean
xbarstar <- apply(xstar, 1, mean)
## the bootsrap estimate of the error of the mean
deltax <- sd(xbarstar)
cat(deltax, "\n")
## for comparison the standard error from the sample
cat(sd(x)/sqrt(N), "\n")
