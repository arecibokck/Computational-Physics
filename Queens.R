S <- function(x) {
  N <- length(x)
  V <- rep(FALSE, times=N)
  counter <- 0
  for(i in c(1:(N-1))) {
    for(j in c((i+1):N)) {
      if((x[i] == x[j]) || (abs(x[i]-x[j]) == j-i)) {
        V[i] <- TRUE
        V[j] <- TRUE
      }
    }
    if(V[i]) counter <- counter + 1 
  }
  if(V[N]) counter <- counter + 1
  return(counter)
}

## alternatively S may be written more in R-style
S2 <- function(x) {
  N <- length(x)
  rows <- c(0:7)
  flag <- rep(FALSE, times=N)
  for(i in c(1:(N-1))) {
    ## equal col
    if(length(ii<-which(x[i] == x)) > 1) {
      flag[ii] <- TRUE
    }
    ## diagonals
    else if(length(ii <- which(abs(x[c(i:8)]-x[i]) == rows[c(1:(8-i+1))] )) > 1) {
      flag[ii + i - 1] <- TRUE
    }
  }
  if(any(flag)) return(length(which(flag))-1)
  return(0)
}


SA <- function(x, T=100, alpha=0.9) {
  Sold <- S(x)
  N <- length(x)
  y <- x
  k <- 0
  while(Sold > 0) {
    for(i in c(1:N)) {
      y[i] <- ceiling(runif(n=1, max=N))
      Snew <- S(y)
      accept <- FALSE
      if(Snew < Sold) accept <- TRUE
      else if(runif(1) < exp((Sold-Snew)/T)) accept <- TRUE
      if(accept) {
        x[i] <- y[i]
        Sold <- Snew
      }
      else {
        y[i] <- x[i]
      }
    }
    T <- alpha*T
    k <- k+1
  }
  cat(k, Sold, "\n")
  
  return(x)
}

SA2 <- function(x, T=100, alpha=0.9) {
  Sold <- S(x)
  N <- length(x)
  k <- 0
  while(Sold > 0) {
    y <- sample.int(N, N)
    Snew <- S(y)
    accept <- FALSE
    if(Snew < Sold) accept <- TRUE
    else if(runif(1) < exp((Sold-Snew)/T)) accept <- TRUE
    if(accept) {
      x <- y
      Sold <- Snew
    }
    T <- alpha*T
    k <- k+1
  }
  cat(k, Sold, "\n")
  
  return(x)
}
