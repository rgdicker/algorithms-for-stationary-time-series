# Here I implement the Durbin-Levinson algorithm for solving Toeplitz systems and
# compare it to my blockPCG implementation to see which is faster and which is 
# more precise.

# It looks like blockPCG is 5 to 10 times faster than Durbin-levinson, and Durbin-Levinson is more accurate.
# However, the blockPCG method results in a solution whose distance from the true solution is less
# than 10^(-18), so increased accuracy is not useful.

# Note: blockPCG is faster because the FFT method of matrix-vector multiplication is used.
# The standard method of blockPCG is much slower than Durbin-Levinson.



## Functions ------
require(ltsa)
require(pracma)
# Durbin-Levinson algorithm
dlev <- function(A,y){
    r <- as.vector(A[1,])[-1]
    d <- 1
    b <- -r[1]
    x <- y[1]
    n <- dim(A)[1]
    # Iterations
    for (k in 2:n) {
        # Compute Scalars
        d <- 1 + t(r[1:(k-1)])%*%b
        e <- r[k] + t(b)%*%r[(k-1):1]
        f <- y[k] - t(r[1:(k-1)])%*%x[(k-1):1]
        # Update Vectors
        b[k] <- -e/d
        x[k] <- f/d
        x[1:(k-1)] <- x[1:(k-1)]+x[k]*b[(k-1):1]
        b[1:(k-1)] <- b[1:(k-1)]+b[k]*b[(k-1):1]
    }
    x
}
# Preconditioned Conjugate Gradient Algorithm
blockPCG <- function(A,b,numits,blocksize) {
    B <- solve(A[1:blocksize,1:blocksize])
    n <- dim(A)[1]
    blockmult <- function(Block,b) {
        xt <- vector(length = n)
        for(j in 1:(n/blocksize)) {
            part <- ((j-1)*blocksize+1):(j*blocksize)
            xt[part] <- Block%*%b[part]
        } 
        xt
    }
    toepmult <- function(A,v) {
        n <- dim(A)[1]
        x <- as.matrix(c(A[1,],0,A[1,][n:2]))
        p <- c(v,rep(0,n))
        f <- fft(p)
        g <- fft(x)
        h <- as.vector(f*g)
        z <- Re(ifft(h)[1:n])
        z
    }
    #initial step
    X <- matrix(0,nrow = n, ncol = numits+1)
    r <- b
    z <- blockmult(B,r)#solve(M,r)
    p <- z
    #iterations
    for (i in 1:numits) {
        Ap <- toepmult(A,p)#A %*% p
        rz <- as.numeric(crossprod(r,z))
        a <- rz/as.numeric(crossprod(p,Ap))
        X[,(i+1)] <- X[,i] + a*p
        r <- r - a*(Ap)
        z <- blockmult(B,r)#solve(M,r)
        c <- as.numeric(crossprod(r,z))/rz
        p <- z + c*p
    }
    X
}

resids <- function(X,A,b) {
    resids <- 0
    x <- TrenchInverse(A)%*%b
    for (i in 1:dim(X)[2]) {
        resids[i] <- sqrt(crossprod(X[,i]-x))
    }
    resids
}


### Make Matrices ----
h <- 0:(4000-1)
#one
one <- (h/4)*besselK(h/4,1)
one[1] <- 1
C1 <- toeplitz(one)


# two
two <- (1+abs(h)/4)*exp(-abs(h)/4)
C2 <- toeplitz(two)


# three
three <- exp(-(h/2)^2)
C3 <- toeplitz(three)


# four 
four <- c(1)
for ( i in 2:length(h) ) {
    four[i] <- ((h[i]-1+0.4)/(h[i]-0.4))*four[i-1]
}
C4 <- toeplitz(four)

# Solve system 1 ----
set.seed(300)
A <- C1
L <- t(chol(A))
z <- rnorm(4000)
b1 <- L%*%z

# Run PCG with 50 iterations
its <- 50
X <- blockPCG(A,b1,its,16)
#Run dlev
dlev_x <- dlev(A,b1)
# Check if PCG and dlev got close enough to the solution 
actual_x <- as.numeric(TrenchInverse(A)%*%b1)
PCG_last_iter1 <- X[,its]
dlev_resid <- crossprod(dlev_x-actual_x)
PCG_resid <- crossprod(PCG_last_iter1-actual_x)
resid_ratio <- PCG_resid/dlev_resid
resids1 <- c(dlev_resid,PCG_resid,resid_ratio)
resids1

require(microbenchmark)
# Check speeds
bench1 <- microbenchmark(dlev(A,b1),blockPCG(A,b1,its,16), times = 10)
summary(bench1)
medians <- summary(bench1)[5]
times_faster <- medians[1,]/medians[2,]
round(times_faster,1)

# When I ran the benchmark, blockPCG was about 10 times faster than dlev, and both functions were able to produce a very accurate solution. 
# The residual norms for both methods were less than 1e-20, but the residual norm from d_lev is about 10 times smaller than for PCG


# Solve system 2 ----
set.seed(300)
A <- C2
L <- t(chol(A))
z <- rnorm(4000)
b2 <- L%*%z

# Run PCG with 50 iterations
its <- 50
X <- blockPCG(A,b2,its,20)
#Run dlev
dlev_x <- dlev(A,b2)
# Check if PCG and dlev got close enough to the solution 
actual_x <- as.numeric(TrenchInverse(A)%*%b2)
PCG_last_iter2 <- X[,its]
dlev_resid <- crossprod(dlev_x-actual_x)
PCG_resid <- crossprod(PCG_last_iter2-actual_x)
resid_ratio <- PCG_resid/dlev_resid
resids2 <- c(dlev_resid,PCG_resid,resid_ratio)
resids2

#require(microbenchmark)
# Check speeds
bench2 <- microbenchmark(dlev(A,b2),blockPCG(A,b2,its,16), times = 10)
summary(bench2)
medians <- summary(bench2)[5]
times_faster <- medians[1,]/medians[2,]
round(times_faster,1)

# When I ran the benchmark, blockPCG was about 11 times faster than dlev, and both functions were able to produce a very accurate solution. 
# The residual norms for both methods were less than 1e-18, but the residual norm from d_lev is about 42 times smaller than for PCG


# Solve system 3 ----
set.seed(300)
A <- C3
L <- t(chol(A))
z <- rnorm(4000)
b3 <- L%*%z

# Run PCG with 50 iterations
its <- 50
X <- blockPCG(A,b3,its,20)
#Run dlev
dlev_x <- dlev(A,b3)
# Check if PCG and dlev got close enough to the solution 
actual_x <- as.numeric(TrenchInverse(A)%*%b3)
PCG_last_iter3 <- X[,its]
dlev_resid <- crossprod(dlev_x-actual_x)
PCG_resid <- crossprod(PCG_last_iter3-actual_x)
resid_ratio <- PCG_resid/dlev_resid
resids3 <- c(dlev_resid,PCG_resid,resid_ratio)
resids3

#require(microbenchmark)
# Check speeds
bench3 <- microbenchmark(dlev(A,b3),blockPCG(A,b3,its,16), times = 10)
summary(bench3)
medians <- summary(bench3)[5]
times_faster <- medians[1,]/medians[2,]
round(times_faster,1)

# When I ran the benchmark, blockPCG was about 10 times faster than dlev, and both functions were able to produce a very accurate solution. 
# The residual norms for both methods were less than 1e-18, but the residual norm from d_lev is about 119 times smaller than for PCG


# Solve system 4 ----
set.seed(300)
A <- C4
L <- t(chol(A))
z <- rnorm(4000)
b4 <- L%*%z

# Run PCG with 50 iterations
its <- 100
X <- blockPCG(A,b4,its,50)
#Run dlev
dlev_x <- dlev(A,b4)
# Check if PCG and dlev got close enough to the solution 
actual_x <- as.numeric(TrenchInverse(A)%*%b4)
PCG_last_iter4 <- X[,its]
dlev_resid <- crossprod(dlev_x-actual_x)
PCG_resid <- crossprod(PCG_last_iter4-actual_x)
resid_ratio <- PCG_resid/dlev_resid
resids4 <- c(dlev_resid,PCG_resid,resid_ratio)
resids4

#require(microbenchmark)
# Check speeds
bench4 <- microbenchmark(dlev(A,b4),blockPCG(A,b4,its,16), times = 10)
summary(bench4)
medians <- summary(bench4)[5]
times_faster <- medians[1,]/medians[2,]
round(times_faster,1)

# When I ran the benchmark, blockPCG was about 5 times faster than dlev, and both functions were able to produce a very accurate solution. 
# The residual norms for both methods were less than 1e-22, and the ratio of the residual norms from PCG and dlev is approximately 1

## Compare to Cholesky Method ----
l2norm <- function(x) { sqrt( sum( abs(x)^2 ) ) }

cholsolve <- function(A,b) {
    actual_x <- TrenchInverse(A)%*%b
    L <- t(chol(A))
    Linv <- solve(L)
    y <- Linv%*%b
    x <- t(Linv)%*%y
    x
}


Chol_Solutions <- cbind(cholsolve(C1,b1), cholsolve(C2,b2) , cholsolve(C3,b3) , cholsolve(C4,b4) )
PCG_solutions <- cbind(PCG_last_iter1, PCG_last_iter2, PCG_last_iter3, PCG_last_iter4)
diff <- PCG_solutions - Chol_Solutions
l2norms <- c()
for (i in 1:4) {
    l2norms[i] <- l2norm(diff[,i])
}
names(l2norms) <- c('system1','system2','system3','system4')
l2norms


# The 2 norm of the difference between the Durbin- Levinson solutions and the Cholesky-based solutions are given in the l2norms vector. 

times <- cbind(summary(bench1)[5],summary(bench2)[5],summary(bench3)[5],summary(bench4)[5])
row.names(times) <- c('dlev','PCG')
names(times) <- names(l2norms)
times
# Median time (in milliseconds) required to solve each system is given in the times vector
