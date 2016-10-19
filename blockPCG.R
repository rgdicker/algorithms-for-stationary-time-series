# This function solves a Toeplitz system Ax=b using the Preconditioned Conjugate Gradient method, 
# where the preconditioner is a block diagonal version of A. It is 40x faster than using solve(A,b).
# Since covariance matrices of stationary time series are Toeplitz, this method can be used to 
# speed up forecasting and maximum likelihood estimation in the time series context.


blockPCG <- function(A,b,numits,blocksize) {
    B <- solve(A[1:blocksize,1:blocksize])
    n <- dim(A)[1]
    blocktoep_solve <- function(Blockinv,b) {
        xt <- vector(length = n)
        for(j in 1:(n/blocksize)) {
            part <- ((j-1)*blocksize+1):(j*blocksize)
            xt[part] <- Blockinv%*%b[part]
        } 
        xt
    }
    #initial step
    X <- matrix(0,nrow = n, ncol = numits+1)
    r <- b
    z <- blocktoep_solve(B,r)#solve(M,r)
    p <- z
    #iterations
    for (i in 1:numits) {
        rz <- as.numeric(crossprod(r,z))
        a <- rz/as.numeric(crossprod(p,A%*%p))
        X[,(i+1)] <- X[,i] + a*p
        r <- r - a*(A%*%p)
        z <- blocktoep_solve(B,r)#solve(M,r)
        c <- as.numeric(crossprod(r,z))/rz
        p <- z + c*p
    }
    X
}


# Test this when A is a Matern covariance matrix and b is a random vector from a normal distribution with covariance A

# make matrix
h <- 0:(4000-1)
one <- (h/4)*besselK(h/4,1)
one[1] <- 1
C1 <- toeplitz(one)
# generate b
L1 <- t(chol(A))
z <- rnorm(4000)
b <- L1%*%z

# To evaluate the performance of the algorithm, I run it using various block sizes and look at the residuals over each iteration

# Create residual calculator

require(ltsa)
resids <- function(X,A,b) {
    resids <- 0
    x <- TrenchInverse(A)%*%b
    for (i in 1:dim(X)[2]) {
        resids[i] <- sqrt(crossprod(X[,i]-x))
    }
    resids
}

# Run PCG and get residuals
num_its <- 50
sizes <- c(1,2,4,8,20)
resids1 <- matrix(nrow = num_its+1, ncol = 5)
for (i in 1:5) {
    X <- blockPCG(A,b,num_its,sizes[i])
    resids1[,i] <- resids(X,A,b)
}


# Plot the residuals

require(ggplot2)
require(RColorBrewer)
col <- brewer.pal(4, "Set1")

p1 <- qplot(1:dim(resids1)[1],resids1[,1], ylab = 'Residual', xlab = 'Iteration') + geom_point(aes(y=resids1[,2]), color = col[1]) + geom_point(aes(y=resids1[,3]), color = col[2]) + geom_point(aes(y=resids1[,4]), color = col[3]) + geom_point(aes(y=resids1[,5]), color = col[4]) + ggtitle(expression(gamma(h)==(abs(h)/4)*K[1](abs(h)/4))) 
p1


# The algorithm appears to get to the solution fairly quickly in terms of the number of iterations required, 
# but the time required to evaluate each iteration makes it quite slow.

# Let's use profvis to figure out what parts are slowing it down
install.packages("profvis")
library(profvis)



profvis({
    blocksize <- 8; numits <- 50
    B <- solve(A[1:blocksize,1:blocksize])
    n <- dim(A)[1]
    blocktoep_solve <- function(Blockinv,b) {
        xt <- vector(length = n)
        for(j in 1:(n/blocksize)) {
            part <- ((j-1)*blocksize+1):(j*blocksize)
            xt[part] <- Blockinv%*%b[part]
        } 
        xt
    }
    #initial step
    X <- matrix(0,nrow = n, ncol = numits+1)
    r <- b
    z <- blocktoep_solve(B,r)#solve(M,r)
    p <- z
    #iterations
    for (i in 1:numits) {
        rz <- as.numeric(crossprod(r,z))
        a <- rz/as.numeric(crossprod(p,A%*%p))
        X[,(i+1)] <- X[,i] + a*p
        r <- r - a*(A%*%p)
        z <- blocktoep_solve(B,r)#solve(M,r)
        c <- as.numeric(crossprod(r,z))/rz
        p <- z + c*p
    }
    X
})


# Looks like the vast majority of time is spent calculating A%*%p. Since A is a Toeplitz matrix, 
# we can put it in the corner of a Circulant matrix and use Fast Fourier Transforms to get the product, 
# rather than doing the multiplications directly.


require(pracma) # for the ifft  (inverse fast fourier transform) function
blockPCG2 <- function(A,b,numits,blocksize) {
    B <- solve(A[1:blocksize,1:blocksize])
    n <- dim(A)[1]
    blocktoep_solve <- function(Blockinv,b) {
        xt <- vector(length = n)
        for(j in 1:(n/blocksize)) {
            part <- ((j-1)*blocksize+1):(j*blocksize)
            xt[part] <- Blockinv%*%b[part]
        } 
        xt
    }
    # New Toeplitz-specific multiplication function
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
    z <- blocktoep_solve(B,r)#solve(M,r)
    p <- z
    #iterations
    for (i in 1:numits) {
        Ap <- toepmult(A,p)#A%*%p
        rz <- as.numeric(crossprod(r,z))
        a <- rz/as.numeric(crossprod(p,Ap))
        X[,(i+1)] <- X[,i] + a*p
        r <- r - a*(Ap)
        z <- blocktoep_solve(B,r)#solve(M,r)
        c <- as.numeric(crossprod(r,z))/rz
        p <- z + c*p
    }
    X
}

# Let's use the microbenchmark to see how much faster the new function is
# Note: Since the microbenchmark function works by running each algorithm multiple times, 
# the function may take a long time (up to 2 mins) to run
require(microbenchmark)
bench <- microbenchmark(blockPCG(A,b,50,8),blockPCG2(A,b,50,8), solve(A,b), times = 10)
medians <- summary(bench)[5]
medians
times_faster <- medians[1,]/medians[2,]
times_faster
times_faster_vs_solve <- medians[3,]/medians[2,]
times_faster_vs_solve
# The outcome of this benchmark will of course vary according to the properties of the machine on which it runs, 
# but it looks like avoiding the matrix-vector product calculation A%*%p makes the algorithm about 20x faster! 
# The original algorithm was already faster than the built-in solve() function, so quickPCG2 is about 40x faster than solve()!


# Let's look at how this algorithm performs for a variety of systems Ax=b

# Create a second covariance matrix and random normal vector

# make matrix
two <- (1+abs(h)/4)*exp(-abs(h)/4)
A <- toeplitz(two)
# generate b
L1 <- t(chol(A))
z <- rnorm(4000)
b <- L1%*%z

# Run PCG and get residuals

sizes <- c(1,2,4,8,20)
resids2 <- matrix(nrow = num_its+1, ncol = 5)
for (i in 1:5) {
    X <- blockPCG2(A,b,num_its,sizes[i])
    resids2[,i] <- resids(X,A,b)
}

# Create a third covariance matrix and random normal vector

# make matrix
three <- exp(-(h/2)^2)
A <- toeplitz(three)
# generate b
L1 <- t(chol(A))
z <- rnorm(4000)
b <- L1%*%z

# Run PCG and get residuals

sizes <- c(1,2,4,8,20)
resids3 <- matrix(nrow = num_its+1, ncol = 5)
for (i in 1:5) {
    X <- blockPCG2(A,b,num_its,sizes[i])
    resids3[,i] <- resids(X,A,b)
}

# Create a fourth covariance matrix and random normal vector

# make matrix
four <- c(1)
for ( i in 2:length(h) ) {
    four[i] <- ((h[i]-1+0.4)/(h[i]-0.4))*four[i-1]
}
A <- toeplitz(four)
# generate b
L1 <- t(chol(A))
z <- rnorm(4000)
b <- L1%*%z

# Run PCG and get residuals

sizes <- c(1,2,4,8,20)
resids4 <- matrix(nrow = num_its+1, ncol = 5)
for (i in 1:5) {
    X <- blockPCG2(A,b,num_its,sizes[i])
    resids4[,i] <- resids(X,A,b)
}




# Create a plot showing the convergence for each of the four systems using various block sizes for the Preconditioner 

require(ggplot2)
require(RColorBrewer)
col <- brewer.pal(4, "Set1")

p1 <- qplot(1:dim(resids1)[1],resids1[,1], ylab = 'Residual', xlab = 'Iteration') + geom_point(aes(y=resids1[,2]), color = col[1]) + geom_point(aes(y=resids1[,3]), color = col[2]) + geom_point(aes(y=resids1[,4]), color = col[3]) + geom_point(aes(y=resids1[,5]), color = col[4]) + ggtitle(expression(gamma(h)==(abs(h)/4)*K[1](abs(h)/4))) 
p2 <- qplot(1:dim(resids2)[1],resids2[,1], ylab = 'Residual', xlab = 'Iteration') + geom_point(aes(y=resids2[,2]), color = col[1]) + geom_point(aes(y=resids2[,3]), color = col[2]) + geom_point(aes(y=resids2[,4]), color = col[3]) + geom_point(aes(y=resids2[,5]), color = col[4]) + ggtitle(expression(gamma(h)==(1+abs(h)/4)*exp(-abs(h)/4)))

p3 <- qplot(1:dim(resids3)[1],resids3[,1], ylab = 'Residual', xlab = 'Iteration') + geom_point(aes(y=resids3[,2]), color = col[1]) + geom_point(aes(y=resids3[,3]), color = col[2]) + geom_point(aes(y=resids3[,4]), color = col[3]) + geom_point(aes(y=resids3[,5]), color = col[4]) + ggtitle(expression(gamma(h)==exp(-(abs(h)/4)^2)))

p4 <- qplot(1:dim(resids4)[1],resids4[,1], ylab = 'Residual', xlab = 'Iteration') + geom_point(aes(y=resids4[,2]), color = col[1]) + geom_point(aes(y=resids4[,3]), color = col[2]) + geom_point(aes(y=resids4[,4]), color = col[3]) + geom_point(aes(y=resids4[,5]), color = col[4]) + ggtitle(expression(gamma(h)==(h-1+d)*gamma(h-1)/(h-d)))

# Make composite plot
require(gridExtra)
require(grid)
require(lattice)
t <- textGrob("Colors indicate size of blocks. Black = 1, Red = 2, Blue = 4, Green = 8, Purple = 20")


PLOT <- grid.arrange(arrangeGrob(p1,p2,p3,p4, ncol = 2), t,nrow = 2,heights = c(10,1))

PLOT
