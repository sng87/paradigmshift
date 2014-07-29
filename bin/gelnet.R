dyn.load( "~/bin/paradigmshift/bin/gelnet.so" )

X = as.matrix(read.table('data.X.matrix', header = TRUE, row.names = 1))
y = read.table('class.y.vector', header = TRUE, row.names = 1)$class
d = read.table('weight.d.vector', header = TRUE, row.names = 1)$weight
P = as.matrix(read.table('penalty.P.matrix', header = TRUE, row.names = 1))
l1 = 0.05 # try sweeping these parameters
l2 = 1.0


## Constructs a logistic regression GELNET
##
## Inputs:
## X - n-by-p matrix of n samples in p dimensions
## y - n-by-1 vector of binary response labels
## d - p-by-1 vector of feature weights
## P - p-by-p feature association penalty matrix
## vL1, vL2 - two k-by-1 vectors of k pairwise L1- and L2-norm penalties
## max.iter - maximum number of iterations
## eps - convergence precision
## w.init, b.init - the initial parameter estimates
##
## Output:
##   w - p-by-k matrix of p model weights for the k regularization points
##   b - k-by-1 vector of bias terms for the k regularization points
gelnet.logreg <- function( X, y, d, P, vL1, vL2,
                          max.iter = 100, eps = 1e-5,
                          w.init = rep(0,p), b.init = 0.5 )
  {
    n <- nrow(X)
    p <- ncol(X)
    k <- length( vL1 )

    ## Verify arguments
    stopifnot( length( unique(y) ) == 2 )
    stopifnot( length(y) == n )
    stopifnot( length(d) == p )
    stopifnot( all( dim(P) == c(p,p) ) )
    if( is.null( colnames(X) ) == FALSE )
      {
        stopifnot( is.null( rownames(P) ) == FALSE )
        stopifnot( is.null( colnames(P) ) == FALSE )
        stopifnot( all( colnames(X) == rownames(P) ) )
        stopifnot( all( colnames(X) == colnames(P) ) )
      }
    stopifnot( length( vL2 ) == k )
    stopifnot( length( w.init ) == p )
    stopifnot( length( b.init ) == 1 )

    ## Convert the labels to {0,1}
    y <- as.integer( y == max(y) )
    
    ## Set the initial parameter estimates
    w <- matrix( w.init, p, k )
    b <- rep( b.init, length=k )
    S <- X %*% w[,1] + b[1]
    Pw <- P %*% w[,1]

    res <- .C( "gelnet_logreg_path",
              as.double(X), as.integer(y), as.double(d), as.double(P),
              as.double(vL1), as.double(vL2), as.integer(k),
              as.double(S), as.double(Pw), as.integer(n), as.integer(p),
              as.integer(max.iter), as.double(eps),
              w = as.double( w ), b = as.double(b) )

    res1 <- list()
    res1$w <- matrix( res$w, p, k )
    res1$b <- res$b

    rownames( res1$w ) <- colnames(X)
    
    res1
  }

gelvec = gelnet.logreg(X, y, d, P, l1, l2)
write.table(gelvec$w[,1], 'gelnet.w.vector', sep = '\t', quote = FALSE, col.names = FALSE)
