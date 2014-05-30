#
# 20130305: Sparse Matrix Version
#


Trace <- function( mat ){ # mat should be a square matrix 
  return( sum( diag( mat ) ) )
}

PSNR <- function( ref, tst ){
    mse <- mean( (tst-ref)**2 )
    rng2 <- (max(ref) - min(ref))**2
    return ( 10 * log(rng2/mse, base=10) )
}

logit <- function( u ){
    return( 0.5 * log(u/(1-u)) )
}

ilogit <- function( x ){
    return( exp(x)/(exp(x)+exp(-x)) )
}

dlogit <- function( u ){
    return( 0.5 * 1/(u*(1-u)) )
}

dilogit <- function( x ){
    return( 2 / ( (exp(x)+exp(-x))**2 ) )
}


NBinom2D <- function( nn, lmdmat, Ksamples=10000 ){
  Ly <- nrow(lmdmat)
  Lx <- ncol(lmdmat)
  LL <- Lx*Ly
  gennum <- rbinom( n = nn * LL, K, lmdmat/K )
  PData <- aperm( array( gennum, dim=c(Lx, Ly, nn) ), perm=c(1,2,3) )
  return( PData )
}


#MyLambda <- function( Lx, Ly ){ # See preliminary result by Lmd.R if you have any questions.
#  sdiag1 <- (rep( c( rep(1, Lx-1), 0 ), Ly ))[1:(Lx*Ly-1)]
#  sdiag2 <- rep( 1, Lx * (Ly-1) )
#  submat1 <- rbind( cbind( rep(0,length(sdiag1)), diag( sdiag1 ) ),
#                   rep( 0, length(sdiag1)+1 ) )
#  submat2 <- rbind( cbind( matrix( 0, nrow=length(sdiag2), ncol=Lx ), diag( sdiag2 ) ),
#                   matrix( 0, nrow=Lx, ncol=Lx*Ly ) )
#  subtmp <- submat1+submat2+t(submat1)+t(submat2)
#  maind <- colSums( subtmp )
#  return( diag( maind ) - subtmp )
#}


MyLambda <- function( Lx, Ly, h=1, v=1 ){
    Ah <- h
    if( length(Ah) == 1 ){
        Ah <- rep( Ah, Ly*(Lx-1) )
    }
    if( length(Ah) != Ly*(Lx-1) ){
        stop( 'num. of horiz. components is wrong' )
    }
    
    Av <- v
    if( length(Av) == 1 ){
        Av <- rep( 1, (Ly-1)*Lx )
    }
    if( length(Av) != (Ly-1)*Lx ){
        stop( 'num. of vert. components is wrong' )
    }
    return( Atilde.csr( Ah, Av, Ly, Lx ) )
}



#
# Generate Adjacent Matrix
# Ah means graph weigths corresponding to the horizontal connection.
# Av means graph weigths corresponding to the vertical connections.
# Prow means horizontal size of image matrix
# Pcol means vertical size of image matrix
#
Atilde <- function( Ah, Av, Prow, Pcol ){
  sdiag1 <- array( rbind( matrix( Av, nrow=Prow-1 ), rep(0, Pcol) ))[1:(Prow*Pcol-1)]
  sdiag2 <- Ah
  submat1 <- rbind( cbind( rep(0,length(sdiag1)), diag( sdiag1 ) ),
                   rep( 0, length(sdiag1)+1 ) )
#  submat2 <- rbind( cbind( matrix( 0, nrow=length(sdiag2), ncol=Pcol ), diag( sdiag2 ) ),
#                   matrix( 0, nrow=Pcol, ncol=Pcol*Prow ) )
  submat2 <- rbind( cbind( matrix( 0, nrow=length(sdiag2), ncol=Prow ), diag( sdiag2 ) ),
                   matrix( 0, nrow=Prow, ncol=Pcol*Prow ) )
  subtmp <- submat1+submat2+t(submat1)+t(submat2)
  maind <- colSums( subtmp )
  return( diag( maind ) - subtmp )
}


Atilde.csr <- function( Ah, Av, Prow, Pcol ){
  sdiag1 <- array( rbind( matrix( Av, nrow=Prow-1 ), rep(0, Pcol) ))[1:(Prow*Pcol-1)]
  sdiag2 <- Ah
  submat1 <- rBind( cBind( rep(0,length(sdiag1)), Diagonal( x=sdiag1 ) ),
                   rep( 0, length(sdiag1)+1 ) )
  submat2 <- rBind(
                   cBind( Matrix( 0, nrow=length(sdiag2), ncol=Prow, sparse=TRUE ), Diagonal( x=sdiag2 ) ),
                   Matrix( 0, nrow=Prow, ncol=Pcol*Prow, sparse=TRUE ) )
  subtmp <- submat1+submat2+t(submat1)+t(submat2)
  maind <- colSums( subtmp )
  return( diag( maind ) - subtmp )
}
