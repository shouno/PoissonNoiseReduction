#
# Loopy Belief Propagation for 2D GMRF
#

library( MASS )


MkMsgLst <- function( Mpars, dim, names = NULL ){ # Make Message Array Mpars * dim
  # if you'd like to 2D Gaussian BP, call MkMsgArray( 2, c(Ly, Lx), c('mu', 'gam') )
  lst <- as.list( seq(1, Mpars) )
  lst <- lapply( lst, function(u) array(0, dim) )
  if( is.null( names ) == FALSE ){
    if( length(names) != Mpars ){
      stop( "names parameter are wrong" )
    }
    names( lst ) <- names
  }
  return( lst )
}


MyLambda2D <- function( Ly, Lx ){ # See preliminary result by Lmd.R if you have any questions.
  # for 1D chain
  #L1 <- rep( 1, L-1 )
  #Maind <- diag( c(L1, 0) ) + diag( c(0, L1) ) 
  #Subd <- rbind( rep( 0, length(L1)+1 ), cbind( diag(L1), rep(0, length(L1) ) ) )
  #return( Maind - Subd - t(Subd) )
#
  sdiag1 <- (rep( c( rep(1, Lx-1), 0 ), Ly ))[1:(Lx*Ly-1)]
  sdiag2 <- rep( 1, Lx * (Ly-1) )
  submat1 <- rbind( cbind( rep(0,length(sdiag1)), diag( sdiag1 ) ),
                   rep( 0, length(sdiag1)+1 ) )
  submat2 <- rbind( cbind( matrix( 0, nrow=length(sdiag2), ncol=Lx ), diag( sdiag2 ) ),
                   matrix( 0, nrow=Lx, ncol=Lx*Ly ) )
  subtmp <- submat1+submat2+t(submat1)+t(submat2)
  maind <- colSums( subtmp )
  return( diag( maind ) - subtmp )
}


InferX2D <- function( Ly, Lx, N, talpha, tbeta, th, alpha, beta, h, x0=0 ){
  M <- Lx*Ly
  Lmd <- MyLambda2D( Ly, Lx )
  Eye <- diag( rep( 1, M ) )
#
                                        # 多次元ガウス分布の真のパラメータ
  tC1 <- talpha * Lmd + th * Eye
  tm1 <- rep(x0, M)
#
  x <- matrix( mvrnorm( N, tm1, solve(tC1) ), nrow=N )
  obsrnd <- matrix( rnorm( n=N*M, mean=0, sd=sqrt( 1/tbeta ) ), nrow=N )
  y <- x + obsrnd
#
  S <- alpha * Lmd + (beta+h) * Eye
  Sinv <- solve(S)
  return( list( x = x, y = y, mu = (beta*y + h*x0) %*% Sinv, S = S ) )
}


clipval <- function( val, cmin, cmax ){ #clip value from [cmin, cmax]
  val <- apply( cbind(cmin, val), 1, max )
  val <- apply( cbind(cmax, val), 1, min )
  return( val )
}


ShiftMat <- function( src, sft ){
  toplft <- c(1,1)  # top, left
  btmrgt <- c(nrow(src), ncol(src) ) # bottom, right

  dst <- matrix( 0, nrow=btmrgt[1], ncol=btmrgt[2] )
                                        #
  # Clipping by Cohen's Algorithm?
  dlt <- toplft+sft # destination pasting position
  dbr <- btmrgt+sft
  if( dlt[1] > btmrgt[1] || dlt[2] > btmrgt[2] )	return( dst )
  if( dbr[1] < toplft[1] || dbr[2] < toplft[2] )	return( dst )
                                        #
  dlt <- clipval( dlt, toplft, btmrgt )
  dbr <- clipval( dbr, toplft, btmrgt )
  slt <- clipval( toplft-sft, toplft, btmrgt )
  sbr <- clipval( btmrgt-sft, toplft, btmrgt )
                                        #
  dst[dlt[1]:dbr[1], dlt[2]:dbr[2]] <- src[slt[1]:sbr[1], slt[2]:sbr[2]]
  return( dst )
}


LoopyBP2D <- function( Lx, Ly, alpha, beta, h, y, x0=0, ITmax=1000, tol=1e-9 ){
    ### Caclulate 2D-GMRF Belief Prop
    # params (Ly, Lx) input matrix size
    # alpha, beta, h, hyper-parameters
    # y: observed matrix array, size (Lx*Ly)
    #
    # Detail is in study memo around 20140327
    #

    # Modified including average x0~ mean(y)

    # shift origin. depricated for averaged model
    # ymean <- mean( y )
    # y <- y - ymean


    dirs <- 4
    M <- Ly*Lx
    bnds <- list( c(0,1), c(-1,0), c(0,-1), c(1,0) ) # bnds[[i]] means i-th direction of message from
    names(bnds) <- c('rt', 'up', 'lt', 'dw' )
    msgs = list( matrix(0, M, dirs), matrix(0, M, dirs) )
    names(msgs) = c( 'mu', 'gam' )
    tmsgs = list( matrix(0, M, dirs), matrix(0, M, dirs) )
    names(tmsgs) = c( 'mu', 'gam' )

    if( length(alpha) == 1 ){
        a0 <- alpha
        alpha <- rep(alpha, 2*M-Lx-Ly)
    }
    ah <- matrix( alpha[1:(Ly*(Lx-1))], Ly, Lx-1 )
    av <- matrix( alpha[(Ly*(Lx-1)+1):(2*M-Lx-Ly)], Ly-1, Lx )
    invart <- as.array( cbind( 1/ah, matrix(0, Ly, 1) ) )
    invalt <- as.array( cbind( matrix(0, Ly, 1), 1/ah ) )
    invaup <- as.array( rbind( matrix(0, 1, Lx), 1/av ) )
    invadw <- as.array( rbind( 1/av, matrix(0, 1, Lx) ) )

    inva <- matrix( c(invart, invaup, invalt, invadw), M, 4 )

#    browser()


    lbpx <- rep(0, M)
    for( t in seq(1:ITmax) ){
        sumMsg1 <- msgs$mu * msgs$gam
        for( d in seq(1:dirs) ){
            od <- ((d-1)+2) %% 4 + 1
            insum1 <- rowSums( sumMsg1[,-od] )
            insum2 <- rowSums( msgs$gam[,-od] )

            tmean <- (beta * y + h * x0 + insum1) / (beta + h + insum2)
            tinvar <- inva[,d] + 1/(beta + h + insum2)
#            tinvar <- 1/alpha + 1/(beta + h + insum2)
            tinvar <- 1/tinvar
            sft <- bnds[[d]]
            tmsgs$mu[,d]  <- as.array( ShiftMat( matrix(tmean, Ly, Lx), sft ) )
            tmsgs$gam[,d] <- as.array( ShiftMat( matrix(tinvar, Ly, Lx), sft ) )
        }
        msgs <- tmsgs

        sumMsg1 <- msgs$mu * msgs$gam
        outsum1 <- rowSums( sumMsg1 )
        outsum2 <- rowSums( msgs$gam )
        sgm2 <- 1/(beta + h + outsum2)
        mu <- (beta * y + h * x0 + outsum1) * sgm2

        err <- mean( abs(lbpx-mu)/sum(abs(mu)) )
        if( err < tol )	break

        lbpx <- mu
    }

    # revert shift origin
    # mu <- mu + ymean
    # msgs$mu <- msgs$mu + ymean

    ret <- list( mu, sgm2, msgs, t )
    names( ret ) <- c( 'mu', 'sigma2', 'message', 't' )
    return( ret )
}



#LoopyBP2D <- function( Lx, Ly, alpha, beta, h, y, ITmax=1000 ){
#    ### Calculate 2D-GMRF Belief Prop
#    # params: (Ly, Lx) input matrix size
#    #         alpha, beta, h, hyper-parameters
#    #         y: observed matrix which dimension = (Ly, Lx)
#
#    # The number of Message passing direction is 4 in the lattice
#    #
#    Mpars <- 2 # num of message parameters: (mu, gamma) 
#    bnds <- list( c(0,1), c(-1,0), c(0,-1), c(1,0) ) # bnds[[i]] means i-th direction of message from
#    msgs <- list()
#    msgs <- lapply( bnds, function(u) MkMsgLst( Mpars,  c( Ly, Lx ), names=c( 'mu', 'gam' ) ) )
#    tmsgs <- list()
#    tmsgs <- lapply( bnds, function(u) MkMsgLst( Mpars,  c( Ly, Lx ), names=c( 'mu', 'gam' ) ) )
#    sumMsg1 <- sapply( msgs, function(u) u$gam * u$mu )  # Components of Message: gamma * mu, which dimensions are (# of units) x (directions)
#    sumMsg2 <- sapply( msgs, function(u) u$gam )         # Components of Message: gamma , which dimensions are (# of units) x (directions)
#
#    lbpx <- matrix( 0, nrow=Ly, ncol=Lx )
#    for( t in seq(1:ITmax) ){
#        for( i in seq( 1, length(bnds) ) ){ # i means incoming bond direction index
##            browser()
#                
#            sft <- bnds[[i]] # sft direction(incoming bond direction)
#            msg <- msgs[[i]] # Message from sft direction
#
#            oi <- ((i-1)+2)%%4 + 1 # Output direction index, it's Oops ugly codes
#
#            # Sum-up message except the output direction
#            insum1 <- matrix( rowSums( sumMsg1[,-oi] ), nrow=Ly, ncol=Lx )
#            insum2 <- matrix( rowSums( sumMsg2[,-oi] ), nrow=Ly, ncol=Lx )
#
#            # Output of message is mean and variance
#            tmean   <- (beta * y + insum1 )/( beta + h + insum2 )
#            tinvvar <- 1/alpha + 1/(beta + h + insum2 )
#            tinvvar <- 1/tinvvar
#            tmsgs[[i]]$mu  <- ShiftMat( tmean, sft )
#            tmsgs[[i]]$gam <- ShiftMat( tinvvar, sft )
#            # print( msgs[[i]] )
#        }
#        msgs <- tmsgs
#        sumMsg1 <- sapply( msgs, function(u) u$gam * u$mu )  # Message of (gamma * mu)=>(# unit) x (directions)
#        sumMsg2 <- sapply( msgs, function(u) u$gam )         # Message of (gamma)ma   =>(# unit) x (directions)
#
#        outsum1 <- matrix( rowSums( sumMsg1 ), nrow=Ly, ncol=Lx )
#        outsum2 <- matrix( rowSums( sumMsg2 ), nrow=Ly, ncol=Lx )
#        sgm2 <- 1/(beta + h + outsum2)
#        mu <- (beta * y + outsum1) * sgm2
#
#        # Check convergence
#        err <- mean(abs(lbpx-mu))/sum( abs(mu) )
#        #  cat( "t =", t, "err =", err, "\n" );
#        if( err < 1.0e-9 ) break
#
#        lbpx <- mu
#    }
#    ret <- list( mu, sgm2, msgs, t )
#    names( ret ) <-  c( 'mu', 'sigma2', 'message', 't' )
#    return( ret )
#}
