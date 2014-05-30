#
# 20130305: コメント付け
# 20130305: スパース行列を使うことを考える．
#           ~/Documents/DOCs/TeX/TVBayes/Rproj/BayesLWD2D6.R を参考に
# 20130712: 前バージョンのが多分うまくいかないだろうという理論計算結果だったので
#           大幅に改定作業．
#           多分，オフセットμを推定パラメータにするとそこに全力で全部喰われる
# 20130713: K がでかすぎる？ K-> 1000 くらいで良い感じか？
#           最適化が大暴れするので Nelder-Mead 法であたりをつけてから CG で，再最適化
# 20130713: ログ取れるように関数化しておく
# 20130714: 実画像で遊べる？


library( Matrix )    # 疎行列パッケージを導入
library( png )

source( 'MyPoisson2D.R' )


Coef1 <- function( x ){
  y <- rep( 0, length( x ) )
  y[x==0] <- 1.0
  y[x!=0] <- tanh( x[x!=0] ) / x[x!=0]
  return( y )
}


NBinom2D <- function( nn, lmdmat, Ksamples=10000 ){
  Ly <- nrow(lmdmat)
  Lx <- ncol(lmdmat)
  LL <- Lx*Ly
  gennum <- rbinom( n = nn * LL, K, lmdmat/K )
  PData <- aperm( array( gennum, dim=c(Lx, Ly, nn) ), perm=c(1,2,3) )
  return( PData )
}

logit <- function( u ){
    return( 0.5 * log(u/(1-u)) )
}

ilogit <- function(x){
    return( exp(x)/(exp(x)+exp(-x)) )
}

PSNR <- function( ref, tst ){
    mse <- mean( (tst-ref)**2 )
    rng2 <- (max(ref) - min(ref))**2
    return ( 10 * log(rng2/mse, base=10) )
}




Trial <- function( filename, lmdmat, K=1000, minoff=2, rng=20, ITmax=100 ){
    Lx <- ncol( lmdmat );
    Ly <- nrow( lmdmat );
    N <- 1

    obsmat <- NBinom2D( N, lmdmat, Ksamples=K ) 
    obsmat <- Matrix( obsmat, nrow=Lx*Ly, ncol=N )
    obsmean <- rowMeans( obsmat )

    rhomat <- lmdmat/K
    xtrue <- logit( rhomat )

    M <- Lx * Ly

    Lmd <- MyLambda( Lx, Ly )
    Eye <- Diagonal( x = rep( 1, M ) )
    zz <- as.vector(2 * obsmean - K)
    Lmdeigs <- eigen( Lmd )
    eigvals <- Lmdeigs$values

    xi <- rep( 0, M )
    alpha <- 1.0
    h <- 1.0
    lnA <- log(alpha)
    lnH <- log(h)
    rhoinit <- mean(obsmean)/K
    for(t in 1:ITmax) {
        Tht <- Diagonal( x = Coef1( xi ) * K )
        C1 <- solve( Tht + alpha * Lmd + h * Eye )
        mnew <- C1 %*% zz

        AA <- sum( diag( Lmd %*% C1 ) )+ as.numeric(t(mnew) %*% Lmd %*% mnew)
        BB <- sum( diag( C1 ) ) + sum(mnew ** 2)

        f1 <- function( params ){
            alpha <- exp( params[1] )
            h <- exp( params[2] )
            return( sum( eigvals/(alpha*eigvals+h) ) - AA )
        }
        f2 <- function( params ){
            alpha <- exp( params[1] )
            h <- exp( params[2] )
            return( sum( 1/(alpha*eigvals+h) ) - BB )
        }
        optf <- function( params ){
            return( f1(params)**2 + f2(params)**2 )
        }
        groptf <- function( params ){
            alpha <- exp( params[1] )
            h <- exp( params[2] )
            tmp2 <- - sum( (eigvals/(alpha*eigvals+h))**2 )
            tmp1 <- - sum( eigvals/((alpha*eigvals+h)**2) )
            tmp0 <- - sum( 1      /((alpha*eigvals+h)**2) )
            df1da <- 2*f1(params)*tmp2
            df2da <- 2*f2(params)*tmp1
            df1dh <- 2*f1(params)*tmp1
            df2dh <- 2*f2(params)*tmp0
            dfda <- ( df1da + df2da ) * alpha
            dfdh <- ( df1dh + df2dh ) * h
            return( c(dfda, dfdh) )
        }
        #  opts <- optim( par=c(alpha,h), fn=optf, gr=groptf, method="L-BFGS-B",
        #                lower=c(0,0), upper=c(Inf, Inf) )
        # opts <- optim( par=c(lnA,lnH), fn=optf, gr=groptf, method="BFGS" )
        # opts <- optim( par=c(lnA,lnH), fn=optf, method="Nelder-Mead" )   
        opts <- optim( par=c(lnA,lnH), fn=optf, method="Nelder-Mead" )   
        lnA <- opts$par[1]
        lnH <- opts$par[2]
        opts <- optim( par=c(lnA,lnH), fn=optf, gr=groptf, method="CG" )   
        lnA <- opts$par[1]
        lnH <- opts$par[2]

        alpha <- exp(lnA)
        h <- exp(lnH)
        xi <- as.array( sqrt( diag( C1 ) + mnew * mnew ) )

        cat( t, minoff, rng, alpha, h, opts$value, '\n' )

#    cat( t, minoff, rng, alpha, h, resmu, resalpha, resh, '\n' )
#
#    if( (resmu < tol1) && (resalpha < tol2) && (resh < tol2) )	break
#    else{
#      alpha <- opts$par[1]
#      if( alpha < 0 ) alpha <- 1e-5
#      h <- opts$par[2]
#      if( h < 0 ) h <- 1e-5
#      mu <- as.numeric(munew)
#    }
        y0 <- as.vector(rhomat)
        y1 <- as.vector(obsmean)/K
        est <- as.vector( ilogit(mnew) )

        cat( '    PSNR: ', PSNR( y0, est ), '\n' )
        yrng <- c(min(c(y0,y1,est)), max(c(y0,y1,est)))
#    yrng <- c(0.0, 0.005)
        plot( y0, type='l', ylim=yrng )
        lines( y1, type='l', col='blue', lwd=1 )
        lines( as.vector(est), type='l', col='red', lwd=1 )
    }
    save( file=filename, rhomat, obsmat, xtrue, xi, alpha, h, y0, y1, est )
}




#tstart <- proc.time()
Lx <- 64
Ly <- 64
N <- 1
K <- 1000

wholeimg <- readPNG( 'cameraman.png' )
xoff <- 72
yoff <- 0
img <- wholeimg[(yoff+1):(yoff+Ly), (xoff+1):(xoff+Lx)]
#img <- 1-img # 黒が 0 に近くなるように



#K0 <- 0.2 * pi
#Kx <- 1.2
#Ky <- .2
#Kx <- 1.0
#Ky <- 1.5
#ix <- matrix( rep( seq(0, Lx-1) * 2*pi / Lx * Kx, Ly ), nrow=Ly, byrow="True" )
#iy <- matrix( rep( seq(0, Ly-1) * 2*pi / Ly * Ky, Lx ), ncol=Lx )

#img <- cos( ix + iy + K0 )
maximg <- max( img )
minimg <- min( img )


rngs <- c( 20, 40, 80, 160 )
#rngs <- c( 20, 80, 160 )
#rngs <- c( 20 )


minoff <- 2

for( rng in rngs ){
    lmdmat <- minoff + rng / (maximg-minimg) * (img-minimg)

    NTrial <- 1
    for(nt in 1:NTrial) {
        fname <- sprintf( "Pois2D13_%02d_%02d_%02d.RData", minoff, rng, nt )
        Trial( file=fname, lmdmat, K=1000, minoff=minoff, rng=rng, ITmax=100 )
    }
}
