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
# 20140325: Loopy BP による高速化
# 20140328: BP の実装変化にともなう変更＋αの取り扱い変更
# 20140329: LBP 使うなら疎行列いらなくね？
# 20141021: 平均値バージョンを追加しておく
#
library( png )
library( fields )
library( gplots )


source( 'MyPoisson2D.R' )
source( 'MyLBP.R' )

Coef1 <- function( x ){
  y <- rep( 0, length( x ) )
  y[x==0] <- 1.0
  y[x!=0] <- tanh( x[x!=0] ) / x[x!=0]
  return( y )
}


Trial <- function( filename, lmdmat, K=1000, minoff=2, rng=20, ITmax=100 ){
    Lx <- ncol( lmdmat );
    Ly <- nrow( lmdmat );
    N <- 1   # number of sampling images

    obsmat <- NBinom2D( N, lmdmat, Ksamples=K ) 
    obsmat <- matrix( obsmat, nrow=Lx*Ly, ncol=N )
    obsmean <- rowMeans( obsmat )
    lmd0 = mean( obsmat )

    rhomat <- lmdmat/K
    xtrue <- logit( rhomat )

    M <- Lx * Ly
    
    zz <- as.vector(2 * obsmean - K)

    xi <- rep( 0, M )
    alpha0 <- 1e-3
    h <- 1.0e-5        # treating as fixed value

    rho0 <- lmd0/K
    muold <- rep( logit(rho0), M )
    muold0 <- rep( logit(rho0), M )
    
#    browser()
    
    for(t in 1:ITmax) {
#        mumat <- matrix( muold0, Ly, Lx )
#        etah <- array( 0.5 * (mumat[,1:(Lx-1)] + mumat[,2:Lx]) )
#        etav <- array( 0.5 * (mumat[1:(Ly-1),] + mumat[2:Ly,]) )
#        etah <- dilogit( etah )
#        etav <- dilogit( etav )
        etah <- rep(1, (Lx-1)*Ly)
        etav <- rep(1, Lx*(Ly-1))
        alph <- alpha0*(etah**2)
        alpv <- alpha0*(etav**2)
        alpha <- c( alph, alpv )

        beta <- K * Coef1( xi )
        y <- zz / beta

        lbp <- LoopyBP2D( Lx, Ly, alpha, beta, h, y, x0=rho0, ITmax=1000 )

        mu <- as.numeric( lbp$mu )
        Sinv <- as.numeric( lbp$sigma2 )
        #20140328 modify implementation
        gam <- lbp$message$gam
        msg1 <- as.numeric( matrix( gam[,1], Ly, Lx )[,2:Lx] )       # right
        msg3 <- as.numeric( matrix( gam[,3], Ly, Lx )[,1:(Lx-1)] )   # left
        msg2 <- as.numeric( matrix( gam[,2], Ly, Lx )[1:(Ly-1),] )   # up
        msg4 <- as.numeric( matrix( gam[,4], Ly, Lx )[2:Ly,] )       # down
        hcor <- (alph - msg1) * (alph - msg3 ) / alph**3
        vcor <- (alpv - msg2) * (alpv - msg4 ) / alpv**3
        corsum <- sum(hcor) + sum(vcor)

        matmu <- matrix( mu, Ly, Lx )
        matsg <- matrix( Sinv, Ly, Lx )
        etah <- matrix( etah, Ly, Lx-1 )
        etav <- matrix( etav, Ly-1, Lx )
        mu2hdif <- (etah * (matmu[,1:(Lx-1)] - matmu[,2:Lx]) )**2
        mu2vdif <- (etav * (matmu[1:(Ly-1),] - matmu[2:Ly,]) )**2
        sg2hadd <- (etah**2) * (matsg[,1:(Lx-1)] + matsg[,2:Lx])
        sg2vadd <- (etav**2) * (matsg[1:(Ly-1),] + matsg[2:Ly,])

        alphainv <- sum( mu2hdif + sg2hadd ) + sum( mu2vdif + sg2vadd ) - 2 * corsum
        alphainv <- as.numeric( alphainv/(M-1) )
        xi <- as.array( sqrt( mu**2 + Sinv ) )

        alpha0 <- 1/alphainv
        err <- (sum(abs(mu-muold))/sum(abs(mu)))
        cat( t, minoff, rng, alpha0, h, err, '\n' )

        # show interractive result
        y0 <- as.vector(rhomat)
        y1 <- as.vector(obsmean)/K
        est <- as.vector( ilogit(mu) )

        cat( '    PSNR: ', PSNR( y0, est ), '<->',  PSNR( y0, y1 ), '\n' )
#        yrng <- c(min(c(y0,y1,est)), max(c(y0,y1,est)))
#        plot( y0[1000:1280], type='l', ylim=yrng )
#        lines( y1[1000:1280], type='l', col='blue', lwd=1 )
#        lines( est[1000:1280], type='l', col='red', lwd=1 )

        if( err < 1.0e-5 )	break
        muold <- mu

    }

    save( file=filename, rhomat, obsmat, xtrue, xi, alpha, h, y0, y1, est )
}




#tstart <- proc.time()
#Lx <- 64
#Ly <- 64
#xoff <- 72
#yoff <- 0
Lx <- 128
Ly <- 128
xoff <- 36
yoff <- 0
N <- 1
K <- 1000

wholeimg <- readPNG( 'cameraman.png' )
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


#rngs <- c( 10, 20, 40, 80, 160 )
#rngs <- c( 20, 80, 160 )
rngs <- c(10, 20, 40, 80, 160)

minoff <- 2

for( rng in rngs ){
    lmdmat <- minoff + rng / (maximg-minimg) * (img-minimg)

    NTrial <- 1
    for(nt in 1:NTrial) {
        fname <- sprintf( "Pois2DLBP3_%02d_%02d_%02d.RData", minoff, rng, nt )
        Trial( file=fname, lmdmat, K=1000, minoff=minoff, rng=rng, ITmax=1000 )
    }
}
