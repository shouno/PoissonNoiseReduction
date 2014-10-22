#
# 20140325: Loopy BP Accelaration
# 20140406: Compare with Exact EM solution
# 20141021: include average model
#

library( Matrix ) # Sparse matrix package for exact solution
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


Trial <- function( lmdmat, K=1000, minoff=2, rng=20, ITmax=100 ){
    Lx <- ncol( lmdmat );
    Ly <- nrow( lmdmat );
    N <- 1

    obsmat <- NBinom2D( N, lmdmat, Ksamples=K ) 
    obsmat <- matrix( obsmat, nrow=Lx*Ly, ncol=N )
    obsmean <- rowMeans( obsmat )
    lmd0 <- mean( obsmat )
    
    rhomat <- lmdmat/K
    xtrue <- logit( rhomat )

    M <- Lx * Ly
    
    zz <- as.vector(2 * obsmean - K)

    lbpxi <- rep( 0, M )
    lbpalpha0 <- 1e-3
    h <- 1.0e-5        # treating as fixed value

    rho0 <- lmd0/K
    lbpmuold <- rep( logit(rho0), M )
    lbpmuold0 <- rep( logit(rho0), M )

    extmu <- rep( logit(rho0), M )
    extmuold <- rep( logit(rho0), M )
    extxi <- rep( 0, M )
    extalpha0 <- 1.0
#    browser()
    extLmd <- MyLambda( Ly, Lx )
    extEye <- Diagonal( x = rep( 1, M ))

    zerotime <- proc.time()
    zerotime <- proc.time() - zerotime
    sumexttime <- zerotime
    sumlbptime <- zerotime

    for(t in 1:ITmax) {
        #
        # Exact solution part
        #
        exttime <- proc.time()
        Tht <- Diagonal( x = Coef1(extxi) * K )
        C1 <- solve( Tht + extalpha0 * extLmd + h * extEye )
        extmu <- C1 %*% (zz + h * rho0)
        extalphainv <- sum( diag(extLmd %*% C1) ) + as.numeric( t(extmu) %*% extLmd %*% extmu )
        extalphainv <- as.numeric( extalphainv/(M-1) )
        extxi <- as.array( sqrt( diag(C1) + extmu**2 ) )
        extalpha0 <- 1/extalphainv
        exttime <- proc.time() - exttime
        sumexttime <- sumexttime + exttime

        #
        # LBP solution part
        #
        lbptime <- proc.time()
        etah <- rep(1, (Lx-1)*Ly)
        etav <- rep(1, Lx*(Ly-1))
        alph <- lbpalpha0*(etah**2)
        alpv <- lbpalpha0*(etav**2)
        alpha <- c( alph, alpv )

        beta <- K * Coef1( lbpxi )
        y <- zz / beta

        lbp <- LoopyBP2D( Lx, Ly, alpha, beta, h, y, x0=logit(rho0), ITmax=1000 )

        lbpmu <- as.numeric( lbp$mu )
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

        matmu <- matrix( lbpmu, Ly, Lx )
        matsg <- matrix( Sinv, Ly, Lx )
        etah <- matrix( etah, Ly, Lx-1 )
        etav <- matrix( etav, Ly-1, Lx )
        mu2hdif <- (etah * (matmu[,1:(Lx-1)] - matmu[,2:Lx]) )**2
        mu2vdif <- (etav * (matmu[1:(Ly-1),] - matmu[2:Ly,]) )**2
        sg2hadd <- (etah**2) * (matsg[,1:(Lx-1)] + matsg[,2:Lx])
        sg2vadd <- (etav**2) * (matsg[1:(Ly-1),] + matsg[2:Ly,])

        lbpalphainv <- sum( mu2hdif + sg2hadd ) + sum( mu2vdif + sg2vadd ) - 2 * corsum
        lbpalphainv <- as.numeric( lbpalphainv/(M-1) )
        lbpxi <- as.array( sqrt( lbpmu**2 + Sinv ) )

        lbpalpha0 <- 1/lbpalphainv
        lbptime <- proc.time() - lbptime
        sumlbptime <- sumlbptime + lbptime

        lbperr <- (sum(abs(lbpmu-lbpmuold))/sum(abs(lbpmu)))
        exterr <- (sum(abs(extmu-extmuold))/sum(abs(extmu)))
        # show interractive result
        y0 <- as.vector(rhomat)
        y1 <- as.vector(obsmean)/K
        lbpest <- as.vector( ilogit(lbpmu) )
        extest <- as.vector( ilogit(extmu) )
        cat( 't: ', t, minoff, rng, '\n' )
        cat( '    Exact:', extalpha0, lbperr, exttime, '\n' )
        cat( '        PSNR(Exact): ', PSNR( y0, extest ), '<->', PSNR( y0, y1 ), '\n' )
        cat( '    LBP:  ', lbpalpha0, exterr, lbptime, '\n' )
        cat( '        PSNR(LBP):   ', PSNR( y0, lbpest ), '<->', PSNR( y0, y1 ), '\n' )

#        yrng <- c(min(c(y0,y1,est)), max(c(y0,y1,est)))
#        plot( y0[1000:1280], type='l', ylim=yrng )
#        lines( y1[1000:1280], type='l', col='blue', lwd=1 )
#        lines( est[1000:1280], type='l', col='red', lwd=1 )

        if( lbperr < 1.0e-5 )	break
        lbpmuold <- lbpmu
        extmuold <- extmu

    }
    ans <- list(t, extalpha0, lbpalpha0, extxi, lbpxi, sumexttime, sumlbptime, 
                y0, y1, extest, lbpest )
    names(ans) <- c('t', 'extalpha', 'lbpalpha', 'extxi', 'lbpxi', 'exttime', 'lbptime', 
                    'y0', 'y1', 'extest', 'lbpest')
    return(ans)
#    save( file=filename, rhomat, obsmat, xtrue, xi, alpha, h, y0, y1, est )
}







#tstart <- proc.time()
#Lx <- 16
#Ly <- 16
#xoff <- 100
#yoff <- 25
#
#Lx <- 32
#Ly <- 32
#xoff <- 92
#yoff <- 20
#
#Lx <- 48
#Ly <- 48
#xoff <- 80
#yoff <- 10
#
Lx <- 64
Ly <- 64
xoff <- 72
yoff <- 0
#
#Lx <- 80
#Ly <- 80
#xoff <- 64
#yoff <- 0

#Lx <- 96
#Ly <- 96
#xoff <- 48
#yoff <- 0


#
#Lx <- 128
#Ly <- 128
#xoff <- 36
#off <- 0
#xoff <- 140
#yoff <- 20
#Lx <- 32
#Ly <- 32

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

    NTrial <- 10
    for(nt in 1:NTrial) {
        ans <- Trial( lmdmat, K=1000, minoff=minoff, rng=rng, ITmax=1000 )
        fname <- sprintf( "Poisson2DCompLBPvsExt_%02d_%02d_%02d.RData", minoff, rng, nt )
        save( file=fname, ans, Lx, Ly, rng )
    }
}


