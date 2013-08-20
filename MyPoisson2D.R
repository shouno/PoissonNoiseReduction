#
# 20130305: Atidlde は ~/Documents/DOCs/TeX/TVBayes/Rproj/MyLWD2D.R から
#


Trace <- function( mat ){ # mat should be a square matrix 
  return( sum( diag( mat ) ) )
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


MyLambda <- function( Lx, Ly ){
  Ah <- rep( 1, Ly*(Lx-1) )
  Av <- rep( 1, (Ly-1)*Lx )
  return( Atilde.csr( Ah, Av, Ly, Lx ) )
}



#
# 隣接要素の行列を生成
# Ah は，横方向のボンドに張り付いた変数
# Av は，縦方向のボンドに張り付いた変数
# Prow は画像の横サイズ
# Pcol は画像の縦サイズ
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
