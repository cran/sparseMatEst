#####################
# Matrix Generation #
#####################

# Tri-Diagonal Matrix
genTriDiag <- function( k, x, d=1 ){
  mat = diag(rep(d,k));
  off = diag(x,k-1);
  off = rbind(off,0);
  off = cbind(0,off);
  mat = mat+off+t(off);
  return(mat);
}

# Auto-Regressive Matrix
genArMat <- function( k, x ){
  indx = matrix(1:k,k,k);
  indx = indx-t(indx);
  mat  = x^abs(indx);
  return(mat);
}

# Banded Matrix
genBandedMat <- function( k, m ){
  mat = matrix(0,k,k);
  for( i in 1:(m-1) ){
    off = diag(rep(1-i/m,k-i));
    off = rbind(off,matrix(0,i,k-i));
    off = cbind(matrix(0,k,i),off);
    mat = mat + off;
  }
  mat = mat + t(mat);
  diag(mat) <- 1;
  return(mat);
}

# Random Sparse Matrix
genRandSpMat <- function( k,p=0.1 ){
  msk = matrix(rbinom(k^2,1,p),k,k);
  rnd = matrix(runif(k^2,0.3,0.8),k,k);
  mat = msk*rnd;
  mat = mat+t(mat);
  eig = min(eigen(mat,only.values=TRUE)$values);
  if(eig>0)
    eig = 0;
  mat = mat + diag(rep(0.01-eig,k));
  return(mat);
}


# Binary Tree Matrix
genBinaryTreeMat <- function( k, p=0.25 ){
  cnt = k;
  tot = k*(k+1)/2;
  mat = cbind(rbind(rep(0,k-1),diag(k-1)),0);
  mat = rbind(mat,c(1,rep(0,k-1)),matrix(0,tot-cnt-1,k));
  for( i in (k-1):2 ){
    tmp = matrix(0,cnt,i);
    cnt = cnt + i;
    tmp = rbind(tmp,  cbind(rbind(rep(0,i-1),diag(i-1)),0) );
    tmp = rbind(tmp,  c(1,rep(0,i-1)),matrix(0,tot-cnt-1,i)  );
    mat = cbind(mat,tmp);
  }
  mat = cbind(mat,0);
  return( diag(tot)+p*(mat+t(mat)) );
}


# Multi-Diagonal Matrix
genMultDiag <- function( k, x, m=1, d=1 ){
  mat = diag(rep(d,k));
  off = matrix(0,k-m,k-m);
  for( i in 1:m ){
    diag(off) <- x;
    off = rbind(off,0);
    off = cbind(0,off);
  }
  mat = mat+off+t(off);
  return(mat);
}

# Block diagonal matrix
genBlockDiag <- function( nBlock, sBlock, x ){
  k   = nBlock*sBlock;
  blk = matrix(x,sBlock,sBlock);
  diag(blk) <- 1;
  row = cbind( blk,matrix(0,sBlock,sBlock*(nBlock-1))  );
  mat = row;
  for( i in 2:nBlock )
    mat = rbind( 
      mat, 
      row[, c( (k-sBlock*(i-1)+1):k,1:(k-sBlock*(i-1)) ) ]
    );
  return(mat);
}