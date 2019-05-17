desparsePrec <- function(
  dat, rho, type = 'glasso'
){
  empCv = cov(dat);
  switch( type,
          glasso = desparsePrecGlasso( empCv, rho ),
          ridge  = desparsePrecRidge( dat,empCv, rho )
  )
}

desparsePrecGlasso <- function(
  emp, rho
){
  glEst = glasso(emp,rho);
  return( 2*glEst$wi - glEst$wi%*%emp%*%glEst$wi );
}

desparsePrecRidge <- function(
  dat, emp, rho
){
  k  = ncol(dat);
  sv = svd(dat);
  glEst = glasso(emp,1)$wi;
  prj= sv$v%*%t(sv$v)%*%glEst;
  diag(prj) <- 0;
  rdEst = solve(emp+diag(rho,k));
  return( rdEst + prj );
}

sparseMat <- function(
  shat, k, alf=0.5, iter = 10, pnrm=Inf, THRSH='hard'
){
  thrLen = length(THRSH);
  if(thrLen==1){
    if(THRSH=='all'){
      THRSH =c("hard","soft","scad","adpt");
      thrLen = length(THRSH);
    }
  }
  toRet  = vector("list",thrLen);
  names(toRet) <- THRSH
  tmp   = array(0,dim=c(k,k,iter+1));
  for( i in 1:thrLen ){
    toRet[[i]] = tmp;
    toRet[[i]][,,1]=shat;
  }

  # Empirical Diagonal and normalize
  origDiag= abs(diag(shat));
  diag(shat) = abs(diag(shat));
  eDiag   = sqrt(abs(diag(diag(shat))));
  invDiag = sqrt(abs(diag(1/diag(shat))));
  shat    = invDiag %*% shat %*% invDiag;
  # Compute eta quantile of off-diagonal entries
  # (see Kashlak & Kong (2018))
  sigmas  = shat[lower.tri(diag(k))];
  if(iter==0){
    #for(i in 1:thrLen)
    #  toRet[[i]][,,2]=origDiag;
    return(toRet);
  }
  if( is.na(sum(abs(sigmas))) )
    warning("NA Error in sigmas");
  sigMedi = c(0,quantile( abs(sigmas),1-alf,na.rm=T));
  # Compute the 0% false positive estimator
  fp00    = diag(diag(shat));
  for( j in 1:thrLen )
    toRet[[j]][,,2] = eDiag%*%thresh( shat,sigMedi[2],THRSH[j] )%*%eDiag;
  if(iter==1)
    return(toRet);

  for( i in 2:iter ){
    # Compute eta % false positive estimator
    fp50    = thresh(shat,sigMedi[i],type='hard');
    if(sum(fp50==diag(k))==k^2){
      for( j in 1:thrLen)
        toRet[[j]][,,i] = origDiag;
      sigMedi = c(sigMedi,sigMedi[i])
      next;
    }
    # Compute distance between above two estimators
    ra50a   = phi( fp50-fp00,pnrm );
    # Compute the confidence ball radius for the alf % estimator
    ra   = ra50a*(sqrt(alf));
    # Search confidence ball for a the alf % estimator
    sigMedi = c(sigMedi,
      spc_CoMZeroBinSearchRev( shat, k, ra, pnrm )
    );
    for( j in 1:thrLen )
      toRet[[j]][,,i+1] = eDiag%*%thresh( shat,sigMedi[i+1],THRSH[j] )%*%eDiag;
  }
  return(toRet);
}

#
#  Function to search confidence ball for sparse estimator
#  Apply a binary search technique
#

spc_CoMZeroBinSearchRev <- function( shat, k, ra, pnrm ){
  maxThr  = 1;
  lmbPrev = 0;
  quant   = 0.5;
  delta   = 0.25;
  spar    = shat;
  # 20 steps of the binary search is more than enough
  # Could probably reduce this
  for( i in 1:10 ){
    lmb = quant;
    if( lmb == lmbPrev )
      break;
    if(lmb > maxThr){
      quant=quant-delta;
      delta=delta/2;
      next;
    }
    snew= thresh( shat,lmb,'hard' );
    if( phi(snew-diag(k),pnrm)<ra ){
      spar = snew;
      quant= quant - delta;
    } else {
      quant= quant + delta;
    }
    lmbPrev = lmb;
    delta   = delta/2;
  }
  return(lmb);
}

#####
#  Some norm and distance functions
#####

# Phi, the distance function
phi <- function( mat,p ){
  return(pschnorm(mat,p))
}

# General p-Schatten norm
pschnorm <- function( mat, p ){
  if(p==2)
    return( hsnorm(mat) );
  ev = abs(eigen(mat,only.values=TRUE,symmetric=TRUE)$values);
  if( is.infinite(p) )
    return(max(ev));
  return( sum(ev^p)^(1/p) );
}

# Hilbert-Schmidt / Frobenius Norm
hsnorm <- function( mat ){
  return( sqrt(sum(mat^2)) );
}

# Matrix Square Root
sqrtMat <-function(A)
{
  svdA= svd( A );
  D   = diag( sqrt(svdA$d) );
  return( svdA$u %*% D %*% t(svdA$v) )
}


#####
#  Threshold functions
#####

# Generic function
thresh <- function( s, lmb, type ){
  switch( type,
          soft = thr_soft(s,lmb),
          hard = thr_hard(s,lmb),
          scad = thr_scad(s,lmb),
          adpt = thr_adpt(s,lmb)
  )
}

# Hard Threshold
thr_hard <- function(s,lmb){
  return( s*(abs(s)>lmb) );
}

# Soft Threshold
thr_soft <- function(s,lmb){
  res = abs(s)-lmb;
  return( sign(s)*res*( res>0 ) );
}

# SCAD Threshold
thr_scad <- function(s,lmb,aa=3.7){
  msk1= s<(2*lmb);
  msk2= (s<(aa*lmb))-msk1;
  msk3= rep(1,length(s))-msk1-msk2;
  scd = (aa-1)/(aa-2)*(s-2*lmb)+lmb;
  res = thr_soft(s*msk1,lmb) + (scd)*msk2 + thr_hard(s*msk3,lmb);
  return(res);
}

# Adaptive LASSO Threshold
thr_adpt <- function(s,lmb,eta=1){
  res = abs(s)-lmb^(eta+1)*abs(s)^(-eta);
  res[which(is.na(res))]=0;
  res[which(is.infinite(res))]=0;
  return( sign(s)*res*( res>0 ) );
}



