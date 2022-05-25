# SVD (With scaling - equivalent to PCA)
svd.impute <- function (x, r) {
  # 1) Initialize data
  mindex <- is.na(x)
  mcol <- colMeans(x, na.rm=T)
  mrow <- rowMeans(x, na.rm=T)
  mtot <- mean(x, na.rm=T)
  last <- 0
  xtemp <- x
  # For single missing values, average col and row means
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    xtemp[i,j] <- (mrow[i] + mcol[j]) / 2
   }
  }
  # If entire row is missing, use column mean
  for (i in 1:nrow(x)) {
    if (is.na(mrow[i])) {
       xtemp[i,] <- mcol
    }
  }
  # If entire column is missing, use row mean
  for (j in 1:ncol(x)) {
    if (is.na(mcol[j])) {
       xtemp[,j] <- mrow
    }
  }
  # If entire row and column are missing, use overall mean
  for (i in 1:nrow(x)) {
   for (j in 1:ncol(x)) {
    if (is.na(mcol[j]) & is.na(mrow[i])) {
       xtemp[i,j] <- mtot
    }
   }
  }
  xtemp[!mindex] <- x[!mindex]

  # Center the data (based on x)
  center <- matrix(rep(apply(x,1,mean,na.rm=T),ncol(x)),nrow=nrow(x))
  center[which(is.na(center))] <- 0
  xtemp <- xtemp - center
  
  # Scale the data
  scalef <- matrix(rep(apply(xtemp,1,sd,na.rm=T),ncol(x)),nrow=nrow(x))
  scalef[which(is.na(scalef))] <- 1
  xtemp <- xtemp / scalef

 ### LOOP 2 and 3 until convergence
 conv <- F
 while(!conv) {
  # 2) Compute rank r SVD
  stemp <- svd(xtemp, nu=r, nv=r)
  
  # 3) Impute missing values using SVD in (2) (EM algorithm)
  if (r>0) {
    xfull <- stemp$u[,1:r] %*% diag(x=stemp$d[1:r], nrow=r) %*% t(stemp$v[,1:r])
  } else if (r==0) {
    xfull <- matrix(0,nrow(x),ncol(x))
  }
  last <- xtemp
  xtemp[mindex] <- xfull[mindex]

  # Check for convergence
  if (norm(xtemp - last, type='f')^2  / norm(xtemp, type='f')^2 < .0001) { conv <- T }
 }
  return(xtemp * scalef + center)
}















