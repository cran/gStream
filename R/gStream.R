gstream = function(distM, L, N0, k, statistics=c("all","o","w","g","m"), n0=0.3*L, n1=0.7*L, ARL=10000,alpha=0.05,skew.corr=TRUE,asymp=FALSE){
  r1 = list()
  n0 = ceiling(n0)
  n1 = floor(n1)
  if(n0<2){
    cat("Note: Starting index has been set to n0 = 2 as the graph-based statistics are not well-defined for t<2. \n")
    n0=2
  }
  if(n1>(L-2)){
    cat("Note: Ending index has been set to n1 =", L-2, " as the graph-based statistics are not well-defined for t>",L-2,". \n")
    n1=L-2
  }
  if(N0<L){
    stop("Warning: Please adjust either N0 or L. The number of historical observations (N0) must be at least L. \n")
  }
  N = dim(distM)[1]
  r1$scanZ = getscanZ(distM,N0,L,N,k,n0,n1,statistics)

  r1$b = getb(distM,ARL,alpha,N0,n0,n1,L,k,statistics,skew.corr,dif=1e-10, nIterMax=100,asymp)

  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    r1$tauhat$ori = which(r1$scanZ$ori>r1$b$ori)
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    r1$tauhat$weighted = which(r1$scanZ$weighted>r1$b$weighted)
  }
  if (length(which(!is.na(match(c("m","max","all"),statistics))))>0){
    r1$tauhat$max.type = which(r1$scanZ$max.type>r1$b$max.type)
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    r1$tauhat$generalized = which(r1$scanZ$generalized>r1$b$generalized)
  }

  return(r1)
}

#obtain test statistics
getZL = function(distM, k = 1){
  L = dim(distM)[1]
  A = matrix(0,L,k)
  for (i in 1:L){
    A[i,] = (sort(distM[i,1:L], index.return=T)$ix)[1:k]
  }
  temp = table(A)
  id = as.numeric(row.names(temp))
  deg = rep(0,L)
  deg[id] = temp
  deg.sumsq = sum(deg^2)
  cn = sum((deg-k)^2)/L/k
  count = 0
  for (i in 1:L){
    ids = A[i,]
    count = count + length(which(A[ids,]==i))
  }
  vn = count/L/k

  ts = 1:(L-1)
  q = (L-ts-1)/(L-2)
  p = (ts-1)/(L-2)
  EX1L = 2*k*(ts)*(ts-1)/(L-1)
  EX2L = 2*k*(L-ts)*(L-ts-1)/(L-1)
  EX = 4*k*ts*(L-ts)/(L-1)

  config1 = (2*k*L + 2*k*L*vn)
  config2 = (3*k^2*L + deg.sumsq -2*k*L -2*k*L*vn)
  config3 = (4*L^2*k^2 + 4*k*L + 4*k*L*vn - 12*k^2*L - 4*deg.sumsq)

  f11 = 2*(ts)*(ts-1)/L/(L-1)
  f21 = 4*(ts)*(ts-1)*(ts-2)/L/(L-1)/(L-2)
  f31 = (ts)*(ts-1)*(ts-2)*(ts-3)/L/(L-1)/(L-2)/(L-3)

  f12 = 2*(L-ts)*(L-ts-1)/L/(L-1)
  f22 = 4*(L-ts)*(L-ts-1)*(L-ts-2)/L/(L-1)/(L-2)
  f32 = (L-ts)*(L-ts-1)*(L-ts-2)*(L-ts-3)/L/(L-1)/(L-2)/(L-3)

  h = 4*(ts-1)*(L-ts-1)/((L-2)*(L-3))
  VX = EX*(h*(1+vn-2*k/(L-1))+(1-h)*cn)

  var1 = config1*f11 + config2*f21 + config3*f31 - EX1L^2
  var2 = config1*f12 + config2*f22 + config3*f32 - EX2L^2
  v12 = config3*((ts)*(ts-1)*(L-ts)*(L-ts-1))/(L*(L-1)*(L-2)*(L-3)) - EX1L*EX2L

  X = X1 = X2 = rep(0,L-1)
  for (t in 1:(L-1)){
    X2[t] = 2*(length(which(A[(t+1):L,]>t)))
    X1[t] = 2*(length(which(A[1:t,]<=t)))
    X[t] = 2*(length(which(A[1:t,]>t))+length(which(A[(t+1):L,]<=t)))
  }

  Rw = q*X1 + p*X2
  ERw = q*EX1L + p*EX2L
  varRw = q^2*var1 + p^2*var2 + 2*p*q*v12
  Zw = (Rw - ERw)/sqrt(varRw)

  Zdiff = ((X1-X2)-(EX1L-EX2L))/sqrt(var1+var2-2*v12)

  S = Zw^2 + Zdiff^2
  M = apply(cbind(abs(Zdiff),Zw),1,max)
  Z = (EX-X)/sqrt(VX)
  list(R=X,R1= X1, R2 = X2, Rw = Rw, Z1 = (X1-EX1L)/sqrt(var1) , Z2 = (X2-EX2L)/sqrt(var2), Zdiff = Zdiff, Zw = Zw, S = S, M =M, Z=Z )
}

getscanZ = function(distM,N0,L,N,k,n0,n1,statistics="all"){
  maxZ = maxZw = maxS = maxM = rep((N0+1):N)
  for (n in (N0+1):N){
    #n0s = n-L+n1
    #n1s = n-n0
    tests = getZL(distM[(n-L+1):n,(n-L+1):n],k)
    maxZ[n-N0] = max(tests$Z[n0:n1])
    maxZw[n-N0] = max(tests$Zw[n0:n1])
    maxS[n-N0] = max(tests$S[n0:n1])
    maxM[n-N0] = max(tests$M[n0:n1])
  }

  scanZ = list()
  if (length(which(!is.na(match(c("o","ori","original","all"),statistics))))>0){
    scanZ$ori = maxZ
  }
  if (length(which(!is.na(match(c("w","weighted","all"),statistics))))>0){
    scanZ$weighted = maxZw
  }
  if (length(which(!is.na(match(c("m","max","g","generalized","all"),statistics))))>0){
    scanZ$max.type = maxM
  }
  if (length(which(!is.na(match(c("g","generalized","all"),statistics))))>0){
    scanZ$generalized = maxS
  }
  return(scanZ)
}

#obtain graph-based quantities
gb_quantities = function(distM,N0,k){
  psum = qsum = psumk = qsumk = psumk1 = qsumk1 = psumk2 = qsumk2 = pLk1 = qLk1 = deg.sum3.n = aaa1.n = aaa2.n = daa.n = dda.n = rep(0,1)
  psumk_hao = qsumk_hao = rep(0,1)
  n = N0
  An = matrix(0,n,k+2)
  for (i in 1:n){
    An[i,] = (sort(distM[i,1:n], index.return=T)$ix)[1:(k+2)]
  }
  temp = table(An[,1:k])
  id = as.numeric(row.names(temp))
  deg = rep(0,n)
  deg[id] = temp
  deg.sumsq = sum(deg^2)
  deg.sum3 = sum(deg^3)
  count = daa = dda = aaa1 = aaa2 = 0
  for (i in 1:n){
    ids = An[i,1:k]
    count = count + length(which(An[ids,1:k]==i))
    daa = daa + deg[i]*length(which(An[ids,1:k]==i))
    dda = dda + deg[i]*sum(deg[ids])
    for (j in ids){
      u = An[j,1:k]
      aaa1 = aaa1 + length(which(An[u,1:k]==i))
      aaa2 = aaa2 + length(which(!is.na(match(ids,u))))
    }
  }
  psum = count/n
  qsum = deg.sumsq/n-k #j,l cannot be the same
  deg.sum3.n = deg.sum3
  aaa1.n = aaa1
  aaa2.n = aaa2
  daa.n = daa
  dda.n = dda

  count1 = count2 = count3 = count4 = count5 = count6 = count7 = count8 =  0
  for (i in 1:n){
    ids = An[i,k]
    count1 = count1 + length(which(An[ids,1:k]==i))
    count2 = count2 + length(which(An[-i,1:k]==ids))

    ids1 = An[i,k+1]
    count3 = count3 + length(which(An[ids1,1:k]==i))
    count4 = count4 + length(which(An[-i,1:k]==ids1))

    count7 = count7 + length(which(An[ids1,k+1]==i))
    count8 = count8 + length(which(An[-i,k+1]==ids1))

    ids2 = An[i,k+2]
    count5 = count5 + length(which(An[ids2,1:k]==i))
    count6 = count6 + length(which(An[-i,1:k]==ids2))
  }
  psumk = count1/n
  qsumk = count2/n
  psumk1 = count3/n
  qsumk1 = count4/n
  psumk2 = count5/n
  qsumk2 = count6/n
  pLk1 = count7/n
  qLk1 = count8/n
  list(psum=psum,qsum=qsum, psumk1=psumk1, qsumk1=qsumk1, psumk2=psumk2, qsumk2=qsumk2, pLk1=pLk1, qLk1=qLk1, psumk = psumk, qsumk=qsumk, deg.sumsq = deg.sumsq, deg.sum3.n= deg.sum3.n, aaa1.n=aaa1.n, aaa2.n=aaa2.n, daa.n=daa.n, dda.n=dda.n)
}

#obtain critical values
Nu = function(x){ # the Nu function
  y = x/2
  (1/y)*(pnorm(y)-0.5)/(y*pnorm(y) + dnorm(y))
}

C1_Z = function(x, L, k, psum,qsum){ # x=n-t
  ((16*(k + 2*psum - psum)*(2*L- 2*x - 1)*(x^2 - x))/(L^3 - 6*L^2 + 11*L - 6) - (16*k^2*x^2*(L- x))/(L - 1)^2 + (16*k^2*x*(L- x)^2)/(L - 1)^2 + (4*x*(3*k^2 + k + 2*qsum - qsum)*(3*L^2 - 10*L*x - 3*L + 8*x^2 + 2*x + 2))/(L^3 - 6*L^2 + 11*L - 6) + (16*L*k^2*x*(3*L*x - L^2 + L - 2*x^2 - 2*x + 1))/((L - 1)*(L - 2)*(L - 3)))/(4*((k^2*(((((-x+ 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k + qsum - k^2))/k + ((-x+ 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L- x)^2)/(L - 1)^2)^(1/2)) - (((16*k^2*(((((-x+ 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x+ 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L- x))/(L - 1)^2 - (16*k^2*(((((-x+ 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x+ 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x*(L- x)^2)/(L - 1)^2 + (64*k*(((((-x+ 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x+ 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))*x^2*(L- 2*x)*(L- x)^2*(psum - qsum - L*psum + L*qsum - L*k^2 + 3*k^2))/((L - 1)^2*(L^3 - 6*L^2 + 11*L - 6)))*(L*((4*x*(L- x))/(L*(L - 1)) - (16*(L- x)*(L - x - 1)*(x^2 - x))/(L*(L - 1)*(L - 2)*(L - 3)))*(3*k^2 + k + 2*qsum - qsum) + (16*(k + 2*psum - psum)*(x^2 - x)*(L^2 - 2*L*x - L + x^2 + x))/(L^3 - 6*L^2 + 11*L - 6) - (16*k^2*x^2*(L- x)^2)/(L - 1)^2 + (16*L*k^2*(L- x)*(L - x - 1)*(x^2 - x))/((L - 1)*(L - 2)*(L - 3))))/(128*((k^2*(((((-x+ 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k + qsum - k^2))/k + ((-x+ 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L- x)^2)/(L - 1)^2)^(3/2))
}

C2_Z = function(x, L, k, psum, qsum, psumk, qsumk){ # x=n-t
  ((((4*x*(L - x))/(L*(L - 1)) - (16*(L - x)*(L - x - 1)*(x^2 - x))/(L*(L - 1)*(L - 2)*(L - 3)))*(2*k + 4*qsum - 2*qsum - 7*L*k - 8*L*qsum + 7*L*qsum - 2*L*qsumk - 9*L*k^2 + 3*L^2*k + 2*L^2*qsum - 3*L^2*qsum + 2*L^2*qsumk + 6*k^2 + 3*L^2*k^2))/(L^2 - 3*L + 2) + L*(3*k^2 + k + 2*qsum - qsum)*((4*(L^2 - 2*L*x - 2*L + x^2 + 2*x))/(L*(L - 1)*(L - 2)) - (4*x*(L - x))/(L^2*(L - 1)) + (4*x*(L - x))/(L*(L - 1)^2*(L - 2)) + (16*(-x + 1)*(L^2 - 3*L*x - L + 2*x^2 + x))/(L*(L - 1)*(L - 2)*(L - 3)) - (16*(4*L^2 - 12*L + 6)*(L - x)*(L - x - 1)*(x^2 - x))/(L^2*(L - 1)^2*(L - 2)^2*(L - 3)^2)) + (16*k^2*x^2*(L - x))/(L - 1)^2 - (16*k^2*x*(L - x)^2)/(L - 1)^2 - (16*(k + 2*psum - psum)*(-x + 1)*(L^6 - 3*L^5*x - 7*L^5 + 2*L^4*x^2 + 23*L^4*x + 17*L^4 - 20*L^3*x^2 - 55*L^3*x - 17*L^3 + 4*L^2*x^3 + 50*L^2*x^2 + 47*L^2*x + 6*L^2 - 12*L*x^3 - 36*L*x^2 - 12*L*x + 6*x^3 + 6*x^2))/(L*(11*L - 6*L^2 + L^3 - 6)^2) + (16*(x^2 - x)*(L^2 - 2*L*x - L + x^2 + x)*(2*k + 4*psum - 2*psum - 7*L*k - 10*L*psum + 7*L*psum - 2*L*psumk + 3*L^2*k + 4*L^2*psum - 3*L^2*psum + 2*L^2*psumk))/(L*(L^2 - 3*L + 2)*(L^3 - 6*L^2 + 11*L - 6)) - (16*L*k^2*(-x + 1)*(L^2 - 3*L*x - L + 2*x^2 + x))/((L - 1)*(L - 2)*(L - 3)) + (16*k^2*(4*L^2 - 12*L + 6)*(L - x)*(L - x - 1)*(x^2 - x))/((L - 1)^2*(L - 2)^2*(L - 3)^2))/(4*((k^2*(((((-x + 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k + qsum - k^2))/k + ((-x + 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L - x)^2)/(L - 1)^2)^(1/2)) + (((16*k^2*(((((-x + 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x + 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L - x))/(L - 1)^2 - (16*k^2*(((((-x + 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x + 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x*(L - x)^2)/(L - 1)^2 + (64*k*(((((-x + 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k - k^2 + qsum))/k + ((-x + 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))*x^2*(L - 2*x)*(L - x)^2*(psum - qsum - L*psum + L*qsum - L*k^2 + 3*k^2))/((L - 1)^2*(L^3 - 6*L^2 + 11*L - 6)))*(L*((4*x*(L - x))/(L*(L - 1)) - (16*(L - x)*(L - x - 1)*(x^2 - x))/(L*(L - 1)*(L - 2)*(L - 3)))*(3*k^2 + k + 2*qsum - qsum) + (16*(k + 2*psum - psum)*(x^2 - x)*(L^2 - 2*L*x - L + x^2 + x))/(L^3 - 6*L^2 + 11*L - 6) - (16*k^2*x^2*(L - x)^2)/(L - 1)^2 + (16*L*k^2*(L - x)*(L - x - 1)*(x^2 - x))/((L - 1)*(L - 2)*(L - 3))))/(128*((k^2*(((((-x + 1)*(4*L - 4*x - 4))/((L - 2)*(L - 3)) + 1)*(k + qsum - k^2))/k + ((-x + 1)*(4*L - 4*x - 4)*(k + psum - L*k - L*psum + 2*k^2))/(k*(L - 1)*(L - 2)*(L - 3)))^2*x^2*(L - x)^2)/(L - 1)^2)^(3/2))
}

C1_w_asy = function(x){
  1/(2*x*(1-x))
}
C2_w_asy = function(x,k,psum,psumk1){
  (x^2-x+1)/(x*(1-x)) - (2*k*psumk1)/(k+psum)
}

C1_d_asy = function(x){
  1/(x*(1-x))
}
C2_d_asy = function(x,k,qsum,qsumk1){
   (10*qsum-4*k*qsumk1-(6*k^2-10*k))/(2*(qsum-k^2+k)) - 1/(2*x*(1-x))
 }

C1_w = function(x,L,k,psum,qsum){
  if(k==1){
    result=-(x^2*(x - 1)^2*(2*x^2 - 2*L*x + L)*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2*(2*psum - 4*L + qsum - 3*L*psum - L*qsum - L*k^2 + L^2*psum + L^2 + 3*k^2 + 3))/(2*(L - 1)^5*(L - 2)^6*(L - 3)^3*((x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
  }
  if(k==5){
    result= -(x^2*(x - 1)^2*(2*x^2 - 2*L*x + L)*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2*(2*psum - 20*L + qsum - 3*L*psum - L*qsum - L*k^2 + L^2*psum + 5*L^2 + 3*k^2 + 15))/(2*(L - 1)^5*(L - 2)^6*(L - 3)^3*((x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
  }
  else{
    result = -(x^2*(x - 1)^2*(2*x^2 - 2*L*x + L)*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^3)/(2*(L - 1)^5*(L - 2)^6*(L - 3)^3*((x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
  }
  return(result)
}

C2_w = function(x,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2){
  if(k==1){
    num=- ((x - L + 1)*(36*L - 72*x + 72*k^2*x^2 + 72*k^2*x^3 + 24*L*psum + 12*L*qsum + 66*L*x - 48*psum*x - 24*qsum*x + 36*L*k^2 - 74*L^2*psum + 85*L^3*psum - 45*L^4*psum + 11*L^5*psum - L^6*psum - 31*L^2*qsum + 27*L^3*qsum - 9*L^4*qsum + L^5*qsum - 312*L*x^2 + 88*L^2*x - 36*L*x^3 - 169*L^3*x + 96*L^4*x - 23*L^5*x + 2*L^6*x - 72*k^2*x + 88*psum*x^2 + 8*psum*x^3 - 48*psumk1*x^2 + 48*psumk1*x^3 - 16*psumk2*x^2 + 16*psumk2*x^3 + 56*qsum*x^2 - 8*qsum*x^3 - 16*qsumk1*x^2 + 16*qsumk1*x^3 - 4*qsumk2*x^2 + 4*qsumk2*x^3 - 105*L^2 + 112*L^3 - 54*L^4 + 12*L^5 - L^6 - 69*L^2*k^2 + 43*L^3*k^2 - 11*L^4*k^2 + L^5*k^2 + 144*x^2 + 313*L^2*x^2 + 33*L^2*x^3 - 153*L^3*x^2 - 10*L^3*x^3 + 35*L^4*x^2 + L^4*x^3 - 3*L^5*x^2 - 210*L*k^2*x^2 + 52*L^2*k^2*x - 66*L*k^2*x^3 - 64*L^3*k^2*x + 20*L^4*k^2*x - 2*L^5*k^2*x + 245*L^2*psum*x^2 + 33*L^2*psum*x^3 - 133*L^3*psum*x^2 - 10*L^3*psum*x^3 + 33*L^4*psum*x^2 + L^4*psum*x^3 - 3*L^5*psum*x^2 + 30*L^2*psumk1*x^2 + 70*L^2*psumk1*x^3 + 14*L^2*psumk2*x^2 - 50*L^3*psumk1*x^2 + 14*L^2*psumk2*x^3 - 20*L^3*psumk1*x^3 - 12*L^3*psumk2*x^2 + 18*L^4*psumk1*x^2 - 2*L^3*psumk2*x^3 + 2*L^4*psumk1*x^3 + 2*L^4*psumk2*x^2 - 2*L^5*psumk1*x^2 + 68*L^2*qsum*x^2 - 20*L^3*qsum*x^2 + 2*L^4*qsum*x^2 + 14*L^2*qsumk1*x^2 + 14*L^2*qsumk1*x^3 + 4*L^2*qsumk2*x^2 - 12*L^3*qsumk1*x^2 + 2*L^2*qsumk2*x^3 - 2*L^3*qsumk1*x^3 - 2*L^3*qsumk2*x^2 + 2*L^4*qsumk1*x^2 + 60*L*psum*x + 48*L*psumk1*x + 16*L*psumk2*x + 6*L*qsum*x + 16*L*qsumk1*x + 4*L*qsumk2*x + 152*L^2*k^2*x^2 + 20*L^2*k^2*x^3 - 42*L^3*k^2*x^2 - 2*L^3*k^2*x^3 + 4*L^4*k^2*x^2 + 66*L*k^2*x - 218*L*psum*x^2 + 40*L^2*psum*x - 38*L*psum*x^3 - 117*L^3*psum*x + 78*L^4*psum*x - 21*L^5*psum*x + 2*L^6*psum*x + 52*L*psumk1*x^2 - 100*L^2*psumk1*x - 100*L*psumk1*x^3 + 12*L*psumk2*x^2 - 28*L^2*psumk2*x + 70*L^3*psumk1*x - 28*L*psumk2*x^3 + 14*L^3*psumk2*x - 20*L^4*psumk1*x - 2*L^4*psumk2*x + 2*L^5*psumk1*x - 94*L*qsum*x^2 + 48*L^2*qsum*x + 2*L*qsum*x^3 - 52*L^3*qsum*x + 18*L^4*qsum*x - 2*L^5*qsum*x + 12*L*qsumk1*x^2 - 28*L^2*qsumk1*x - 28*L*qsumk1*x^3 + 2*L*qsumk2*x^2 - 6*L^2*qsumk2*x + 14*L^3*qsumk1*x - 6*L*qsumk2*x^3 + 2*L^3*qsumk2*x - 2*L^4*qsumk1*x))/((L - 2)^3*(L - 4)*(L^2 - 4*L + 3)^2*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(1/2)) - (x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(2*L^2*x - L^2 - 6*L*x^2 + 2*L*x + L + 4*x^3 - 2*x)*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2*(2*psum - 4*L + qsum - 3*L*psum - L*qsum - L*k^2 + L^2*psum + L^2 + 3*k^2 + 3))/(2*(L - 3)^3*(L^2 - 3*L + 2)^6*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
    den=1
    part2=0
  }
  if(k==5){
    num=((x - L + 1)*(75600*L - 151200*x + 30240*k^2*x^2 + 30240*k^2*x^3 + 10080*L*psum + 5040*L*qsum + 196740*L*x - 20160*psum*x - 10080*qsum*x + 15120*L*k^2 - 34956*L^2*psum + 48188*L^3*psum - 34305*L^4*psum + 13857*L^5*psum - 3282*L^6*psum + 450*L^7*psum - 33*L^8*psum + L^9*psum - 14958*L^2*qsum + 16615*L^3*qsum - 8845*L^4*qsum + 2506*L^5*qsum - 388*L^6*qsum + 31*L^7*qsum - L^8*qsum - 771480*L*x^2 + 123450*L^2*x - 75600*L*x^3 - 418250*L^3*x + 347605*L^4*x - 145120*L^5*x + 34290*L^6*x - 4640*L^7*x + 335*L^8*x - 10*L^9*x - 30240*k^2*x + 36960*psum*x^2 + 3360*psum*x^3 - 47040*psumk1*x^2 + 47040*psumk1*x^3 - 14400*psumk2*x^2 + 14400*psumk2*x^3 + 23520*qsum*x^2 - 3360*qsum*x^3 - 6720*qsumk1*x^2 + 6720*qsumk1*x^3 - 1800*qsumk2*x^2 + 1800*qsumk2*x^3 - 249570*L^2 + 324015*L^3 - 215750*L^4 + 81815*L^5 - 18350*L^6 + 2405*L^7 - 170*L^8 + 5*L^9 - 34794*L^2*k^2 + 30009*L^3*k^2 - 13141*L^4*k^2 + 3222*L^5*k^2 - 448*L^6*k^2 + 33*L^7*k^2 - L^8*k^2 + 302400*x^2 + 925350*L^2*x^2 + 98370*L^2*x^3 - 609605*L^3*x^2 - 51675*L^3*x^3 + 233495*L^4*x^2 + 14030*L^4*x^3 - 53130*L^5*x^2 - 2080*L^5*x^3 + 7060*L^6*x^2 + 160*L^6*x^3 - 505*L^7*x^2 - 5*L^7*x^3 + 15*L^8*x^2 - 99828*L*k^2*x^2 + 9570*L^2*k^2*x - 39348*L*k^2*x^3 - 33736*L^3*k^2*x + 19838*L^4*k^2*x - 5548*L^5*k^2*x + 830*L^6*k^2*x - 64*L^7*k^2*x + 2*L^8*k^2*x + 140076*L^2*psum*x^2 + 20176*L^2*psum*x^3 - 100385*L^3*psum*x^2 - 10387*L^3*psum*x^3 + 41021*L^4*psum*x^2 + 2808*L^4*psum*x^3 - 9792*L^5*psum*x^2 - 416*L^5*psum*x^3 + 1348*L^6*psum*x^2 + 32*L^6*psum*x^3 - 99*L^7*psum*x^2 - L^7*psum*x^3 + 3*L^8*psum*x^2 - 8800*L^2*psumk1*x^2 + 137880*L^2*psumk1*x^3 - 600*L^2*psumk2*x^2 - 63210*L^3*psumk1*x^2 + 38880*L^2*psumk2*x^3 - 74670*L^3*psumk1*x^3 - 19470*L^3*psumk2*x^2 + 52530*L^4*psumk1*x^2 - 19410*L^3*psumk2*x^3 + 22140*L^4*psumk1*x^3 + 14400*L^4*psumk2*x^2 - 18540*L^5*psumk1*x^2 + 5010*L^4*psumk2*x^3 - 3600*L^5*psumk1*x^3 - 4380*L^5*psumk2*x^2 + 3300*L^6*psumk1*x^2 - 630*L^5*psumk2*x^3 + 300*L^6*psumk1*x^3 + 600*L^6*psumk2*x^2 - 290*L^7*psumk1*x^2 + 30*L^6*psumk2*x^3 - 10*L^7*psumk1*x^3 - 30*L^7*psumk2*x^2 + 10*L^8*psumk1*x^2 + 44994*L^2*qsum*x^2 - 502*L^2*qsum*x^3 - 21536*L^3*qsum*x^2 + 52*L^3*qsum*x^3 + 5678*L^4*qsum*x^2 - 2*L^4*qsum*x^3 - 834*L^5*qsum*x^2 + 64*L^6*qsum*x^2 - 2*L^7*qsum*x^2 + 280*L^2*qsumk1*x^2 + 17200*L^2*qsumk1*x^3 + 270*L^2*qsumk2*x^2 - 8990*L^3*qsumk1*x^2 + 4290*L^2*qsumk2*x^3 - 8210*L^3*qsumk1*x^3 - 2400*L^3*qsumk2*x^2 + 6220*L^4*qsumk1*x^2 - 1890*L^3*qsumk2*x^3 + 1990*L^4*qsumk1*x^3 + 1500*L^4*qsumk2*x^2 - 1760*L^5*qsumk1*x^2 + 390*L^4*qsumk2*x^3 - 230*L^5*qsumk1*x^3 - 360*L^5*qsumk2*x^2 + 220*L^6*qsumk1*x^2 - 30*L^5*qsumk2*x^3 + 10*L^6*qsumk1*x^3 + 30*L^6*qsumk2*x^2 - 10*L^7*qsumk1*x^2 + 32952*L*psum*x + 47040*L*psumk1*x + 14400*L*psumk2*x + 6396*L*qsum*x + 6720*L*qsumk1*x + 1800*L*qsumk2*x + 99366*L^2*k^2*x^2 + 20670*L^2*k^2*x^3 - 46952*L^3*k^2*x^2 - 5612*L^3*k^2*x^3 + 12056*L^4*k^2*x^2 + 832*L^4*k^2*x^3 - 1728*L^5*k^2*x^2 - 64*L^5*k^2*x^3 + 130*L^6*k^2*x^2 + 2*L^6*k^2*x^3 - 4*L^7*k^2*x^2 + 39348*L*k^2*x - 105772*L*psum*x^2 + 6036*L^2*psum*x - 17252*L*psum*x^3 - 54214*L^3*psum*x + 52495*L^4*psum*x - 24070*L^5*psum*x + 6084*L^6*psum*x - 866*L^7*psum*x + 65*L^8*psum*x - 2*L^9*psum*x + 82040*L*psumk1*x^2 - 129080*L^2*psumk1*x - 129080*L*psumk1*x^3 + 23880*L*psumk2*x^2 - 38280*L^2*psumk2*x + 137880*L^3*psumk1*x - 38280*L*psumk2*x^3 + 38880*L^3*psumk2*x - 74670*L^4*psumk1*x - 19410*L^4*psumk2*x + 22140*L^5*psumk1*x + 5010*L^5*psumk2*x - 3600*L^6*psumk1*x - 630*L^6*psumk2*x + 300*L^7*psumk1*x + 30*L^7*psumk2*x - 10*L^8*psumk1*x - 48524*L*qsum*x^2 + 18654*L^2*qsum*x + 2132*L*qsum*x^3 - 29436*L^3*qsum*x + 17026*L^4*qsum*x - 4954*L^5*qsum*x + 774*L^6*qsum*x - 62*L^7*qsum*x + 2*L^8*qsum*x + 10760*L*qsumk1*x^2 - 17480*L^2*qsumk1*x - 17480*L*qsumk1*x^3 + 2760*L*qsumk2*x^2 - 4560*L^2*qsumk2*x + 17200*L^3*qsumk1*x - 4560*L*qsumk2*x^3 + 4290*L^3*qsumk2*x - 8210*L^4*qsumk1*x - 1890*L^4*qsumk2*x + 1990*L^5*qsumk1*x + 390*L^5*qsumk2*x - 230*L^6*qsumk1*x - 30*L^6*qsumk2*x + 10*L^7*qsumk1*x))
    den=((L - 2)^3*(L^2 - 4*L + 3)^2*(L^4 - 26*L^3 + 251*L^2 - 1066*L + 1680)*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(1/2))
    part2 =  (x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(2*L^2*x - L^2 - 6*L*x^2 + 2*L*x + L + 4*x^3 - 2*x)*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2*(2*psum - 20*L + qsum - 3*L*psum - L*qsum - L*k^2 + L^2*psum + 5*L^2 + 3*k^2 + 15))/(2*(L - 3)^3*(L^2 - 3*L + 2)^6*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
  }else{
    num = - ((x - L + 1)*(18*k*x - 18*k^2*x^3 - 9*L*k - 6*L*psum - 3*L*qsum - 18*k^2*x^2 + 12*psum*x + 6*qsum*x - 9*L*k^2 + 24*L^2*k - 22*L^3*k + 8*L^4*k - L^5*k + 17*L^2*psum - 17*L^3*psum + 7*L^4*psum - L^5*psum + 7*L^2*qsum - 5*L^3*qsum + L^4*qsum - 36*k*x^2 + 18*k^2*x + 2*psum*x^2 - 26*psum*x^3 - 12*psumk*x^2 + 12*psumk*x^3 + 4*qsum*x^2 - 16*qsum*x^3 - 6*qsumk*x^2 + 6*qsumk*x^3 + 15*L^2*k^2 - 7*L^3*k^2 + L^4*k^2 + 48*L*k^2*x^2 - 61*L^2*k*x^2 - 16*L^2*k^2*x + 12*L*k^2*x^3 - 6*L^2*k*x^3 + 23*L^3*k*x^2 + 12*L^3*k^2*x + L^3*k*x^3 - 3*L^4*k*x^2 - 2*L^4*k^2*x - 67*L^2*psum*x^2 - 20*L^2*psum*x^3 + 33*L^3*psum*x^2 + 3*L^3*psum*x^3 - 5*L^4*psum*x^2 + 10*L^2*psumk*x^2 + 12*L^2*psumk*x^3 - 10*L^3*psumk*x^2 - 2*L^3*psumk*x^3 + 2*L^4*psumk*x^2 - 26*L^2*qsum*x^2 - 4*L^2*qsum*x^3 + 6*L^3*qsum*x^2 + 6*L^2*qsumk*x^2 + 2*L^2*qsumk*x^3 - 2*L^3*qsumk*x^2 - 12*L*k*x - 36*L*psum*x + 12*L*psumk*x - 18*L*qsum*x + 6*L*qsumk*x - 26*L^2*k^2*x^2 - 2*L^2*k^2*x^3 + 4*L^3*k^2*x^2 + 69*L*k*x^2 - 12*L*k^2*x - 25*L^2*k*x + 9*L*k*x^3 + 36*L^3*k*x - 15*L^4*k*x + 2*L^5*k*x + 41*L*psum*x^2 + 19*L^2*psum*x + 41*L*psum*x^3 + 12*L^3*psum*x - 11*L^4*psum*x + 2*L^5*psum*x + 10*L*psumk*x^2 - 22*L^2*psumk*x - 22*L*psumk*x^3 + 12*L^3*psumk*x - 2*L^4*psumk*x + 20*L*qsum*x^2 + 6*L^2*qsum*x + 18*L*qsum*x^3 + 6*L^3*qsum*x - 2*L^4*qsum*x + 2*L*qsumk*x^2 - 8*L^2*qsumk*x - 8*L*qsumk*x^3 + 2*L^3*qsumk*x))
    den =((L - 2)^3*(L^2 - 4*L + 3)^2*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(1/2))
    part2 = (x^2*(x - 1)^2*(L^2 - 2*L*x - L + x^2 + x)^2*(2*L^2*x - L^2 - 6*L*x^2 + 2*L*x + L + 4*x^3 - 2*x)*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^3)/(2*(L - 3)^3*(L^2 - 3*L + 2)^6*((x^2*(x - 1)^2*(x - L - 2*L*x + L^2 + x^2)^2*(3*k + 2*psum + qsum - 4*L*k - 3*L*psum - L*qsum - L*k^2 + L^2*k + L^2*psum + 3*k^2)^2)/((L - 3)^2*(L^2 - 3*L + 2)^4))^(3/2))
  }
  return(num/den - part2)
}

C1_d = function(x,L,k,qsum){
  if(k==1){
    result = (L*x^2*(L - x)^2*(- k^2 + k + qsum)^2*(- k^2 + qsum + 1))/(2*(L - 1)^3*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2))
  }
  if(k==5){
    result = (L*x^2*(L - x)^2*(- k^2 + k + qsum)^2*(- k^2 + qsum + 5))/(2*(L - 1)^3*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2))
  }
  else{
    result = (L*x^2*(L - x)^2*(- k^2 + k + qsum)^3)/(2*(L - 1)^3*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2))
  }
  return(result)
}

C2_d = function(x,L,k,psum,qsum,qsumk,psumk1,qsumk1,psumk2,qsumk2){
  if(k==1){
    result=-(x^2*(L - x)^2*(- k^2 + k + qsum)^2*(48*k^2*x^2 + 144*L*x - 12*L^2*qsum + 19*L^3*qsum - 8*L^4*qsum + L^5*qsum + 204*L*x^2 - 204*L^2*x + 82*L^3*x - 10*L^4*x - 144*qsum*x^2 + 32*qsumk1*x^2 + 8*qsumk2*x^2 - 12*L^2 + 19*L^3 - 8*L^4 + L^5 + 12*L^2*k^2 - 19*L^3*k^2 + 8*L^4*k^2 - L^5*k^2 - 144*x^2 - 82*L^2*x^2 + 10*L^3*x^2 - 100*L*k^2*x^2 + 100*L^2*k^2*x - 46*L^3*k^2*x + 6*L^4*k^2*x - 82*L^2*qsum*x^2 + 10*L^3*qsum*x^2 + 28*L^2*qsumk1*x^2 + 4*L^2*qsumk2*x^2 - 4*L^3*qsumk1*x^2 + 144*L*qsum*x - 32*L*qsumk1*x - 8*L*qsumk2*x + 46*L^2*k^2*x^2 - 6*L^3*k^2*x^2 - 48*L*k^2*x + 204*L*qsum*x^2 - 204*L^2*qsum*x + 82*L^3*qsum*x - 10*L^4*qsum*x - 56*L*qsumk1*x^2 + 56*L^2*qsumk1*x - 12*L*qsumk2*x^2 + 12*L^2*qsumk2*x - 28*L^3*qsumk1*x - 4*L^3*qsumk2*x + 4*L^4*qsumk1*x))/(2*(L - 1)^4*(L^3 - 9*L^2 + 26*L - 24)*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2))
  }
  if(k==5){
    result=-(x^2*(L - x)^2*(- k^2 + k + qsum)^2*(6720*k^2*x^2 + 100800*L*x - 1680*L^2*qsum + 2746*L^3*qsum - 1317*L^4*qsum + 277*L^5*qsum - 27*L^6*qsum + L^7*qsum + 147960*L*x^2 - 147960*L^2*x + 68360*L^3*x - 14110*L^4*x + 1360*L^5*x - 50*L^6*x - 20160*qsum*x^2 + 4480*qsumk1*x^2 + 1200*qsumk2*x^2 - 8400*L^2 + 13730*L^3 - 6585*L^4 + 1385*L^5 - 135*L^6 + 5*L^7 + 1680*L^2*k^2 - 2746*L^3*k^2 + 1317*L^4*k^2 - 277*L^5*k^2 + 27*L^6*k^2 - L^7*k^2 - 100800*x^2 - 68360*L^2*x^2 + 14110*L^3*x^2 - 1360*L^4*x^2 + 50*L^5*x^2 - 14344*L*k^2*x^2 + 14344*L^2*k^2*x - 7400*L^3*k^2*x + 1610*L^4*k^2*x - 160*L^5*k^2*x + 6*L^6*k^2*x - 13672*L^2*qsum*x^2 + 2822*L^3*qsum*x^2 - 272*L^4*qsum*x^2 + 10*L^5*qsum*x^2 + 8080*L^2*qsumk1*x^2 + 1980*L^2*qsumk2*x^2 - 2780*L^3*qsumk1*x^2 - 600*L^3*qsumk2*x^2 + 400*L^4*qsumk1*x^2 + 60*L^4*qsumk2*x^2 - 20*L^5*qsumk1*x^2 + 20160*L*qsum*x - 4480*L*qsumk1*x - 1200*L*qsumk2*x + 7400*L^2*k^2*x^2 - 1610*L^3*k^2*x^2 + 160*L^4*k^2*x^2 - 6*L^5*k^2*x^2 - 6720*L*k^2*x + 29592*L*qsum*x^2 - 29592*L^2*qsum*x + 13672*L^3*qsum*x - 2822*L^4*qsum*x + 272*L^5*qsum*x - 10*L^6*qsum*x - 10160*L*qsumk1*x^2 + 10160*L^2*qsumk1*x - 2640*L*qsumk2*x^2 + 2640*L^2*qsumk2*x - 8080*L^3*qsumk1*x - 1980*L^3*qsumk2*x + 2780*L^4*qsumk1*x + 600*L^4*qsumk2*x - 400*L^5*qsumk1*x - 60*L^5*qsumk2*x + 20*L^6*qsumk1*x))/(2*(L - 1)^4*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2)*(L^5 - 28*L^4 + 303*L^3 - 1568*L^2 + 3812*L - 3360))
  }
  else{
    result = -(x^2*(L - x)^2*(- k^2 + k + qsum)^2*(4*k^2*x^2 - L^2*k + L^3*k - L^2*qsum + L^3*qsum - 12*k*x^2 - 4*qsumk*x^2 + L^2*k^2 - L^3*k^2 - 6*L*k^2*x^2 + 6*L^2*k^2*x + 12*L*k*x + 4*L*qsumk*x + 10*L*k*x^2 - 4*L*k^2*x - 10*L^2*k*x + 2*L*qsum*x^2 - 2*L^2*qsum*x + 4*L*qsumk*x^2 - 4*L^2*qsumk*x))/(2*(L - 1)^4*(L - 2)*((x^2*(L - x)^2*(- k^2 + k + qsum)^2)/(L - 1)^2)^(3/2))
  }
  return(result)
}

T3.lambda = function(b, n0, n1, L, k, psum, qsum, psumk, qsumk){
  C3 = C1_Z(n0:n1,L,k,psum,qsum)
  C4 = C2_Z(n0:n1,L,k,psum,qsum,psumk,qsumk)
  dnorm(b)*b^3*sum(C3*C4*Nu(sqrt(2*C3*b^2))*Nu(sqrt(2*C4*b^2)))
}

T3.lambdaZw = function(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp){
  if(asymp==TRUE){
     C1 = C1_w_asy((n0:n1)/L)
     C2 = C2_w_asy((n0:n1)/L,k,psum,psumk1)
     C2[C2<0] = 0.00000001
  }
  if(asymp==FALSE){
     C1  = C1_w(n0:n1,L,k,psum,qsum)
     C2 = C2_w(n0:n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2,qsumk2)

     C2[C2<0] = 0.00000001
  }

  dnorm(b)*b^3*sum(C1*C2*Nu(sqrt(2*C1*b^2))*Nu(sqrt(2*C2*b^2)))
}

T3.lambdaZdiff = function(b,n0,n1,L,k,psum,qsum,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp){
  if(asymp==TRUE){
    C1 = C1_d_asy((n0:n1)/L)
    C2 = C2_d_asy((n0:n1)/L,k,qsum,qsumk1)
    C2[C2<0] = 0.00000001
  }
  if(asymp==FALSE){
    C1 =C1_d(n0:n1,L,k,qsum)
    C2 =C2_d(n0:n1,L,k,psum,qsum,qsumk,psumk1,qsumk1,psumk2,qsumk2)

    C2[C2<0] = 0.00000001

  }

  nu1 = Nu(sqrt(2*C1*b^2))
  nu2 = Nu(sqrt(2*C2*b^2))
  nu1[is.na(nu1)]=0
  nu2[is.na(nu2)]=0

  dnorm(b)*b^3*sum(C1*C2*Nu(sqrt(2*C1*b^2))*Nu(sqrt(2*C2*b^2)))
}

T3.lambdaM = function(D,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp){
  pval_Zd = T3.lambdaZdiff(b,n0,n1,L,k,psum,qsum,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
  pval_Zw = T3.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)

  return(1-(1-D*2*pval_Zd)*(1-D*pval_Zw))
}

T3.lambdaS = function(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp){
  integrand = function(t,w){
    if(asymp==TRUE){
      C1w  = C1_w_asy(t/L)
      C2w = C2_w_asy(t/L,k,psum,psumk1)
      C1d  =C1_d_asy(t/L)
      C2d = C2_d_asy(t/L,k,qsum,qsumk1)
      C1d[C1d<0]=0
      C2d[C2d<0]=0
      C1w[C1w<0]=0
      C2w[C2w<0]=0
     }
   if(asymp==FALSE){
      C1w  = C1_w(t,L,k,psum,qsum)
      C2w = C2_w(t,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2)
      C1d  =C1_d(t,L,k,qsum)
      C2d = C2_d(t,L,k,psum,qsum,qsumk,psumk1,qsumk1, psumk2 ,qsumk2)
      C1d[C1d<0]=0
      C2d[C2d<0]=0
      C1w[C1w<0]=0
      C2w[C2w<0]=0
    }

    nu1 = Nu(sqrt(2*b*(C1d*cos(w)^2+C1w*sin(w)^2)))
    nu2 = Nu(sqrt(2*b*(C2d*cos(w)^2+C2w*sin(w)^2)))
    (4*(C1d*cos(w)^2+C1w*sin(w)^2)*(C2d*cos(w)^2+C2w*sin(w)^2)*b^2*nu1*nu2)/(2*pi)
  }
  integrand0 = function(t) {integrate(integrand,0,2*pi,t=t,subdivisions=3000, stop.on.error=FALSE)$value}
  result=Vectorize(integrand0)
  dchisq(b,2)*integrate(result, n0, n1, subdivisions=3000, stop.on.error=FALSE)$value
}

EX3_new.f = function(t, n, k, psum, deg.sumsq, deg.sum3, daa, dda, aaa1, aaa2){
  x1 = 2*k*n+6*k*n*psum/k
  x2 = 3*k^2*n+deg.sumsq+2*k*n*psum+2*daa-x1
  x3 = 4*k^2*n^2*(1+psum/k)-4*(3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa)+2*(2*k*n+6*k*n*psum/k)
  x4 = 4*k^3*n+3*k*deg.sumsq+deg.sum3-3*(3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa)+2*(2*k*n+6*k*n*psum/k)
  x5 = 4*k^3*n+2*k*deg.sumsq+2*dda-2*(3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa)+x1-2*aaa1-6*aaa2
  x6 = 2*aaa1+6*aaa2
  x7 = 2*k*n*(3*k^2*n+deg.sumsq)-4*k^2*n^2*(1+psum/k) - 4*(4*k^3*n+2*k*deg.sumsq+2*dda-(3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa))-2*(4*k^3*n+3*k*deg.sumsq+deg.sum3-(3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa))+ 4*((3*k^2*n+deg.sumsq+2*k^2*n*psum/k+2*daa)-(2*k*n+6*k*n*psum/k))+2*x6
  x8 = 8*k^3*n^3-4*x1-24*x2-6*x3-8*x4-24*x5-8*x6-12*x7
  p1 = t*(t-1)/n/(n-1)
  p2 = t*(t-1)*(t-2)/n/(n-1)/(n-2)
  p3 = t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))
  p4 = p5 = t*(t-1)*(t-2)*(t-3)/(n*(n-1)*(n-2)*(n-3))
  p6 = p2
  p7 = t*(t-1)*(t-2)*(t-3)*(t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4))
  p8 = t*(t-1)*(t-2)*(t-3)*(t-4)*(t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

  q1 = (n-t)*(n-t-1)/n/(n-1)
  q2 = (n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)
  q3 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)/(n*(n-1)*(n-2)*(n-3))
  q4 = q5 = q3
  q6 = q2
  q7 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)/(n*(n-1)*(n-2)*(n-3)*(n-4))
  q8 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(n-t-4)*(n-t-5)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))

  A11 = 4*x1*p1 + 24*x2*p2 + 6*x3*p3 + 8*x4*p4 + 24*x5*p5 + 8*x6*p6 + 12*x7*p7 + x8*p8
  A12 = 2*x3*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)) + 4*x7*t*(t-1)*(t-2)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) + x8*t*(t-1)*(t-2)*(t-3)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  A21 = 2*x3*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3)) + 4*x7*(n-t)*(n-t-1)*(n-t-2)*(t)*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)) + x8*(n-t)*(n-t-1)*(n-t-2)*(n-t-3)*(t)*(t-1)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  A22 = 4*x1*q1 + 24*x2*q2 + 6*x3*q3 + 8*x4*q4 + 24*x5*q5 + 8*x6*q6 + 12*x7*q7 + x8*q8

  q = (n-t-1)/(n-2)
  p = (t-1)/(n-2)
  ERw3 = q^3*A11 + 3*q^2*p*A12 + 3*q*p^2*A21 + p^3*A22
  ERd3 =  A11 - 3*A12 + 3*A21 - A22
  list(A11=A11, A12=A12, A21=A21, A22=A22, ERw3=ERw3, ERd3=ERd3)
}

EX3.f = function(n, t, k, vn, deg.sumsq, deg.sum3, daa, dda, aaa1, aaa2){
  x1 = 2*k*n+6*k*n*vn
  x2 = 3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa-x1
  x3 = 4*k^2*n^2*(1+vn)-4*(3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa)+2*(2*k*n+6*k*n*vn)
  x4 = 4*k^3*n+3*k*deg.sumsq+deg.sum3-3*(3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa)+2*(2*k*n+6*k*n*vn)
  x5 = 4*k^3*n+2*k*deg.sumsq+2*dda-2*(3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa)+x1-2*aaa1-6*aaa2
  x6 = 2*aaa1+6*aaa2
  x7 = 2*k*n*(3*k^2*n+deg.sumsq)-4*k^2*n^2*(1+vn) - 4*(4*k^3*n+2*k*deg.sumsq+2*dda-(3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa))-2*(4*k^3*n+3*k*deg.sumsq+deg.sum3-(3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa))+ 4*((3*k^2*n+deg.sumsq+2*k^2*n*vn+2*daa)-(2*k*n+6*k*n*vn))+2*x6
  x8 = 8*k^3*n^3-4*x1-24*x2-6*x3-8*x4-24*x5-8*x6-12*x7
  p1 = 2*t*(n-t)/n/(n-1)
  p2 = p1/2
  p3 = 4*t*(t-1)*(n-t)*(n-t-1)/(n*(n-1)*(n-2)*(n-3))
  p4 = t*(n-t)*((n-t-1)*(n-t-2)+(t-1)*(t-2))/(n*(n-1)*(n-2)*(n-3))
  p5 = p7 = p3/2
  p8 = 8*t*(t-1)*(t-2)*(n-t)*(n-t-1)*(n-t-2)/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5))
  4*x1*p1 + 24*x2*p2 + 6*x3*p3 + 8*x4*p4 + 24*x5*p5 + 12*x7*p7 + x8*p8
}

ERw = function(n,t,k){
  q=(n-t-1)/(n-2)
  p=(t-1)/(n-2)
  ER1 = 2*k*t*(t-1)/(n-1)
  ER2 = 2*k*(n-t)*(n-t-1)/(n-1)
  result = q*ER1 + p*ER2
  return(result)
}

ERd = function(n,t,k){
  result = 2*k*t*(t-1)/(n-1)- 2*k*(n-t)*(n-t-1)/(n-1)
  return(result)
}

varRw = function(n,t,k,psum,deg.sumsq){
  vn= psum/k
  config1 = (2*k*n + 2*k*n*vn)
  config2 = (3*k^2*n + deg.sumsq -2*k*n -2*k*n*vn)
  config3 = (4*n^2*k^2 + 4*k*n + 4*k*n*vn - 12*k^2*n - 4*deg.sumsq)
  f11 = 2*(t)*(t-1)/n/(n-1)
  f21 = 4*(t)*(t-1)*(t-2)/n/(n-1)/(n-2)
  f31 = (t)*(t-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
  V1 = config1*f11 + config2*f21  + config3*f31 - (2*k*t*(t-1)/(n-1))^2

  f12 = 2*(n-t)*(n-t-1)/n/(n-1)
  f22 = 4*(n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)
  f32 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)/n/(n-1)/(n-2)/(n-3)
  V2 = config1*f12 + config2*f22  + config3*f32 - (2*k*(n-t)*(n-t-1)/(n-1))^2

  P3 = (t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3))
  V12 = config3*P3 - (2*k*t*(t-1)/(n-1))*(2*k*(n-t)*(n-t-1)/(n-1))

  q=(n-t-1)/(n-2)
  p=(t-1)/(n-2)
  result = q^2*V1  + p^2*V2 + 2*p*q*V12
  return(result)
}

varRd = function(n,t,k,psum,deg.sumsq){
  vn= psum/k
  config1 = (2*k*n + 2*k*n*vn)
  config2 = (3*k^2*n + deg.sumsq -2*k*n -2*k*n*vn)
  config3 = (4*n^2*k^2 + 4*k*n + 4*k*n*vn - 12*k^2*n - 4*deg.sumsq)
  f11 = 2*(t)*(t-1)/n/(n-1)
  f21 = 4*(t)*(t-1)*(t-2)/n/(n-1)/(n-2)
  f31 = (t)*(t-1)*(t-2)*(t-3)/n/(n-1)/(n-2)/(n-3)
  V1 = config1*f11 + config2*f21  + config3*f31 - (2*k*t*(t-1)/(n-1))^2

  f12 = 2*(n-t)*(n-t-1)/n/(n-1)
  f22 = 4*(n-t)*(n-t-1)*(n-t-2)/n/(n-1)/(n-2)
  f32 = (n-t)*(n-t-1)*(n-t-2)*(n-t-3)/n/(n-1)/(n-2)/(n-3)
  V2 = config1*f12 + config2*f22  + config3*f32 - (2*k*(n-t)*(n-t-1)/(n-1))^2

  P3 = (t*(t-1)*(n-t)*(n-t-1))/(n*(n-1)*(n-2)*(n-3))
  V12 = config3*P3 - (2*k*t*(t-1)/(n-1))*(2*k*(n-t)*(n-t-1)/(n-1))

  result = V1  + V2 - 2*V12
  return(result)
}

varR = function(n,t,k,psum,qsum){
  vn=psum/k
  cn=(qsum+k-k^2)/k
  h = 4*(t-1)*(n-t-1)/(n-2)/(n-3)
  4*k*t*(n-t)/(n-1)*(h*(1+vn-2*k/(n-1))+(1-h)*cn)
}

T3.skewed.lambda = function(b, n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1, aaa2, daa, dda){
  C3 = C1_Z(n0:n1,L,k,psum,qsum)
  C4 = C2_Z(n0:n1,L,k,psum,qsum,psumk,qsumk)

  n = L
  ts = 1:(n-1)
  vn=psum/k
  EX = 4*k*ts*(n-ts)/(n-1)
  EX2 = 4*k*(1+vn)*2*ts*(n-ts)/(n-1) + 4*(3*k^2*n+deg.sumsq-2*k*n*(1+vn))*ts*(n-ts)/n/(n-1) + (4*k^2*n^2-4*(3*k^2*n+deg.sumsq)+4*k*n*(1+vn))*4*ts*(ts-1)*(n-ts)*(n-ts-1)/(n*(n-1)*(n-2)*(n-3))
  EX3 = EX3.f(n,ts,k,vn, deg.sumsq,deg.sum3,daa,dda,aaa1,aaa2)
  VX = EX2-EX^2
  gamma = -(EX3-3*EX*VX-EX^3)/(VX^(3/2))
  theta = rep(0,n-1)
  pos = which(1+2*gamma*b>0)
  theta[pos] = (sqrt((1+2*gamma*b)[pos])-1)/gamma[pos]
  S = (1+gamma*theta)^(-1/2)*exp((b-theta)^2/2 + gamma*theta^3/6)
  nn = n-length(pos)
  if (nn>0.75*n){
    print("Not enough points for extrapolation!")
    return(0)
  }
  if (nn>=2*n0){
    neg = which(1+2*gamma*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    if (nn<n){
      id2 = id1 + ceiling(0.02*n)  # use this
      id3 = id2 + ceiling(0.02*n)
      inc = (S[id3]-S[id2])/(id3-id2)
      S[id2:1] = S[id2+1]-inc*(1:id2)
      S[(n-id2):(n-1)] = S[id2:1]
    }else{
      ymax = S[ceiling(n/2)]
      id = id1+0.05*n
      a = (ymax-S[id])/(id-n/2)^2
      S[1:id] = ymax-a*((1:id)-n/2)^2
      S[(n-id):(n-1)] = S[id:1]
    }
    neg2 = which(S<0)
    S[neg2] = 0
  }
  dnorm(b)*b^3*sum(S[(n-n0):(n-n1)]*C3*C4*Nu(sqrt(2*C3*b^2))*Nu(sqrt(2*C4*b^2)))
}

T3.skewed.lambdaZw = function(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda){
  C1  = C1_w(n0:n1,L,k,psum,qsum)
  C2 = C2_w(n0:n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2,qsumk2 )
  n = L
  ts = 1:(n-1)

  C2[C2<0] = 0.00000001
  EX = ERw(L,ts,k)
  EX3 = EX3_new.f(ts, L, k, psum, deg.sumsq, deg.sum3, daa, dda, aaa1, aaa2)$ERw3
  VX = varRw(L,ts,k,psum,deg.sumsq)
  gamma = (EX3-3*EX*VX-EX^3)/(VX^(3/2))
  theta = rep(0,n-1)
  pos = which(1+2*gamma*b>0)
  theta[pos] = (sqrt((1+2*gamma*b)[pos])-1)/gamma[pos]
  S = (1+gamma*theta)^(-1/2)*exp((b-theta)^2/2 + gamma*theta^3/6)
  nn = n-length(pos)
  if (nn>0.75*n){
    print("Not enough points for extrapolation!")
    return(0)
  }
  if (nn>=(n0-1)+(n-n0)){
    #print("extrapolate")
    neg = which(1+2*gamma*b<=0)
    dif = neg[2:nn]-neg[1:(nn-1)]
    id1 = which.max(dif)
    id2 = id1 + ceiling(0.03*n)
    id3 = id2 + ceiling(0.09*n)
    inc = (S[id3]-S[id2])/ceiling(0.09*n)
    S[id2:1] = S[id2+1]-inc*(1:id2)
    S[(n/2+1):n] = S[(n/2):1]
    neg2 = which(S<0)
    S[neg2] = 0
  }
  dnorm(b)*b^3*sum(S[(n-n0):(n-n1)]*C1*C2*Nu(sqrt(2*C1*b^2))*Nu(sqrt(2*C2*b^2)))
}

T3.skewed.lambdaZd = function(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3,aaa1,aaa2,daa,dda){
  C1  = C1_d(n0:n1,L,k,qsum)
  C2 = C2_d(n0:n1,L,k,psum,qsum,qsumk,psumk1,qsumk1,psumk2,qsumk2)

  C1[C1<0] = 0.0001
  C2[C2<0] = 0.0001
  n = L
  ts = 1:(n-1)

  EX = ERd(L,ts,k)
  EX3 = EX3_new.f(ts, L, k, psum, deg.sumsq, deg.sum3, daa, dda, aaa1, aaa2)$ERd3
  VX = varRd(L,ts,k,psum,deg.sumsq)
  gamma = (EX3-3*EX*VX-EX^3)/(VX^(3/2))
  theta = rep(0,n-1)
  pos = which(1+2*gamma*b>0)
  theta[pos] = (sqrt((1+2*gamma*b)[pos])-1)/gamma[pos]
  S = (1+gamma*theta)^(-1/2)*exp((b-theta)^2/2 + gamma*theta^3/6)
  S[n/2]=S[n/2-1]
  nn = n-length(pos)

  nn.l = ceiling(n/2)-length(which(1+2*gamma[1:ceiling(n/2)]*b>0))
  nn.r = ceiling(n/2)-length(which(1+2*gamma[ceiling(n/2+1):n]*b>0))
  # print(nn.l)
  # print(nn.r)
  if (nn>0.75*n){
    print("Not enough points for extrapolation!")
    return(0)
  }
  if (nn.r>=(n-n1)){
    #print('extrapolate')
    neg = which(1+2*gamma[ceiling(n/2+1):(n-1)]*b<=0)
    id1 = neg[1]+ceiling(n/2)-1
    id2 = id1 - ceiling(0.06*n)
    id3 = id2 - ceiling(0.06*n)
    inc = (S[id3]-S[id2])/(id3-id2)
    S[id2:n] = S[id2]+inc*((id2:n)-id2)+inc^2*((id2:n)-id2)+inc^3*((id2:n)-id2)

    if(n0<=0.15*200){
      S[S<0]=0
    }
  }
  dnorm(b)*b^3*sum(S[(n-n0):(n-n1)]*C1*C2*Nu(sqrt(2*C1*b^2))*Nu(sqrt(2*C2*b^2)))
}

T3.skewed.lambdaM = function(D,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda){
  pval_Zd = T3.skewed.lambdaZd(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
  pval_Zw = T3.skewed.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)

  return(1-(1-D*2*pval_Zd)*(1-D*pval_Zw))
}

getbZ = function(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=5,skew.corr=FALSE,dif=1e-10, nIterMax=100){
  m0=ARL*alpha
  if(skew.corr==FALSE){
    pm = T3.lambda(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.lambda(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
    }
    pM = T3.lambda(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.lambda(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
    }
    b = (bmin+bmax)/2
    p = T3.lambda(b,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.lambda(b,n0,n1,L,k,psum,qsum,psumk,qsumk)*m0
      nIter = nIter + 1
    }
    return(b)
  }else{
    pm = T3.skewed.lambda(bmin,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1,aaa2, daa, dda)*m0
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.skewed.lambda(bmin,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1, aaa2, daa, dda)*m0
    }
    pM = T3.skewed.lambda(bmax,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1,aaa2, daa, dda)*m0
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.skewed.lambda(bmax,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1,aaa2, daa, dda)*m0
    }
    b = (bmin+bmax)/2
    p = T3.skewed.lambda(b,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1,aaa2, daa, dda)*m0
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.skewed.lambda(b,n0, n1, L, k, psum, qsum, psumk, qsumk, deg.sumsq, deg.sum3, aaa1,aaa2, daa, dda)*m0
      nIter = nIter + 1
    }
    return(b)
  }
}

getbZw = function(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda, bmin=3,bmax=5,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp){
  m0=ARL*alpha
  if(skew.corr==FALSE){

    pm = T3.lambdaZw(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.lambdaZw(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    }
    pM = T3.lambdaZw(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.lambdaZw(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    }
    b = (bmin+bmax)/2
    p = T3.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
      nIter = nIter + 1
    }
    return(b)
  }else{
    pm = T3.skewed.lambdaZw(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.skewed.lambdaZw(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
    }
    pM = T3.skewed.lambdaZw(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.skewed.lambdaZw(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
    }
    b = (bmin+bmax)/2
    p = T3.skewed.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk1,psumk,qsumk,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.skewed.lambdaZw(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1, psumk2, qsumk2, deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)*m0
      nIter = nIter + 1
    }
    return(b)
  }

}

getbS = function(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=5,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp){
  m0=ARL*alpha
  pm = T3.lambdaS(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
  while (pm<alpha){
    bmin = bmin-1
    pm = T3.lambdaS(bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
  }
  pM = T3.lambdaS(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
  while (pM>alpha){
    bmax = bmax+1
    pM = T3.lambdaS(bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
  }
  b = (bmin+bmax)/2
  p = T3.lambdaS(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
  nIter = 1
  while (abs(p-alpha)>dif && nIter<nIterMax){
    if (p<alpha){
      bmax = b
    }else{
      bmin = b
    }
    b = (bmin+bmax)/2
    p = T3.lambdaS(b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,asymp)*m0
    nIter = nIter + 1
  }
  b
}

getbM = function(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=21,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp){
  m0=ARL*alpha
  if(skew.corr==FALSE){
    pm = T3.lambdaM(m0,bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.lambdaM(m0,bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
    }
    pM = T3.lambdaM(m0,bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.lambdaM(m0,bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
    }
    b = (bmin+bmax)/2
    p = T3.lambdaM(m0,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.lambdaM(m0,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,asymp)
      nIter = nIter + 1
    }
    return(b)
  }else{
    pm = T3.skewed.lambdaM(m0,bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
    while (pm<alpha){
      bmin = bmin-1
      pm = T3.skewed.lambdaM(m0,bmin,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
    }
    pM = T3.skewed.lambdaM(m0,bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
    while (pM>alpha){
      bmax = bmax+1
      pM = T3.skewed.lambdaM(m0,bmax,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
    }
    b = (bmin+bmax)/2
    p = T3.skewed.lambdaM(m0,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
    nIter = 1
    while (abs(p-alpha)>dif && nIter<nIterMax){
      if (p<alpha){
        bmax = b
      }else{
        bmin = b
      }
      b = (bmin+bmax)/2
      p = T3.skewed.lambdaM(m0,b,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2, qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda)
      nIter = nIter + 1
    }
    return(b)
  }
}
#need user to specify ARL and probability of early stop error (for example, ARL = 10,000 and alpha = 0.01 implies D = 100)

getb = function(distM,ARL,alpha,N0,n0,n1,L,k,statistics,skew.corr,dif=1e-10, nIterMax=100,asymp){

    quantities = gb_quantities(distM,N0,k)
    psum = quantities$psum
    qsum = quantities$qsum
    psumk = quantities$psumk
    qsumk = quantities$qsumk
    psumk1 = quantities$psumk1
    qsumk1 = quantities$qsumk1
    psumk2 = quantities$psumk2
    qsumk2 = quantities$qsumk2
    deg.sumsq = quantities$deg.sumsq
    deg.sum3 = quantities$deg.sum3.n
    aaa1 = quantities$aaa1.n
    aaa2 = quantities$aaa2.n
    dda = quantities$dda.n
    daa = quantities$daa.n

  output=list()
  if (skew.corr==FALSE){
    if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
      output$ori = getbZ(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=5,skew.corr=FALSE,dif=1e-10, nIterMax=100)
    }
    if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
      output$weighted = getbZw(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda, bmin=3,bmax=5,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp)

    }
    if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
      output$max.type = getbM(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=7,bmax=21,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp)
    }
    if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
      output$generalized = getbS(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=20,bmax=30,skew.corr=FALSE,dif=1e-10, nIterMax=100,asymp)
    }
  return(output)
  }

  #skewness correction
  if (length(which(!is.na(match(c("o","ori","original","all"), statistics))))>0){
    output$ori = getbZ(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=5,skew.corr=TRUE)
  }
  if (length(which(!is.na(match(c("w","weighted","all"), statistics))))>0){
    output$weighted = getbZw(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=3,bmax=5,skew.corr=TRUE,dif=1e-10, nIterMax=100,asymp)
  }
  if (length(which(!is.na(match(c("m","max","all"), statistics))))>0){
    output$max.type = getbM(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=8,bmax=20,skew.corr=TRUE,dif=1e-10, nIterMax=100,asymp)
  }
  if (length(which(!is.na(match(c("g","generalized","all"), statistics))))>0){
    output$generalized = getbS(ARL,alpha,n0,n1,L,k,psum,qsum,psumk,qsumk,psumk1,qsumk1,psumk2,qsumk2,deg.sumsq,deg.sum3, aaa1,aaa2,daa,dda,bmin=20,bmax=30,skew.corr=TRUE,dif=1e-10, nIterMax=100,asymp)
  }
  return(output)


}





