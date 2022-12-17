

## Pointwise version
GCT.tstat = function(X, khat, lrv.ba, type=1){
  p = ncol(X);
  X1 = X[(1:khat),]; X2 =X[-(1:khat),];
  if(missing(lrv.ba)){
#    Xtilde = hdbinseg::clean.cp(t(X), type="sbs")$x;
#    lrv.ba = apply(Xtilde, 1, longrun_bartlett);
#    lrv.ba = sqrt(lrv.ba);
    lrv.ba = LRV(X1, X2, "bartlett");
    lrv.ba = sqrt(lrv.ba);
  }

  xtilde = t(X)/lrv.ba; # xtilde is dim*length
  tmp = hdbinseg:::func_dc(xtilde, .5);
  aa = tmp$acs;
  bb = aa[,khat]; bb2 = bb*bb;

  if(type ==1 ){    Gn = max(apply(aa^2, 2, max));}
  if(type ==2 ){    Gn = max(apply(aa, 2, max));}
  if(type == 3){    Gn = max(apply(aa^2, 2, sum));} ## Sum binary
  if(type == 4){    Gn = max(tmp$res);} ## DC
  return(Gn)
}

## Binary Segmentation version
GCT.break = function(X, c=1, wd, B=299, bsize, hsize, spatial = spatial, do.parallel, n.cl, type=1){

  length = n = nrow(X); p = ncol(X);
  khat0  = lse.break(X, trim=wd, type=type);

  if(missing(bsize)){
    K.d = apply(X, 2, blocklength.andrew);
#    bsize = min(1.5*median(K.d), log(n*p));
    bsize = c*median(K.d);
    bsize = floor(bsize);
  }
  if(missing(hsize)){
    L.d = apply(X, 1, blocklength.andrew);
##    hsize = min(1.5*median(L.d), log(n*p));
    hsize = c*median(L.d);
    hsize = floor(hsize);
  }

  if(!spatial){hsize = p};

  GCT.novote2 = function(X, khat, B, bsize, hsize, do.parallel, n.cl, type){
    ## Implement several GCT-based tests
    ## Two LRV are usied; sd and bartlett
    ## do not use voting

    n <- nrow(X);    p <- ncol(X);
    tstat = GCT.tstat(X, khat=khat, type=type);

    blength = floor(n/bsize)+1;
    hlength = floor(p/hsize)+1;

    X1 = Demean(X, colMeans(X));

    if(do.parallel){
    library(doParallel)
    cl <- parallel::makeCluster(n.cl, setup_strategy = "sequential")
    registerDoParallel(cl)
    Tstat.WB = foreach(b=1:B, .combine='c', .packages=c('factorcpt', 'hdbinseg'), .export=c("LRV", "longrun_bartlett", "Demean", "GCT.tstat") ) %dopar% {

      one = matrix(1, bsize, hsize);
      two = matrix(rnorm(blength*hlength), blength, hlength)
      three = kronecker(two, one);
      four = three[1:n, 1:p];
      X.WB = X1*four;

      tstat.WB = GCT.tstat(X.WB, khat=khat, type=type);
    }
    stopCluster(cl)
    } else{
      Tstat.WB = numeric(B);
      for(b in 1:B){

        one = matrix(1, bsize, hsize);
        two = matrix(rnorm(blength*hlength), blength, hlength)
        three = kronecker(two, one);
        four = three[1:n, 1:p];
        X.WB = X1*four;

        Tstat.WB[b] = GCT.tstat(X.WB, khat=khat, type=type);
    }

    }
    out = list();
    out$tstat = tstat;
    names(out$tstat) = c("Max2", "Max", "Sum", "DC")[type];
    out$pval = mean(tstat < Tstat.WB);
    out$q = quantile(Tstat.WB, .95);
    return(out)
  }


  result = GCT.novote2(X=X, khat = khat0, B=B, bsize=bsize, hsize=hsize, do.parallel=do.parallel, n.cl=n.cl, type=type);
  out = list();
  out$khat = khat0;
  out$pval = result$pval;
  out$tstat = result$tstat;
  out$q = result$q;
  return(out)
}


GCT.binary <- function(X, type, do.parallel=TRUE, n.cl=7, maxK=3, spatial=TRUE, critTF=FALSE){

  if(missing(type)){ type=1; }

  n = length = nrow(X); p = ncol(X);

  out1 = GCT.break(X, wd=5, do.parallel=do.parallel, spatial = spatial, n.cl=n.cl, type=type);
  result = list(pvalhist = out1$pval, Brhist = out1$khat, tstathist = out1$tstat);
  br = c(0, n);
  nmax = 1;

  if(critTF){
    crit = out1$q;
    tstat = out1$tstat;
    khat = out1$khat;
    # Step1 : 1 break
    if (tstat > crit) {
      br <- c(0, khat, n)
      st <- sum(abs(tstat) > crit)

      # Step1 : 1 break
      spindex <- c(1, 1)

      while (st > 0 && nmax < maxK) {
        nmax <- nmax + 1
        lbr <- length(br)
        Ts <- NULL
        sst <- Br <- NULL
        brindex <- seq(from = 1, to = lbr - 1, by = 1)
        for (j in brindex) {
          if (spindex[j]) {
            id <- seq(from = br[j] + 1, to = br[j + 1], by = 1)
            dat = X[id,];
            if(nrow(dat) > 20 ){ # set the minimum segment length as 11
  #            out2 = GCT.break(dat, wd=5, do.parallel=do.parallel, n.cl=n.cl, type=type);
              khat2  = lse.break(dat, trim=5, type=type);
              tstat2 = GCT.tstat(dat, khat2, type=type);
              Br <- c(Br, br[j] + khat2)
              Ts <- c(Ts, tstat2)
              idf <- 1 * (tstat2 > crit)
              sst <- c(sst, rep(idf, idf + 1))
            } else {
              sst <- c(sst, 0)
            }
          } else {
            sst <- c(sst, 0)
          }
        }

        st <- sum(sst)
        if (!is.null(Ts)) {
          newbr <- (abs(Ts) > crit) * Br
          br <- unique(c(br, newbr))
        }
        br <- sort(br)
        spindex <- sst
        result$crit = crit;
        result$tstathist <- c(result$tstathist, Ts)
        result$Brhist <- c(result$Brhist, Br)
      } }
    } else {
  # Step1 : 1 break
  if(out1$pval < .05){
    khat = out1$khat;
    br = c(0, khat, n);
    st = sum( abs(out1$pval) < .05);

  # Step1 : 1 break
  spindex = c(1, 1);

  while(st > 0 && nmax < maxK){
    nmax = nmax+1;
    lbr = length(br); P = NULL; sst = Br = NULL;
    brindex = seq(from=1, to=lbr-1, by=1);
    for(j in brindex){
      if(spindex[j]){
            id = seq(from=br[j]+1, to= br[j+1], by=1);
            dat = X[id,];
            if(nrow(dat) > 20 ){ # set the minimum segment length as 11
              out2 = GCT.break(dat, wd=5, spatial = spatial, do.parallel=do.parallel, n.cl=n.cl, type=type);
              Br = c(Br, br[j]+out2$khat);
              P = c(P, out2$pval);
      	      idf = 1*(out2$pval < .05);
              sst = c(sst, rep(idf, idf+1));
      	  } else{ sst = c(sst, 0);}
      } else{ sst = c(sst, 0);}
    }

     st = sum(sst);
     if(!is.null(P)){ newbr = ( abs(P) < .05)*Br;
     br = unique(c(br, newbr)); }
     br = sort(br);
     spindex = sst;

    result$pvalhist = c(result$pvalhist, P);
    result$Brhist = c(result$Brhist, Br);
      }
    } }

    result$khat = br;
    result$iter = length(br)-1;
  return(result);
}


## Sliding version

## For MOSUM version
GCT.mosum = function(X, khat, lrv.ba){
  p = ncol(X);
  X1 = X[(1:khat),]; X2 =X[-(1:khat),];
  if(missing(lrv.ba)){
#    Xtilde = hdbinseg::clean.cp(t(X), type="sbs")$x;
#    lrv.ba = apply(Xtilde, 1, longrun_bartlett);
#    lrv.ba = sqrt(lrv.ba);
    lrv.ba = LRV(X1, X2, "bartlett");
    lrv.ba = sqrt(lrv.ba);
  }

  xtilde = t(X)/lrv.ba; # xtilde is dim*length
  tmp = hdbinseg:::func_dc(xtilde, .5);
  aa = tmp$acs;
  bb = aa[,khat]; bb2 = bb*bb;

  #  tmp1 = apply(xtilde, 1, cusum.stat, khat, adjTF=TRUE, isMax=FALSE) ## gives length*dim
  #  bb = tmp1[khat, ]; bb2 = bb*bb;

  Gn1 = max(bb2);
  Gn2 = max(bb);
  Gn3 = max(apply(aa^2, 2, sum)); ## Sum binary
  Gn4 = sum(bb2); ## Sum for mosum
  return(c(Gn1, Gn2, Gn3=0, Gn4))
}

GCT.slide3 = function(X, G, B, bsize, dwd, n.cl=50, spatial=TRUE, do.parallel=TRUE){

  length = n = nrow(X); p = ncol(X);
  if(missing(dwd)){ dwd = round(1.5*log(n)); }

  if(missing(G)){
#    Xtilde = clean.cp(t(X), type="sbs")$x;
    K.d = apply(X, 2, blocklength.andrew);
    G = floor(2*max(K.d));
    G = min(G, floor(.2*n));
    # wd = max(10, floor(log(n)));
  }
  if(missing(bsize)){ bsize = floor(sqrt(2*G));}
  time.seq = seq(from = G, to = length-G, by=1);

  if(missing(B)){ B = 21; }
#  bsize = floor(log(2*G));
  bsize = floor(sqrt(2*G));
  blength = floor((2*G)/bsize)+1;

  if(!spatial){hsize = p;  hlength = 2;
  } else{
    #  hsize = floor(log(p));
    L.d = apply(X, 1, blocklength.andrew);
    hsize = min(1.5*median(L.d), log(n*p)); hsize = floor(hsize);
    hsize =floor(sqrt(p));
    hlength = floor(p/hsize)+1;
  }

  if(do.parallel){
  library(doParallel)
  cl <- parallel::makeCluster(n.cl, setup_strategy = "sequential")
  registerDoParallel(cl)
  BQ = foreach(b=1:length(time.seq), .combine = rbind, .export=c("LRV", "longrun_bartlett", "Demean", "GCT.mosum"), .packages=c("factorcpt") ) %dopar% {
    i = time.seq[b];
    id1 = seq(from=i-G+1, to = i, by=1);
    id2 = seq(from=i+1, to = i+G, by=1);
    subY1 = X[id1,];
    subY2 = X[id2,];
    subY = rbind(subY1, subY2);
    m <- colMeans(subY)
    centeredY <- rbind(Demean(subY1, m), Demean(subY2, m))
    out = GCT.mosum(subY, G);
    r1 = matrix(0, 4, B);
    for( j in 1:B){
      one = matrix(1, bsize, hsize);
      two = matrix(rnorm(blength*hlength), blength, hlength)
      three = kronecker(two, one);
      four = three[1:(2*G), 1:p];
      X.WB <- centeredY * four;
      r1[,j] = GCT.mosum(X.WB, G);
    }

    cbind(out, r1);
  }

  stopCluster(cl)
  } else{
    BQ = matrix(0, 4, B+1);
    for (b in 1:length(time.seq)){
      i = time.seq[b];
      id1 = seq(from=i-G+1, to = i, by=1);
      id2 = seq(from=i+1, to = i+G, by=1);
      subY1 = X[id1,];
      subY2 = X[id2,];
      subY = rbind(subY1, subY2);
      m <- colMeans(subY)
      centeredY <- rbind(Demean(subY1, m), Demean(subY2, m))
      out = GCT.mosum(subY, G);
      r1 = matrix(0, 4, B);
      for( j in 1:B){
        one = matrix(1, bsize, hsize);
        two = matrix(rnorm(blength*hlength), blength, hlength)
        three = kronecker(two, one);
        four = three[1:(2*G), 1:p];
        X.WB <- centeredY * four;
        r1[,j] = GCT.mosum(X.WB, G);
      }
      BQ = rbind(BQ, cbind(out, r1));
    }
  }

  Ts = matrix(BQ[,1], nrow=4);
  tt = matrix(BQ[,-1], nrow=4)
  aq = apply(tt, 1, quantile, .95) # critical value
  br1 = locatecp.th(Ts[1,], aq[1], dwd);
  br2 = locatecp.th(Ts[2,], aq[2], dwd);
  br3 = locatecp.th(Ts[3,], aq[3], dwd);
  br4 = locatecp.th(Ts[4,], aq[4], dwd);

  ret = list();
  ret$Ts = Ts;
  ret$time = time.seq;
  ret$critical = aq;
  ret$G = G;
  ret$br.GCT = c(0, time.seq[br1], n);
  ret$br.max = c(0, time.seq[br2], n);
  ret$br.maxsum = c(0, time.seq[br3], n);
  ret$br.sum = c(0, time.seq[br4], n);
  return(ret)
}



##### Auxilary functions

lse.break = function(X, trim, type=1){

    n = nrow(X)
    if(missing(trim)){ trim= 0; }


    if(type < 4){
    cusum.stat = function(data, khat, adjTF=TRUE, isMax=FALSE){
      n = length(data)
      gmean = mean(data);
      x1 = cumsum(data); x1 = x1[-n];
      j=seq(from=1, to=(n-1), by=1);
      x2 = abs(x1 - j*gmean);
      if(adjTF){
        den = sqrt(j*(n-j)/n);
        lse = x2/den;
      } else{
        lse = x2/sqrt(n);  }

      if(isMax){
        lse = lse[khat]; }
      return(lse)
    }

    cc =apply(X, 2, cusum.stat, khat=10, adjTF=FALSE, isMax=FALSE)
    lse = rowSums(cc);
    trimid = seq(from=trim+1, to=(n-trim-1), by=1);
    lsetrim = lse[trimid];
    khat = trim + which.max(lsetrim);
    } else{
      DC.stat = function(X, sn){
        ## Scaled version
        if(missing(sn)){ sn = apply(X, 2, mad); }
        X1 = t(X)/sn; # X is a length*dim vector
        Y = hdbinseg:::func_dc(X1, .5)$res;
        trim = 5; n= dim(X)[1];
        trimid = seq(from=trim+1, to=(n-trim-1), by=1);
        Y = Y[trimid];
        out = list();
        out$stat = max(Y);
        out$khat = trim+which.max(Y);
        return(out)
      }
      khat = DC.stat(X)$khat;
    }
    return(khat)
}


longrun_bartlett <-function(data, q){
  n = length(data);
  if(missing(q)){
    rhohat =  ar.ols(data, FALSE, order=1);
    rhohat = rhohat$ar[1];
    q = 1.1447*(n*4*rhohat^2/(1-rhohat^2)^2)^(1/3);
  } # Data dependent bw

  q = min(floor(q+1), floor(n/2-2));
  wq = c(seq(from=1, to=q, by=1), q+1, seq(from=q, to=1, by=-1));
  xcov = acf(data, lag.max = q, type = "covariance", plot = FALSE, na.action = na.fail, demean = TRUE)
  sn2 = xcov$acf;
  sn2 = as.numeric(sn2);
  id2 = seq(from=q+1, to=2, by=-1)
  sn2 = c(sn2[id2], sn2)
  sn2 = sum(wq*sn2)/(q+1);
  return(sn2);
}

LRV = function(X1, X2, type){
  m1 <- colMeans(X1)
  m2 <- colMeans(X2);
  X <- rbind(X1, X2);

  devX1 <- Demean(X1, m1)
  devX2 <- Demean(X2, m2)

  if(type == "bartlett"){
    se = sqrt(apply(rbind(devX1, devX2), 2, longrun_bartlett));
  } else if (type == "sd") {
    #  se <- sqrt(colMeans(devX1^2) + colMeans(devX2^2))
    se = apply(rbind(devX1, devX2), 2, sd)
  } else {
    se = sqrt(apply(rbind(devX1, devX2), 2, longrun_bartlett, q=type));
  }
  return(se)
}


###################################################
## Bartlett long-run variance calculation
###################################################

Demean <- function(X, m)
{
  X - matrix(m, nrow = nrow(X), ncol = ncol(X), byrow = TRUE)
}

blocklength.andrew = function(data){
  n = length(data);
  rhohat =  ar.ols(data, FALSE, order=1);
  rhohat = rhohat$ar[1];
  q = 1.1447*(n*4*rhohat^2/(1-rhohat^2)^2)^(1/3);
  return(round(q))
}

# ###
# HH.stat = function(X1, X2, sn){
#
#   X = rbind(X1, X2);
#   if(missing(sn)){ sn = apply(X, 2, mad); }
#   X = t(X)/sn; # X1 is a dim*length vector
#
#   d = nrow(X); n = ncol(X); khat = nrow(X1);
#   Y = apply(X, 1, cusum.stat, adjTF=TRUE, isMax=FALSE); ## apply function gives length*dim as output
#   HH = rowSums(Y^2);
# #  HH = wt1*(rowSums(Y^2)-p)/p;
#   out = list();
#   out$HH = max(HH);
# #  sst = rowSums(Y^2);
# #  j = seq(from=1, to=n-1, by=1);
# #  wt1 =   j*(n-j)/n^2; wt1 = wt1/sqrt(d);
# #  HH = (sst-d)*wt1;
# #  EH = sst/sqrt(2*d);
# #  out = list();
# #  out$HH = HH;
# #  out$EH = EH;
#   return(out)
# }

#
#
# DC.stat = function(X, sn){
#   ## Scaled version
#   if(missing(sn)){ sn = apply(X, 2, mad); }
#   X1 = t(X)/sn; # X is a length*dim vector
#   Y = func_dc(X1)$res;
#   trim = 5; n= dim(X)[1];
#   trimid = seq(from=trim+1, to=(n-trim-1), by=1);
#   Y = Y[trimid];
#   out = list();
#   out$stat = max(Y);
#   out$khat = trim+which.max(Y);
#   return(out)
# }


locatecp.th = function(ts1, cq, dwd){
  if(missing(dwd)){ dwd = 3}

  if(sum(1*(ts1 > cq)) > .5*length(ts1)){ cq = quantile(ts1, .5) };

  id = 1*(ts1 >= cq); id = c(0, id);
  idplus = which( diff(id) == 1);
  idminus = which( diff(id) == -1);

  br = NULL;
  nid = min(length(idplus), length(idminus));
  if(nid >0){
    for(j in 1:nid){
      s = idplus[j]; e = idminus[j];
      if( e-s > dwd){
        subD = ts1[seq(from=s, to= e, by=1)];
        br = c(br, s+ which.max(subD) -1);
      }
    }
  } else{
    br= NULL;
  }
  return(br)
}


#  Restore the mean function from the breaks
## X is length*dim matrix so be careful here
## return mf-hat, which is dim*length
Mf.hat = function(X, br){
 nbr = sort(br);
  M=NULL;
  for(k in 1:(length(nbr)-1) ){
    ind = seq(from = nbr[k]+1, to = nbr[k+1], by=1)
    mfhat = replicate(nbr[k+1]-nbr[k], colMeans(X[ind,]));
    M = cbind(M, mfhat)
  }
  return(M)
}


###########
# Error model
############
GenError <- function(n, p, sigma, rho)
  # To generate AR(1) errors such that
  #           e_t = rho*e_{t-1} + a_t, a_t: iid N(0, sigma^2)
  #    n: sample size
  #    p: dimension
  #  rho: ar coefficient for the error process
{
  if ( rho == 0) {
    e <- matrix(rnorm(n*p, sd = sigma), n, p)
  } else {
    e <- matrix(NA, n, p)
    for ( j in 1:p ) {
      e[,j] <- arima.sim(n = n, list(ar = rho), sd = sigma)
    }
  }

  return(e)
}

##############
# ARMA(2,2)
##############
ARMA22 = function(n, p, rho){
  id = seq(from=0, to = 99, by=1)
  rhoi = rho/(id+1);
  sig1 = .2;
  burn = 50;

  E = matrix(rnorm((n+burn)*(p+n+burn), mean=0, sd=sig1), ncol=(burn+n))
  E1 = apply(E, 2, filter, filter=rhoi, method=c("convolution"), sides=1)
  E1 = E1[-(1:(burn+n)),];

  innov = arima.sim(n = p+burn , model=list(ma = rhoi), sd = sig1);
  innov = innov[-(1:burn)];
  X = matrix(0, nrow(E1), ncol(E1))
  for(j in 1:nrow(E1)){
    X[j, ] = arima.sim(n = n+burn, list(order=c(2,0,1), ar = c(.2, -0.3), ma = c(.2)), innov = E1[j,])
  }

  X = X[,-(1:burn)]

  return(X)
}


###
## Factor error
##

Factor.err = function(n, p, rho){
  id = seq(from=0, to = 99, by=1)
  rhoi = rho/(id+1);
  sig1 = .2*sqrt(1-rho^2);
  burn = 50;

  E = matrix(rnorm((n+burn+2)*(p+100+burn), mean=0, sd=sig1), ncol=(burn+n+2))
  E1 = apply(E, 2, filter, filter=rhoi, method=c("convolution"), sides=1)
  E1 = E1[-(1:(burn+100)),];

  h = rnorm((n+burn), mean=0, sd=.1);
  X = matrix(0, nrow(E1), ncol(E1))

  U = matrix(0, nrow(E1), ncol(E1));
  for(i in 1:p){
    nid = n+burn;
    for(j in 1:nid){
      U[i,j+2] = rho*h[j] + .2*U[i,j+1] - .3*U[i,j] + E1[i,j+1] + .2*E1[i, j];
    }
  }

  U = U[,-(1:(burn+2))];
  return(U)
}


####
# Generate Mean Function
####

GenMean = function(n, p, location, eta, delta){
  ## locations will be [0, k1, k2, n] etc..
  khat = c(0, location, n);
  nk = length(location);
  M = matrix(0, p, khat[2]-khat[1]);
  for(i in 1:nk){
    psize = floor(eta[i]*p);
    dim.jump = sample(1:p, psize, replace=FALSE);
    jump = runif(psize, min = delta[i]*.75, max = delta[i]*1.25);
    sjump = 2*rbinom(psize, 1, prob=1/2)-1;
    jump = jump*sjump
    mf = matrix(0, p, 1);
    mf[dim.jump, 1] = jump;
    M = cbind(M, mf %x% matrix(1, 1, khat[i+2] - khat[i+1]));
  }
  return(M)
}
