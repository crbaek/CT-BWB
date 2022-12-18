
rm(list=ls()) 
library(hdbinseg)
source("CT-BWB-library.R")

## This code illustrate CTBWB with one change-point
## E(2) error is used
## dimension =100, TS length=100
## Change-point is at 40 with jump size .4
## Sparse signal subject to change only 4 dimensions (out of 100)
## For E1 error in Ducker et al. (2023), take DGP=1
## For E3 error in Ducker et al. (2023), take DGP=3

DGP=2;
rho=.8; n=100; p=100;
location = c(floor(.4*n));
delta = .4 # jump size
eta = c(.04); # sparsity level
set.seed(12345)

if(DGP == 1){
      e <- GenError(n = n, p = p, rho = rho, sigma = sqrt(1 - rho^2)*.2);
      e = t(e);
} else if (DGP == 2){
      e <- ARMA21(n = n, p = p, rho = rho);
} else{
      e <- Factor.err(n = n, p = p, rho = rho);
}
    
mf = GenMean(n, p, location, eta, delta);
X = mf + e; ## Generated seris of dim*length
Y= t(X); ## Be careful here that Generated series should be length*dim for the below.

## If spatial=TRUE, it implements CTBWB
## BS-max2(type=1), BS-max(type=2), BS-sum(type=3), BS-DC(type=4)              
out1 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=1, spatial=TRUE);
out2 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=2, spatial=TRUE);
out3 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=3, spatial=TRUE);
out4 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=4, spatial=TRUE);

out1$br
out2$br
out3$br
out4$br

## CTBWB.mosum with spatial=TRUE implement MS-max2, MS-max, MS-sum
out5 = CTBWB.mosum(Y, spatial=TRUE, do.parallel = TRUE);

out5$ms.max2
out5$ms.max
out5$ms.sum


## If you do not want to consider spatial correlations in BWB, take spatial=FALSE.
## Then, it gives BWB of Jirak(2015) only considering temporal correlations
out11 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=1, spatial=FALSE);
out12 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=2, spatial=FALSE);
out13 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=3, spatial=FALSE);
out14 = CTBWB.binary(Y, do.parallel = TRUE, n.cl=10, type=4, spatial=FALSE);
out15 = CTBWB.mosum(Y, spatial=FALSE, n.cl=10, do.parallel =TRUE);

out11$br
out12$br
out13$br
out14$br
out15$ms.max2
out15$ms.max
out15$ms.sum

## DC of Cho (2016)
out21 = dcbs.alg(X, cp.type = 1, phi = 0.5, diag = FALSE, B = 500, q = 0.05, do.parallel = 0)
out21$ecp
