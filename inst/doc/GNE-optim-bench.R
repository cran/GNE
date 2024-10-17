## ----load, message=FALSE, echo=TRUE, eval=TRUE--------------------------------
library(GNE)

## ----argphi, fig=FALSE, echo=TRUE, eval=TRUE----------------------------------
myarg <- list(C=c(2, 3), D=c(4,0))
dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	dij <- 1*(i == j)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[i] - arg$C[i])*(x[other] - arg$D[i])^4*dij 
	res + 4*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^3*(1-dij) 
}

dimlam <- c(1, 1)
#g_i(x)
g <- function(x, i)
	ifelse(i == 1, sum(x[1:2]) - 1, 2*x[1]+x[2]-2)
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	ifelse(i == 1, 1, 1 + 1*(i == j))

## ----argJacF, fig=FALSE, echo=TRUE, eval=TRUE---------------------------------
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
	dij <- 1*(i == j)
	dik <- 1*(i == k)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[other] - arg$D[i])^4*dij*dik 
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*dij*(1-dik)
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*(1-dij)*dik
	res + 12*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^2*(1-dij)*(1-dik)
}
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k) 0

## ----testGNE, fig=FALSE, echo=TRUE, eval=TRUE---------------------------------
set.seed(1234)
z0 <- rexp(sum(dimx)+sum(dimlam))
GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=0))

## ----fig=FALSE, echo=FALSE, eval=TRUE-----------------------------------------
#list of true GNEs
trueGNE <- rbind(c(2, -2, 0, 5*2^5),
	c(-2, 3, 8, 0),
	c(0, 1, 4*3^4, 0),
	c(1, 0, 2^9, 6))
colnames(trueGNE) <- c("x1", "x2", "lam1", "lam2")
rownames(trueGNE) <- 1:4

## ----fig=FALSE, echo=FALSE, eval=TRUE-----------------------------------------
print(trueGNE)

## ----bench, fig=FALSE, echo=TRUE, eval=FALSE----------------------------------
#  wholebench <- function(z0)
#  {
#    #min function
#    resMin <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMin),
#                             argjac=list(gphia= GrAphiMin, gphib= GrBphiMin), echo=FALSE)
#  
#    #FB function
#    resFB <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiFB),
#                            argjac=list(gphia= GrAphiFB, gphib= GrBphiFB), echo=FALSE)
#  
#    #Mangasarian function
#    resMan <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMan, f=function(t) t^3),
#                          argjac=list(gphia= GrAphiMan, gphib= GrBphiMan, fprime=function(t) 3*t^2),
#                          echo=FALSE, control=list(maxit=200))
#  
#    #LT function
#    resLT <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiLT, q=4),
#                            argjac=list(gphia= GrAphiLT, gphib= GrBphiLT, q=4))
#  
#    #KK function
#    resKK <- bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiKK, lambda=3/2),
#                            argjac=list(gphia= GrAphiKK, gphib= GrBphiKK, lambda=3/2))
#  
#    list(resMin=resMin, resFB=resFB, resMan=resMan, resLT=resLT, resKK=resKK)
#  }

## ----benchcall, fig=FALSE, echo=TRUE, eval=FALSE------------------------------
#  initialpt <- cbind(c(4, -4), c(-4, 4), c(3, 0), c(0, 3), c(-1, -1), c(0, 0))
#  mytablelist <- list()
#  for(i in 1: NCOL(initialpt))
#  {
#  	z0 <- c(initialpt[, i], 1, 1)
#  	mybench <- wholebench(z0)
#  
#  	cat("z0", z0, "\n")	
#  
#  	mytable12 <- data.frame(method=mybench[[1]]$compres[, 1],
#  	round(
#  		cbind(mybench[[1]]$compres[,c(-1, -4)], mybench[[2]]$compres[,c(-1, -4)])
#  		, 3) )
#  
#  	mytable35 <- data.frame(method=mybench[[1]]$compres[, 1],
#  	round(
#  		cbind(mybench[[3]]$compres[,c(-1, -4)], mybench[[5]]$compres[,c(-1, -4)])
#  		, 3) )
#  
#  	mytablelist <- c(mytablelist, z0=list(z0), MINFB=list(mytable12), MANKK=list(mytable35))
#  }

## ----benchessai, fig=FALSE, echo=TRUE, eval=FALSE-----------------------------
#  z0 <- c(-4, 4, 1, 1)
#  bench.GNE.nseq(z0, F, JacF, argPhi=list(phi=phiMin),
#                 argjac=list(gphia= GrAphiMin, gphib= GrBphiMin), echo=FALSE)$compres

## ----singjac, fig=FALSE, echo=TRUE, eval=TRUE---------------------------------
z0 <- c(0, 0, 1, 1)
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiMin, gcomplb=GrBphiMin)

## ----singjac2, fig=FALSE, echo=TRUE, eval=TRUE--------------------------------
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiKK, gcomplb=GrBphiKK, argcompl=3/2)

## ----testGNEceq, fig=FALSE, echo=TRUE, eval=TRUE------------------------------
z0 <- 1+rexp(sum(dimx)+2*sum(dimlam))
GNE.ceq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	method="PR", control=list(trace=0))

## ----testNI, fig=FALSE, echo=TRUE, eval=TRUE----------------------------------
#O_i(x)
obj <- function(x, i, arg)
  (x[i] - arg$C[i])^2*(x[-i] - arg$D[i])^4
#g(x)
gtot <- function(x)
  sum(x[1:2]) - 1
#Gr_x_j g(x)
jacgtot <- function(x)
	cbind(1, 1)

z0 <- rexp(sum(dimx))

GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, 
         jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), 
         problem="NIR", merit="NI")


GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, 
         jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), 
         problem="VIR", merit="VI")

