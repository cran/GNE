\name{SSR}
\alias{SSR}
\alias{funSSR}
\alias{jacSSR}
\title{SemiSmooth Reformulation}
\description{
		functions of the SemiSmooth Reformulation of the GNEP
	}
\usage{
funSSR(z, dimx, dimlam, grobj, arggrobj, constr, argconstr,  grconstr, arggrconstr, 
	compl, argcompl, dimmu, joint, argjoint, grjoint, arggrjoint, echo=FALSE)
jacSSR(z, dimx, dimlam, heobj, argheobj, constr, argconstr,  grconstr, arggrconstr, 
	heconstr, argheconstr, gcompla, gcomplb, argcompl, dimmu, joint, argjoint,
	grjoint, arggrjoint, hejoint, arghejoint, echo=FALSE)
}
\arguments{
  \item{z}{a numeric vector \code{z} containing \eqn{(x, lambda, mu)} values.}
  \item{dimx}{a vector of dimension for \eqn{x}.}
  \item{dimlam}{a vector of dimension for \eqn{lambda}.}
  \item{grobj}{gradient of the objective function, see details.}
  \item{arggrobj}{a list of additional arguments of the objective gradient.}
  \item{constr}{constraint function, see details.}
  \item{argconstr}{a list of additional arguments of the constraint function.}
  \item{grconstr}{gradient of the constraint function, see details.}
  \item{arggrconstr}{a list of additional arguments of the constraint gradient.}
  \item{compl}{the complementarity function with (at least) two arguments: \code{compl(a,b)}.}
  \item{argcompl}{list of possible additional arguments for \code{compl}.}
  
  \item{dimmu}{a vector of dimension for \eqn{mu}.}
  \item{joint}{joint function, see details.}
  \item{argjoint}{a list of additional arguments of the joint function.}
  \item{grjoint}{gradient of the joint function, see details.}
  \item{arggrjoint}{a list of additional arguments of the joint gradient.}

  \item{heobj}{Hessian of the objective function, see details.}
  \item{argheobj}{a list of additional arguments of the objective Hessian.}	
  \item{heconstr}{Hessian of the constraint function, see details.}
  \item{argheconstr}{a list of additional arguments of the constraint Hessian.}
  \item{gcompla}{derivative of the complementarity function w.r.t. the first argument.}
  \item{gcomplb}{derivative of the complementarity function w.r.t. the second argument.}
  
  \item{hejoint}{Hessian of the joint function, see details.}
  \item{arghejoint}{a list of additional arguments of the joint Hessian.}
  \item{echo}{a logical to show some traces.}    
}
\details{
Compute the SemiSmooth Reformulation of the GNEP: the Generalized Nash equilibrium problem is defined
by objective functions \eqn{Obj} with player variables \eqn{x} defined in \code{dimx} and 
may have player-dependent constraint functions \eqn{g} of dimension \code{dimlam} 
and/or a common shared joint function \eqn{h} of dimension \code{dimmu},
where the Lagrange multiplier are \eqn{lambda} and \eqn{mu}, respectively,
see F. Facchinei et al.(2009) where there is no joint function.
\describe{
    \item{\bold{Arguments of the Phi function}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{grobj}}{
			The gradient \eqn{Grad Obj} of an objective function \eqn{Obj} (to be minimized) must have 3 arguments for \eqn{Grad Obj(z, playnum, ideriv)}: vector \code{z}, player number, derivative index
			, and optionnally additional arguments in \code{arggrobj}.
			}
		\item{\code{constr}}{
			The constraint function \eqn{g} must have 2 arguments: vector \code{z}, player number,
			such that \eqn{g(z, playnum) <= 0}. Optionnally, \eqn{g} may have additional arguments in \code{argconstr}.
			}
		\item{\code{grconstr}}{
			The gradient of the constraint function \eqn{g}  must have 3 arguments: vector \code{z}, player number, derivative index,
			and optionnally additional arguments in \code{arggrconstr}.
			}
		\item{\code{compl}}{It must have two arguments and optionnally additional arguments in \code{argcompl}. 
			A typical example is the minimum function.}
		\item{\code{joint}}{
			The constraint function \eqn{h} must have 1 argument: vector \code{z},
			such that \eqn{h(z) <= 0}. Optionnally, \eqn{h} may have additional arguments in \code{argjoint}.
			}
		\item{\code{grjoint}}{
			The gradient of the constraint function \eqn{h}  must have 2 arguments: vector \code{z}, derivative index,
			and optionnally additional arguments in \code{arggrjoint}.
			}		
	}
	}
	\item{\bold{Arguments of the Jacobian of Phi}}{	
The arguments which are functions must respect the following features
	\describe{
		\item{\code{heobj}}{It must have 4 arguments: vector \code{z}, player number, two derivative indexes and optionnally additional arguments in \code{argheobj}.}
		\item{\code{heconstr}}{It must have 4 arguments: vector \code{z}, player number, two derivative indexes and optionnally additional arguments in \code{argheconstr}.}
		\item{\code{gcompla},\code{gcomplb}}{It must have two arguments and optionnally additional arguments in \code{argcompl}.}
		\item{\code{hejoint}}{It must have 3 arguments: vector \code{z}, two derivative indexes and optionnally additional arguments in \code{arghejoint}.}
	}
	}
}
See the example	below.
	
}
\value{
A vector for \code{funSSR} or a matrix for \code{jacSSR}.
}
\references{
 F. Facchinei, A. Fischer and V. Piccialli (2009), 
 \emph{Generalized Nash equilibrium problems and Newton methods},
 Math. Program.
} 
\seealso{
  See also \code{\link{GNE.nseq}}.
}
\author{
 Christophe Dutang
}
\examples{

# (1) associated objective functions
#

dimx <- c(2, 2, 3)

#Gr_x_j O_i(x)
grfullob <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- 3*(x - 1:7)^2
	}
	if(i == 2)
	{
		grad <- 1:7*(x - 1:7)^(0:6)
	}
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2 - 5	
		grad <- c(1, 0, 1, 0, 4*x[5]*s, 4*x[6]*s, 4*x[7]*s)
			
	}
	grad[j]	
}


#Gr_x_k Gr_x_j O_i(x)
hefullob <- function(x, i, j, k)
{
	x <- x[1:7]
	if(i == 1)
	{
		he <- diag( 6*(x - 1:7) )
	}
	if(i == 2)
	{
		he <- diag( c(0, 2, 6, 12, 20, 30, 42)*(x - 1:7)^c(0, 0:5) )
	}
	if(i == 3)
	{
		s <- x[5]^2 + x[6]^2 + x[7]^2	
		
		he <- rbind(rep(0, 7), rep(0, 7), rep(0, 7), rep(0, 7),
			c(0, 0, 0, 0, 4*s+8*x[5]^2, 8*x[5]*x[6], 8*x[5]*x[7]),
			c(0, 0, 0, 0, 8*x[5]*x[6], 4*s+8*x[6]^2, 8*x[6]*x[7]),
			c(0, 0, 0, 0,  8*x[5]*x[7], 8*x[6]*x[7], 4*s+8*x[7]^2) )
	}
	he[j,k]	
}



# (2) constraint linked functions
#

dimlam <- c(1, 2, 2)

#constraint function g_i(x)
g <- function(x, i)
{
	x <- x[1:7]
	if(i == 1)
		res <- sum( x^(1:7) ) -7
	if(i == 2)
		res <- c(sum(x) + prod(x) - 14, 20 - sum(x))
	if(i == 3)
		res <- c(sum(x^2) + 1, 100 - sum(x))
	res
}


#Gr_x_j g_i(x)
grfullg <- function(x, i, j)
{
	x <- x[1:7]	
	if(i == 1)
	{
		grad <- (1:7) * x ^ (0:6)
	}
	if(i == 2)
	{
		grad <- 1 + sapply(1:7, function(i) prod(x[-i]))
		grad <- cbind(grad, -1)
	}
	if(i == 3)
	{
		grad <- cbind(2*x, -1)
	}


	if(i == 1)
		res <- grad[j]	
	if(i != 1)
		res <- grad[j,]	
	as.numeric(res)
}



#Gr_x_k Gr_x_j g_i(x)
hefullg <- function(x, i, j, k)
{
	x <- x[1:7]
	if(i == 1)
	{
		he1 <- diag( c(0, 2, 6, 12, 20, 30, 42) * x ^ c(0, 0, 1:5) )
	}
	if(i == 2)
	{
		he1 <- matrix(0, 7, 7)
		he1[1, -1] <- sapply(2:7, function(i) prod(x[-c(1, i)]))
		he1[2, -2] <- sapply(c(1, 3:7), function(i) prod(x[-c(2, i)]))
		he1[3, -3] <- sapply(c(1:2, 4:7), function(i) prod(x[-c(3, i)]))
		he1[4, -4] <- sapply(c(1:3, 5:7), function(i) prod(x[-c(4, i)]))
		he1[5, -5] <- sapply(c(1:4, 6:7), function(i) prod(x[-c(5, i)]))
		he1[6, -6] <- sapply(c(1:5, 7:7), function(i) prod(x[-c(6, i)]))
		he1[7, -7] <- sapply(1:6, function(i) prod(x[-c(7, i)]))
						
						
		he2 <- matrix(0, 7, 7)
		
	}
	if(i == 3)
	{
		he1 <- diag(rep(2, 7))
		he2 <- matrix(0, 7, 7)
	}
	if(i != 1)
		return( c(he1[j, k], he2[j, k])	)
	else				
		return( he1[j, k] )
}


# (3) compute Phi
#

z <- rexp(sum(dimx) + sum(dimlam))

n <- sum(dimx)
m <- sum(dimlam)
x <- z[1:n]
lam <- z[(n+1):(n+m)]

resphi <- funSSR(z, dimx, dimlam, grobj=grfullob, constr=g, grconstr=grfullg, compl=phiFB)


check <- c(grfullob(x, 1, 1) + lam[1] * grfullg(x, 1, 1), 
	grfullob(x, 1, 2) + lam[1] * grfullg(x, 1, 2), 
	grfullob(x, 2, 3) + lam[2:3] \%*\% grfullg(x, 2, 3),
	grfullob(x, 2, 4) + lam[2:3] \%*\% grfullg(x, 2, 4), 	
	grfullob(x, 3, 5) + lam[4:5] \%*\% grfullg(x, 3, 5),
	grfullob(x, 3, 6) + lam[4:5] \%*\% grfullg(x, 3, 6),
	grfullob(x, 3, 7) + lam[4:5] \%*\% grfullg(x, 3, 7),
	phiFB( -g(x, 1), lam[1]), 
	phiFB( -g(x, 2)[1], lam[2]), 
	phiFB( -g(x, 2)[2], lam[3]), 
	phiFB( -g(x, 3)[1], lam[4]), 
	phiFB( -g(x, 3)[2], lam[5]))
	
	

#check
cat("\n\n________________________________________\n\n")

#part A
print(cbind(check, res=as.numeric(resphi))[1:n, ])
#part B
print(cbind(check, res=as.numeric(resphi))[(n+1):(n+m), ])
	
# (4) compute Jac Phi
#
	
resjacphi <- jacSSR(z, dimx, dimlam, heobj=hefullob, constr=g, grconstr=grfullg, 
	heconstr=hefullg, gcompla=GrAphiFB, gcomplb=GrBphiFB)

	
#check
cat("\n\n________________________________________\n\n")


cat("\n\npart A\n\n")	


checkA <- 
rbind(
c(hefullob(x, 1, 1, 1) + lam[1]*hefullg(x, 1, 1, 1), 
hefullob(x, 1, 1, 2) + lam[1]*hefullg(x, 1, 1, 2),
hefullob(x, 1, 1, 3) + lam[1]*hefullg(x, 1, 1, 3),
hefullob(x, 1, 1, 4) + lam[1]*hefullg(x, 1, 1, 4),
hefullob(x, 1, 1, 5) + lam[1]*hefullg(x, 1, 1, 5),
hefullob(x, 1, 1, 6) + lam[1]*hefullg(x, 1, 1, 6),
hefullob(x, 1, 1, 7) + lam[1]*hefullg(x, 1, 1, 7)
),
c(hefullob(x, 1, 2, 1) + lam[1]*hefullg(x, 1, 2, 1), 
hefullob(x, 1, 2, 2) + lam[1]*hefullg(x, 1, 2, 2),
hefullob(x, 1, 2, 3) + lam[1]*hefullg(x, 1, 2, 3),
hefullob(x, 1, 2, 4) + lam[1]*hefullg(x, 1, 2, 4),
hefullob(x, 1, 2, 5) + lam[1]*hefullg(x, 1, 2, 5),
hefullob(x, 1, 2, 6) + lam[1]*hefullg(x, 1, 2, 6),
hefullob(x, 1, 2, 7) + lam[1]*hefullg(x, 1, 2, 7)
),
c(hefullob(x, 2, 3, 1) + lam[2:3] \%*\% hefullg(x, 2, 3, 1), 
hefullob(x, 2, 3, 2) + lam[2:3] \%*\% hefullg(x, 2, 3, 2),
hefullob(x, 2, 3, 3) + lam[2:3] \%*\% hefullg(x, 2, 3, 3),
hefullob(x, 2, 3, 4) + lam[2:3] \%*\% hefullg(x, 2, 3, 4),
hefullob(x, 2, 3, 5) + lam[2:3] \%*\% hefullg(x, 2, 3, 5),
hefullob(x, 2, 3, 6) + lam[2:3] \%*\% hefullg(x, 2, 3, 6),
hefullob(x, 2, 3, 7) + lam[2:3] \%*\% hefullg(x, 2, 3, 7)
),
c(hefullob(x, 2, 4, 1) + lam[2:3] \%*\% hefullg(x, 2, 4, 1), 
hefullob(x, 2, 4, 2) + lam[2:3] \%*\% hefullg(x, 2, 4, 2), 
hefullob(x, 2, 4, 3) + lam[2:3] \%*\% hefullg(x, 2, 4, 3), 
hefullob(x, 2, 4, 4) + lam[2:3] \%*\% hefullg(x, 2, 4, 4), 
hefullob(x, 2, 4, 5) + lam[2:3] \%*\% hefullg(x, 2, 4, 5), 
hefullob(x, 2, 4, 6) + lam[2:3] \%*\% hefullg(x, 2, 4, 6), 
hefullob(x, 2, 4, 7) + lam[2:3] \%*\% hefullg(x, 2, 4, 7)
),
c(hefullob(x, 3, 5, 1) + lam[4:5] \%*\% hefullg(x, 3, 5, 1),  
hefullob(x, 3, 5, 2) + lam[4:5] \%*\% hefullg(x, 3, 5, 2),  
hefullob(x, 3, 5, 3) + lam[4:5] \%*\% hefullg(x, 3, 5, 3),  
hefullob(x, 3, 5, 4) + lam[4:5] \%*\% hefullg(x, 3, 5, 4),  
hefullob(x, 3, 5, 5) + lam[4:5] \%*\% hefullg(x, 3, 5, 5),  
hefullob(x, 3, 5, 6) + lam[4:5] \%*\% hefullg(x, 3, 5, 6),  
hefullob(x, 3, 5, 7) + lam[4:5] \%*\% hefullg(x, 3, 5, 7)
),
c(hefullob(x, 3, 6, 1) + lam[4:5] \%*\% hefullg(x, 3, 6, 1),   
hefullob(x, 3, 6, 2) + lam[4:5] \%*\% hefullg(x, 3, 6, 2),  
hefullob(x, 3, 6, 3) + lam[4:5] \%*\% hefullg(x, 3, 6, 3),  
hefullob(x, 3, 6, 4) + lam[4:5] \%*\% hefullg(x, 3, 6, 4),  
hefullob(x, 3, 6, 5) + lam[4:5] \%*\% hefullg(x, 3, 6, 5),  
hefullob(x, 3, 6, 6) + lam[4:5] \%*\% hefullg(x, 3, 6, 6),  
hefullob(x, 3, 6, 7) + lam[4:5] \%*\% hefullg(x, 3, 6, 7)
),
c(hefullob(x, 3, 7, 1) + lam[4:5] \%*\% hefullg(x, 3, 7, 1),   
hefullob(x, 3, 7, 2) + lam[4:5] \%*\% hefullg(x, 3, 7, 2),  
hefullob(x, 3, 7, 3) + lam[4:5] \%*\% hefullg(x, 3, 7, 3),  
hefullob(x, 3, 7, 4) + lam[4:5] \%*\% hefullg(x, 3, 7, 4),  
hefullob(x, 3, 7, 5) + lam[4:5] \%*\% hefullg(x, 3, 7, 5),  
hefullob(x, 3, 7, 6) + lam[4:5] \%*\% hefullg(x, 3, 7, 6),  
hefullob(x, 3, 7, 7) + lam[4:5] \%*\% hefullg(x, 3, 7, 7)
)
)


print(resjacphi[1:n, 1:n] - checkA)


cat("\n\n________________________________________\n\n")


cat("\n\npart B\n\n")	


checkB <- 
rbind(
cbind(c(grfullg(x, 1, 1), grfullg(x, 1, 2)), c(0, 0), c(0, 0), c(0, 0), c(0, 0)),
cbind(c(0, 0), rbind(grfullg(x, 2, 3), grfullg(x, 2, 4)), c(0, 0), c(0, 0)),
cbind(c(0, 0, 0), c(0, 0, 0), c(0, 0, 0), 
 rbind(grfullg(x, 3, 5), grfullg(x, 3, 6), grfullg(x, 3, 7)))
)


print(resjacphi[1:n, (n+1):(n+m)] - checkB)	


cat("\n\n________________________________________\n\n")
cat("\n\npart C\n\n")	


gx <- c(g(x,1), g(x,2), g(x,3))

checkC <- 
- t(
cbind(
rbind(
grfullg(x, 1, 1) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 2) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 3) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 4) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 5) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 6) * GrAphiFB(-gx, lam)[1],
grfullg(x, 1, 7) * GrAphiFB(-gx, lam)[1]
),
rbind(
grfullg(x, 2, 1) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 2) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 3) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 4) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 5) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 6) * GrAphiFB(-gx, lam)[2:3],
grfullg(x, 2, 7) * GrAphiFB(-gx, lam)[2:3]
),
rbind(
grfullg(x, 3, 1) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 2) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 3) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 4) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 5) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 6) * GrAphiFB(-gx, lam)[4:5],
grfullg(x, 3, 7) * GrAphiFB(-gx, lam)[4:5]
)
)
)



print(resjacphi[(n+1):(n+m), 1:n] - checkC)


cat("\n\n________________________________________\n\n")

cat("\n\npart D\n\n")	


checkD <- diag(GrBphiFB(-gx, lam)) 

print(resjacphi[(n+1):(n+m), (n+1):(n+m)] - checkD)

}
\keyword{math}
\keyword{optimize}

