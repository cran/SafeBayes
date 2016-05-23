\name{rinvGauss}
\alias{rinvGauss}
\title{The inverse Gaussian and Wald distributions}
\description{
  Random generator for the inverse Gaussian and Wald distributions.
}
\usage{
	rinvGauss(n, nu, lambda)
}
\arguments{
	\item{n}{vector of numbers of observations}
	\item{nu}{vector real and non-negative parameter -- the Wald distribution results when nu=1}
	\item{lambda}{vector real and non-negative parameter}
	}
\value{
  The output values conform to the output from other such functions in R. \code{rinvGauss()} generates random numbers.
  }
\details{
This function is copied from the SuppDists package by Bob Wheeler. I have copied this function here, because the SuppDists pacakge is no longer maintained, so that I can maintain the rinvGauss function for use in the functions in this package.

Probability functions:

\deqn{f(x,\nu,\lambda)=\sqrt{\frac{\lambda}{2\pi x^3}}\exp\left[-\lambda\frac{(x-\nu)^2}{2\nu^2 x}\right] \mbox{-- the density}}{f(x,nu,lambda)=sqrt[lambda/(2 pi x^3)]exp[-lambda(x-nu)^2/(2 x nu^2)] -- the density
}

\deqn{F(x,\nu,\lambda)=\Phi\left[\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\nu}-1\right)\right]+e^{2\lambda/\nu}\Phi\left[\sqrt{\frac{\lambda}{x}}\left(\frac{x}{\nu}+1\right)\right] \mbox{-- the distribution function}}{F(x,nu,lambda)=Phi[sqrt(lambda/x)((x/nu)-1)]+exp(2*lambda/nu)*Phi[sqrt(lambda/x)((x/nu)+1)] -- the distribution function}

where \eqn{\Phi[]}{Phi[]} is the standard normal distribution function.


The calculations are accurate to at least seven significant figures over
an extended range - much larger than that of any existing tables. We
have tested them for \eqn{\lambda / \nu = 10e-20}{(lambda/nu)=10e-20},
and \eqn{\lambda/\nu=10e4}{lambda / nu = 10e4}. Accessible tables are
those of Wassan and Roy (1969), which unfortunately, are sometimes good
to only two significant digits. Much better tables are available in an
expensive CRC Handbook (1989), which are accurate to at least 7 significant digits for \eqn{\lambda/\nu \ge 0.02}{lambda / nu >= 0.02} to \eqn{\lambda/\nu \le 4000}{lambda / nu <= 4000}. 



These are first passage time distributions of Brownian motion with positive drift. See Whitmore and Seshadri (1987) for a heuristic derivation. The Wald (1947) form represents the average sample number in sequential analysis. The distribution has a non-monotonic failure rate, and is of considerable interest in lifetime studies: Ckhhikara and Folks (1977). A general reference is Seshadri (1993). 



This is an extremely difficult distribution to treat numerically, and it would not have been possible without some extraordinary contributions. An elegant derivation of the distribution function is to be found in Shuster (1968). The first such derivation seems to be that of Zigangirov (1962), which because of its inaccessibility, the author has not read. The method of generating random numbers is due to Michael, Schucany, and Haas (1976). The approximation of Whitmore and Yalovsky (1978) makes it possible to find starting values for inverting the distribution. All three papers are short, elegant, and non- trivial. 
}

\references{
Ckhhikara, R.S. and Folks, J.L. (1977) The inverse Gaussian distribution as a lifetime model. \emph{Technometrics.} \bold{19-4.} 461-468. 

CRC Handbook. (1989). \emph{Percentile points of the inverse Gaussian distribution.} J.A. Koziol (ed.) Boca Raton, FL. 

Michael, J.R., Schucany, W.R. and Haas, R.W. (1976). Generating random variates using transformations with multiple roots. \emph{American Statistician.} \bold{30-2.} 88-90. 

Seshadri, V. (1993). \emph{The inverse Gaussian distribution.} Clarendon, Oxford 

Shuster, J. (1968). On the inverse Gaussian distribution function. Jour. \emph{Am. Stat. Assoc.} \bold{63.} 1514-1516. 

Wasan, M.T. and Roy, L.K. (1969). Tables of inverse Gaussian percentage points. \emph{Technometrics.} \bold{11-3.} 591-604. 

Wald, A. (1947). \emph{Sequential analysis.} Wiley, NY 

Whitmore, G.A. and Seshadri, V. (1987). A heuristic derivation of the inverse Gaussian distribution. \emph{American Statistician.} \bold{41-4.} 280-281. 

Whitmore, G.A. and Yalovsky, M. (1978). A normalizing logarithmic transformation for inverse Gaussian random variables. \emph{Technometrics.} \bold{20-2.} 207-208. 

Zigangirov, K.S. (1962). Expression for the Wald distribution in terms of normal distribution. \emph{Radiotech.Electron.} \bold{7.} 164-166. 
}

\author{
	Bob Wheeler \email{bwheelerg@gmail.com}
}
\examples{
rinvGauss(1, 1, 16)
}
\keyword{distribution}

