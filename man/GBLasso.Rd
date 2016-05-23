\name{GBLasso}
\alias{GBLasso}
\title{Generalized Bayesian Lasso
}
\description{The function GBLasso (Generalized Bayesian Lasso) provides a Gibbs sampler to sample the posterior of generalized Bayesian lasso regression models with learning rate \eqn{\eta}{eta}.
}
\usage{
GBLasso(y, X = NULL, eta = 1, prior = NULL, nIter = 1100, burnIn = 100, 
	thin = 10, minAbsBeta = 1e-09, weights = NULL, pIter = TRUE)
}
\arguments{
  \item{y}{
Vector of outcome variables, numeric, NA allowed, length \eqn{n}{n}.
}
  \item{X}{
Design matrix, numeric, dimension \eqn{n \times p}{n x p}, \eqn{n \geq 2}{n >= 2}.

}
  \item{eta}{
Learning rate \eqn{\eta}{eta}, numeric, \eqn{0 \leq \eta \leq 1}{0 <= eta <= 1}. Default 1.

}
  \item{prior}{
List containing the following elements
	\itemize{
	\item prior$varE: prior for the variance parameter \eqn{\sigma^2}{\sigma^2} with parameters $df and $S for respectively degrees of freedom and scale parameters for an inverse-chi-square distribution. Default \eqn{(0,0)}.
	\item prior$lambda: prior for the penalty parameter \eqn{\lambda}{lambda} with three items. 
	\itemize{
	\item $value Initial value for \eqn{\lambda}{lambda}. Default \eqn{50}{50}
	\item $type Can be `fixed': initial value is used as fixed penalty parameter or `random', in which case a prior for \eqn{\lambda}{lambda} is specified. Default `random'.
	\item For a Gamma prior on \eqn{\lambda^2}{lambda^2}: $shape for shape parameter and $rate for the rate parameter; for a Beta prior on \eqn{\lambda}{lambda}: $shape1, $shape2 and $max for \eqn{\lambda}{lambda} proportional to Beta(\eqn{\lambda}{lambda}/max, shape1, shape2). Default: Gamma\eqn{(0,0)}.
}}}
  \item{nIter}{
Number of iterations, integer. Default 1100.

}
  \item{burnIn}{
Number of iterations for burn-in, integer. Default 100.

}
  \item{thin}{
Number of iterations for thinning, integer. Default 10.

}
  \item{minAbsBeta}{
Minimum absolute value of sampled coefficients beta to avoid numerical problems, numeric. Default \eqn{10^{-9}}{10^-9}.

}
  \item{weights}{
Vector of weights, numeric, length \eqn{n}{n}. Default NULL, in which case all weights are set to 1.
}
  \item{pIter}{
Print iterations, logical. Default TRUE.
}
}
\details{
Details on the generalized Bayesian lasso can be found in Chapter 2 of (de Heide, 2016). The implementation is heavily based on the BLR package of (de los Campos \emph{et al.}, 2009).



Several authors have brought forward the idea of equipping Bayesian updating with a learning rate \eqn{\eta}{eta}, resulting in an \eqn{\eta}{eta}-generalized posterior (Vovk (1990), McAllester (2003), Seeger (2002), Catoni (2007), Audibert (2004), Zhang (2004)). Gruenwald (2012) suggested its use as a method to deal with model misspecification. In the \eqn{\eta}{eta}-generalized posterior, the likelihood is raised to the power \eqn{\eta}{eta} in order to trade the relative weight of the likelihood and the prior, where \eqn{\eta = 1}{eta = 1} corresponds to standard Bayes.
}
\value{
\item{$y}{
	Vector of original outcome variables.}
\item{$weights}{
	Vector of weights.}
\item{$mu}{
	Posterior mean of the intercept.}
\item{$varE}{
	Posterior mean of of the variance.}
\item{$yHat}{
	Posterior mean of mu + X*beta + epsilon.}
\item{$SD.yHat}{
	Corresponding standard deviation.}
\item{$whichNa}{
	Vector with indices of missing values of y.}
\item{$fit$pD}{
	Estimated number of effective parameters.}
\item{$fit$DIC}{
	Deviance Information Criterion.}
\item{$lambda}{
	Posterior mean of \eqn{\lambda}{lambda}.}
\item{$bL}{
	Posterior mean of \eqn{\beta}{beta}.}
\item{$SD.bL}{
	Corresponding standard deviation.}
\item{$tau2}{
	Posterior mean of \eqn{\tau^2}{tau^2}.}
\item{$prior}{
	List containing the priors used.}
\item{$nIter}{
	Number of iterations.}
\item{$burnIn}{
	Number of iterations for burn-in.}
\item{$thin}{
	Number of iterations for thinning.}
\item{$eta}{ 
	Learning rate \eqn{\eta}{eta}.}
}
\references{
de Heide, R. 2016.
The Safe-Bayesian Lasso. \emph{Master Thesis, Leiden University}.

de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

Audibert, J.Y. 2013.
Bayesian generalized double pareto shrinkage. \emph{Statistica Sinica} \bold{23} 119-143.

Catoni, O. 2007.
PAC-Bayesian Supervised Classification. Lecture Notes - Monograph Series. IMS.

Gruenwald, P.D. 2012.
chapter The Safe Bayesian. \emph{Algorithmic Learning Theory: 23rd International Conference, ALT 2012, Lyon, France, October 29-31, 2012. Proceedings}. 169-183. Springer Berlin Heidelberg

McAllester, D. 2003.
PAC-Bayesian stochastic model selection. \emph{Machine Learning} \bold{51(1)} 5-21

Vovk, V.G. 1990.
Aggregating strategies. In \emph{Proc. COLT 90} 371-383

Zhang, T. 2004.
Learning bounds for a generalized family of Bayesian posterior distributions. \emph{Advances of Neural Information Processing Systems 16}, Thrun, L.S. and Schoelkopf, B. eds., MIT Press 1149-1156

}
\author{
R. de Heide
}

\examples{
rm(list=ls())
library(SafeBayes)

# Simulate data
x <- runif(100, -1, 1) # 100 random uniform x's between -1 and 1
y <- NULL

# for each x, a y that is 0 + Gaussian noise
for (i in 1:100) {
  y[i] <-  0 + rnorm(1, mean=0, sd=1/4) 
  }
  
# Now sample 100 zero's and ones (coin toss)
cointoss <- sample(0:1, 100, replace=TRUE)
# indices of the ones
indices <- which(cointoss==1)

# we replace x and y with (0,0) for the indices the cointoss 
# landed tail (1)
x[indices] <- 0
y[indices] <- 0

plot(x,y)

# Determine the generalized posterior for eta = 0.25
obj <- GBLasso(y, x, eta=0.25)

# posterior means of the coefficients beta and intercept mu
betafour <- obj$bL
mufour <- obj$mu
}
\keyword{SafeBayes}
