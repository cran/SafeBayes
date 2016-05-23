\name{SBRidgeRSq}
\alias{SBRidgeRSq}
\title{
R-square Safe-Bayesian Ridge Regression}
\description{
The function SBRidgeRSq (R-square-Safe-Bayesian Ridge Regression) provides a Gibbs sampler together with the R-square-Safe-Bayesian algorithm for Ridge regression models with fixed variance.}
\usage{
SBRidgeRSq(y, X = NULL, sigma2 = NULL, etaseq = 1, prior = NULL, nIter = 1100, 
	burnIn = 100, thin = 10, minAbsBeta = 1e-09, pIter = TRUE)
}
\arguments{
  \item{y}{
Vector of outcome variables, numeric, NA allowed, length \eqn{n}{n}.
}
  \item{X}{
Design matrix, numeric, dimension \eqn{n \times p}{n x p}, \eqn{n \geq 2}{n >= 2}. 
}
  \item{sigma2}{
Fixed variance parameter \eqn{\sigma^2}{sigma^2}, numeric. Default NULL, in which case the variance will be estimated from the data \emph{per addition of new data point} in the Safe-Bayesian algorithm.
}
  \item{etaseq}{
Vector of learning rates \eqn{\eta}{eta}, numeric, \eqn{0 \leq \eta \leq 1}{0 <= eta <= 1}. Default 1.
}
  \item{prior}{
List containing the following elements
	\itemize{
	\item prior$varE: prior for the variance parameter \eqn{\sigma^2}{\sigma^2} with parameters $df and $S for respectively degrees of freedom and scale parameters for an inverse-chi-square distribution. Default \eqn{(0,0)}.
	\item prior$varBR: prior for the variance of the Gaussian prior for the coefficients \eqn{\beta}{beta}, with parameters $df and $S for respectively degrees of freedom and scale parameters for an inverse-chi-square distribution. Default \eqn{(0,0)}.}}
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
  \item{pIter}{
Print iterations, logical. Default TRUE.
}
}
\details{
Details on generalized Bayesian regression can be found in (de Heide, 2016). The implementation of the Gibbs sampler is based on the BLR package of (de los Campos \emph{et al.}, 2009).

The Safe-Bayesian algorithm was proposed by Gruenwald (2012) as a method to learn the learning rate for the generalized posterior to deal with model misspecification.}
\value{
\item{$y}{
	Vector of original outcome variables.}
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
\item{$bR}{
	Posterior mean of \eqn{\beta}{beta}.}
\item{$SD.bR}{
	Corresponding standard deviation.}
\item{$prior}{
	List containing the priors used.}
\item{$nIter}{
	Number of iterations.}
\item{$burnIn}{
	Number of iterations for burn-in.}
\item{$thin}{
	Number of iterations for thinning.}
\item{$CMRSEallen}{
	List of cumulative posterior-expected posterior-randomized square-loss per \eqn{\eta}{eta}.}
\item{$eta.min}{
	Learning rate \eqn{\eta}{eta} minimizing the cumulative posterior-expected posterior-randomized square-loss.}
}
\references{
de Heide, R. 2016.
The Safe-Bayesian Lasso. \emph{Master Thesis, Leiden University}.

de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

Gruenwald, P.D. 2012.
chapter The Safe Bayesian. \emph{Algorithmic Learning Theory: 23rd International Conference, ALT 2012, Lyon, France, October 29-31, 2012. Proceedings}. 169-183. Springer Berlin Heidelberg}
\author{
R. de Heide}
\examples{
rm(list=ls())
# Simulate data
x <- runif(10, -1, 1) # 10 random uniform x's between -1 and 1
y <- NULL

# for each x, an y that is 0 + Gaussian noise
for (i in 1:10) {
  y[i] <-  0 + rnorm(1, mean=0, sd=1/4) 
  }

plot(x,y)

# Let R-square-SafeBayes learn the learning rate 
sbobj <- SBRidgeRSq(y, x, etaseq=c(1, 0.5, 0.25))

# eta 
sbobj$eta.min
}
\keyword{SafeBayes}