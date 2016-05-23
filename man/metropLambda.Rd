\name{metropLambda}
\alias{metropLambda}
\title{
Metropolis-Hastings algorithm to sample lambda with a Beta prior for the Bayesian Lasso}
\description{
Metropolis-Hastings algorithm to sample lambda with a Beta prior from (de los Campos \emph{et al.}, 2009) for the Bayesian Lasso regression model.}
\usage{
metropLambda(tau2, lambda, shape1 = 1.2, shape2 = 1.2, max = 200, ncp = 0)
}
\arguments{
  \item{tau2}{
Latent parameter tau-squared to form the Laplace prior on the coefficients of the Lasso from a normal-mixture.
}
  \item{lambda}{
Initial value for lambda.}
  \item{shape1}{
First shape parameter for the Beta distribution.}
  \item{shape2}{
Second shape parameter for the Beta distribution.
}
  \item{max}{
Maximum value of lambda.}
  \item{ncp}{
Dummy parameter.}
}
\details{
Metropolis-Hastings algorithm to sample lambda with a Beta prior from (de los Campos \emph{et al.}, 2009) for the Bayesian Lasso regression model.}
\value{
Returns a value for lambda to use in the Gibbs samplers of the functions in the SafeBayes package.
}
\references{
de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.}
\author{
Copied from (de los Campos \emph{et al.}, 2009).}
\examples{
rm(list=ls())
library(SafeBayes)
tau2 <- 1/4
lambda <- 50

metropLambda(tau2=tau2, lambda=lambda)
}
\keyword{SafeBayes}