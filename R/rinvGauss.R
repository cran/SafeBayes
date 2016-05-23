rinvGauss <-
function (n, nu, lambda) 
{
    n <- if (length(n) > 1) 
        length(n)
    else n
    N <- max(length(nu), length(lambda))
    nu <- rep(nu, length.out = N)
    lambda <- rep(lambda, length.out = N)
    .C("rinvGaussR", as.double(nu), as.double(lambda), as.integer(n), 
        as.integer(N), value = double(n),PACKAGE="SafeBayes")$value
}