SBLassoRSq <-
function (y, X = NULL, sigma2=NULL, etaseq = 1, prior = NULL, nIter = 1100, burnIn = 100, 		thin = 10, minAbsBeta = 1e-09, pIter = TRUE) {

    	y <- as.numeric(y) 
    	n <- length(y)
    	XL <- as.matrix(X)
    	
   	MRSError <- NULL
    CMRSEallen <- NULL
   
   	ytemp <- y
   	XLorigineel <- XL
    XLtemp <- as.matrix(XL)
    
    for (eta in etaseq) {
    	
    for (specfold in 2:(n-1)) { 
    
    		ytest <- ytemp[specfold+1]
    		y <- ytemp[1:specfold]
    		XLtrain <- XLtemp[1:specfold,]
    		XLtest <- XLtemp[specfold+1,]
    		XL <- as.matrix(XLtrain)
    		
	n <- length(y)
    	
    	eta <- as.numeric(eta)
    	if (any(eta > 1) | any(eta < 0)) {
    		stop("One or more of the eta's provided is larger than 1 or smaller than 0 ")
    	}

    if (!is.null(XL)) {
        if (any(is.na(XL)) | nrow(XL) != n) 
            stop("The number of rows in X does not correspond with that of y or it contains missing values")
    }

    if (is.null(prior)) {
        cat("No prior was provided, improper priors are used.\n")
        prior = list(lambda = list(shape = 0, rate = 0, type = "random", value = 50))
    }
    
    nSums <- 0
    whichNa <- which(is.na(y))
    nNa <- sum(is.na(y))
    weights <- rep(1, n)
    sumW2 <- sum(weights^2)
    mu <- weighted.mean(x = y, w = weights, na.rm = TRUE)
    yStar <- y * weights
    yStar[whichNa] <- mu * weights[whichNa]
    e <- (yStar - weights * mu)
   
    if (is.null(sigma2)) {
    		varE <- var(e, na.rm= TRUE)/2 
    		cat("Sigma^2 is not provided and will be estimated from the data")
    		}
    else {
    		varE <- sigma2
    		}
    		
    post_mu <- 0
    post_varE <- 0
    post_logLik <- 0
    post_yHat <- rep(0, n)
    post_yHat2 <- rep(0, n)
    
    if (is.null(prior$lambda)) {
            cat("No prior was provided for lambda, improper prior is used.\n")
            prior$lambda <- list(shape = 0, rate = 0, value = 50, 
                type = "random")
        }
        
        pL <- ncol(XL)
        xL2 <- colSums(XL*XL)
        bL <- rep(0, pL)
        namesBL <- colnames(XL)
		tmp <- 1/2/sum(xL2/n)
        tau2 <- rep(tmp, pL)
        lambda <- prior$lambda$value
        lambda2 <- lambda^2
        post_lambda <- 0
        post_bL <- rep(0, pL)
        post_bL2 <- post_bL
        post_tau2 <- rep(0, pL)
        
        XLstacked <- as.vector(XL)
        rm(XL)
    
	time <- proc.time()[3]

# Sampling procedure 

    bLsave <- NULL
    esave <- NULL

    for (i in 1:nIter) {
            varBj <- tau2 * varE
            ans <- .Call("safe_sample_beta", n, pL, XLstacked, xL2, 
                bL, e, varBj, varE, minAbsBeta, eta)
           
            bL <- ans[[1]]
            e <- ans[[2]]
            
            bLsave <- rbind(bLsave, bL) 
            esave <- rbind(esave, e)
            
            nu <- sqrt(varE) * lambda/abs(bL)
            tmp <- NULL
            
            try(tmp <- rinvGauss(n = pL, nu = nu, lambda = lambda2))
            if(!is.null(tmp))
            {
               if(!any(is.na(sqrt(tmp))))
               { 
                  tau2 <- 1/tmp
               }else{
                  cat("WARNING: tau2 was not updated due to numeric problems with beta\n");
               }
            }else{
               cat("WARNING: tau2 was not updated due to numeric problems with beta\n");
            }
            
            if (prior$lambda$type == "random") {
                if (is.null(prior$lambda$rate)) {
                  lambda <- metropLambda(tau2 = tau2, lambda = lambda, 
                    shape1 = prior$lambda$shape1, shape2 = prior$lambda$shape2, 
                    max = prior$lambda$max)
                  lambda2 <- lambda^2
                }
                
                else {
                  rate <- sum(tau2)/2 + prior$lambda$rate
                  shape <- pL + prior$lambda$shape
                  lambda2 <- rgamma(rate = rate, shape = shape,n = 1)
                  if(!is.na(lambda2))
                  {
                       lambda <- sqrt(lambda2)
                  }else{
                     cat("WARNING: lambda was not updated due to numeric problems with beta\n");   
                  }
                }
            }
            
        e <- e + weights * mu
        rhs <- sum(weights * e)/varE
        C <- sumW2/varE
        sol <- rhs/C
        mu <- rnorm(n = 1, sd = sqrt(1/C)) + sol
        e <- e - weights * mu
		
		sdE <- sqrt(varE)
        yHat <- yStar - e
        if (nNa > 0) {
            e[whichNa] <- rnorm(n = nNa, sd = sdE)
            yStar[whichNa] <- yHat[whichNa] + e[whichNa]
        }
        
        if ((i%%thin == 0)) {
            if (i >= burnIn) {
                nSums <- nSums + 1
                k <- (nSums - 1)/(nSums)
                tmpE <- e/weights
                tmpSD <- sqrt(varE)/weights
                if (nNa > 0) {
                  tmpE <- tmpE[-whichNa]
                  tmpSD <- tmpSD[-whichNa]
                }
                
                logLik <- sum(dnorm(tmpE, sd = tmpSD, log = TRUE))
                post_logLik <- post_logLik * k + logLik/nSums
                post_mu <- post_mu * k + mu/nSums
                post_varE <- post_varE * k + varE/nSums
                post_yHat <- post_yHat * k + yHat/nSums
                post_yHat2 <- post_yHat2 * k + (yHat^2)/nSums
                post_lambda <- post_lambda * k + lambda/nSums
                post_bL <- post_bL * k + bL/nSums
                post_bL2 <- post_bL2 * k + (bL^2)/nSums
                post_tau2 <- post_tau2 * k + tau2/nSums
            }
        }
        
        tmp <- proc.time()[3]
        
        if (pIter == TRUE) {
         	cat(paste(c("Iter: ", "time/iter: ", "varE: ", "lambda: "), 
             	c(i, round(tmp - time, 3), round(varE, 3), round(lambda, 
                 	3))))
         	cat("\n")
         	cat(paste("------------------------------------------------------------"))
         	cat("\n")
        }
        time <- tmp
    }
    
        tempbetas <- as.matrix(bLsave)
        rind <- sample((burnIn+1):nIter, 100, replace=TRUE)

        bLRand <- as.matrix(tempbetas[rind,])
        blrlang <- dim(bLRand)[1]
        
        RSEtemp <- NULL
        for (k in 1:blrlang) {
        
        		bL2 <- bLRand[k,]
        		ypred <- NULL
			ynew <- NULL
			
			for (i in 1:length(XLtest)) {
		  		ynew[i] <- bL2[i]*(XLtest[i])
				}		
			ypred <- post_mu + sum(ynew)
			
        		RSEtemp[k] <- (ytest - ypred)^2
        }
        
        MRSError[specfold-1] <- mean(RSEtemp)
        #print("MRSError")
        #print(MRSError)
        }
        
        CMRSE <- sum(MRSError)
        CMRSEallen <- cbind(CMRSEallen, CMRSE)
        }
        
       	indeta <- which.min(CMRSEallen)
       	#print(indeta)
       	eta.min <- etaseq[indeta]

    tmp <- sqrt(post_yHat2 - (post_yHat^2))
    out <- list(y = y, weights = weights, mu = post_mu, varE = post_varE, 
        yHat = I(post_yHat/weights), SD.yHat = I(tmp/weights), 
        whichNa = whichNa)
    names(out$yHat) <- names(y)
    names(out$SD.yHat) <- names(y)
    tmpE <- (yStar - post_yHat)/weights
    tmpSD <- sqrt(post_varE)/weights
    if (nNa > 0) {
        tmpE <- tmpE[-whichNa]
        tmpSD <- tmpSD[-whichNa]
    }
    
    out$fit <- list()
    logLikAtPostMean <- sum(dnorm(tmpE, sd = tmpSD, log = TRUE))
    out$fit$pD <- -2 * (post_logLik - logLikAtPostMean)
    out$fit$DIC <- out$fit$pD - 2 * post_logLik
    
    out$lambda <- post_lambda
    out$bL <- as.vector(post_bL)
    tmp <- as.vector(sqrt(post_bL2 - (post_bL^2)))
    out$SD.bL <- tmp
    out$tau2 <- post_tau2
    names(out$bL) <- namesBL
    names(out$SD.bL) <- namesBL
    
    out$prior <- prior
    out$nIter <- nIter
    out$burnIn <- burnIn
    out$thin <- thin
    
    out$CMRSEallen <- CMRSEallen
    out$eta.min <- eta.min
    
    return(out)
}
