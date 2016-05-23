SBRidgeRlog <-
function (y, X = NULL, etaseq = 1, prior = NULL, nIter = 1100, burnIn = 100, thin = 10, minAbsBeta = 1e-09, pIter=TRUE) {

    	y <- as.numeric(y) 
    	n <- length(y)
    	XL <- as.matrix(X)

   		MRlogError <- NULL
    	CMRlogEallen <- NULL
   
   		ytemp <- y
   		XLorigineel <- XL
    	XLtemp <- XL
    
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
        prior = list(varE = list(S = 0, df = 0), varBR = list(S = 0,df = 0))
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
    varE <- var(e, na.rm = TRUE)/2
    
    if (is.null(prior$varE)) {
        cat("No prior was provided for residual variance, an improper prior is used.\n")
        prior$varE <- list(df = 0, S = 0)       
    }

    post_mu <- 0
    post_varE <- 0
    post_logLik <- 0
    post_yHat <- rep(0, n)
    post_yHat2 <- rep(0, n)

   
    if (is.null(prior$varBR)) {
			cat("No prior was provided for varBR, improper prior is used.\n")
            prior$varBR <- list(df = 0, S = 0)
        }
        
        pR <- ncol(XL)
        xR2 <- colSums(XL * XL)
		bR <- rep(0, pR)
        namesBR <- colnames(XL)
        varBR <- varE/2/sum(xR2/n)
        
        post_bR <- rep(0, pR)
        post_bR2 <- post_bR
        post_varBR <- 0
        
        XLstacked <- as.vector(XL)
        rm(XL)
        
		time <- proc.time()[3]

# Sampling procedure

    bLsave <- NULL
    esave <- NULL
    varEsave <- NULL
    
    for (i in 1:nIter) {
    	pL <- pR
    	xL2 <- xR2
    	bL <- bR
    	varBj <- rep(varBR, pR)
    	
       	ans <- .Call("safe_sample_beta", n, pL, XLstacked, xL2, 
                bL, e, varBj, varE, minAbsBeta, eta)
            
       	bR <- bL <- ans[[1]]
       	e <- ans[[2]]
       	
       	bLsave <- rbind(bLsave, bL) 
        esave <- rbind(esave, e)
                        
        SS <- crossprod(bR) + prior$varBR$S
        df <- pR + prior$varBR$df
        varBR <- SS/rchisq(df = df, n = 1)
    	
    	e <- e + weights * mu
        rhs <- sum(weights * e)/varE
        C <- sumW2/varE
        sol <- rhs/C
        mu <- rnorm(n = 1, sd = sqrt(1/C)) + sol
        e <- e - weights * mu
        SS <- eta*crossprod(e) + prior$varE$S
        df <- eta*n + prior$varE$df
        
        varE <- as.numeric(SS)/rchisq(n = 1, df = df)
        sdE <- sqrt(varE)
        yHat <- yStar - e
        if (nNa > 0) {
            e[whichNa] <- rnorm(n = nNa, sd = sdE)
            yStar[whichNa] <- yHat[whichNa] + e[whichNa]
        }
        
        varEsave <- rbind(varEsave, varE)
                
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
                post_bR <- post_bR * k + bR/nSums
                post_bR2 <- post_bR2 * k + (bR^2)/nSums
                post_varBR <- post_varBR * k + varBR/nSums
             }
        }

        tmp <- proc.time()[3]
        
        if (pIter == TRUE) {
        cat(paste(c("Iter: ", "time/iter: ", "varE: "), 
            c(i, round(tmp - time, 3), round(varE, 3))))
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
        
        varERand <- varEsave[rind]

        RlogEtemp <- NULL
        for (k in 1:blrlang) {
        
        	bL2 <- bLRand[k,]
        	varE2 <- varERand[k]
        	# Let op, matrix wordt hier gedropt, wordt numeric
        	ypred <- NULL
			ynew <- NULL
			
			for (i in 1:length(XLtest)) {
		  		ynew[i] <- bL2[i]*(XLtest[i])
				}		
			ypred <- post_mu + sum(ynew)
			
        	RlogEtemp[k] <- (ytest - ypred)^2/(2*varE2) + log(2*pi*varE2)/2
        }
        
        MRlogError[specfold-1] <- mean(RlogEtemp)
        print("MRlogError")
        print(MRlogError)
        }
        
        CMRlogE <- sum(MRlogError)
        
        CMRlogEallen <- cbind(CMRlogEallen, CMRlogE)
      	print("allen")
      	print(CMRlogEallen)
        }
        
       	indeta <- which.min(CMRlogEallen)
       	print(indeta)
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
    
    out$bR <- as.vector(post_bR)
    tmp <- as.vector(sqrt(post_bR2 - (post_bR^2)))
    out$SD.bR <- tmp
    out$varBR <- post_varBR
    names(out$bR) <- namesBR
    names(out$SD.bR) <- namesBR
    
    out$prior <- prior
    out$nIter <- nIter
    out$burnIn <- burnIn
    out$thin <- thin
    
    out$CMRlogEallen <- CMRlogEallen
    out$eta.min <- eta.min
    
    return(out)
}
