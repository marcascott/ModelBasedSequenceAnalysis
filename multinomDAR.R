require(flexmix)

#Extension of FLXMRmultinom class to incorporate DAR(1) errors
#and implicitly allows for missing data as per Scott et al. (2020) article
#Additional modifications were made to allow for imbalance across dimensions, but these are still in progress.
#

setClass("FLXMRmultinomDAR",
         slots = c(b.na = "logical",ids="integer"),
         contains = "FLXMRmultinom")

expit <- function(x) exp(x)/(1+exp(x))

FLXMRmultinomDAR <- function(formula=.~., DAR=TRUE, ids=NULL, times=NULL, maxit=NULL, subset=NULL, ...)
{
  z <- new("FLXMRmultinomDAR", weighted=TRUE, formula=formula,
           family = "multinom", name=paste("FLXMRglm", "multinom with DAR effects", sep=":"))
  
  z@preproc.y <- function(x){
    #cat("in preproc.y\n")
    x <- as.matrix(x)
    # 01July19: modified to allow NAs in outcome (via na.omit & na.rm param)
    x <- as.integer(factor(x))  #will not use NAs as factor level by default.
    if (min(x,na.rm=T) < 1 | length(unique(na.omit(x))) != max(x,na.rm=T))
      stop("x needs to be coercible to an integer vector containing all numbers from 1 to max(x)")
    y <- matrix(0, nrow = length(x), ncol = max(x,na.rm=T))
    y[cbind(seq_along(x), x)] <- 1  #NAs will leave the row as zeros
    y
  }
  
  #z@preproc.x <- function(x){x}
    # under the current scheme there is nothing to do to x

  
  z@defineComponent <- expression({
    
    #this is typically "called" from the "with" function, and the context is the fit function, 
    # which has access to x, y, w, component, but not the 'whole' model object, which is why
    # we can't "see" some of the info in the slots we filled.  
    
    predict <- function(x) {
      #the zero rows of y (accessible in this internal method) mostly handle the NA issue
      p <- tcrossprod(x, coef)
      eta <- cbind(1, exp(p))
      pUnadj <- eta/rowSums(eta)
      if (!DAR) {
          return(pUnadj)  #currently, it seems NAs are removed by here.
      } else {
        rho <- expit(logitRho)  #DAR effect.  Need "lagged" y - use id (simple approach first)
        #correct for missing y here:
        #b.na isn't "available" from the model object calling this.  
        b.na <- apply(y,1,sum)==0 #these rows are NA on outcome
        lastRow <- sum(!b.na)
        lagged <- rbind(0,y[!b.na,,drop=F][-lastRow,,drop=F]) #drop last row
        #if we have 'time', we can compute the gap:
        tlag <- diff(c(0,times[!b.na]))  #amount of gap (wrong at start points)
        idx <- cumsum(c(1,rle(ids[!b.na])$lengths))  #to get positions of first in each id
        idx <- idx[-length(idx)] #don't need last posn
        #correct start points for tlag:
        tlag[idx] <- NA # will be overwritten anyway
        pAdj <- (rho^tlag)*lagged + ((1-(rho^tlag))*pUnadj[!b.na,,drop=F])
        pAdj[idx,] <- pUnadj[!b.na,,drop=F][idx,,drop=F] #these just get the marginals (so corrects the 'initial' NAs)
        #now to put the new values back in:
        pUnadj[!b.na,] <- pAdj
        #browser()
        return(pUnadj)
      }
    }
    logLik <- function(x, y) {
      #original
      #ll <- log(predict(x))[cbind(seq_len(nrow(y)), max.col(y, "first"))]
      #less efficient, probably, but handles the 0 rows elegantly.
      ll <- apply(log(predict(x))*y,1,sum)
      ll
    }
    new("FLXcomponent",
        parameters=list(coef=coef,logitRho=logitRho), logLik=logLik, predict=predict,
        df=df+1*(DAR))
  })
  
  z@fit <- function(x, y, w, component){
    wrapLL <- function(par,x,y,w) {
      #wrapper for optim call
      newCoef <- coef #borrow structure
      newCoef[] <- par[-1]
      para<-list(coef=newCoef,logitRho=par[1], df = length(coef))
      sum(w*with(para,eval(z@defineComponent))@logLik(x,y))
    }

    #test if first time through..  if so, initialize
    if (length(component)==0 || !DAR) {
      r <- ncol(x) ##subsetting not needed here
      p <- ncol(y)
      if (p < 2) stop("Multinom requires at least two components.")
      #run nnet if !DAR; use it as init value if DAR
      mask <- c(rep(0, r + 1), rep(c(0, rep(1, r)), p - 1))
      fit <- nnet::nnet.default(x, y, w, mask = mask, size = 0,
                          skip = TRUE, softmax = TRUE, censored = FALSE,
                          rang = 0, trace=FALSE, ...)
      fit$coefnames <- colnames(x) 
      fit$weights <- w
      fit$vcoefnames <- fit$coefnames[seq_len(ncol(x))] 
      fit$lab <- seq_len(ncol(y))
      class(fit) <- c("multinom", "nnet")
      coef <- coef(fit)
      logitRho <- NA
    } 
    if (length(component)!=0 && DAR) { #take prior fit
      coef <- component$coef
      logitRho <- component$logitRho
    } 
    #take fit from MN (or prior run) and initialize call to optim 
    if (DAR) {
      #cat("in DAR code...\n")
      if (is.null(maxit)) maxit <- 500
      if (is.na(logitRho)) logitRho <- -3.5 # good intial value near 0 prob.
      par <- c(logitRho,coef) #initial value should be small
      opt <- optim(par,wrapLL,x=x,y=y,w=w, control = list(maxit = maxit, fnscale = -1, reltol=5e-4)) #changed to subset x
      coef[] <- opt$par[-1] #maintains structure
      logitRho <- opt$par[1]
    }
    with(list(coef = coef, logitRho=logitRho, df = length(coef)+1*(DAR)),
         eval(z@defineComponent))
  }
  z
}

setMethod("existGradient", signature(object = "FLXMRmultinomDAR"),
          function(object) FALSE)

#this function is out of scope of the next one, so I'm copying it in here:
.FLXgetGroupingVar <- function(x)
{
  lf <- length(x)
  while (lf > 1) {
    x <- x[[lf]]
    lf <- length(x)
  }
  x
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRmultinomDAR"),
          function(model, data, formula, lhs=TRUE, ...)
          {
            model <- callNextMethod()
            #callNextMethod(model, data, formula, lhs=TRUE, ...)  #should call the flexmix built in version
            #cat("in FLXgetModelmatrix - DAR\n")
            #these two lines have little impact on the runs, but prove that one can store in the model object.
            model@b.na <- apply(model@y,1,sum)==0
            model@ids <- data[,as.character(.FLXgetGroupingVar(formula))]
            model
          })


          
