#' The contingency periodogram for periodicity in categorical time series
#' 
#' A function to estimate the contingency periodogram to test for periodicity in
#' categorical time series.
#' 
#' This is the contingency periodogram of Pierre Legedre and Pierre Dutielle to
#' test for periodicity in categorical time series. I have coded the function
#' so as to provide both the Fisher exact test and the asymptotic chi-square
#' test.
#' 
#' @param x A vector representing the categorical time series.
#' @param maxper the maximum lag (period) considered.
#' @param exact If TRUE the FISHER exact test is calculated
#' @return An object of class "contingency.periodogram" is returned consisting
#' of a matrix with a row for each period considered.  The columns are:
#' \item{exact.p}{the Fisher exact test at each lag (if exact=TRUE).}
#' \item{chi2}{the asymptotic chi-square value.} \item{df}{the chi-square
#' degrees-of-freedom.} \item{asympt.p}{the chi-squared asymptotic p-value.}
#' @references Legendre et al. (1981) The contingency periodogram: A method for
#' identifying rhytms in series of nonmetric ecological data. Journal of
#' Ecology, 69, 965-979. https://doi.org/10.2307/2259648
#' @keywords ts
#' @examples
#'     data(plodia)
#'     data<-as.factor((scale(plodia) > 0))
#'     fit <- contingency.periodogram(data, maxper = 9) 
#'     \dontrun{plot(fit)}
#' @export
#' @importFrom graphics lines plot 
#' @importFrom stats acf ar ar.mle chisq.test fisher.test lm pchisq pf predict qchisq quantile residuals rnorm spec.ar var
contingency.periodogram <- function(x, maxper = 6, exact=FALSE){
    t1 <- as.factor(x)
    n <- length(x)
    kji <- matrix(NA, nrow = maxper, ncol = 4)
    dimnames(kji) <- list(c(1:maxper), c("exact.p", "chi2", "df", "asympt.p"))
    for(i in 2:maxper) {
        t2 <- t(1:i)
        kast <- (as.integer(n/i)) * i
        t3 <- cbind(as.factor(t1[1:kast]), as.factor(t2))
        if(exact) kji[i, 1] <- as.vector(fisher.test(table(t3[, 1], t3[, 2]))[1])$p.value
        t4 <- chisq.test(table(t3[, 1], t3[, 2]))
        kji[i, 2] <- as.vector(t4$statistic)
        kji[i, 3] <- as.vector(t4$parameter)
        kji[i, 4] <- as.vector(t4$p.value)
    }
    res <- as.matrix(kji[-1,])
    class(res) <- "contingency.periodogram"
    res
}

##############################################################################################
#' Plot contingency periodograms
#' 
#' `plot' method for "contingency.periodogram" class object.
#' 
#' 
#' @param x an object of class "contingency.periodogram", usually, as a result
#' of a call to \code{contingency.periodogram}.
#' @param ... generic plot arguments.
#' @return A contingency periodogram is plotted. The line represents the
#' critical value based on the chi-squared test (95\%).
#' @seealso \code{\link{contingency.periodogram}}
#' @keywords ts
#' @export
plot.contingency.periodogram<-function(x, ...){
##############################################################################################
#this is the generic plot function for contingency.periodogram objects
##############################################################################################
 args.default <- list(xlab = "Period", ylab = "Chi-squared value", 
                  type="b")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  do.call(plot, c(list(x = c(2:(nrow(x)+1)), y = x[,"chi2"]), args))
    n <- nrow(x)
    lines(2:n, qchisq(0.95, 1:(n - 1)))
}

#' The order of a time series using cross-validation of the linear
#' autoregressive model (conditional least-squares).
#' 
#' A function to estimate the order of a time series using cross-validation of
#' the linear autoregressive model. Coefficients are estimated using
#' conditional least-squares. I coded this functions to estimate the order of ecological time series.
#' Bjornstad et al. (1998, 2001)
#' 
#' The time series is normalized prior to cross-validation.
#' 
#' Note that if the dynamics is highly nonlinear, the nonparametric
#' order-estimators (\code{\link{ll.order}}) may be more appropriate.  (I coded
#' this function to use for comparison with the nonparametric methods, because
#' these also uses (nonlinear) conditional least-squares.)
#' 
#' @param x A time series without missing values
#' @param order The candidate orders. The default is 1:5
#' @param n.cond The number of observation to condition on.  The default is 5
#' (must be >= max(order))
#' @param echo if TRUE a counter for the data points and the orders is produced
#' to monitor progress.
#' @return An object of class "lin.order" is returned consisting of the
#' following components: \item{order}{the grid of orders considered.}
#' \item{CVd}{the cross-validation errors across the grid of orders.}
#' @references Bjornstad, O.N., Begon, M., Stenseth, N. C., Falck, W., Sait, S. M. and Thompson, D. J. 1998. Population dynamics of the Indian meal moth: demographic stochasticity and delayed regulatory mechanisms. Journal of Animal Ecology 67:110-126. https://doi.org/10.1046/j.1365-2656.1998.00168.x
#' Bjornstad, O.N., Sait, S.M., Stenseth, N.C., Thompson, D.J. & Begon, M. 2001. Coupling and the impact of specialised enemies on the dimensionality of prey dynamics. Nature 401: 1001-1006. https://doi.org/10.1038/35059003
#' @seealso \code{\link{ll.order}}
#' @keywords ts
#' @examples
#' 
#'     data(plodia)
#'     fit <- lin.order.cls(sqrt(plodia), order=1:5)
#'     \dontrun{plot(fit)}
#'     summary(fit)
#' @export
lin.order.cls <- function(x, order = 1:5, n.cond = 5, echo = TRUE){
##############################################################################################
    ans <- as.data.frame(matrix(NA, ncol = 2, nrow = length(order)))
    names(ans) <- c("order", "CVd")
    ans[, 1] <- order
    n <- length(x)
    p <- length(order)
    cvseries <- (x - mean(x[(n.cond + 1):n]))/sqrt(var(x[(n.cond + 1):n]))
    xmat <- matrix(0, (n - n.cond), (p + 1))

    for(i in 1:length(order)){
        xmat[, i] <- cvseries[(n.cond + 1 - order[i]):(n - order[i])]
    }

    xmat[, (p + 1)] <- cvseries[(n.cond + 1):n]

    for(j in 1:p) {
        cv <- 0
        conv <- 0
        for(i in 1:(n - n.cond)) {
            dati <- xmat[ - i,  ]
            coef <- solve(t(dati[, (1:j)]) %*% dati[, (1:j)]) %*% t(dati[, (1:j)]) %*% dati[, (p + 1)]
            tpred <- t(xmat[i, 1:j]) %*% coef
            cv <- cv + (tpred - xmat[i, (p + 1)])^2
            }
        if(echo){cat("\n order ",j, " of ",p, "\r")}
        ans[[2]][j] <- cv/(n - n.cond)
    }
    class(ans) <- "lin.order"
    ans
}


##############################################################################################
#' Summarize linear cross-validation for time-series order
#' 
#' `summary' method for class "lin.order".
#' 
#' 
#' @param object an object of class "lin.order", usually, as a result of a call
#' to \code{lin.order.cls} or \code{lin.order.mle}.
#' @param \dots no other arguments currently allowed
#' @return A slightly prettyfied version of the object is printed.
#' @seealso \code{\link{lin.order.cls}}
#' @keywords ts
#' @export
summary.lin.order <- function(object, ...){
##############################################################################################
#this is the generic summary function for lin.order objects
##############################################################################################
    max <- object$order[order(object$CVd)][1]
    cat(paste("The estimated order is ", max, " with a cross-validation error of ", round(
        object$CVd[order(object$CVd)[1]], 3), "\n\n", sep = ""))
    out <- data.frame(order = object[[1]], CVd = object[[2]])
    out
}


#' Plot linear cross-validation for time-series order
#' 
#' `plot' method for class "lin.order".
#' 
#' 
#' @param x an object of class "lin.order", usually, as a result of a call to
#' \code{lin.order.cls} or \code{lin.order.mle}.
#' @param ... generic plot arguments.
#' @return A xy-plot of order against cross-validation error is produced.
#' @seealso \code{\link{lin.order.cls}}
#' @keywords ts
#' @export
plot.lin.order <- function(x, ...){
##############################################################################################
#this is the generic plot function for lin.order objects
##############################################################################################
args.default <- list(ylab = "Cross-validation error", xlab = "Order", 
                       type="b")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  do.call(plot, c(list(x = x$order, y = x$CVd), args))
}

#' Consistent nonlinear estimate of the order using local polynomial
#' regression.
#' 
#' A function to estimate the order of a time series using the nonparametric
#' order selection method of Cheng and Tong (1992, 1994) as modified by Yao &
#' Tong (1994; see also Fan, Yao & Tong 1996).  The method uses leave-one-out
#' cross-validation of the locally linear regression against lagged-abundances.
#' 
#' The time series is normalized prior to cross-validation.
#' 
#' A Gaussian kernel is used for the locally linear regression.
#' 
#' The bandwidth is optimized using cross-validation. If a single bandwidth is
#' provided, no cross validation of bandwidth will be carried out. Highly
#' nonlinear data will require more narrow bandwidths. If NA is returned it may
#' be because the min bandwidth considered is too small relative to the density
#' of data.
#' 
#' Missing values are NOT permitted.
#' 
#' If \code{deg} is set to 0, the order is estimated on the basis of the
#' Nadaraya-Watson (locally constant) estimator of the conditional expectation
#' against lagged-abundances (Cheng and Tong 1992, 1994). 
#' 
#' @param x A time series without missing values.
#' @param order The candidate orders. The default is 1:5.
#' @param step The time step for prediction.
#' @param deg The degree of the local polynomial.
#' @param bandwidth The candidate bandwidths to be considered.
#' @param cv if TRUE leave-one-out cross-validation will be performed.
#' @param echo if TRUE a counter shows the progress
#' @return An object of class "ll.order" is returned consisting of the
#' following components: \item{grid}{the grid of orders, bandwidths, and CV's.}
#' \item{grid$order}{the orders.} \item{grid$CV}{the cross-validation score
#' across the grid of orders and bandwidths. (If \code{cv = TRUE}).}
#' \item{grid$GCV}{the generalized cross-validation score.}
#' \item{grid$bandwidth}{the bandwidths.} \item{grid$df}{the degrees of freedom
#' of the fitted model.}
#' 
#' \item{order}{the vector of orders considered.} \item{deg}{The degree of the
#' local polynomial.}
#' @references Cheng, B. & Tong, H. (1992) On consistent nonparametric order
#' determination and chaos. Journal of Royal Statistical Society B, 54,
#' 427-449.
#' 
#' Cheng, B. & Tong, H. (1994) Orthogonal projection, embedding dimension and
#' sample size in chaotic time series from a statistical perspective.
#' Philosophical Transactions of the Royal Society London, A. , 348, 325-341. https://doi.org/10.1098/rsta.1994.0094
#' 
#' Fan, J., Yao, Q., & Tong, H. (1996) Estimation of conditional densities and
#' sensitivity measures in nonlinear dynamical systems. Biometrika, 83,
#' 189-206. ttps://doi.org/10.1093/biomet/83.1.189
#' 
#' Yao, Q. & Tong, H. (1994) Quantifying the influence of initial values on
#' non-linear prediction.  Journal of Royal Statistical Society B, 56, 701-725.
#' 
#' Bjornstad, O.N., Sait, S.M., Stenseth, N.C., Thompson, D.J., & Begon, M.
#' (2001) Coupling and the impact of specialised enemies on the dimensionality
#' of prey dynamics. Nature, 409, 1001-1006. https://doi.org/10.1038/35059003
#' 
#' Loader, C. (1999) Local Regression and Likelihood.  Springer, New York. https://doi.org/10.1007/b98858
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#'    fit <- ll.order(sqrt(plodia), order=1:3, bandwidth
#'                = seq(0.5, 1.5, by = 0.5)) 
#' 
#'     \dontrun{plot(fit)}
#' 
#'     summary(fit)
#' 
#' @export
#' @importFrom locfit locfit.raw dat
ll.order<- function(x, order = 1:5, step=1, deg = 2, bandwidth = c(seq(0.3, 1.5, by = 0.1), 2:10), cv=TRUE, echo=TRUE){
##############################################################################################
    res<-as.data.frame(matrix(NA, ncol = 6, nrow = length(order)*length(bandwidth)))
    names(res) <- c("order", "CV", "GCV", "bandwidth", "df", "GCV.df")

    bogrid<-expand.grid(bandwidth, order)
    res$order<-bogrid[,2]
    res$bandwidth<-bogrid[,1]

    T <- length(x)

    cvseries <- (x - mean(x[(max(order) + 1):T]))/sqrt(var(x[(max(order) + 1):T]))
    ldata<-mkx(cvseries, order+step-1)

    n<-dim(ldata)[1]

    for(i in 1:dim(bogrid)[1]){
    tmp<-NULL
    if(cv == TRUE){
        tmp<-locfit.raw(lpx(ldata[,1:(bogrid[i,2])], deg=deg, h=bogrid[i,1]),
                y=ldata[,(length(order)+1)], kern='gauss', ev=dat(cv=TRUE))
        res$CV[i]<--2*tmp$dp["lk"]/n
        if(res$CV[i]==0) res$CV[i]<-NA
        res$df[i]<-tmp$dp["df1"]
        }

        tmp<-locfit.raw(lpx(ldata[,1:(bogrid[i,2])], deg=deg, h=bogrid[i,1]),
                y=ldata[,(length(order)+1)], kern='gauss', ev=dat(cv=FALSE))

    #GCV (p 50 in Loader 1999)
    res$GCV[i]<--2*tmp$dp["lk"]*n/(n-tmp$dp["df1"])^2
    if(res$GCV[i]==0) res$GCV[i]<-NA
    #df
    res$GCV.df[i]<-tmp$dp["df1"]

    if(echo == TRUE){
            cat(i, " / ", dim(bogrid)[1], "\n")
    }
    }

    out<-list(grid=res, order=order, deg=deg, step=step, call=deparse(match.call()), cv=cv, x=x)
    class(out) <- "ll.order"
    out
}

#' Utility function
#' 
#' hack to make ll.order work with locfit >1.5. not to be called by the user.
#' 
#' not to be called by the user.
#' 
#' @param x \dots{}
#' @param nn \dots{}
#' @param h \dots{}
#' @param adpen \dots{}
#' @param deg \dots{}
#' @param acri \dots{}
#' @param scale \dots{}
#' @param style \dots{}
#' @keywords misc
#' @author Catherine Loader
#' @export
lpx<-function (x, nn = 0, h = 0, adpen = 0, deg = 2, acri = "none",
    scale = FALSE, style = "none"){
##############################################################################################
#locfit hack to make ll.order work with locfit >1.5
##############################################################################################
    x <- cbind(x)
    z <- as.list(match.call())
    z[[1]] <- z$nn <- z$h <- z$adpen <- z$deg <- z$acri <- z$scale <- z$style <- NULL
    #dimnames(x) <- list(NULL, z)
    if (missing(nn) & missing(h) & missing(adpen))
        nn <- 0.7
    attr(x, "alpha") <- c(nn, h, adpen)
    attr(x, "deg") <- deg
    attr(x, "acri") <- acri
    attr(x, "style") <- style
    attr(x, "scale") <- scale
    class(x) <- "lp"
    x
}

#' Summarize nonparametric cross-validation for time-series order
#' 
#' `summary' method for class "ll.order".
#' 
#' See \code{\link{ll.order}} for details.
#' 
#' @param object an object of class "ll.order", usually, as a result of a call
#' to \code{ll.order}.
#' @param GCV if TRUE (or if cross-validation was not done), uses GCV values.
#' @param \dots no other arguments currently allowed
#' @return A matrix summarizing the minimum cross-validation error (cv.min) and
#' the associated Gaussian-kernel bandwidth (bandwidth.opt) and model
#' degrees-of-freedom for each order considered.
#' @seealso \code{\link{ll.order}}
#' @keywords ts
#' @export
summary.ll.order <- function(object, GCV=FALSE, ...){
##############################################################################################
#this is the generic summary function for ll.order objects
##############################################################################################
    ans <- as.data.frame(matrix(NA, ncol = 4, nrow = length(object$order)))
    names(ans) <- c("order", "cv.min", "bandwidth.opt", "df")
    if(object$cv==FALSE) GCV = TRUE

    ans$GCV.min<-NA
    ans$GCV.bandwidth.opt <- NA
    ans$GCV.df <- NA

    ans[, 1] <- object$order
    for(i in 1:length(object$order)) {
       if(object$cv==TRUE){
        ans[i, 2] <- min(object$grid$CV[object$grid$order == i])
        wh<- which(object$grid$CV[object$grid$order == i]==ans[i,2])
        ans[i, 3] <- object$grid$bandwidth[object$grid$order == i][wh]
        ans[i, 4] <- object$grid$df[object$grid$order == i][wh]
        }

        ans$GCV.min[i] <- min(object$grid$GCV[object$grid$order == i], na.rm=TRUE)
        wh <- which(object$grid$GCV[object$grid$order == i]==ans$GCV.min[i])
        ans$GCV.bandwidth.opt[i] <- object$grid$bandwidth[object$grid$order == i][wh]
        ans$GCV.df[i] <- object$grid$GCV.df[object$grid$order == i][wh]
    }

    max <- ans$order[order(ans$cv.min)][1]
    cat(paste("The estimated order is ", max, " with a cross-validation error of ", round(ans$
        cv.min[order(ans$cv.min)[1]], 3), "\nand Gaussian bandwidth ", round(as.numeric(
        ans$bandwidth.opt[order(ans$cv.min)][1]), 3), ". (using a local polynomial with ", object$deg,
         " degrees).\n\n", sep = ""))
    ans
}

##############################################################################################
#' Plot nonparametric cross-validation for time-series order
#' 
#' `plot' method for class "ll.order".
#' 
#' See \code{\link{ll.order}} for details.
#' 
#' @param x an object of class "ll.order", usually, as a result of a call to
#' \code{ll.order}.
#' @param ... generic plot arguments.
#' @return A xy-plot of minimum cross-validation error against order is
#' produced.
#' @seealso \code{\link{ll.order}}
#' @keywords ts
#' @export
plot.ll.order <- function(x, ...){
##############################################################################################
#this is the generic plot function for ll.order objects
##############################################################################################
 args.default <- list(xlab = "Cross-validation error", ylab = "order", 
                       type="b")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
   ans <- as.data.frame(matrix(NA, ncol = 2, nrow = length(x$order)))
    names(ans) <- c("order", "cv.min")
    order<-x$order
    ans[, 1] <- order
    for(i in 1:length(order)) {
        ans[i, 2] <- min(x$grid$CV[x$grid$order == i])
    }
    do.call(plot, c(list(x = ans$order, y = ans$cv.min), args))
}

##############################################################################################
#' Print nonparametric cross-validation for time-series order
#' 
#' `print' method for class "ll.order".
#' 
#' See \code{\link{ll.order}} for details.
#' 
#' @param x an object of class "ll.order", usually, as a result of a call to
#' \code{ll.order}.
#' @param verbose if TRUE provides a raw-printing of the object.
#' @param \dots no other arguments currently allowed
#' @return A matrix summarizing the minimum cross-validation error (cv.min) and
#' the associated Gaussian-kernel bandwidth (bandwidth.opt) and model
#' degrees-of-freedom for each order considered.
#' @seealso \code{\link{ll.order}}
#' @keywords ts
#' @export
print.ll.order <- function(x, verbose = FALSE, ...){
##############################################################################################
#this is the generic print function for ll.order objects
#
#ARGUMENTS
#verbose   if FALSE, summary is used. If TRUE, the raw list is echoed
##############################################################################################
    x
    if(!verbose) {
    out <- summary(x)
    print(out)
    cat("\n\nFor a raw listing use print(x, verbose=TRUE)\n")
    }
    if(verbose) {
        print.default(x)
    }
}

##############################################################################################
#' Predict values from ll.order object.
#' 
#' Calculates the leave-one-out predicted values for the optimal ll.order
#' object
#' 
#' See \code{\link{ll.order}} for details.
#' 
#' @param object an object of class "ll.order", usually, as a result of a call
#' to \code{ll.order}.
#' @param \dots no other arguments currently allowed
#' @return A data frame with observed and predicted values for the optimal
#' ll-model is returned.
#' @seealso \code{\link{ll.order}}
#' @keywords ts
#' @export
predict.ll.order<- function(object, ...){
##############################################################################################
    x2=object$x
    ans <- as.data.frame(matrix(NA, ncol = 4, nrow = length(object$order)))
    names(ans) <- c("order", "cv.min", "bandwidth.opt", "df")
    ans[, 1] <- object$order
    for (i in 1:length(object$order)) {
        if (object$cv == TRUE) {
            ans[i, 2] <- min(object$grid$CV[object$grid$order ==
                i])
            wh <- which(object$grid$CV[object$grid$order == i] ==
                ans[i, 2])
            ans[i, 3] <- object$grid$bandwidth[object$grid$order ==
                i][wh]
            ans[i, 4] <- object$grid$df[object$grid$order ==
                i][wh]
        }
        if (object$cv == FALSE) {
        ans[i, 2] <- min(object$grid$GCV[object$grid$order ==
            i], na.rm = TRUE)
        wh <- which(object$grid$GCV[object$grid$order == i] ==
            ans$GCV.min[i])
        ans[i, 3] <- object$grid$bandwidth[object$grid$order ==
            i][wh]
        ans[i, 4] <- object$grid$GCV.df[object$grid$order ==
            i][wh]
    }
    }

    ord = ans$order[order(ans$cv.min)][1]
    bw=as.numeric(ans$bandwidth.opt[order(ans$cv.min)][1])
    deg=object$deg
    step=object$step

    res<-data.frame(obs=x2, pred=rep(NA, length(x2)))

    bogrid<-expand.grid(bw, ord)

    T <- length(x2)

    tmu=mean(x2[(max(ord) + 1):T])
    tsd=sqrt(var(x2[(max(ord) + 1):T]))
    cvseries <- (x2 - tmu)/tsd
    ldata<-mkx(cvseries, step:(ord+step-1))

    n<-dim(ldata)[1]

    for(k in 1:(T-ord)){
        tmp<-NULL
        tdata=ldata[-k,]
        tmp<-locfit.raw(lpx(tdata[,1:(bogrid[,2])], deg=deg, h=bogrid[,1]),
                y=tdata[,ord+1], kern='gauss', ev=ldata[k,1:(bogrid[,2])])
        res$pred[k+ord]=tsd*predict(tmp)+tmu
        }
      res
}

##############################################################################################
#' Nonlinear forecasting at varying lags using local polynomial regression.
#' 
#' A wrapper function around \code{ll.order} to calculate prediction profiles
#' (a la Sugihara and May 1990 and Yao and Tong 1994). The method uses
#' leave-one-out cross-validation of the local regression (with CV optimized
#' bandwidth) against lagged-abundances at various lags.
#' 
#' see \code{\link{ll.order}} for details.
#' 
#' @aliases prediction.profile.ll print.ppll
#' @param x A time series without missing values.
#' @param step The vector of time steps for forward prediction.
#' @param order The candidate orders. The default is 1:5.
#' @param deg The degree of the local polynomial.
#' @param bandwidth The candidate bandwidths to be considered.
#' @return An object of class "ppll" consisting of a list with the following
#' components: \item{step}{the prediction steps considered.} \item{CV}{the
#' cross-validation error.} \item{order}{the optimal order for each step.}
#' \item{bandwidth}{the optimal bandwidth for each step.} \item{df}{the degrees
#' of freedom for each step.}
#' @seealso \code{\link{ll.order}}
#' @references Sugihara, G., and May, R.M. (1990) Nonlinear forecasting as a
#' way of distinguishing chaos from measurement error in time series. Nature
#' 344, 734-741. https://doi.org/10.1038/344734a0
#' 
#' Yao, Q. and Tong, H. (1994) Quantifying the influence of initial values on
#' non-linear prediction.  Journal of Royal Statistical Society B, 56, 701-725.
#' 
#' Fan, J., Yao, Q., and Tong, H. (1996) Estimation of conditional densities
#' and sensitivity measures in nonlinear dynamical systems. Biometrika, 83,
#' 189-206. https://doi.org/10.1093/biomet/83.1.189
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#'      fit <- prediction.profile.ll(sqrt(plodia), step=1:3, order=1:3,
#'           bandwidth = seq(0.5, 1.5, by = 0.5))
#' 
#'     \dontrun{plot(fit)}
#' @export
prediction.profile.ll<- function(x, step=1:10,order = 1:5, deg = 2, bandwidth = c(seq(0.3, 1.5, by = 0.1), 2:10)){
##############################################################################################
    res<-as.data.frame(matrix(NA, ncol=5, nrow=length(step)))
    names(res)<- c("step", "CV", "order", "bandwidth", "df")
    res$step<-step
    for(k in 1:length(step)){
          tmp<-ll.order(x, order = order, step=step[k], deg = deg, bandwidth = bandwidth, cv=TRUE, echo=FALSE)
          wh<-which(tmp$grid$CV==min(tmp$grid$CV))
          res[k,"CV"]<-tmp$grid$CV[wh]
          res[k,"bandwidth"]<-tmp$grid$bandwidth[wh]
          res[k,"order"]<-tmp$grid$order[wh]
          res[k,"df"]<-tmp$grid$df[wh]
          cat(k, '\n')
          }
    res2<-list(ppll=res)
    class(res2) <- "ppll"
    return(res2)
}

#################################I############################################################
#' Plot function for prediction profile objects
#' 
#' `plot' method for class "ppll".
#' 
#' See \code{\link{prediction.profile.ll}} for details.
#' 
#' @param x an object of class "ppll", usually, as a result of a call to
#' \code{prediction.profile.ll}.
#' @param ... generic plot arguments.
#' @return A xy-plot of one minus the cross-validation error (i.e. the
#' prediction accuracy against prediction time step.
#' @seealso \code{\link{prediction.profile.ll}}
#' @keywords ts
#' @export
plot.ppll<-function(x, ...){
##############################################################################################
args.default <- list(ylab = "Predictability", xlab = "Time lag", 
                       type="b", ylim=c(min(c(0,1-x$ppll$CV)), 1))
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  do.call(plot, c(list(x = x$ppll$step, y = 1-x$ppll$CV), args))
}

##############################################################################################
#' Nonlinear forecasting of local polynomial `empirical dynamic model'.
#' 
#' A function to forcaste a local polynomial `empirical dynamic model'.
#' 
#' The function produces a nonlinear (nonparametric) forecast using the
#' conditional mean method of Fan et al (1996). A Gaussian kernel is used for
#' the local polynomial autoregression.
#' 
#' The bandwidth and order is best estimated with the
#' \code{\link{ll.order}}-function.
#' 
#' Missing values are NOT permitted.
#' 
#' If \code{deg} is set to 0, the forecast uses the Nadaraya-Watson (locally
#' constant) estimator of the conditional expectation against lagged-abundances.
#' 
#' @param x A time series without missing values.
#' @param order The order for the nonparametric (local polynomial)
#' autoregression.
#' @param bandwidth The bandwidth for the nonparametric (local polynomial)
#' autoregression.
#' @param len The length of the predicted time-series. If NA the length of the
#' training time series will be used.
#' @param deg The degree of the local polynomial.
#' @return A time series with the nonlinear (nonparametric) forecast is
#' returned
#' @seealso \code{\link{ll.order}}
#' @references Fan, J., Yao, Q., & Tong, H. (1996) Estimation of conditional
#' densities and sensitivity measures in nonlinear dynamical systems.
#' Biometrika, 83, 189-206. https://doi.org/10.1093/biomet/83.1.189
#' 
#' Loader, C. (1999) Local Regression and Likelihood.  Springer, New York. https://doi.org/10.2307/1270956
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#'    sim1 <- ll.edm(sqrt(plodia), order=2, bandwidth = 1.5) 
#' @export
ll.edm=function (x, order, bandwidth, len=NA, deg = 2){
##############################################################################################
    T <- length(x)
    if(is.na(len)){len=T}
    mu=mean(x[(order + 1):T])
    sig=sqrt(var(x[(order + 1):T]))
    cvseries <- (x - mu)/sig
    ldata <- data.frame(mkx(cvseries, 1:order))
    names(ldata)[order+1]="Y"
    xnam=names(ldata)[1:order]

    tmp <- eval(parse(text=paste("locfit(Y~lp(", paste(xnam, collapse= ','), ", deg = deg, h = bandwidth), kern = 'gauss', data=ldata)")))

    pdata=ldata[1,1:order]
    ptmp=predict(tmp, newdata=as.matrix(pdata))

    sim=NA
    sim[1]=ptmp
    for(i in 2:len){
    pdata=cbind(sim[i-1],pdata[1,1:(order-1)]) 
    if(any(!is.finite(as.matrix(pdata)))){cat("Inf produced \n"); break}  
    ptmp=predict(tmp, newdata=as.matrix(pdata))
    sim[i]=ptmp
    }
    sim=sim*sig+mu
return(sim)
}


##############################################################################################
#' The Lomb periodogram for unevenly sampled data
#' 
#' The function to estimate the Lomb periodogram for a spectral analysis of
#' unevenly sampled data.
#' 
#' This is the Lomb periodogram to test for periodicity in time series of
#' unevenly sampled data.
#' 
#' Missing values should be deleted in both x and y before execution.
#' 
#' @param y vector of length n representing the unevenly sampled time series.
#' @param x the a vector (of length n) representing the times of observation.
#' @param freq the frequencies at which the periodogram is to be calculated. If
#' NULL the canonical frequencies (the Fourier frequencies) are used.
#' @return An object of class "lomb" is returned consisting of the following
#' components: \item{freq}{the frequencies as supplied.} \item{spec}{the
#' estimated amplitudes at the different frequencies.} \item{f.max}{the
#' frequency of maximum amplitude.} \item{per.max}{the corresponding period of
#' maximum amplitude.} \item{p}{the level of significance associated with the
#' max period.}
#' @references Lomb, N.R. (1976) Least-squares frequency-analysis of unequally
#' spaced data. Astrophysics and Space Science 39, 447-462.
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#'     y <- sqrt(plodia)
#'     x <- 1:length(y) 
#' 
#'     #make some missing values
#'     y[10:19] <- NA; x[10:19] <- NA 
#'     #omit NAs
#'     y <- na.omit(y); x <- na.omit(x) 
#' 
#'     #the lomb p'gram
#'     fit <- spec.lomb(y, x) 
#'     summary(fit)
#'     \dontrun{plot(fit)}
#' @export
spec.lomb <- function (y=stop("no data arg"), x=stop("no time arg"), freq=NULL){
##############################################################################################
  if(is.null(freq)){
    nyear <- max(x)-min(x)+1
    f <- seq(0,.5,length=nyear/2)
  }
  else{
    f <- freq
  }

  # Check arguments
  if (length(y) != length(x)) stop("y and x different lengths");
  if (min(f) < 0 || max(f) > 1) stop("freq must be between 0 and 1");
  if (min(f) == 0 ) f <- f[f>0];    # Get rid of zeros

  nt <- length(x);          # Number of datapoints
  nf <- length(f);          # Number of frequencies
  ones.t <- rep(1,nt);          # Useful unit vectors
  ones.f <- rep(1,nf);

  ## Convert to angular frequencies
  omega <- 2 * pi * f;

  ## Stats of the time series
  hbar <- mean(y);
  hvar <- var(y);
  hdev <- y - hbar;

  ## Calculate the vector of taus
  two.omega.t <- 2 * omega %*% t(x);
  sum.sin <- sin(two.omega.t) %*% ones.t;
  sum.cos <- cos(two.omega.t) %*% ones.t;
  tau <- atan(sum.sin/sum.cos) / (2*omega);

  ## Calculate the trig functions that go into the main expression
  t.m.tau <- (ones.f %*% t(x)) - (tau %*% t(ones.t));
  omega.ttau <- (omega %*% t(ones.t)) * t.m.tau;
  sin.ott <- sin(omega.ttau);
  cos.ott <- cos(omega.ttau);
  z <- ((cos.ott %*% hdev)^2 / ((cos.ott^2) %*% ones.t) +
    (sin.ott %*% hdev)^2 / ((sin.ott^2) %*% ones.t)) / (2 * hvar);

  max <- z == max(z,na.rm=TRUE)
  max <- max[is.na(max)==FALSE]
  P <- 1 - (1-exp(-z[max]))^(length(x))

  res <- list(spec=z[,1], freq=f, f.max=f[max], per.max=1/f[max], p = P)
  class(res) <- "lomb"
  res
}

##############################################################################################
#' Plot Lomb periodograms
#' 
#' `plot' method for objects of class "lomb".
#' 
#' 
#' @param x an object of class "lomb", usually, as a result of a call to
#' \code{spec.lomb}.
#' @param ... generic plot arguments.
#' @return A Lomb periodogram is composed of a xy-plot of amplitude against
#' frequency.
#' @seealso \code{\link{spec.lomb}}
#' @keywords ts
#' @export
plot.lomb <- function(x, ...){
##############################################################################################
args.default <- list(xlab = "Frequency", ylab = "Amplitude", 
                       type="l")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  
  do.call(plot, c(list(x = x$freq, y = x$spec), args))
}

##############################################################################################
#' Summarizes Lomb periodograms
#' 
#' `summary' method for objects of class "lomb".
#' 
#' 
#' @param object an object of class "lomb", usually, as a result of a call to
#' \code{spec.lomb}.
#' @param ... generic plot arguments.
#' @return A list summarizing the analysis is printed: \item{period}{the
#' dominant period.} \item{p.val}{the p.value.}
#' @seealso \code{\link{spec.lomb}}
#' @keywords ts
#' @export
summary.lomb <- function(object, ...){
##############################################################################################
list(period=object$per.max,p.val=object$p)
}


##############################################################################################
#' Lagrange multiplier test for additivity in a timeseries
#' 
#' add.test is a function to test the permissibility of the additive
#' autoregressive model:
#' 
#' N(t) = f1(N(t-1)) + f2(N(t-2)) + ... + fd(N(t-d)) + e(t )
#' 
#' against the alternative:
#' 
#' N(t) = F(N(t-1), N(t-2), ..., N(t-d)) + e(t)
#' 
#' This is the Lagrange multiplier test for additivity developed by Chen et al.
#' (1995: test II).
#' 
#' @param x A time series (vector without missing values).
#' @param order a scalar representing the order to be considered.
#' @param n.cond The number of observation to condition on.  The default is
#' \code{order} (must be >= \code{order})
#' @return a vector is returned consisting of the asymtpotic chi-square value,
#' the associated d.f.  and asymptotic p.val for the test of additivity.
#' @references Chen, R., Liu, J.S. & Tsay, R.S. (1995) Additivity tests for
#' nonlinear autoregression. Biometrika, 82, 369-383. https://doi.org/10.1093/biomet/82.2.369
#' 
#' Bjornstad, O.N., Begon, M., Stenseth, N.C., Falck, W., Sait, S.M., &
#' Thompson, D.J. (1998) Population dynamics of the Indian meal moth:
#' demographic stochasticity and delayed regulatory mechanisms. Journal of
#' Animal Ecology, 67, 110-126. https://doi.org/10.1046/j.1365-2656.1998.00168.x
#' @keywords ts
#' @examples
#' 
#'      data(plodia)
#'      add.test(sqrt(plodia), order = 3)
#' @export
#' @importFrom acepack ace
add.test<-function(x, order, n.cond = FALSE){
##############################################################################################
    resid.ace <- function(aceobj){
    aceobj$ty - apply(aceobj$tx, 1, sum)
    }
    if(!n.cond){
        n.cond <- order
        }
    nx <- length(x)
    tmp.mkx <- matrix(0, (nx - n.cond), (n.cond + 1))
    for(i in 1:n.cond)
        tmp.mkx[, i] <- x[(n.cond + 1 - i):(nx - i)]
    tmp.mkx[, (n.cond + 1)] <- x[(n.cond + 1):nx]

    tmp.ace <- ace(tmp.mkx[, 1:order], tmp.mkx[, (n.cond + 1)], lin = 0)
    tmp.resid1 <- resid.ace(tmp.ace)
    h <- 0
    K <- ((order - 1) * order * (order + 7))/6
    tmp.resid2 <- matrix(NA, ncol = K, nrow = dim(tmp.mkx)[1])
    for(i in 1:order)
        for(j in i:order)
            if(i != j) {
                tmp <- apply(tmp.mkx[, c(i, j)], 1, prod)
                h <- h + 1
                tmp.resid2[, h] <- resid.ace(ace(tmp.mkx[, 1:order], tmp, lin = 0))
            }
    for(i in 1:order)
        for(j in i:order)
            for(k in j:order)
                if(i != j | i != k | j != k) {
                  tmp <- apply(tmp.mkx[, c(i, j, k)], 1, prod)
                  h <- h + 1
                  tmp.resid2[, h] <- resid.ace(ace(tmp.mkx[, 1:order], tmp, lin = 0))
                }
    resid.lm <- lm(tmp.resid1 ~ tmp.resid2)
    unlist(list(chisq = round(dim(tmp.mkx)[1] * summary(resid.lm)$r.squared,4), df=K,
               p.val = round(1 - pchisq(dim(tmp.mkx)[1] * summary(resid.lm)$r.squared, K),4)))
}


##############################################################################################
#' A Tukey one-degree-of-freedom test for linearity in time series.
#' 
#' a function to test the permissibility of the linear autoregressive model:
#' 
#' N(t) = a0 + a1N(t-1) + a2N(t-2) + ... + adN(t-d) + e(t )
#' 
#' against the alternative:
#' 
#' Nt = F(N(t-1), N(t-2), ..., N(t-d)) + e(t)
#' 
#' This is the Tukey one-degree-of-freedom test of linearity developed by Tsay
#' (1986). Orders up to 5 is permissible. [although the code is easily
#' extended].
#' 
#' @param x A time series (vector without missing values).
#' @param order a scalar representing the order to be considered.
#' @return A vector is returned consisting of the asymtpotic F-value, the
#' associated numerator and denominator d.f.'s and asymptotic p.val for the
#' test of linearity
#' @references Tsay, R.S. (1986) Nonlinearity tests for time series.
#' Biometrika, 73, 461-466. https://doi.org/10.1093/biomet/73.2.461
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#'    lin.test(sqrt(plodia), order = 3)
#' @export
lin.test <- function(x, order){
##############################################################################################
    nx <- length(x)
    Y <- matrix(0, (nx - order), (order + 1))
    for(i in 1:order)
        Y[, i] <- x[(order + 1 - i):(nx - i)]
    Y[, (order + 1)] <- x[(order + 1):nx]
    D <- dim(Y)
    ur0 <- residuals(switch(as.character(D[2] - 1),
        "1" = lm(Y[, D[2]] ~ Y[, 1]),
        "2" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2]),
        "3" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3]),
        "4" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3] + Y[, 4]),
        "5" = lm(Y[, D[2]] ~ Y[, 1] + Y[, 2] + Y[, 3] + Y[, 4] + Y[, 5])
        ))
    ur1 <- residuals(switch(as.character(D[2] - 1),
        "1" = lm(ur0 ~ poly(Y[, 1], degree= 2)),
        "2" = lm(ur0 ~ poly(Y[, 1], Y[, 2], degree= 2)),
        "3" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], degree= 2)),
        "4" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], Y[, 4], degree= 2)),
        "5" = lm(ur0 ~ poly(Y[, 1], Y[, 2], Y[, 3], Y[, 4], Y[, 5], degree= 2))
        ))
    m <- switch(as.character(D[2] - 1),
        "1" = 1,
        "2" = 3,
        "3" = 6,
        "4" = 10,
        "5" = 15)
    Fval <- ((sum(ur0^2) - sum(ur1^2))/m)/(sum(ur1^2)/(D[1] - m - 1))
    pval <- 1 - pf(Fval, m, D[1] - m - 1)
    unlist(list(order = order, F = round(Fval, 4), df1 = m, df2 = D[1] - m - 1, p = round(pval, 4)))
}

##############################################################################################
#'  Ljung-Box test for whiteness in a time series.
#' 
#' portman.Q uses the cummulative ACF to test for whiteness of a time series.
#' 
#' This is the Ljung-Box version of the the Portemanteau test for whiteness
#' (Tong 1990). It may in particular be usefull to test for whiteness in the
#' residuals from time series models.
#' 
#' @param x A time series (vector without missing values).
#' @param K the maximum lag of the ACF to be used in the test.
#' @return A vector is returned consisting of the asymtpotic chi-square value,
#' the associated d.f. and asymptotic p.val for the test of whiteness.
#' @references Tong, H. (1990) Non-linear time series : a dynamical system
#' approach. Clarendon Press, Oxford.
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#'    portman.Q(sqrt(plodia), K = 10) 
#' 
#'    fit <- ar(sqrt(plodia)) 
#'    portman.Q(na.omit(fit$resid), K = 10) 
#' @export
portman.Q <- function(x, K){
##############################################################################################
    Q <- 0
    n <- length(x)
    p <- acf(x, plot = FALSE, lag.max = K)$acf[2:(K + 1)]
    for(k in 1:K)
        Q <- Q + p[k]^2/(n - k)
    Q <- n * (n + 2) * Q
    res <- list(chisq = round(Q,4), df = K, p.val = round(1 - pchisq(Q, K),4))
    unlist(res)
}

##############################################################################################
#' Confidence interval for the ar-spectrum and the dominant period.
#' 
#' A function to estimate a "confidence interval" for the power spectrum and in
#' particular a confidence interval for the dominant period. The function uses
#' resampling of the autoregressive parameters to attain the estimate.
#' 
#' A "confidence interval" for the periodogram is obtained by resampling the
#' ar-coefficients using the variance-covariance matrix from the ar.mle object.
#' 
#' If a zero'th order process is chosen by using the AIC criterion, a first
#' order process will be used.
#' 
#' If the dynamics is highly nonlinear, the parametric estimate of the power
#' spectrum may be inappropriate.
#' 
#' @param x A time series without missing values.
#' @param order a scalar representing the order to be considered. If
#' \code{"aic"} the orderis be selected automatically using the AIC criterion.
#' @param resamp the number of resamples of the ar-coefficients from the
#' covariance matrix.
#' @param nfreq the number of points at which to save the value for the power
#' spectrum (and confidence envelope).
#' @param echo If \code{TRUE}, a counter for each nrun shows the progress.
#' @return An object of class "specar.ci" is returned consisting of the
#' following components: \item{order}{the ar-order.} \item{spectrum$freq}{the
#' spectral frequencies.} \item{spectrum$spec}{the estimated power-spectrum of
#' the data.} \item{resamp$spectrum}{gives the quantile summary for the
#' resampling distribution of the spectral powers.} \item{resamp$maxfreq}{the
#' full vector of output for the resampled max.frequencies.}
#' @seealso \code{\link{plot.specar.ci}} \code{\link{summary.specar.ci}}
#' @keywords ts
#' @examples
#' 
#'    data(plodia)
#' 
#' 
#'     fit <- specar.ci(sqrt(plodia), order=3, resamp=10) 
#' 
#'     \dontrun{plot(fit, period=FALSE)}
#' 
#'     summary(fit)
#' @export
specar.ci <- function(x, order, resamp = 500, nfreq = 100, echo = TRUE){
##############################################################################################
    if(order == "aic") {
        s.ar <- ar(x, aic = TRUE)
        if(s.ar$order == 0) {
            s.ar <- ar.mle(x, order.max = 1, aic = FALSE)
            order <- s.ar$order
        }
        else {
            order <- s.ar$order
        }
    }
    else {
        s.ar <- ar.mle(x, order.max = order, aic = FALSE)
    }
    real <- spec.ar(s.ar, n.freq = nfreq, plot = FALSE)
    trekk <- matrix(NA, ncol = nfreq, nrow = resamp)
    maxfreq <- 1:resamp

    s.ar2<-s.ar

    for(i in 1:resamp) {
        if(order > 1) {
            vs <- svd(s.ar$asy.var.coef)
            vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
            ans <- matrix(rnorm(order), nrow = 1) %*% vsqrt
            ans <- sweep(ans, 2, s.ar$ar, "+")
        }
        if(order == 1) {
            ans <- rnorm(1, s.ar$ar, sqrt(s.ar$var.coef))
        }
        s.ar2$ar <- as.vector(ans)
        s.ar.mle3 <- spec.ar(s.ar2, n.freq = nfreq, plot = FALSE)
        trekk[i,  ] <- s.ar.mle3$spec
        maxfreq[i] <- s.ar.mle3$freq[match(max(s.ar.mle3$spec), s.ar.mle3$spec)]
        if(echo) {
            cat(i, "\r")
        }
    }
    trekk <- apply(trekk, 2, quantile, probs = c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9,
        0.95, 0.975, 1))
    dimnames(trekk) <- list(c(0, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 1),
        NULL)
    res <- list(spectrum = list(freq = real$freq, spec = real$spec, maxfreq = real$freq[match(
        max(real$spec), real$spec)]), order = order, resamp = list(spectrum = trekk,
        maxfreq = maxfreq))
    class(res) <- "specar.ci"
    res
}

##############################################################################################
#' Summarize ar-spectra with CI's
#' 
#' `summary' method for objects of class "specar.ci".
#' 
#' 
#' @param object an object of class "specar.ci", usually, as a result of a call
#' to \code{specar.ci}.
#' @param period If TRUE the summary is given in terms of the period, if false
#' it is in terms of the frequency
#' @param ... generic plot arguments.
#' @return A list summarizing the analysis is printed: \item{period}{the
#' dominant period.} \item{p.val}{the p.value.}
#' @seealso \code{\link{specar.ci}}
#' @keywords ts
#' @export
summary.specar.ci <- function(object, period = TRUE, ...){
##############################################################################################
#this is the generic summary function for specar.ci objects
#
#ARGUMENTS
#period    if T, the summary is in terms of period (=1/freq) rather than frequency
##############################################################################################
    if(period == TRUE) {
        ans <- list(period = unlist(list(period = 1/object$spectrum$maxfreq,
            ci.lower = as.vector(1/quantile(object$resamp$maxfreq, probs = c(
            0.975))), ci.upper = as.vector(1/quantile(object$resamp$maxfreq,
            probs = c(0.025))))), order = as.vector(object$order),
            resamp.summary = summary(1/object$resamp$maxfreq))
    }
    if(period == FALSE) {
        ans <- list(frequency = unlist(list(freq = object$spectrum$maxfreq,
            ci.lower = as.vector(quantile(object$resamp$maxfreq, probs = c(
            0.025))), ci.upper = as.vector(quantile(object$resamp$maxfreq,
            probs = c(0.975))))), order = as.vector(object$order),
            resamp.summary = summary(object$resamp$maxfreq))
    }
    ans
}

##############################################################################################
#' Plot ar-spectra with CI's
#' 
#' `plot' method for class "specar.ci".
#' 
#' 
#' @param x an object of class "specar.ci", usually, as a result of a call to
#' \code{\link{specar.ci}}.
#' @param period if TRUE x-axis is period, if FALSE frequency.
#' @param ... generic plot arguments.
#' @return A xy-plot of amplitude against period (or frequency).
#' @seealso \code{\link{specar.ci}}
#' @keywords ts
#' @importFrom grDevices gray
#' @importFrom graphics polygon
#' @export
plot.specar.ci <- function(x, period = TRUE, ...){
##############################################################################################
#this is the generic plot function for specar.ci objects
#
#ARGUMENTS
#period    if T, the summary is in terms of period (=1/freq) rather than frequency
##############################################################################################
           n <- length(x$spectrum$freq)
 
    if(period == TRUE) {
args.default <- list(xlab = "Period", ylab = "Amplitude", ylim = range(x$resamp$spectrum[c(2, 10), 2:n]),
                type = "l")
  args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  x2=1/x$spectrum$freq[2:n]
}  

    if(period == FALSE) {
args.default <- list(xlab = "Frequency", ylab= "Amplitude", ylim = range(x$resamp$spectrum[c(2, 10), 2:n]),
                 type = "l")
   args.input <- list(...)
  args <- c(args.default[!names(args.default) %in% names(args.input)], args.input)
  x2=x$spectrum$freq[2:n]
}  

do.call(plot, c(list(x = x2, y = x$spectrum$spec[2:n]), args))
  if(period==TRUE){
      polygon(c(1/x$spectrum$freq[2:n], rev(1/x$spectrum$freq[2:n])), 
            c(x$resamp$spectrum[
                "0.025", 2:n], 
              rev(x$resamp$spectrum[
                "0.975", 2:n])), col = gray(0.8), 
            lty = 0)
     lines(x2, x$spectrum$spec[2:n]) 
   }
  if(period==FALSE){
     polygon(c(x$spectrum$freq[2:n], rev(x$spectrum$freq[2:n])), 
            c(x$resamp$spectrum[
                "0.025", 2:n], 
              rev(x$resamp$spectrum[
                "0.975", 2:n])), col = gray(0.8), 
            lty = 0)
     lines(x2, x$spectrum$spec[2:n]) 
  }
}

##############################################################################################
#' Utility function
#' 
#' A function to create matrix of lagged time series.  Called by various
#' functions.
#' 
#' If lags is \code{c(1,4)}, say, then the function returns a matrix that
#' consist of columns x(t-1), x(t-4), x(t).
#' 
#' @param x A univariate time series.
#' @param lags The vector of time lags.
#' @return A matrix of lagged abundances. The last column is the current
#' @author Upmanu Lall
#' @references Lall, U. & Sharma, A. (1996) A nearest neighbor
#' bootstrap for time series resampling. Water Resources Research, 32, 679-693. https://doi.org/10.1029/95wr02966
#' @keywords misc
#' @export
mkx<-function(x, lags){
##############################################################################################
# U. Lall and A. Sharma - Lall, U. & Sharma, A. (1996) A nearest neighbor
#bootstrap for time series resampling. Water Resources Research, 32, 679-693.
#
#function to create matrix of lagged time series.
#x is the univariate time series
#lags is the vector of lags. If lags contains 1 and 4 (say) then
#x1 (output) would consist of xt-1, xt-4, xt.
    nx <- length(x)
    nl <- length(lags)
    ml <- max(lags)
    x1 <- matrix(0, (nx - ml), (nl + 1))
    for(i in 1:nl)
        x1[, i] <- x[(ml + 1 - lags[i]):(nx - lags[
            i])]
    x1[, (nl + 1)] <- x[(ml + 1):nx]
    x1
}



