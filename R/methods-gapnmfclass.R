gapnmfclass <- function(x, alpha, a, b, K, smoothness=100) {
    M <- dim(x)[1]
    N <- dim(x)[2]

    rhow <- 10000 * matrix(rgamma(M*K, smoothness, smoothness), nrow=M, ncol=K)
    tauw <- 10000 * matrix(rgamma(M*K, smoothness, smoothness), nrow=M, ncol=K)
    rhoh <- 10000 * matrix(rgamma(M*K, smoothness, smoothness), nrow=K, ncol=M)
    tauh <- 10000 * matrix(rgamma(M*K, smoothness, smoothness), nrow=K, ncol=M)
    rhot <- K * 10000 * matrix(rgamma(M*K, smoothness, smoothness),
                               nrow=K, ncol=1)
    taut <- 1 / K * 10000 * matrix(rgamma(M*K, smoothness, smoothness),
                                   nrow=K, ncol=1)

    Ew <- computegigexpectations(a, rhow, tauw)
    Eh <- computegigexpectations(b, rhoh, tauh)
    Et <- computegigexpectations(alpha/K, rhot, taut)

    obj <- new("gapnmfclass",
               x=x / mean(x),
               alpha=alpha,
               a=a,
               b=b,
               K=as.integer(K),
               M=M,
               N=N,
               rhow=rhow,
               tauw=tauw,
               rhoh=rhoh,
               tauh=tauh,
               rhot=rhot,
               taut=taut,
               Ew=Ew$Ex,
               Ewinv=Ew$Exinv,
               Ewinvinv=Ew$Exinv ^ (-1),
               Eh=Eh$Ex,
               Ehinv=Eh$Exinv,
               Ehinvinv=Eh$Exinv ^ (-1),
               Et=Et$Ex,
               Etinv=Et$Exinv,
               Etinvinv=Et$Exinv ^ (-1))

}

setMethod("goodk", "gapnmfclass",
    function(obj, cutoff) {
        if (missing(cutoff)) {
            cutoff <- 1e-10 * max(obj@x)
        }

        # NOT QUITE RIGHT
        powers <- obj@Et * cbind(apply(obj@Ew, 2, max)) * 
                  rbind(apply(obj@Eh, 1, max))
        sorted <- order(powers, descending=TRUE)
        temp <- powers[sorted]
        tmpk <- which(temp / max(temp) > cutoff)
        tmpk <- tmpk[length(tmpk)]
        goodk <- sorted[1:tmpk]

        if (powers[goodk[length(goodk)]] < cutoff) {
            goodk[length(goodk)] <- NULL
        }
    }
)

setMethod("recomputeexpectations", "gapnmfclass",
    function(obj) {
        w <- computegigexpectations(obj@a, obj@rhow, obj@tauw)
        h <- computegigexpectations(obj@b, obj@rhoh, obj@tauh)
        t <- computegigexpectations(obj@alpha / obj@K, obj@rhot, obj@taut)

        obj@Ew <- w$Ex
        obj@Ewinv <- w$Exinv
        obj@Ewinvinv <- obj@Ewinv ^ (-1)
        obj@Eh <- h$Ex
        obj@Ehinv <- h$Exinv
        obj@Ehinvinv <- obj@Ehinv ^ (-1)
        obj@Et <- t$Ex
        obj@Etinv <- t$Exinv
        obj@Etinvinv <- obj@Etinv ^ (-1)
    }
)

setMethod("bound", "gapnmfclass",
    function(obj, varargin) {
        verbose <- 0
        if (varargin == "verbose") {
            verbose <- 1
            lastscore <- 0
        }

        score <- 0

        # FINISH THIS!!!!
    }
)

setMethod("xbar", "gapnmfclass",
    function (obj, goodk) {
        if (missing(goodk)) {
            goodk <- 1:obj@K
        }
        return(obj@Ew[, goodk] %*% diag(obj@Et[goodk]) %*% obj@Eh[goodk, ])
    }
)

setMethod("xtwid", "gapnmfclass",
    function(obj, goodk) {
        if (missing(goodk)) {
            goodk <- 1:obj@K
        }

        return(obj@Ewinvinv[, goodk] %*% diag(obj@Etinvinv[goodk]) %*% 
               obj@Ehinvinv[goodk, ])
    }
)
