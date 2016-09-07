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
