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

    Ew <- computegigexpecations(a, rhow, tauw)
    Eh <- computegigexpecations(b, rhoh, tauh)
    Et <- computegigexpecations(alpha/K, rhot, taut)

    obj <- new("gapnmfclass",
               x=x / mean(x),
               alpha=alpha,
               a=a,
               b=b,
               K=K,
               M=M,
               N=N,
               rhow=rhow,
               tauw=tauw,
               rhoh=rhoh,
               tauh=tauh,
               rhot=rhot,
               taut=taut,
               Ew=Ew$Ew,
               Ewinv=Ew$Ewinv,
               Ewinvinv=Ew$Ewinv ^ (-1),
               Eh=Eh$Ew,
               Ehinv=Eh$Ewinv,
               Ehinvinv=Eh$Ewinv ^ (-1),
               Et=Et$Ew,
               Etinv=Et$Ewinv,
               Etinvinv=Et$Ewinv ^ (-1))

}
