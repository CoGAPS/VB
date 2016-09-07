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
               taut=taut)

#     [obj.Ew obj.Ewinv] = computegigexpectations(a, obj.rhow, obj.tauw);
#     obj.Ewinvinv = obj.Ewinv.^-1;
#     [obj.Eh obj.Ehinv] = computegigexpectations(b, obj.rhoh, obj.tauh);
#     obj.Ehinvinv = obj.Ehinv.^-1;
#     [obj.Et obj.Etinv] = computegigexpectations(alpha/K, obj.rhot, obj.taut);
#     obj.Etinvinv = obj.Etinv.^-1;
}
