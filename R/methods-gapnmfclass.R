gapnmfclass <- function(x, alpha, a, b, K, smoothness=100) {
    M <- dim(x)[1]
    N <- dim(x)[2]

    obj <- new("gapnmfclass",
               x=x / mean(x),
               alpha=alpha,
               a=a,
               b=b,
               K=K,
               M=M,
               N=N)

#     obj.rhow = 10000*gamrnd(smoothness, 1/smoothness, M, K);
#     obj.tauw = 10000*gamrnd(smoothness, 1/smoothness, M, K);
#     obj.rhoh = 10000*gamrnd(smoothness, 1/smoothness, K, N);
#     obj.tauh = 10000*gamrnd(smoothness, 1/smoothness, K, N);
#     obj.rhot = K*10000*gamrnd(smoothness, 1/smoothness, K, 1);
#     obj.taut = 1/K*10000*gamrnd(smoothness, 1/smoothness, K, 1);
#     
#     [obj.Ew obj.Ewinv] = computegigexpectations(a, obj.rhow, obj.tauw);
#     obj.Ewinvinv = obj.Ewinv.^-1;
#     [obj.Eh obj.Ehinv] = computegigexpectations(b, obj.rhoh, obj.tauh);
#     obj.Ehinvinv = obj.Ehinv.^-1;
#     [obj.Et obj.Etinv] = computegigexpectations(alpha/K, obj.rhot, obj.taut);
#     obj.Etinvinv = obj.Etinv.^-1;
}
