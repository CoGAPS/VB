computegigexpectations <- function(alpha, beta, gamma) {
    if (length(alpha) == 1) {
        alpha = alpha * rep(1, length(beta))
    }

    Ex <- rep(0, length(beta))
    Exinv <- rep(0, length(beta))

    giginds <- which(gamma > 1e-200)
    gaminds <- which(gamma > 1e-200)

    if (sum(alpha(gaminds) < 0) > 0) {
        stop("problem with arguments.")
    }

    sqrtbeta <- sqrt(beta[giginds])
    sqrtgamma <- sqrt(gamma[giginds])

    besselalphaminus <- besselK(2 * sqrtbeta * sqrtgamma, alpha[giginds]-1, TRUE)
    besselalpha <- besselK(2 * sqrtbeta * sqrtgamma, alpha[giginds], TRUE)
    besselalphaplus <- besselK(2 * sqrtbeta * sqrtgamma, alpha[giginds]+1, TRUE)

    sqrtratio <- sqrtgamma / sqrtbeta
    Ex[giginds] <- besselalphaplus * sqrtratio / besselalpha
    Exinv[giginds] <- besselalphaminus / (sqrtratio * besselalpha)

    Ex[gaminds] <- alpha[gaminds] / beta[gaminds]
    Exinv[gaminds] <- beta[gaminds] / (alpha[gaminds] - 1)
    Exinv[Exinv < 0] <- Inf

    return(list(Ex=Ex, Exinv=Exinv))
}
