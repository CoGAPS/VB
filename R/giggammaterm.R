giggammaterm <- function(Ex, Exinv, rho, tau, a, b) {
    score <- 0
    cutoff <- 1e-200
    zerotau <- which(tau <= cutoff)
    nonzerotau <- find(tau > cutoff)
    score <- score + length(Ex)*(a*log(b) - lgamma(a))
    score <- score - sum((b - rho) * Ex)

    score <- score - size(nonzerotau) * log(0.5)
    score <- score + sum(tau(nonzerotau) * Exinv[nonzerotau])
    score <- score - 0.5 * a * sum(log(rho[nonzerotau]) - log(tau[nonzerotau]))
    score <- score + sum(log(besselK(2*sqrt(rho[nonzerotau] * tau[nonzerotau]), a, TRUE)) - 2*sqrt(rho[nonzerotau]*tau[nonzerotau]))

    score <- score + sum(-a*log(rho[zerotau]) + lgamma(a))

    return(score)
}
