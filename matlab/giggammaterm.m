function score = giggammaterm(Ex, Exinv, rho, tau, a, b)
% function score = giggammaterm(Ex, Exinv, rho, tau, a, b)

score = 0;
cutoff = 1e-200;
zerotau = find(tau(:) <= cutoff);
nonzerotau = find(tau(:) > cutoff);
score = score + numel(Ex)*(a*log(b) - gammaln(a));
score = score - sum((b - rho(:)) .* Ex(:));

score = score - numel(nonzerotau) * log(0.5);
score = score + sum(tau(nonzerotau) .* Exinv(nonzerotau));
score = score - 0.5 * a * sum(log(rho(nonzerotau)) - log(tau(nonzerotau)));
% It's numerically safer to use scaled version of besselk.
score = score + sum(log(besselk(a, 2*sqrt(rho(nonzerotau) .* tau(nonzerotau)), 1)) - 2*sqrt(rho(nonzerotau).*tau(nonzerotau)));

score = score + sum(-a*log(rho(zerotau)) + gammaln(a));