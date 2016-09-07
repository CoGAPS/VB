function [Ex Exinv] = computegigexpectations(alpha, beta, gamma)
% function [Ex Exinv] = computegigexpectations(alpha, beta, gamma)

if (numel(alpha) == 1)
    alpha = alpha * ones(size(beta));
end

Ex = zeros(size(beta));
Exinv = zeros(size(beta));

% For very small values of gamma and positive values of alpha, the GIG
% distribution becomes a gamma distribution, and its expectations are both
% cheaper and more stable to compute that way.
giginds = find(gamma(:) > 1e-200);
gaminds = find(gamma(:) <= 1e-200);
if (sum(alpha(gaminds) < 0) > 0)
    error('problem with arguments.');
end

% Compute expectations for GIG distribution.
sqrtbeta = sqrt(beta(giginds));
sqrtgamma = sqrt(gamma(giginds));
% besselparam = sqrtbeta .* sqrtgamma;
% % Note that we're using the *scaled* version here, since we're just
% % computing ratios and it's more stable.
besselalphaminus = besselk(alpha(giginds)-1, 2*sqrtbeta.*sqrtgamma, 1);
besselalpha = besselk(alpha(giginds), 2*sqrtbeta.*sqrtgamma, 1);
besselalphaplus = besselk(alpha(giginds)+1, 2*sqrtbeta.*sqrtgamma, 1);
sqrtratio = sqrtgamma ./ sqrtbeta;
Ex(giginds) = besselalphaplus .* sqrtratio ./ besselalpha;
Exinv(giginds) = besselalphaminus ./ (sqrtratio .* besselalpha);

% Compute expectations for gamma distribution where we can get away with 
% it.
Ex(gaminds) = alpha(gaminds) ./ beta(gaminds);
Exinv(gaminds) = beta(gaminds) ./ (alpha(gaminds) - 1);
Exinv(Exinv < 0) = inf;
