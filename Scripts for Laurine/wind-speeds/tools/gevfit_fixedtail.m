%
% Parameter estimates and confidence intervals for GEV data with fixed shape.


function [parmhat,parmci] = gevfit_fixedtail(x,tail,alpha,options)

if ~isvector(x)
    error(message('stats:gevfit:VectorRequired'));
end

if nargin < 2 || isempty(tail)
    error('The shape parameter tail must be specified.');
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
else
    alpha = cast(alpha, "like", []);
end

if nargin < 4 || isempty(options)
    options = statset('gevfit');
else
    options = statset(statset('gevfit'), options);
end

classX = underlyingType(x);
if strcmp(classX, 'single')
    x = double(x);
end

n = length(x);
x = sort(x(:));
xmin = x(1);
xmax = x(end);
rangex = range(x);

% Can't make a fit.
if n == 0 || ~isfinite(rangex)
    parmhat = NaN(1, 2, "like", x);
    parmci = NaN(2, 2, "like", x);
    return
elseif rangex < realmin(classX)
    % When all observations are equal, try to return something reasonable.
    parmhat = [0, x(1)];
    if n == 1
        parmci = cast([-Inf -Inf; Inf Inf], "like", x);
    else
        parmci = [parmhat; parmhat];
    end
    return
end

% Initial guesses for sigma and mu.
evparms = evfit(x);
sigma0 = evparms(2);
mu0 = evparms(1);

parmhat0 = [log(sigma0), mu0]; % Initial guess for [lnsigma, mu]

% Maximize the log-likelihood with respect to lnsigma and mu, with k fixed.
[parmhat, ~, err, output] = fminsearch(@negloglike_fixed_shape, parmhat0, options, x, tail);
parmhat(1) = exp(parmhat(1)); % Transform log(sigma) back to sigma

if err == 0
    if output.funcCount >= options.MaxFunEvals
        warning(message('stats:gevfit:EvalLimit'));
    else
        warning(message('stats:gevfit:IterLimit'));
    end
elseif err < 0
    error(message('stats:gevfit:NoSolution'));
end

if nargout > 1
    probs = [alpha/2; 1-alpha/2];
    [~, acov] = gevlike([tail, parmhat], x);
    se = sqrt(diag(acov))';

    % Compute the CI for sigma using a normal approximation for
    % log(sigmahat), and transform back to the original scale.
    lnsigci = norminv(probs, log(parmhat(1)), se(1) ./ parmhat(1));

    % Compute the CI for mu using a normal distribution for muhat.
    muci = norminv(probs, parmhat(2), se(2));

    parmci = [exp(lnsigci), muci];
end

if strcmp(classX, 'single')
    parmhat = single(parmhat);
    if nargout > 1
        parmci = single(parmci);
    end
end

function nll = negloglike_fixed_shape(parms, data, k)
% Negative log-likelihood for the GEV with fixed shape parameter.
lnsigma = parms(1);
sigma = exp(lnsigma);
mu = parms(2);

n = numel(data);
z = (data - mu) ./ sigma;

if abs(k) > eps
    u = 1 + k .* z;
    if min(u) > 0
        lnu = log1p(k .* z); % log(1 + k .* z)
        t = exp(-(1 / k) * lnu); % (1 + k .* z).^(-1/k)
        nll = n * lnsigma + sum(t) + (1 + 1 / k) * sum(lnu);
    else
        nll = Inf; % The support of the GEV is 1+k*z > 0, or x > mu - sigma/k.
    end
else % limiting extreme value distribution as k -> 0
    nll = n * lnsigma + sum(exp(-z) + z);
end
