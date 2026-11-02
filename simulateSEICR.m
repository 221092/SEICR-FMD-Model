function out = simulateSEICR(p, x0, days, dt, M, method, seed, opts)
%SIMULATESEICR Monte Carlo ensemble of the stochastic SEICR model.
%   out.t        - stored times (days)
%   out.X        - 5-by-numel(t)-by-min(M,opts.keep) stored paths
%   out.Xend     - 5-by-M final states
%   out.tExt     - 1-by-M fade-out times (Inf if no fade-out)
%   out.negative - 1-by-M flag: path produced a negative component
%   opts.absorb  - fade-out level for E+I+C in animals (default 1; [] = off)
%   opts.every   - store every k-th step (default 10); opts.keep (default 60)
if nargin < 8, opts = struct; end
if ~isfield(opts,'absorb'), opts.absorb = 1;  end
if ~isfield(opts,'every'),  opts.every  = 10; end
if ~isfield(opts,'keep'),   opts.keep   = 60; end
rng(seed, 'twister');
n = round(days/dt); dt = days/n;
X = repmat(x0(:), 1, M);
tExt = inf(1,M); negative = false(1,M);
idx = 0:opts.every:n; K = min(M,opts.keep);
out.t = idx*dt; out.X = zeros(5,numel(idx),K); out.X(:,1,:) = X(:,1:K); s = 1;
for k = 1:n
    dW = sqrt(dt)*randn(5,M);
    X = seicrStep(method, (k-1)*dt, X, dt, dW, p);
    negative = negative | any(X < 0, 1);
    if ~isempty(opts.absorb)
        hit = sum(X(2:4,:),1) < opts.absorb & isinf(tExt);
        tExt(hit) = k*dt;
        X(2:4, ~isinf(tExt)) = 0;      % fade-out is absorbing: no fractional reseeding
    end
    if mod(k,opts.every) == 0
        s = s + 1; out.X(:,s,:) = X(:,1:K);
    end
end
out.Xend = X; out.tExt = tExt; out.negative = negative;
end
