function [dts, err, order] = seicrConvergence(scenario)
%SEICRCONVERGENCE Strong-error test along common Brownian paths.
%   Error measure: E|ln I_K - ln I_ref(T)| at T = 30 days, reference = Milstein, dt = 2^-10.
%   Rows of err and entries of order: 'em','milstein','srk','nsfd-naive','nsfd'.
p = seicrParams(scenario); x0 = [p.Lambda/p.mu - 10, 5, 3, 2, 0]';
T = 30; M = 400; kref = 10; ks = 1:6;
rng(5,'twister'); dtf = 2^-kref; nf = round(T/dtf);
dWf = sqrt(dtf)*randn(5, M, nf);
X = repmat(x0,1,M);
for k = 1:nf, X = seicrStep('milstein', (k-1)*dtf, X, dtf, dWf(:,:,k), p); end
Iref = X(3,:);
methods = {'em','milstein','srk','nsfd-naive','nsfd'};
dts = 2.^-ks; err = zeros(numel(methods),numel(ks)); order = zeros(1,numel(methods));
for m = 1:numel(methods)
    for q = 1:numel(ks)
        r = 2^(kref-ks(q)); dt = dts(q);
        dW = reshape(sum(reshape(dWf, 5, M, r, nf/r), 3), 5, M, nf/r);  % aggregate fine increments
        X = repmat(x0,1,M);
        for k = 1:nf/r, X = seicrStep(methods{m}, (k-1)*dt, X, dt, dW(:,:,k), p); end
        err(m,q) = mean(abs(log(max(X(3,:),1e-12)) - log(Iref)));
    end
    c = polyfit(log(dts), log(err(m,:)), 1); order(m) = c(1);
    fprintf('%s %-10s observed strong order %.2f\n', scenario, methods{m}, order(m));
end
end
