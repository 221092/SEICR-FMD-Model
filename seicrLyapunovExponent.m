function [lam, se] = seicrLyapunovExponent(p, years, dt, M, seed)
%SEICRLYAPUNOVEXPONENT Almost-sure growth exponent of the linearised infected
%   subsystem dY = (F(t)-V) Y dt + diag(sigma_2..4) Y dW at S/N = 1.
%   lam < 0 <=> invasion fails almost surely (the exact stochastic threshold).
rng(seed,'twister');
s = p.sigma(2:4).'; [~, V] = seicrNGM(p, 0);
Y = ones(3,M)/3; logn = zeros(1,M); n = round(years*365/dt);
for k = 1:n
    t = (k-1)*dt; [F,~] = seicrNGM(p, seicrBeta(t,p));
    dW = sqrt(dt)*randn(3,M);
    Y = (Y + (F - V)*Y*dt) .* exp(s.*dW - 0.5*s.^2*dt);
    if mod(k,50) == 0
        nr = sum(Y,1); logn = logn + log(nr); Y = Y./nr;
    end
end
logn = logn + log(sum(Y,1));
g   = logn/(n*dt);
lam = mean(g); se = std(g)/sqrt(M);
end
