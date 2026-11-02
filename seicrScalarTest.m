function seicrScalarTest()
%SEICRSCALARTEST Bias of the two NSFD constructions on dX = sigma X dW, X(0) = 1.
%   Exact Ito values: E X(T) = 1 and E ln X(T) = -sigma^2 T/2.
sigma = 0.5; dt = 0.1; T = 10; M = 2e5; n = round(T/dt);
rng(7,'twister');
Xa = ones(1,M); Xb = ones(1,M);
for k = 1:n
    dW = sqrt(dt)*randn(1,M);
    Xa = (Xa + sigma*Xa.*max(dW,0)) ./ (1 - sigma*min(dW,0));   % indicator split
    Xb = Xb .* exp(sigma*dW - 0.5*sigma^2*dt);                  % splitting scheme
end
fprintf('exact:          E X = %.3f, E ln X = %.3f\n', 1, -0.5*sigma^2*T);
fprintf('indicator split: E X = %.3f, E ln X = %.3f\n', mean(Xa), mean(log(Xa)));
fprintf('splitting:       E X = %.3f, E ln X = %.3f\n', mean(Xb), mean(log(Xb)));
end
