function b = seicrBeta(t, p)
%SEICRBETA Periodic transmission coefficient beta(t) = beta0(1+eps cos(wt+phi)).
b = p.beta0 * (1 + p.eps * cos(2*pi*t/p.T + p.phi));
end
