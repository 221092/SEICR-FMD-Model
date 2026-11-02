function th = seicrThresholds(p)
%SEICRTHRESHOLDS Deterministic and noise-corrected threshold quantities.
%   th.R0       deterministic R0 (period-mean transmission beta0), Subsection 5.8
%   th.Rs_mean  heuristic noise-corrected index (period mean), Subsection 5.9
%   th.Rs_peak  the same index with beta0 replaced by beta0(1+eps)
%   th.bound    rigorous extinction bound lambda* of Theorem 11.2:
%               s(Fhat - V) - 1/(2*sum_i 1/sigma_i^2), i = 2,3,4; < 0 => extinction
s = p.sigma(:);
th.R0      = indexR(p, p.beta0,           zeros(5,1));
th.Rs_mean = indexR(p, p.beta0,           s);
th.Rs_peak = indexR(p, p.beta0*(1+p.eps), s);
[F, V] = seicrNGM(p, p.beta0*(1+p.eps));
sInf = s(2:4);
if all(sInf > 0), noiseTerm = 1/(2*sum(1./sInf.^2)); else, noiseTerm = 0; end
th.bound = max(real(eig(F - V))) - noiseTerm;
end

function R = indexR(p, B, s)
a = p.mu + p.alpha + 0.5*s(2)^2;
b = p.mu + p.gamma + p.theta + 0.5*s(3)^2;
d = p.mu + p.delta + 0.5*s(4)^2;
R = p.alpha*B/(a*b) * (1 + p.kappa*p.p*p.gamma/d);
end
