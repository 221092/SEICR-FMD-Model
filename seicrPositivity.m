function frac = seicrPositivity(scenario, dts, M, days)
%SEICRPOSITIVITY Fraction of paths that produce a negative component, per scheme and step.
%   frac(m,q): method m in {'em','milstein','nsfd'}, step dts(q); horizon 'days'.
p = seicrParams(scenario); x0 = [p.Lambda/p.mu - 10, 5, 3, 2, 0];
methods = {'em','milstein','nsfd'}; frac = zeros(numel(methods), numel(dts));
for m = 1:numel(methods)
    for q = 1:numel(dts)
        o = simulateSEICR(p, x0, days, dts(q), M, methods{m}, 21, ...
                          struct('absorb',[],'every',1e9,'keep',1));
        frac(m,q) = mean(o.negative);
        fprintf('%s %-8s dt=%-5g negative in %.1f%% of paths\n', ...
                scenario, methods{m}, dts(q), 100*frac(m,q));
    end
end
end
