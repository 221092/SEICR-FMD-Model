function [names, U] = seicrSensitivity()
%SEICRSENSITIVITY Normalised forward sensitivity indices d ln R / d ln q at baseline.
%   Columns of U: R0, heuristic index (low noise), heuristic index (high noise).
%   Central differences in ln q; noise intensities have no effect on R0 (NaN).
names = {'beta0','gamma','kappa','p','delta','alpha','mu','theta','sigma2','sigma3','sigma4'};
P = {seicrParams('deterministic'), seicrParams('low-noise'), seicrParams('high-infected-noise')};
h = 1e-6; U = nan(numel(names), 3);
for j = 1:3
    for i = 1:numel(names)
        if startsWith(names{i},'sigma') && j == 1, continue; end
        up = P{j}; dn = P{j};
        [up, dn] = perturb(up, dn, names{i}, h);
        U(i,j) = (log(Rof(up, j)) - log(Rof(dn, j))) / (2*h);
    end
end
fprintf('%-8s %9s %9s %9s\n', 'q', 'R0', 'Rs low', 'Rs high');
for i = 1:numel(names)
    fprintf('%-8s %9.3f %9.3f %9.3f\n', names{i}, U(i,:));
end
end

function [up, dn] = perturb(up, dn, name, h)
if startsWith(name,'sigma')
    k = str2double(name(end));
    up.sigma(k) = up.sigma(k)*exp(h); dn.sigma(k) = dn.sigma(k)*exp(-h);
else
    up.(name) = up.(name)*exp(h); dn.(name) = dn.(name)*exp(-h);
end
end

function R = Rof(p, j)
th = seicrThresholds(p);
if j == 1, R = th.R0; else, R = th.Rs_mean; end
end
