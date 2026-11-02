function p = seicrParams(scenario)
%SEICRPARAMS Baseline parameters (rates per day, states in animals).
%   scenario: 'deterministic' | 'low-noise' | 'high-infected-noise'
p.Lambda = 8;            % recruitment (animals/day); Lambda/mu = 23,360 head
p.mu     = 1/(8*365);    % natural mortality, 8-year productive life
p.beta0  = 0.30;         % mean transmission rate
p.eps    = 0.40;         % seasonal amplitude (0 <= eps < 1)
p.phi    = 0;            % seasonal phase (radians)
p.T      = 365;          % forcing period (days)
p.kappa  = 0.08;         % relative infectiousness of carriers
p.alpha  = 1/4;          % 1/latent period
p.gamma  = 1/6;          % 1/acute infectious period
p.theta  = 0.002;        % disease-induced mortality
p.p      = 0.25;         % probability of entering the carrier state
p.delta  = 1/180;        % 1/carrier duration
switch lower(scenario)
    case 'deterministic',       p.sigma = zeros(1,5);
    case 'low-noise',           p.sigma = [0.05 0.05 0.05 0.05 0.05];
    case 'high-infected-noise', p.sigma = [0.05 0.55 0.60 0.45 0.05];
    otherwise, error('seicrParams:scenario','Unknown scenario "%s".',scenario);
end
end
