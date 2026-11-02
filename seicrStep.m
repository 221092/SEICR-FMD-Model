function Xn = seicrStep(method, t, X, dt, dW, p)
%SEICRSTEP One step of the chosen integrator for all M paths at once.
%   X, dW : 5-by-M. method: 'em' | 'milstein' | 'srk' | 'nsfd' | 'nsfd-naive'
%   'nsfd' is the positivity-preserving splitting scheme of Section 15.4:
%   an NSFD drift step followed by the exact Ito factor exp(sigma dW - sigma^2 dt/2).
%   'nsfd-naive' is the indicator-split NSFD scheme of Section 15.4, kept only
%   for comparison: it is consistent with the Stratonovich, not the Ito, equation.
sig = p.sigma(:);                    % 5-by-1, broadcasts over columns
switch lower(method)
    case 'em'
        Xn = X + seicrDrift(t,X,p)*dt + sig.*X.*dW;
    case 'milstein'
        Xn = X + seicrDrift(t,X,p)*dt + sig.*X.*dW + 0.5*(sig.^2).*X.*(dW.^2 - dt);
    case 'srk'   % derivative-free order-1.0 scheme (Kloeden-Platen), Section 15.3
        f  = seicrDrift(t,X,p);
        Xh = X + f*dt + sig.*X*sqrt(dt);                   % supporting value
        Xn = X + f*dt + sig.*X.*dW + (sig.*Xh - sig.*X).*(dW.^2 - dt)/(2*sqrt(dt));
    case {'nsfd','nsfd-naive'}
        [gain, loss, ph] = gainLoss(t, X, dt, p);
        if strcmpi(method,'nsfd')
            Y  = (X + ph*gain) ./ (1 + ph*loss);           % positive for any dt
            Xn = Y .* exp(sig.*dW - 0.5*(sig.^2)*dt);      % exact multiplicative factor
        else
            Xn = (X + ph*gain + sig.*X.*max(dW,0)) ./ ...
                 (1 + ph*loss - sig.*min(dW,0));            % indicator split
        end
    otherwise
        error('seicrStep:method','Unknown method "%s".',method);
end
end

function [gain, loss, ph] = gainLoss(t, X, dt, p)
%GAINLOSS Nonnegative gains g_i and per-capita losses l_i of the drift, f_i = g_i - l_i X_i.
mu = p.mu; ph = -expm1(-mu*dt)/mu;                         % denominator function psi(dt)
S = X(1,:); E = X(2,:); I = X(3,:); C = X(4,:);
N   = max(sum(X,1), realmin);
lam = seicrBeta(t,p) .* (I + p.kappa*C) ./ N;              % force of infection
o   = ones(size(S));
gain = [ p.Lambda*o; lam.*S; p.alpha*E; p.p*p.gamma*I; (1-p.p)*p.gamma*I + p.delta*C ];
loss = [ mu + lam; (mu+p.alpha)*o; (mu+p.gamma+p.theta)*o; (mu+p.delta)*o; mu*o ];
end
