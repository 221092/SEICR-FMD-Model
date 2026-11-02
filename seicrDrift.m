function f = seicrDrift(t, X, p)
%SEICRDRIFT Drift of the SEICR Ito system. X is 5-by-M (one column per path).
S = X(1,:); E = X(2,:); I = X(3,:); C = X(4,:); R = X(5,:);
N   = max(sum(X,1), realmin);
inc = seicrBeta(t,p) .* S .* (I + p.kappa*C) ./ N;
f = [ p.Lambda - inc - p.mu*S;
      inc - (p.mu + p.alpha)*E;
      p.alpha*E - (p.mu + p.gamma + p.theta)*I;
      p.p*p.gamma*I - (p.mu + p.delta)*C;
      (1 - p.p)*p.gamma*I + p.delta*C - p.mu*R ];
end
