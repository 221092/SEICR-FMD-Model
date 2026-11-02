function out = seicrPeriodicLNA(p)
%SEICRPERIODICLNA Periodic orbit x*(t) of the skeleton and periodic solution
%   Sigma(t) of the differential Lyapunov equation (Section 14).
f   = @(t,x) seicrDrift(t, x, p);
opt = odeset('RelTol',1e-9,'AbsTol',1e-9);
x0  = [p.Lambda/p.mu - 10; 5; 3; 2; 0];
[~, X] = ode15s(f, [0 300*p.T], x0, opt); xs = X(end,:)';    % relax onto attractor
for it = 1:3                                                % refine Poincare fixed point
    [~, X] = ode15s(f, [0 p.T], xs, opt); xs = X(end,:)';
end
z   = [xs; zeros(25,1)];
rhs = @(t,z) lnaRHS(t, z, p);
for per = 1:60
    [~, Z] = ode15s(rhs, [0 p.T], z, odeset('RelTol',1e-8,'AbsTol',1e-6)); z = Z(end,:)';
end
tt = linspace(0, p.T, 73);
[~, Z] = ode15s(rhs, tt, z, odeset('RelTol',1e-8,'AbsTol',1e-6));
out.t   = tt; out.x = Z(:,1:5)';
out.var = Z(:, 5 + [1 7 13 19 25])';                % diagonal of Sigma(t)
out.cv  = sqrt(out.var)./out.x;
mx = mean(out.x,2); mv = mean(out.var,2);
fprintf('Mean periodic orbit (S,E,I,C,R): %s\n', sprintf('%.1f ', mx));
share = out.x(4,:)./sum(out.x(2:4,:),1);                % instantaneous carrier share
fprintf('Carrier share of infected animals on the orbit: %.1f%% (time average), %.1f%% (ratio of means)\n', ...
        100*mean(share), 100*mx(4)/sum(mx(2:4)));
fprintf('Variance share, whole herd (S,E,I,C,R): %s(%%)\n', sprintf('%.1f ', 100*mv/sum(mv)));
fprintf('Variance share, infected (E,I,C): %s(%%)\n', sprintf('%.1f ', 100*mv(2:4)/sum(mv(2:4))));
fprintf('Mean LNA CV (S,E,I,C,R): %s\n', sprintf('%.2f ', mean(out.cv,2)));
end

function dz = lnaRHS(t, z, p)
x = z(1:5); P = reshape(z(6:end), 5, 5);
J = numJac(t, x, p); G = diag((p.sigma(:).*x).^2);
dz = [seicrDrift(t, x, p); reshape(J*P + P*J' + G, [], 1)];
end

function J = numJac(t, x, p)
J = zeros(5); f0 = seicrDrift(t, x, p);
for j = 1:5
    h = 1e-6*max(1, abs(x(j))); e = zeros(5,1); e(j) = h;
    J(:,j) = (seicrDrift(t, x + e, p) - f0)/h;
end
end
