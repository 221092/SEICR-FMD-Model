function [F, V] = seicrNGM(p, B)
%SEICRNGM New-infection and transition matrices of the (E,I,C) subsystem at S/N = 1.
F = [0, B, p.kappa*B; 0 0 0; 0 0 0];
V = [ p.mu+p.alpha, 0, 0;
     -p.alpha, p.mu+p.gamma+p.theta, 0;
      0, -p.p*p.gamma, p.mu+p.delta ];
end
