%RUN_EXPERIMENTS Reproduces Experiments 1-7 of Section 19: every table and figure.
%   Figures are written to the folder 'fig' (fan.pdf, fadeout.pdf, herd.pdf,
%   IC_trajectory.png, map.pdf, conv.pdf, lna.pdf); table entries are printed.
clear; close all;
outDir = 'fig'; if ~exist(outDir,'dir'), mkdir(outDir); end
scen = {'low-noise','high-infected-noise'};

%% Experiment 1 - threshold quantities (Table: threshold quantities at baseline)
for sc = {'deterministic','low-noise','high-infected-noise'}
    p  = seicrParams(sc{1}); th = seicrThresholds(p);
    [lam, se] = seicrLyapunovExponent(p, 30, 0.02, 100, 1);
    fprintf('%-20s R0=%.3f Rs(mean)=%.3f Rs(peak)=%.3f bound=%+.4f lambda_s=%+.4f (se %.4f)\n', ...
        sc{1}, th.R0, th.Rs_mean, th.Rs_peak, th.bound, lam, se);
end

%% Experiment 2 - skeleton, ensembles and fade-out (fan.pdf, fadeout.pdf, herd.pdf)
days = 5*365; dt = 0.05; M = 1000;
p  = seicrParams('deterministic'); x0 = [p.Lambda/p.mu - 10, 5, 3, 2, 0];
[td, Xd] = ode15s(@(t,x) seicrDrift(t,x,p), linspace(0,days,5*365+1), x0(:), ...
                  odeset('RelTol',1e-9,'AbsTol',1e-9));
Xd = Xd'; infd = sum(Xd(2:4,:),1);
[~, ipk] = max(Xd(3,:)); [mn, imn] = min(infd(td > 365));  tAfter = td(td > 365);
fprintf('Skeleton: peak I = %.0f at day %.0f; post-epidemic minimum of E+I+C = %.1f at day %.0f\n', ...
    Xd(3,ipk), td(ipk), mn, tAfter(imn));
fprintf('Skeleton: total herd stays within %.2f%% of Lambda/mu\n', ...
    100*max(abs(sum(Xd,1)/(p.Lambda/p.mu) - 1)));

res = cell(1,2);
for j = 1:2
    q  = seicrParams(scen{j});
    res{j} = simulateSEICR(q, x0, days, dt, M, 'milstein', 11, struct('every',20,'keep',M));
    e  = res{j}.tExt; pk = squeeze(max(res{j}.X(3,:,:), [], 2));
    alt = simulateSEICR(q, x0, days, dt, M, 'nsfd', 12, struct('keep',1));   % independent noise
    fprintf(['%s: P(tau<=1y)=%.3f P(tau<=5y)=%.3f median tau=%.0f d median peak I=%.0f ', ...
             'peak I range [%.0f, %.0f]; splitting scheme P(tau<=5y)=%.3f; ', ...
             'Milstein paths with a negative component: %d\n'], ...
        scen{j}, mean(e <= 365), mean(isfinite(e)), median(e(isfinite(e))), median(pk), ...
        min(pk), max(pk), mean(isfinite(alt.tExt)), sum(res{j}.negative));
end

% Fan chart of I and C (fan.pdf)
names = 'SEICR'; titles = {'Low noise','High infected-class noise'};
figure('Position',[100 100 900 600]); tiledlayout(2,2,'TileSpacing','compact');
for k = [3 4]
    for j = 1:2
        nexttile; Y = max(squeeze(res{j}.X(k,:,:)), 1e-2);
        q = prctile(Y, [5 50 95], 2); t = res{j}.t/365;
        fill([t, fliplr(t)], [q(:,1)', fliplr(q(:,3)')], [0.8 0.88 1], 'EdgeColor','none'); hold on;
        plot(t, q(:,2), 'b-', 'LineWidth', 1.4);
        plot(td/365, max(Xd(k,:),1e-2), 'k--', 'LineWidth', 1.1);
        set(gca,'YScale','log'); ylim([1e-2 1e5]); yline(1,':'); grid on;
        xlabel('Time (years)'); ylabel(sprintf('%s (animals)', names(k)));
        title(sprintf('%s - %s', titles{j}, names(k)));
        if j == 1 && k == 3
            legend('5-95% band','ensemble median','deterministic skeleton','Location','southwest');
        end
    end
end
exportgraphics(gcf, fullfile(outDir,'fan.pdf'));

% Probability of fade-out by time t (fadeout.pdf)
figure; hold on; tt = linspace(0, days, 500);
for j = 1:2, plot(tt/365, mean(res{j}.tExt(:) <= tt, 1), 'LineWidth', 1.5); end
plot(tt/365, zeros(size(tt)), 'k--', 'LineWidth', 1.2);
xlabel('Time (years)'); ylabel('P(fade-out by t)'); ylim([0 1]); grid on;
legend('Low noise','High infected-class noise','Deterministic skeleton','Location','east');
exportgraphics(gcf, fullfile(outDir,'fadeout.pdf'));

% Total herd size, low noise (herd.pdf)
Nlow = squeeze(sum(res{1}.X,1));                          % numel(t)-by-M
qN = prctile(Nlow, [5 50 95], 2); t = res{1}.t/365;
fprintf('Low noise, year 5: median N = %.0f, 5-95%% = [%.0f, %.0f], mean N = %.0f\n', ...
    qN(end,2), qN(end,1), qN(end,3), mean(Nlow(end,:)));
figure; fill([t, fliplr(t)], [qN(:,1)', fliplr(qN(:,3)')], [0.8 0.88 1], 'EdgeColor','none'); hold on;
plot(t, qN(:,2), 'b-', 'LineWidth', 1.4);
plot(td/365, sum(Xd,1), 'k--', 'LineWidth', 1.2);
yline(p.Lambda/p.mu, ':', 'LineWidth', 1.2);
set(gca,'YScale','log'); grid on; xlabel('Time (years)'); ylabel('Total herd N(t)');
legend('5-95% band','median','deterministic N(t)','\Lambda/\mu','Location','southwest');
exportgraphics(gcf, fullfile(outDir,'herd.pdf'));

% Single path without the fade-out convention (IC_trajectory.png):
% the first seed whose infected classes fall below one animal and later regrow.
q = seicrParams('low-noise');
for seed = 1:500
    o = simulateSEICR(q, x0, days, 0.1, 1, 'milstein', seed, struct('absorb',[],'every',1,'keep',1));
    inf1 = squeeze(sum(o.X(2:4,:,1),1)); kLow = find(inf1 < 1, 1);
    found = ~isempty(kLow) && any(squeeze(o.X(3,kLow:end,1)) > 10);
    if found, break; end
end
if found
    fprintf('Resurgence path: seed %d, E+I+C first below one animal at day %.0f\n', seed, o.t(kLow));
else
    warning('No resurgence path found in seeds 1-500; plotting seed 500.');
end
figure; plot(o.t/365, squeeze(o.X(3,:,1)), 'LineWidth', 1); hold on;
plot(o.t/365, squeeze(o.X(4,:,1)), 'LineWidth', 1); grid on;
xlabel('Time (years)'); ylabel('Animals'); legend('I','C');
title('Acute infection and carrier dynamics');
exportgraphics(gcf, fullfile(outDir,'IC_trajectory.png'), 'Resolution', 300);

%% Experiment 3 - is the heuristic index a threshold? (Table: index against exact exponent)
p = seicrParams('high-infected-noise'); base = p.sigma(2:4);
for c = [0.4 0.6 0.8 1.0 1.2 1.4]
    q = p; q.sigma(2:4) = c*base; th = seicrThresholds(q);
    [lam, se] = seicrLyapunovExponent(q, 30, 0.02, 100, 1);
    fprintf('c=%.1f sigma=(%.2f,%.2f,%.2f) Rs=%.3f bound=%+.4f lambda_s=%+.4f (se %.4f)\n', ...
        c, q.sigma(2:4), th.Rs_mean, th.bound, lam, se);
end
% noise scale c at which the rigorous bound lambda* crosses zero:
% lambda*(c) = s(Fhat - V) - c^2 n0 with n0 = 1/(2 sum_i base_i^-2)
n0 = 1/(2*sum(1./base.^2)); th = seicrThresholds(p);
cstar = sqrt((th.bound + n0)/n0);
fprintf('Rigorous bound certifies extinction for c >= %.2f\n', cstar);

%% Experiment 4 - extinction map in the (beta0, sigma3) plane (map.pdf)
B = linspace(0.10,0.60,11); S3 = linspace(0,1.2,11);
Pext = zeros(numel(S3),numel(B)); LAM = Pext; RS = Pext; R0 = Pext;
for i = 1:numel(S3)
    for j = 1:numel(B)
        q = seicrParams('low-noise'); q.beta0 = B(j); q.sigma(3) = S3(i);
        th = seicrThresholds(q); RS(i,j) = th.Rs_mean; R0(i,j) = th.R0;
        LAM(i,j) = seicrLyapunovExponent(q, 10, 0.05, 60, 100*i+j);
        o = simulateSEICR(q, x0, 2*365, 0.1, 300, 'milstein', 100*i+j, struct('keep',1));
        Pext(i,j) = mean(isfinite(o.tExt));
    end
end
figure; imagesc(B, S3, Pext); set(gca,'YDir','normal'); caxis([0 1]); cb = colorbar;
ylabel(cb, 'P(fade-out within 2 years)'); hold on;
contour(B, S3, RS,  [1 1], 'k--', 'LineWidth', 1.2);
contour(B, S3, LAM, [0 0], 'k-',  'LineWidth', 1.5);
contour(B, S3, R0,  [1 1], 'k:',  'LineWidth', 1.2);
xlabel('\beta_0 (day^{-1})'); ylabel('\sigma_3 (day^{-1/2})');
legend('R_0^\sigma = 1 (heuristic index)','\lambda_s = 0 (exact)','R_0 = 1','Location','northeast');
exportgraphics(gcf, fullfile(outDir,'map.pdf'));

%% Experiment 5 - strong convergence and positivity (conv.pdf)
lbl = {'Euler-Maruyama','Milstein','SRK','NSFD, indicator split','NSFD, splitting'};
figure('Position',[100 100 900 380]); tiledlayout(1,2,'TileSpacing','compact');
for j = 1:2
    [dts, err] = seicrConvergence(scen{j});
    nexttile; loglog(dts, err, 'o-', 'LineWidth', 1.3); hold on;
    ref = err(2,end);
    loglog(dts, ref*(dts/dts(end)).^0.5, ':', 'Color', [0.5 0.5 0.5]);
    loglog(dts, ref*(dts/dts(end)),      '--', 'Color', [0.5 0.5 0.5]);
    grid on; xlabel('\Delta t (days)'); ylabel('E|ln I_K - ln I(30)|'); title(titles{j});
    if j == 2, legend([lbl, {'slope 1/2','slope 1'}], 'Location','southeast'); end
end
exportgraphics(gcf, fullfile(outDir,'conv.pdf'));
seicrScalarTest();
seicrPositivity('high-infected-noise', [0.05 0.25 0.5 1 2], 1000, 365);

%% Experiment 6 - periodic orbit and linear-noise covariance (lna.pdf)
out = seicrPeriodicLNA(seicrParams('low-noise'));
[~, ipk] = max(out.x(3,:));
fprintf('Acute prevalence on the orbit: %.0f to %.0f animals, peak near day %.0f\n', ...
    min(out.x(3,:)), max(out.x(3,:)), out.t(ipk));
figure('Position',[100 100 900 360]); tiledlayout(1,2,'TileSpacing','compact');
nexttile; plot(out.t, out.x(3,:), out.t, out.x(4,:), 'LineWidth', 1.4); grid on;
xlabel('Day of year'); ylabel('Animals'); legend('I: periodic mean x^*(t)','C: periodic mean x^*(t)');
nexttile; plot(out.t, out.cv', 'LineWidth', 1.3); hold on; yline(1,':'); grid on;
xlabel('Day of year'); ylabel('LNA coefficient of variation'); legend('S','E','I','C','R');
exportgraphics(gcf, fullfile(outDir,'lna.pdf'));

%% Experiment 7 - sensitivity indices (Table: normalised sensitivity indices)
seicrSensitivity();
