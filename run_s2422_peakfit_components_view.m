function run_s2422_peakfit_components_view()
% run_s2422_peakfit_components_view
% Interactive plot viewer for component-wise peak fits on S2422 DRTs
% using residual-min lambda and L-curve lambda.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

data = dlmread(input_path, '\t');
freq_vec = data(:,1);
mag = data(:,2);
phase_deg = data(:,3);

Z_exp = mag .* exp(1i * phase_deg * pi / 180);
valid = freq_vec > 0;
freq_vec = freq_vec(valid);
Z_exp = Z_exp(valid);
[freq_vec, idx] = sort(freq_vec);
Z_exp = Z_exp(idx);
tau = 1 ./ (2 * pi * freq_vec);
xlog = log(tau);

% Method 1: minimum residual lambda
lambda_candidates = [1e-4, 1e-3, 1e-2, 1e-1, 1e0].';
nc = numel(lambda_candidates);
mean_resid = zeros(nc,1);
for i = 1:nc
    [g_i, r_i] = tr_drt_local(freq_vec, Z_exp, lambda_candidates(i));
    z_i = calc_eis_local(freq_vec, g_i, r_i);
    mean_resid(i) = mean(abs(Z_exp - z_i));
end
[~, idx_min] = min(mean_resid);
lambda_resid = lambda_candidates(idx_min);

% Method 2: L-curve corner lambda
lambda_scan = logspace(-6, 1, 40).';
ns = numel(lambda_scan);
res_norm = zeros(ns,1);
sol_norm = zeros(ns,1);
for i = 1:ns
    [g_i, r_i] = tr_drt_local(freq_vec, Z_exp, lambda_scan(i));
    z_i = calc_eis_local(freq_vec, g_i, r_i);
    dz = Z_exp - z_i;
    res_norm(i) = norm([real(dz); imag(dz)], 2);
    sol_norm(i) = norm(g_i, 2);
end
curv = lcurve_curvature_local(lambda_scan, res_norm, sol_norm);
[~, idx_corner] = max(curv);
lambda_lcurve = lambda_scan(idx_corner);

[gamma_resid, ~] = tr_drt_local(freq_vec, Z_exp, lambda_resid);
[gamma_lcurve, ~] = tr_drt_local(freq_vec, Z_exp, lambda_lcurve);

[idx_resid, hidden_resid] = detect_peaks_multiscale_hidden(gamma_resid);
[idx_lcurve, hidden_lcurve] = detect_peaks_multiscale_hidden(gamma_lcurve);

fit_resid = fit_peaks_gaussian_local(tau, gamma_resid, idx_resid, hidden_resid);
fit_lcurve = fit_peaks_gaussian_local(tau, gamma_lcurve, idx_lcurve, hidden_lcurve);

sum_resid = build_gaussian_sum(xlog, fit_resid);
sum_lcurve = build_gaussian_sum(xlog, fit_lcurve);
comps_resid = build_gaussian_components(xlog, fit_resid);
comps_lcurve = build_gaussian_components(xlog, fit_lcurve);

fprintf('Residual-min lambda: %.6e | peaks: %d\n', lambda_resid, numel(fit_resid));
fprintf('L-curve lambda     : %.6e | peaks: %d\n', lambda_lcurve, numel(fit_lcurve));

% Plot 1: residual-min method with components
figure('Color', 'w', 'Name', 'S2422 components - residual-min method');
semilogx(tau, gamma_resid, 'b-', 'LineWidth', 1.2); hold on;
semilogx(tau, sum_resid, 'k--', 'LineWidth', 1.2);
plot_components(tau, comps_resid, fit_resid);
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('Residual-min lambda = %.3e', lambda_resid));
grid on;
legend('DRT', 'Gaussian sum', 'Components', 'Location', 'best');

% Plot 2: L-curve method with components
figure('Color', 'w', 'Name', 'S2422 components - L-curve method');
semilogx(tau, gamma_lcurve, 'r-', 'LineWidth', 1.2); hold on;
semilogx(tau, sum_lcurve, 'k--', 'LineWidth', 1.2);
plot_components(tau, comps_lcurve, fit_lcurve);
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('L-curve lambda = %.3e', lambda_lcurve));
grid on;
legend('DRT', 'Gaussian sum', 'Components', 'Location', 'best');

% Plot 3: both DRTs in one figure
figure('Color', 'w', 'Name', 'S2422 two-method DRT overview');
semilogx(tau, gamma_resid, 'b-', 'LineWidth', 1.2); hold on;
semilogx(tau, gamma_lcurve, 'r-', 'LineWidth', 1.2);
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title('DRT comparison: residual-min vs L-curve lambda');
grid on;
legend(sprintf('Residual-min %.3e', lambda_resid), sprintf('L-curve %.3e', lambda_lcurve), 'Location', 'best');

h = msgbox('Peak-fit plots are open. Inspect them and click OK when done.','S2422 Peak Viewer','modal');
uiwait(h);
end

function plot_components(tau, comps, fit_struct)
if isempty(fit_struct)
    return;
end
cc = lines(size(comps,2));
for i = 1:size(comps,2)
    plot(tau, comps(:,i), '-', 'Color', cc(i,:), 'LineWidth', 1.0);
    if fit_struct(i).is_hidden > 0
        plot(fit_struct(i).tau, fit_struct(i).amplitude, 'mo', 'MarkerSize', 7, 'LineWidth', 1.2);
    else
        plot(fit_struct(i).tau, fit_struct(i).amplitude, 'go', 'MarkerSize', 7, 'LineWidth', 1.2);
    end
end
end

function comps = build_gaussian_components(xlog, fit_struct)
if isempty(fit_struct)
    comps = zeros(numel(xlog), 0);
    return;
end
comps = zeros(numel(xlog), numel(fit_struct));
for i = 1:numel(fit_struct)
    a = fit_struct(i).amplitude;
    m = log(fit_struct(i).tau);
    s = max(fit_struct(i).sigma, 1e-6);
    comps(:,i) = a * exp(-0.5 * ((xlog - m) / s).^2);
end
end

function ysum = build_gaussian_sum(xlog, fit_struct)
ysum = zeros(size(xlog));
for i = 1:numel(fit_struct)
    a = fit_struct(i).amplitude;
    m = log(fit_struct(i).tau);
    s = max(fit_struct(i).sigma, 1e-6);
    ysum = ysum + a * exp(-0.5 * ((xlog - m) / s).^2);
end
end

function [peak_idx, hidden_flag] = detect_peaks_multiscale_hidden(y)
y = y(:);
n = numel(y);
ys3 = smooth_mavg(y, 3);
ys7 = smooth_mavg(y, 7);
ys11 = smooth_mavg(y, 11);

base_thr = max(0.03 * max(y), eps);
[p3, ~] = find_peaks_simple_local(ys3, base_thr, 6);
[p7, ~] = find_peaks_simple_local(ys7, base_thr, 6);
[p11, ~] = find_peaks_simple_local(ys11, base_thr, 6);
base_peaks = unique([p3(:); p7(:); p11(:)]);

baseline = smooth_mavg(y, 25);
resid = y - baseline;
resid_thr = max(0.01 * max(y), 1.5 * std(resid));
[p_hidden, ~] = find_peaks_simple_local(resid, resid_thr, 4);

cand = unique([base_peaks(:); p_hidden(:)]);
if isempty(cand)
    peak_idx = [];
    hidden_flag = [];
    return;
end

cand = sort(cand(:));
merged = cand(1);
for i = 2:numel(cand)
    if abs(cand(i) - merged(end)) <= 3
        if y(cand(i)) > y(merged(end))
            merged(end) = cand(i);
        end
    else
        merged(end+1,1) = cand(i); %#ok<AGROW>
    end
end

hidden_flag = zeros(size(merged));
for i = 1:numel(merged)
    ii = merged(i);
    near_hidden = any(abs(p_hidden - ii) <= 2);
    near_base = any(abs(base_peaks - ii) <= 2);
    if near_hidden && ~near_base
        hidden_flag(i) = 1;
    end
end

peak_idx = merged;
valid = peak_idx > 2 & peak_idx < n-1;
peak_idx = peak_idx(valid);
hidden_flag = hidden_flag(valid);
end

function ysm = smooth_mavg(y, w)
y = y(:);
n = numel(y);
if w <= 1
    ysm = y;
    return;
end
half = floor(w/2);
ysm = zeros(n,1);
for i = 1:n
    i1 = max(1, i-half);
    i2 = min(n, i+half);
    ysm(i) = mean(y(i1:i2));
end
end

function [idx, vals] = find_peaks_simple_local(y, min_height, min_distance)
y = y(:);
n = numel(y);

cand = [];
for i = 2:n-1
    if y(i) >= y(i-1) && y(i) > y(i+1) && y(i) >= min_height
        cand(end+1,1) = i; %#ok<AGROW>
    end
end

if isempty(cand)
    idx = [];
    vals = [];
    return;
end

v = y(cand);
[~, ord] = sort(v, 'descend');
sel = [];
for k = 1:numel(ord)
    ii = cand(ord(k));
    if isempty(sel) || all(abs(sel - ii) >= min_distance)
        sel(end+1,1) = ii; %#ok<AGROW>
    end
end

sel = sort(sel);
idx = sel;
vals = y(sel);
end

function fit_struct = fit_peaks_gaussian_local(tau, gamma, peak_idx, hidden_flag)
if isempty(peak_idx)
    fit_struct = [];
    return;
end

tau = tau(:);
gamma = gamma(:);
xall = log(tau);
n = numel(gamma);

fit_struct = repmat(struct('peak_id',0,'tau',0,'amplitude',0,'sigma',0,'is_hidden',0), numel(peak_idx), 1);

for i = 1:numel(peak_idx)
    c = peak_idx(i);
    w = max(8, round(n/40));
    i1 = max(1, c-w);
    i2 = min(n, c+w);

    x = xall(i1:i2);
    y = gamma(i1:i2);

    a0 = max(gamma(c), eps);
    m0 = xall(c);
    s0 = max(0.15, 0.20 * (max(x)-min(x) + eps));
    b0 = max(min(y), 0);

    p0 = [log(a0); m0; log(s0); b0];
    obj = @(p) sum((exp(p(1)) * exp(-0.5*((x-p(2))/exp(p(3))).^2) + p(4) - y).^2);
    opts = optimset('Display','off','MaxFunEvals',4000,'MaxIter',2000);
    p = fminsearch(obj, p0, opts);

    fit_struct(i).peak_id = i;
    fit_struct(i).tau = exp(p(2));
    fit_struct(i).amplitude = max(exp(p(1)), 0);
    fit_struct(i).sigma = max(exp(p(3)), 1e-6);
    fit_struct(i).is_hidden = hidden_flag(i);
end

all_tau = [fit_struct.tau].';
[~, ord] = sort(all_tau, 'descend');
fit_struct = fit_struct(ord);
for i = 1:numel(fit_struct)
    fit_struct(i).peak_id = i;
end
end

function curv = lcurve_curvature_local(lambda_vals, res_norm, sol_norm)
n = numel(lambda_vals);
curv = zeros(n,1);
x = log10(max(res_norm, eps));
y = log10(max(sol_norm, eps));
t = log10(lambda_vals(:));

for i = 2:n-1
    dt1 = t(i) - t(i-1);
    dt2 = t(i+1) - t(i);
    dt = t(i+1) - t(i-1);
    if dt1 <= 0 || dt2 <= 0 || dt <= 0
        continue;
    end

    dx = (x(i+1) - x(i-1)) / dt;
    dy = (y(i+1) - y(i-1)) / dt;

    d2x = 2 * ((x(i+1)-x(i))/dt2 - (x(i)-x(i-1))/dt1) / dt;
    d2y = 2 * ((y(i+1)-y(i))/dt2 - (y(i)-y(i-1))/dt1) / dt;

    den = (dx*dx + dy*dy)^(3/2);
    if den > 0
        curv(i) = abs(dx * d2y - dy * d2x) / den;
    end
end
end

function [gamma, R_inf] = tr_drt_local(freq_vec, Z_exp, lam)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);

Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(lam/2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];

x = lsqnonneg(M, b);
gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calc_eis_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function A_re = calc_A_re_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);
A_re = zeros(n, n);
for p = 1:n
    for q = 1:n
        if q == 1
            log_term = log(tau(q+1) / tau(q));
        elseif q == n
            log_term = log(tau(q) / tau(q-1));
        else
            log_term = log(tau(q+1) / tau(q-1));
        end
        A_re(p, q) = -0.5 / (1 + (omega(p) * tau(q))^2) * log_term;
    end
end
end

function A_im = calc_A_im_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);
A_im = zeros(n, n);
for p = 1:n
    for q = 1:n
        if q == 1
            log_term = log(tau(q+1) / tau(q));
        elseif q == n
            log_term = log(tau(q) / tau(q-1));
        else
            log_term = log(tau(q+1) / tau(q-1));
        end
        A_im(p, q) = 0.5 * (omega(p) * tau(q)) / (1 + (omega(p) * tau(q))^2) * log_term;
    end
end
end
