function run_s2022_random10_paper_peak_refine()
% run_s2022_random10_paper_peak_refine
% MATLAB 2011-compatible workflow for random 10 S2022 temperatures.
% Uses both lambda-selection styles for context, then refines peaks on the
% paper-style DRT and searches for shoulder features.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_png = fullfile(repo_root, 's2022_random10_paper_peak_refine.png');
out_csv = fullfile(repo_root, 's2022_random10_paper_peak_refine.csv');
out_txt = fullfile(repo_root, 's2022_random10_paper_peak_refine_summary.txt');

rng_seed = 20260710;
n_select = 10;
lambda_values = logspace(-4, -1, 11).';  % Range [1e-4, 1e-1] constrained
paper_lambda_override = 0;

% Shoulder classifier tuning (sensitive setting).
shoulder_amp_ratio = 0.75;
shoulder_overlap_factor = 2.8;
shoulder_curvature_min = 0.12;

% Peak threshold tuning for robustness across lambda choices.
peak_seed_rel_height = 0.03;
peak_component_min_ratio = 0.10;
feature_visibility_min_ratio = 0.03;

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text] = xlsread(input_path);
temperature = parse_temperature_labels_local(text, size(data, 2));

rand('twister', rng_seed);
randn('state', rng_seed + 1);
selected_idx = randperm(numel(temperature));
selected_idx = selected_idx(1:n_select);
[selected_temperature, order] = sort(temperature(selected_idx));
selected_idx = selected_idx(order);

fig = figure('Color', 'w', 'Name', 'S2022 random 10 paper-method DRT peak refinement');

fid_csv = fopen(out_csv, 'w');
if fid_csv < 0
    error('Failed to open CSV file for writing: %s', out_csv);
end
fprintf(fid_csv, 'TemperatureK,PaperLambda,ResidualLambda,FeatureType,FeatureIndex,Tau_s,Freq_Hz,Gamma,PeakRank\n');

fid_txt = fopen(out_txt, 'w');
if fid_txt < 0
    fclose(fid_csv);
    error('Failed to open summary file for writing: %s', out_txt);
end

fprintf(fid_txt, 'S2022 random 10 paper-method DRT peak refinement\n');
fprintf(fid_txt, 'Generated: %s\n', datestr(now));
fprintf(fid_txt, 'Random seed: %d\n\n', rng_seed);

for t = 1:n_select
    src_cols = (3 * selected_idx(t) - 2):(3 * selected_idx(t));
    sel = data(:, src_cols);

    freq_vec = sel(:,1);
    mag = sel(:,2);
    phase_deg = sel(:,3);
    valid = isfinite(freq_vec) & isfinite(mag) & isfinite(phase_deg) & (freq_vec > 0);
    freq_vec = freq_vec(valid);
    mag = mag(valid);
    phase_deg = phase_deg(valid);

    Z_exp = mag .* exp(1i * phase_deg * pi / 180);
    [freq_vec, sort_order] = sort(freq_vec);
    Z_exp = Z_exp(sort_order);

    [idx_residual_best, idx_paper_best] = choose_lambda_indices_local(freq_vec, Z_exp, lambda_values);

    lambda_res = lambda_values(idx_residual_best);
    lambda_paper = lambda_values(idx_paper_best);
    if paper_lambda_override > 0
        lambda_paper = paper_lambda_override;
    end

    [gamma_paper, ~] = TR_DRT_local(freq_vec, Z_exp, lambda_paper);

    tau = 1 ./ (2 * pi * freq_vec);
    [peak_idx, shoulder_idx, peak_rank] = refine_peaks_and_shoulders_local( ...
        tau, gamma_paper, shoulder_amp_ratio, shoulder_overlap_factor, shoulder_curvature_min, ...
        peak_seed_rel_height, peak_component_min_ratio, feature_visibility_min_ratio);

    subplot(5, 2, t);
    h_paper = semilogx(tau, gamma_paper, 'b-', 'LineWidth', 1.1); hold on;
    h_peak = [];
    h_shoulder = [];

    if ~isempty(peak_idx)
        h_peak = semilogx(tau(peak_idx), gamma_paper(peak_idx), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'c');
    end
    if ~isempty(shoulder_idx)
        h_shoulder = semilogx(tau(shoulder_idx), gamma_paper(shoulder_idx), 'mv', 'MarkerSize', 5, 'MarkerFaceColor', 'm');
    end

    set(gca, 'XDir', 'reverse');
    grid on;
    title(sprintf('T=%.2f K | paper %.1e', selected_temperature(t), lambda_paper), 'FontSize', 8);
    if mod(t, 2) == 1
        ylabel('gamma(tau)');
    end
    if t > 8
        xlabel('tau (s)');
    end
    if t == 1
        hlist = h_paper;
        llist = {'Paper-style DRT'};
        if ~isempty(h_peak)
            hlist(end+1) = h_peak(1); %#ok<AGROW>
            llist{end+1} = 'Refined peak'; %#ok<AGROW>
        end
        if ~isempty(h_shoulder)
            hlist(end+1) = h_shoulder(1); %#ok<AGROW>
            llist{end+1} = 'Shoulder'; %#ok<AGROW>
        end
        legend(hlist, llist, 'Location', 'Best');
    end

    fprintf(fid_txt, 'Temperature %.6f K\n', selected_temperature(t));
    fprintf(fid_txt, '  Residual lambda: %.8e\n', lambda_res);
    fprintf(fid_txt, '  Paper lambda   : %.8e\n', lambda_paper);
    fprintf(fid_txt, '  Refined peaks  : %d\n', numel(peak_idx));
    fprintf(fid_txt, '  Shoulders      : %d\n', numel(shoulder_idx));

    for k = 1:numel(peak_idx)
        idx = peak_idx(k);
        f_hz = 1 / (2 * pi * tau(idx));
        fprintf(fid_csv, '%.6f,%.8e,%.8e,Peak,%d,%.8e,%.8f,%.8f,%d\n', ...
            selected_temperature(t), lambda_paper, lambda_res, idx, tau(idx), f_hz, gamma_paper(idx), peak_rank(k));
        fprintf(fid_txt, '    Peak %d: tau=%.4e s, gamma=%.4f\n', k, tau(idx), gamma_paper(idx));
    end

    for k = 1:numel(shoulder_idx)
        idx = shoulder_idx(k);
        f_hz = 1 / (2 * pi * tau(idx));
        fprintf(fid_csv, '%.6f,%.8e,%.8e,Shoulder,%d,%.8e,%.8f,%.8f,0\n', ...
            selected_temperature(t), lambda_paper, lambda_res, idx, tau(idx), f_hz, gamma_paper(idx));
        fprintf(fid_txt, '    Shoulder %d: tau=%.4e s, gamma=%.4f\n', k, tau(idx), gamma_paper(idx));
    end
    fprintf(fid_txt, '\n');
end

try
    saveas(fig, out_png);
catch
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, '-dpng', '-r200', out_png);
end
close(fig);

fclose(fid_csv);
fclose(fid_txt);

fprintf('Wrote plot   : %s\n', out_png);
fprintf('Wrote CSV    : %s\n', out_csv);
fprintf('Wrote summary: %s\n', out_txt);
end

function [idx_residual_best, idx_paper_best] = choose_lambda_indices_local(freq_vec, Z_exp, lambda_values)
n_lambda = numel(lambda_values);
mean_abs_residual = zeros(n_lambda, 1);
rid = zeros(n_lambda, 1);
cv_total = zeros(n_lambda, 1);
var_metric = zeros(n_lambda, 1);

for i = 1:n_lambda
    lam = lambda_values(i);
    [gamma_i, R_inf_i] = TR_DRT_local(freq_vec, Z_exp, lam);
    Z_fit = calculate_EIS_local(freq_vec, gamma_i, R_inf_i);
    mean_abs_residual(i) = mean(abs(Z_exp - Z_fit));

    [g_re, R_re] = TR_DRT_component_local(freq_vec, Z_exp, lam, 'real');
    [g_im, ~] = TR_DRT_component_local(freq_vec, Z_exp, lam, 'imag');
    Z_from_re = calculate_EIS_local(freq_vec, g_re, R_re);
    R_from_im = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, g_im);
    Z_from_im = calculate_EIS_local(freq_vec, g_im, R_from_im);

    rid(i) = mean((g_re - g_im).^2);
    cv_total(i) = mean((imag(Z_exp) - imag(Z_from_re)).^2) + mean((real(Z_exp) - real(Z_from_im)).^2);
    var_metric(i) = estimate_resample_variance_local(freq_vec, Z_exp, lam, 8);
end

[~, idx_residual_best] = min(mean_abs_residual);
idx_paper_best = select_lambda_idx_local(rid, cv_total, var_metric);
end

function [peak_idx, shoulder_idx, peak_rank] = refine_peaks_and_shoulders_local(x, y, shoulder_amp_ratio, shoulder_overlap_factor, shoulder_curvature_min, peak_seed_rel_height, peak_component_min_ratio, feature_visibility_min_ratio)
% Peak fitting pipeline:
% 1) Local maxima search on smoothed signal
% 2) Second-derivative amplification to expose shoulders
% 3) Gaussian mixture fit via Levenberg-Marquardt
% 4) Shoulder labeling from fitted component overlap and relative amplitude

y = y(:);
x = x(:);
n = numel(y);

if n < 7
    peak_idx = [];
    shoulder_idx = [];
    peak_rank = [];
    return;
end

[x_sorted, ord] = sort(x);
y_sorted = y(ord);
x_log = log(x_sorted);

w = max(5, 2 * floor(0.02 * n) + 1);
y_sm = moving_mean_local(y_sorted, w);

d1 = gradient_local(y_sm, x_log);
d2 = gradient_local(d1, x_log);
curv = -d2;
curv(curv < 0) = 0;
if max(curv) > 0
    curv = curv / max(curv);
end

% Amplify shoulders using curvature channel before maxima detection.
amp_signal = y_sm + 0.35 * max(y_sm) * curv;

[cand_idx, cand_rank] = find_maxima_rank_local(x_log, amp_signal);
if isempty(cand_idx)
    [~, imax] = max(amp_signal);
    cand_idx = imax;
    cand_rank = 1;
end

keep = amp_signal(cand_idx) >= peak_seed_rel_height * max(amp_signal);
cand_idx = cand_idx(keep);
cand_rank = cand_rank(keep);
if isempty(cand_idx)
    [~, imax] = max(amp_signal);
    cand_idx = imax;
    cand_rank = 1;
end

% Rank-select peaks while enforcing minimum spacing.
[~, ord_rank] = sort(cand_rank, 'descend');
min_dist = max(4, round(0.025 * n));
max_peaks = min(8, max(3, floor(n / 24)));
seed_idx = [];
for k = 1:numel(ord_rank)
    idx = cand_idx(ord_rank(k));
    if isempty(seed_idx) || all(abs(seed_idx - idx) >= min_dist)
        seed_idx(end+1,1) = idx; %#ok<AGROW>
        if numel(seed_idx) >= max_peaks
            break;
        end
    end
end
seed_idx = sort(seed_idx);
m = numel(seed_idx);

baseline0 = max(min(y_sm), eps);
A0 = max(y_sm(seed_idx) - baseline0, 0.02 * max(y_sm));
mu0 = x_log(seed_idx);
sigma0 = zeros(m,1);
for k = 1:m
    if m == 1
        span = max(x_log) - min(x_log);
    elseif k == 1
        span = abs(mu0(2) - mu0(1));
    elseif k == m
        span = abs(mu0(m) - mu0(m-1));
    else
        span = min(abs(mu0(k) - mu0(k-1)), abs(mu0(k+1) - mu0(k)));
    end
    sigma0(k) = min(max(0.05, 0.40 * span), 0.80);
end

theta0 = [log(A0(:)); mu0(:); log(sigma0(:)); log(baseline0)];
theta = fit_gaussian_lm_local(theta0, x_log, y_sorted, m);
[y_fit, A, mu, sigma] = gaussian_sum_from_theta_local(theta, x_log, m);

[mu, ord_mu] = sort(mu);
A = A(ord_mu);
sigma = sigma(ord_mu);

if isempty(A)
    peak_idx = [];
    shoulder_idx = [];
    peak_rank = [];
    return;
end

maxA = max(A);
curv_at_mu = zeros(m, 1);
for k = 1:m
    [~, i_near_curv] = min(abs(x_log - mu(k)));
    curv_at_mu(k) = curv(i_near_curv);
end

% Drop very weak components unless curvature strongly supports them.
valid_comp = (A >= peak_component_min_ratio * maxA) | (curv_at_mu >= shoulder_curvature_min);
if any(valid_comp)
    A = A(valid_comp);
    mu = mu(valid_comp);
    sigma = sigma(valid_comp);
    curv_at_mu = curv_at_mu(valid_comp);
    m = numel(A);
    maxA = max(A);
end

is_shoulder = false(m,1);
for k = 1:m
    if A(k) < shoulder_amp_ratio * maxA
        curv_here = curv_at_mu(k);
        for j = 1:m
            if j == k
                continue;
            end
            close_to_strong = abs(mu(k) - mu(j)) <= shoulder_overlap_factor * max(sigma(k), sigma(j));
            stronger = A(j) > A(k);
            curvature_supported = curv_here >= shoulder_curvature_min;
            if close_to_strong && stronger && curvature_supported
                is_shoulder(k) = true;
                break;
            end
        end
    end
end

if all(is_shoulder)
    [~, imax] = max(A);
    is_shoulder(imax) = false;
end

peak_idx = zeros(sum(~is_shoulder),1);
shoulder_idx = zeros(sum(is_shoulder),1);

ip = 0;
is = 0;
for k = 1:m
    tau_k = exp(mu(k));
    [~, idx_orig] = min(abs(x - tau_k));
    if is_shoulder(k)
        is = is + 1;
        shoulder_idx(is) = idx_orig;
    else
        ip = ip + 1;
        peak_idx(ip) = idx_orig;
    end
end

peak_idx = unique(peak_idx);
shoulder_idx = setdiff(unique(shoulder_idx), peak_idx);

% Final visibility filter: keep only features visible on the plotted DRT.
ymax = max(y);
if ~isempty(peak_idx)
    keep_peak = y(peak_idx) >= feature_visibility_min_ratio * ymax;
    if any(keep_peak)
        peak_idx = peak_idx(keep_peak);
    end
end
if ~isempty(shoulder_idx)
    keep_shoulder = y(shoulder_idx) >= feature_visibility_min_ratio * ymax;
    shoulder_idx = shoulder_idx(keep_shoulder);
end

if isempty(peak_idx)
    [~, idx_best] = max(y_fit);
    peak_idx = ord(idx_best);
end

% Peak rank as integer order by descending gamma value at detected peak positions.
peak_rank = zeros(numel(peak_idx), 1);
if ~isempty(peak_idx)
    [~, ord_pk] = sort(y(peak_idx), 'descend');
    for k = 1:numel(ord_pk)
        peak_rank(ord_pk(k)) = k;
    end
end
end

function theta = fit_gaussian_lm_local(theta0, x, y, m)
if exist('lsqnonlin', 'file') == 2
    opts = optimset('Display', 'off', 'Algorithm', 'levenberg-marquardt', ...
        'MaxFunEvals', 6000, 'MaxIter', 2500, 'TolFun', 1e-8, 'TolX', 1e-8);
    obj = @(th) gaussian_sum_from_theta_local(th, x, m) - y;
    theta = lsqnonlin(obj, theta0, [], [], opts);
else
    % Fallback when Optimization Toolbox is unavailable.
    obj = @(th) sum((gaussian_sum_from_theta_local(th, x, m) - y).^2);
    opts = optimset('Display', 'off', 'MaxFunEvals', 8000, 'MaxIter', 4000);
    theta = fminsearch(obj, theta0, opts);
end
end

function [yhat, A, mu, sigma, b0] = gaussian_sum_from_theta_local(theta, x, m)
A = exp(theta(1:m));
mu = theta(m+1:2*m);
sigma = exp(theta(2*m+1:3*m));
b0 = exp(theta(end));

Phi = zeros(numel(x), m);
for k = 1:m
    Phi(:, k) = exp(-0.5 * ((x - mu(k)) ./ sigma(k)).^2);
end
yhat = Phi * A + b0;
end

function [idx, rank] = find_maxima_rank_local(x, y)
N = numel(x);
if N < 4
    idx = [];
    rank = [];
    return;
end

d = diff(y) ./ diff(x);
idx = [];
rank = [];
for i = 1:(N - 2)
    if d(i) > 0 && d(i + 1) < 0
        idx(end + 1, 1) = i + 1; %#ok<AGROW>
        rank(end + 1, 1) = mean(abs(d(i:i + 1))); %#ok<AGROW>
    end
end
end

function ysm = moving_mean_local(y, w)
n = numel(y);
w = max(1, round(w));
if w > n
    w = n;
end
half = floor(w / 2);
ysm = zeros(size(y));
for i = 1:n
    i1 = max(1, i - half);
    i2 = min(n, i + half);
    ysm(i) = mean(y(i1:i2));
end
end

function g = gradient_local(y, x)
n = numel(y);
g = zeros(size(y));
if n < 2
    return;
end
g(1) = (y(2) - y(1)) / (x(2) - x(1));
g(n) = (y(n) - y(n - 1)) / (x(n) - x(n - 1));
for i = 2:(n - 1)
    g(i) = (y(i + 1) - y(i - 1)) / (x(i + 1) - x(i - 1));
end
end

function temperature = parse_temperature_labels_local(text, n_cols)
if iscell(text)
    header = text(1, :);
else
    header = cellstr(char(text(1, :)));
end

temperature = [];
for k = 1:numel(header)
    tok = header{k};
    if ischar(tok)
        tok(tok == 'K') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end + 1, 1) = val; %#ok<AGROW>
        end
    end
end

n_sets = floor(n_cols / 3);
temperature = temperature(1:min(numel(temperature), n_sets));
end

function idx_best = select_lambda_idx_local(rid_vec, cv_vec, var_vec)
% Paper Method v2: Score-based selection (FIXED Issue #3)
% Combines three metrics equally (normalized):
% - RID: Re/Im divergence (fit consistency)
% - CV: Cross-validation error (stability)
% - Var: Bootstrap resampling variance (robustness)
%
% Lambda selected: argmin(norm(RID) + norm(CV) + norm(Var))

norm_rid = normalize_metric_local(rid_vec);
norm_cv = normalize_metric_local(cv_vec);
norm_var = normalize_metric_local(var_vec);

score = norm_rid + norm_cv + norm_var;
[~, idx_best] = min(score);
end

function [gamma, R_inf] = TR_DRT_local(freq_vec, Z_exp, el)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];

ub = max(abs(Z_exp)) * ones(n + 1, 1);
if exist('lsqlin', 'file') == 2
    lb = zeros(n + 1, 1);
    options = optimset('Display', 'off');
    x = lsqlin(M, b, [], [], [], [], lb, ub, [], options);
else
    x = lsqnonneg(M, b);
    x = min(x, ub);
end

gamma = x(1:end-1);
R_inf = x(end);
end

function [gamma, R_inf] = TR_DRT_component_local(freq_vec, Z_exp, el, component)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

switch lower(component)
    case 'real'
        M = [A_re, ones(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
        b = [Z_re; zeros(n,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
        b = [Z_im; zeros(n,1)];
    otherwise
        error('Unknown component mode: %s', component);
end

ub = max(abs(Z_exp)) * ones(n + 1, 1);
if exist('lsqlin', 'file') == 2
    lb = zeros(n + 1, 1);
    options = optimset('Display', 'off');
    x = lsqlin(M, b, [], [], [], [], lb, ub, [], options);
else
    x = lsqnonneg(M, b);
    x = min(x, ub);
end

gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calculate_EIS_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function R_inf = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, gamma)
A_re = calc_A_re_local(freq_vec);
R_inf = mean(real(Z_exp(:)) - A_re * gamma(:));
if ~(R_inf > 0)
    R_inf = 0;
end
end

function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, el, n_boot)
n = numel(freq_vec);
if nargin < 4 || isempty(n_boot)
    n_boot = 8;
end

sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);
sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);

gamma_samples = zeros(n, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n,1) + 1i * sigma_im * randn(n,1);
    Z_boot = Z_exp(:) + noise;
    [gk, ~] = TR_DRT_local(freq_vec, Z_boot, el);
    gamma_samples(:, k) = gk(:);
end

mean_var = mean(var(gamma_samples, 0, 2));
end

function nvec = normalize_metric_local(vec)
vec = vec(:);
vmin = min(vec);
vmax = max(vec);
if vmax > vmin
    nvec = (vec - vmin) / (vmax - vmin);
else
    nvec = zeros(size(vec));
end
end

function A_re = calc_A_re_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);

A_re = zeros(n, n);
for p = 1:n
    for q = 1:n
        if q == 1
            log_term = log(tau(2) / tau(1));
        elseif q == n
            log_term = log(tau(n) / tau(n-1));
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
            log_term = log(tau(2) / tau(1));
        elseif q == n
            log_term = log(tau(n) / tau(n-1));
        else
            log_term = log(tau(q+1) / tau(q-1));
        end
        A_im(p, q) = 0.5 * (omega(p) * tau(q)) / (1 + (omega(p) * tau(q))^2) * log_term;
    end
end
end