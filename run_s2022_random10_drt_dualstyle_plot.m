function run_s2022_random10_drt_dualstyle_plot()
% run_s2022_random10_drt_dualstyle_plot
% MATLAB 2011-compatible DRT comparison plot for random 10 S2022 temperatures.
% For each selected temperature, plots gamma(tau) for:
%   1) residual-only best lambda
%   2) paper-style best lambda (RID + CV + variance)

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_png = fullfile(repo_root, 's2022_random10_lambda_drt_compare.png');
out_txt = fullfile(repo_root, 's2022_random10_lambda_drt_compare_selection.txt');

rng_seed = 20260710;
n_select = 10;
lambda_values = logspace(-7, -2, 11).';

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

residual_best_idx = zeros(n_select, 1);
paper_best_idx = zeros(n_select, 1);

fig = figure('Color', 'w', 'Name', 'S2022 random 10 DRT compare (dual style)');

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

    [~, residual_best_idx(t)] = min(mean_abs_residual);
    paper_best_idx(t) = select_lambda_idx_local(rid, cv_total, var_metric);

    [gamma_res, ~] = TR_DRT_local(freq_vec, Z_exp, lambda_values(residual_best_idx(t)));
    [gamma_paper, ~] = TR_DRT_local(freq_vec, Z_exp, lambda_values(paper_best_idx(t)));
    tau = 1 ./ (2 * pi * freq_vec);

    subplot(5, 2, t);
    semilogx(tau, gamma_res, 'r-', 'LineWidth', 1.0); hold on;
    semilogx(tau, gamma_paper, 'b-', 'LineWidth', 1.0);
    set(gca, 'XDir', 'reverse');
    grid on;
    title(sprintf('T=%.2fK, res %.1e, paper %.1e', selected_temperature(t), ...
        lambda_values(residual_best_idx(t)), lambda_values(paper_best_idx(t))), 'FontSize', 8);
    if mod(t, 2) == 1
        ylabel('gamma(tau)');
    end
    if t > 8
        xlabel('tau (s)');
    end
    if t == 1
        legend('Residual best', 'Paper best', 'Location', 'Best');
    end
end

saveas(fig, out_png);
close(fig);

fid = fopen(out_txt, 'w');
if fid > 0
    fprintf(fid, 'S2022 random 10 DRT dual-style selection\n');
    fprintf(fid, 'Random seed: %d\n', rng_seed);
    fprintf(fid, 'TemperatureK,ResidualBestLambda,PaperBestLambda\n');
    for t = 1:n_select
        fprintf(fid, '%.6f,%.8e,%.8e\n', selected_temperature(t), ...
            lambda_values(residual_best_idx(t)), lambda_values(paper_best_idx(t)));
    end
    fclose(fid);
end

fprintf('Wrote DRT comparison plot: %s\n', out_png);
fprintf('Wrote selection table   : %s\n', out_txt);
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
[~, idx_cv_min] = min(cv_vec);
[~, idx_var_min] = min(var_vec);
idx_low = min(idx_cv_min, idx_var_min);
idx_high = max(idx_cv_min, idx_var_min);

if idx_low == idx_high
    idx_candidates = (1:numel(rid_vec)).';
else
    idx_candidates = (idx_low:idx_high).';
end

[~, rid_rel_idx] = min(rid_vec(idx_candidates));
idx_best = idx_candidates(rid_rel_idx);

score = normalize_metric_local(rid_vec) + normalize_metric_local(cv_vec) + normalize_metric_local(var_vec);
[~, idx_score_best] = min(score);
if idx_score_best ~= idx_best
    idx_best = idx_score_best;
end
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