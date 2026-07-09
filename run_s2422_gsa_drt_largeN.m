function run_s2422_gsa_drt_largeN()
% run_s2422_gsa_drt_largeN
% Larger-N GSA-DRT with repeated Sobol runs for stability statistics.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 's2422_gsa_drt_largeN.xlsx');
out_png = fullfile(repo_root, 's2422_gsa_drt_largeN.png');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

% ------------------ Load data ------------------
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

% ------------------ Parameter ranges ------------------
lam_log_min = -5;
lam_log_max = -1;
noise_min = 0.00;
noise_max = 0.03;

% Larger sample and repeated runs
N = 128;
d = 2;
R = 5;

fprintf('Large-N GSA-DRT on S2422: N=%d, repeats=%d\n', N, R);

% Baseline residual scales
lambda_ref = 1e-3;
[gamma_ref, R_ref] = tr_drt_local(freq_vec, Z_exp, lambda_ref);
Z_ref = calc_eis_local(freq_vec, gamma_ref, R_ref);
res_ref = Z_exp - Z_ref;
sig_re = std(real(res_ref));
sig_im = std(imag(res_ref));
if ~(sig_re > 0)
    sig_re = 0.005 * max(max(abs(real(Z_exp))), 1);
end
if ~(sig_im > 0)
    sig_im = 0.005 * max(max(abs(imag(Z_exp))), 1);
end

S1_runs = zeros(R, d, 3);
ST_runs = zeros(R, d, 3);

for r = 1:R
    seed_base = 9000 + 100 * r;
    rand('twister', seed_base);
    randn('state', seed_base + 1);

    eta_re = randn(numel(freq_vec), 1);
    eta_im = randn(numel(freq_vec), 1);
    noise_template = sig_re * eta_re + 1i * sig_im * eta_im;

    if exist('sobolset', 'file') == 2
        p = sobolset(2 * d);
        p = scramble(p, 'MatousekAffineOwen');
        U = net(p, 2 * N);
        used_sobol = 1;
    else
        U = rand(2 * N, 2 * d);
        used_sobol = 0;
    end

    A = U(1:N, 1:d);
    B = U(N+1:2*N, d+1:2*d);

    XA = map_to_params(A, lam_log_min, lam_log_max, noise_min, noise_max);
    XB = map_to_params(B, lam_log_min, lam_log_max, noise_min, noise_max);

    XAB1 = XA; XAB1(:,1) = XB(:,1);
    XAB2 = XA; XAB2(:,2) = XB(:,2);

    YA = zeros(N, 3);
    YB = zeros(N, 3);
    YAB1 = zeros(N, 3);
    YAB2 = zeros(N, 3);

    for k = 1:N
        YA(k,:) = model_outputs(freq_vec, Z_exp, tau, XA(k,1), XA(k,2), noise_template);
        YB(k,:) = model_outputs(freq_vec, Z_exp, tau, XB(k,1), XB(k,2), noise_template);
        YAB1(k,:) = model_outputs(freq_vec, Z_exp, tau, XAB1(k,1), XAB1(k,2), noise_template);
        YAB2(k,:) = model_outputs(freq_vec, Z_exp, tau, XAB2(k,1), XAB2(k,2), noise_template);
    end

    for j = 1:3
        VY = var([YA(:,j); YB(:,j)], 1);
        if VY <= eps
            continue;
        end

        S1_runs(r,1,j) = 1 - mean((YB(:,j) - YAB1(:,j)).^2) / (2 * VY);
        ST_runs(r,1,j) = mean((YA(:,j) - YAB1(:,j)).^2) / (2 * VY);

        S1_runs(r,2,j) = 1 - mean((YB(:,j) - YAB2(:,j)).^2) / (2 * VY);
        ST_runs(r,2,j) = mean((YA(:,j) - YAB2(:,j)).^2) / (2 * VY);
    end

    fprintf('Repeat %d done\n', r);
end

S1_mean = squeeze(mean(S1_runs, 1));
ST_mean = squeeze(mean(ST_runs, 1));
S1_std = squeeze(std(S1_runs, 0, 1));
ST_std = squeeze(std(ST_runs, 0, 1));

param_names = {'lambda','noise_scale'};
out_names = {'mean_residual','log10_tau_peak','peak_amplitude'};

fprintf('\nMean Sobol S1 (large-N):\n');
for i = 1:d
    fprintf('  %-11s: [%8.4f  %8.4f  %8.4f]\n', param_names{i}, S1_mean(i,1), S1_mean(i,2), S1_mean(i,3));
end

fprintf('\nMean Sobol ST (large-N):\n');
for i = 1:d
    fprintf('  %-11s: [%8.4f  %8.4f  %8.4f]\n', param_names{i}, ST_mean(i,1), ST_mean(i,2), ST_mean(i,3));
end

% ------------------ Plot with error bars ------------------
figure('Color', 'w', 'Name', 'S2422 GSA-DRT large-N');
for j = 1:3
    subplot(3,1,j);
    x = [1; 2];
    hold on;
    errorbar(x - 0.08, S1_mean(:,j), S1_std(:,j), 'bo', 'LineWidth', 1.2);
    errorbar(x + 0.08, ST_mean(:,j), ST_std(:,j), 'rs', 'LineWidth', 1.2);
    set(gca, 'XTick', x);
    set(gca, 'XTickLabel', param_names);
    ylabel('Index value');
    title(sprintf('Output: %s', out_names{j}));
    legend('S1 mean \pm std', 'ST mean \pm std', 'Location', 'best');
    grid on;
end
saveas(gcf, out_png);

% ------------------ Export ------------------
xlswrite(out_xlsx, {'parameter_code','S1_mean_residual_mean','S1_log10_tau_mean','S1_peak_amp_mean', ...
                    'S1_mean_residual_std','S1_log10_tau_std','S1_peak_amp_std', ...
                    'ST_mean_residual_mean','ST_log10_tau_mean','ST_peak_amp_mean', ...
                    'ST_mean_residual_std','ST_log10_tau_std','ST_peak_amp_std'}, 'summary', 'A1');

summary_table = [1, S1_mean(1,1), S1_mean(1,2), S1_mean(1,3), S1_std(1,1), S1_std(1,2), S1_std(1,3), ST_mean(1,1), ST_mean(1,2), ST_mean(1,3), ST_std(1,1), ST_std(1,2), ST_std(1,3); ...
                 2, S1_mean(2,1), S1_mean(2,2), S1_mean(2,3), S1_std(2,1), S1_std(2,2), S1_std(2,3), ST_mean(2,1), ST_mean(2,2), ST_mean(2,3), ST_std(2,1), ST_std(2,2), ST_std(2,3)];
xlswrite(out_xlsx, summary_table, 'summary', 'A2');
xlswrite(out_xlsx, {'1=lambda, 2=noise_scale'}, 'summary', 'A5');

% Flat per-repeat export
repeat_rows = R * d;
rep_export = zeros(repeat_rows, 8);
row = 1;
for r = 1:R
    for i = 1:d
        rep_export(row,:) = [r, i, S1_runs(r,i,1), S1_runs(r,i,2), S1_runs(r,i,3), ST_runs(r,i,1), ST_runs(r,i,2), ST_runs(r,i,3)];
        row = row + 1;
    end
end

xlswrite(out_xlsx, {'repeat_id','parameter_code','S1_mean_residual','S1_log10_tau_peak','S1_peak_amp','ST_mean_residual','ST_log10_tau_peak','ST_peak_amp'}, 'per_repeat', 'A1');
xlswrite(out_xlsx, rep_export, 'per_repeat', 'A2');

xlswrite(out_xlsx, {'N','R','used_sobol','lambda_ref','sig_re','sig_im','lambda_min','lambda_max','noise_min','noise_max'}, 'meta', 'A1');
xlswrite(out_xlsx, [N, R, used_sobol, lambda_ref, sig_re, sig_im, 10^lam_log_min, 10^lam_log_max, noise_min, noise_max], 'meta', 'A2');

fprintf('\nSaved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_png);
end

function X = map_to_params(U, lam_log_min, lam_log_max, noise_min, noise_max)
X = zeros(size(U,1), 2);
X(:,1) = 10 .^ (lam_log_min + (lam_log_max - lam_log_min) * U(:,1));
X(:,2) = noise_min + (noise_max - noise_min) * U(:,2);
end

function y = model_outputs(freq_vec, Z_exp, tau, lambda, noise_scale, noise_template)
Z_model = Z_exp + noise_scale * noise_template;
[gamma, R_inf] = tr_drt_local(freq_vec, Z_model, lambda);
Z_fit = calc_eis_local(freq_vec, gamma, R_inf);
res = abs(Z_model - Z_fit);

mean_resid = mean(res);
[amp, idx] = max(gamma);
if isempty(idx) || ~(idx >= 1)
    tau_peak = tau(round(numel(tau)/2));
    amp = 0;
else
    tau_peak = tau(idx);
end

y = [mean_resid, log10(max(tau_peak, eps)), amp];
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
