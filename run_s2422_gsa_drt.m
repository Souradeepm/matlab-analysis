function run_s2422_gsa_drt()
% run_s2422_gsa_drt
% Global Sensitivity Analysis (GSA) for DRT on S2422 dataset.
% Inputs (uncertain):
%   x1 = lambda (log-uniform)
%   x2 = noise scale (uniform)
% Outputs:
%   y1 = mean absolute residual
%   y2 = log10(dominant peak tau)
%   y3 = dominant peak amplitude
% Sobol indices are computed using Jansen estimators.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 's2422_gsa_drt.xlsx');
out_png = fullfile(repo_root, 's2422_gsa_drt.png');

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
% x1: lambda in [1e-5, 1e-1] (log-uniform)
lam_log_min = -5;
lam_log_max = -1;
% x2: additive-noise scale in [0, 0.03]
noise_min = 0.00;
noise_max = 0.03;

% Sobol sample size (total model runs ~ (d+2)*N = 4N for d=2)
N = 48;
d = 2;

fprintf('GSA-DRT on S2422 with N=%d, d=%d\n', N, d);
fprintf('lambda range: [1e%d, 1e%d]\n', lam_log_min, lam_log_max);
fprintf('noise scale range: [%.3f, %.3f]\n', noise_min, noise_max);

% ------------------ Baseline for deterministic noise template ------------------
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

rand('twister', 12345);
randn('state', 23456);
eta_re = randn(numel(freq_vec), 1);
eta_im = randn(numel(freq_vec), 1);
noise_template = sig_re * eta_re + 1i * sig_im * eta_im;

% ------------------ Sobol sampling matrices ------------------
if exist('sobolset', 'file') == 2
    p = sobolset(2*d);
    U = net(p, 2*N);
    used_sobol = 1;
else
    U = rand(2*N, 2*d);
    used_sobol = 0;
end

A = U(1:N, 1:d);
B = U(N+1:2*N, d+1:2*d);

% Map unit-cube to physical parameters.
XA = map_to_params(A, lam_log_min, lam_log_max, noise_min, noise_max);
XB = map_to_params(B, lam_log_min, lam_log_max, noise_min, noise_max);

% Build A_Bi matrices.
XAB1 = XA; XAB1(:,1) = XB(:,1);
XAB2 = XA; XAB2(:,2) = XB(:,2);

% ------------------ Model evaluations ------------------
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

% ------------------ Sobol indices (Jansen) ------------------
S1 = zeros(d, 3);
ST = zeros(d, 3);

for j = 1:3
    VY = var([YA(:,j); YB(:,j)], 1);
    if VY <= eps
        continue;
    end

    % Parameter 1
    S1(1,j) = 1 - mean((YB(:,j) - YAB1(:,j)).^2) / (2 * VY);
    ST(1,j) = mean((YA(:,j) - YAB1(:,j)).^2) / (2 * VY);

    % Parameter 2
    S1(2,j) = 1 - mean((YB(:,j) - YAB2(:,j)).^2) / (2 * VY);
    ST(2,j) = mean((YA(:,j) - YAB2(:,j)).^2) / (2 * VY);
end

param_names = {'lambda','noise_scale'};
out_names = {'mean_residual','log10_tau_peak','peak_amplitude'};

fprintf('\nSobol S1 indices:\n');
for i = 1:d
    fprintf('  %-11s: [%8.4f  %8.4f  %8.4f]\n', param_names{i}, S1(i,1), S1(i,2), S1(i,3));
end

fprintf('\nSobol ST indices:\n');
for i = 1:d
    fprintf('  %-11s: [%8.4f  %8.4f  %8.4f]\n', param_names{i}, ST(i,1), ST(i,2), ST(i,3));
end

% ------------------ Visualization ------------------
figure('Color', 'w', 'Name', 'S2422 GSA-DRT Sobol indices');
for j = 1:3
    subplot(3,1,j);
    M = [S1(:,j), ST(:,j)];
    bar(M);
    set(gca, 'XTickLabel', param_names);
    ylabel('Index value');
    title(sprintf('Output: %s', out_names{j}));
    legend('S1','ST','Location','best');
    grid on;
end
saveas(gcf, out_png);

% ------------------ Export ------------------
xlswrite(out_xlsx, {'lambda','noise_scale','y_mean_residual','y_log10_tau_peak','y_peak_amplitude'}, 'A_samples', 'A1');
xlswrite(out_xlsx, [XA, YA], 'A_samples', 'A2');

xlswrite(out_xlsx, {'lambda','noise_scale','y_mean_residual','y_log10_tau_peak','y_peak_amplitude'}, 'B_samples', 'A1');
xlswrite(out_xlsx, [XB, YB], 'B_samples', 'A2');

xlswrite(out_xlsx, {'parameter','S1_mean_residual','S1_log10_tau_peak','S1_peak_amplitude','ST_mean_residual','ST_log10_tau_peak','ST_peak_amplitude'}, 'sobol_indices', 'A1');
idx_table = [1, S1(1,1), S1(1,2), S1(1,3), ST(1,1), ST(1,2), ST(1,3); ...
             2, S1(2,1), S1(2,2), S1(2,3), ST(2,1), ST(2,2), ST(2,3)];
xlswrite(out_xlsx, idx_table, 'sobol_indices', 'A2');
xlswrite(out_xlsx, {'1=lambda, 2=noise_scale'}, 'sobol_indices', 'A5');

xlswrite(out_xlsx, {'N','used_sobol','lambda_ref','sig_re','sig_im','lambda_min','lambda_max','noise_min','noise_max'}, 'meta', 'A1');
xlswrite(out_xlsx, [N, used_sobol, lambda_ref, sig_re, sig_im, 10^lam_log_min, 10^lam_log_max, noise_min, noise_max], 'meta', 'A2');

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
