function run_s2422_plot_two_drt()
% run_s2422_plot_two_drt - Plot two DRT curves for lambdas chosen by
% residual-min and L-curve methods on S2422 dataset.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
plot_path = fullfile(repo_root, 's2422_two_drt_compare.png');
out_path = fullfile(repo_root, 's2422_two_drt_compare.xlsx');

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

% Lambda by residual minimum (same logic as workflow)
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

% Lambda by L-curve corner
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

% Compute two DRT curves
[gamma_resid, R_resid] = tr_drt_local(freq_vec, Z_exp, lambda_resid);
[gamma_lcurve, R_lcurve] = tr_drt_local(freq_vec, Z_exp, lambda_lcurve);
tau = 1 ./ (2 * pi * freq_vec);

figure('Color', 'w', 'Name', 'S2422 two DRT comparison');
semilogx(tau, gamma_resid, 'b-', 'LineWidth', 1.5); hold on;
semilogx(tau, gamma_lcurve, 'r--', 'LineWidth', 1.5);
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title('S2422 DRT comparison for two lambda-selection methods');
legend(sprintf('Residual-min lambda = %.3e', lambda_resid), ...
       sprintf('L-curve lambda = %.3e', lambda_lcurve), ...
       'Location', 'best');
grid on;

saveas(gcf, plot_path);

xlswrite(out_path, {'tau_s','gamma_residual_method','gamma_lcurve_method'}, 'drt', 'A1');
xlswrite(out_path, [tau, gamma_resid, gamma_lcurve], 'drt', 'A2');
xlswrite(out_path, {'lambda_residual_method','lambda_lcurve_method','Rinf_residual_method','Rinf_lcurve_method'}, 'summary', 'A1');
xlswrite(out_path, [lambda_resid, lambda_lcurve, R_resid, R_lcurve], 'summary', 'A2');

fprintf('Residual-min lambda: %.6e\n', lambda_resid);
fprintf('L-curve lambda     : %.6e\n', lambda_lcurve);
fprintf('Saved plot          : %s\n', plot_path);
fprintf('Saved DRT table     : %s\n', out_path);
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

    denom = (dx*dx + dy*dy)^(3/2);
    if denom > 0
        curv(i) = abs(dx * d2y - dy * d2x) / denom;
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
