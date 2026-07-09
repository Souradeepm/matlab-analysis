function run_s2422_lambda_method_compare()
% Compare lambda sensitivity for two selection methods on S2422 dataset.
% Method 1: minimum mean residual
% Method 2: L-curve corner (maximum curvature)

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_file = fullfile(repo_root, 'lambda_method_sensitivity_s2422.xlsx');

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

% ---------------- Method 1: residual-min lambda ----------------
lambda_candidates = [1e-4, 1e-3, 1e-2, 1e-1, 1e0].';
nc = numel(lambda_candidates);
mean_resid_candidates = zeros(nc,1);

for i = 1:nc
    [gamma_i, R_i] = tr_drt_local(freq_vec, Z_exp, lambda_candidates(i));
    Z_fit_i = calc_eis_local(freq_vec, gamma_i, R_i);
    mean_resid_candidates(i) = mean(abs(Z_exp - Z_fit_i));
end

[~, idx_min] = min(mean_resid_candidates);
lambda_resid = lambda_candidates(idx_min);

% ---------------- Method 2: L-curve corner lambda ----------------
lambda_scan = logspace(-6, 1, 40).';
ns = numel(lambda_scan);
res_norm = zeros(ns,1);
sol_norm = zeros(ns,1);

for i = 1:ns
    [gamma_i, R_i] = tr_drt_local(freq_vec, Z_exp, lambda_scan(i));
    Z_fit_i = calc_eis_local(freq_vec, gamma_i, R_i);
    dz = Z_exp - Z_fit_i;
    res_norm(i) = norm([real(dz); imag(dz)], 2);
    sol_norm(i) = norm(gamma_i, 2);
end

curv = lcurve_curvature_local(lambda_scan, res_norm, sol_norm);
[~, idx_corner] = max(curv);
lambda_lcurve = lambda_scan(idx_corner);

% ---------------- Common perturbation sensitivity ----------------
factors = [0.3; 0.5; 0.7; 1.0; 1.5; 2.0; 3.0];
resid_result = perturbation_result_local(freq_vec, Z_exp, lambda_resid, factors);
lcurve_result = perturbation_result_local(freq_vec, Z_exp, lambda_lcurve, factors);

fprintf('Residual-min method lambda: %.6e\n', lambda_resid);
fprintf('  Min delta(%%): %.4f | Max delta(%%): %.4f | MaxAbs delta(%%): %.4f\n', ...
    resid_result.min_delta_pct, resid_result.max_delta_pct, resid_result.max_abs_delta_pct);

fprintf('L-curve method lambda: %.6e\n', lambda_lcurve);
fprintf('  Min delta(%%): %.4f | Max delta(%%): %.4f | MaxAbs delta(%%): %.4f\n', ...
    lcurve_result.min_delta_pct, lcurve_result.max_delta_pct, lcurve_result.max_abs_delta_pct);

figure('Color', 'w', 'Name', 'S2422 lambda sensitivity comparison');
subplot(2,1,1);
semilogx(resid_result.lambda_values, resid_result.delta_pct, 'bo-'); hold on;
plot(lambda_resid, 0, 'bs', 'MarkerSize', 8, 'LineWidth', 1.2);
semilogx(lcurve_result.lambda_values, lcurve_result.delta_pct, 'ro-');
plot(lambda_lcurve, 0, 'rs', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('lambda');
ylabel('Residual change (%)');
title('Lambda perturbation sensitivity (S2422)');
legend('Residual-min perturbation', 'Residual-min reference', ...
       'L-curve perturbation', 'L-curve reference', 'Location', 'best');
grid on;

subplot(2,1,2);
bar([resid_result.max_abs_delta_pct, lcurve_result.max_abs_delta_pct]);
set(gca, 'XTickLabel', {'Residual-min', 'L-curve'});
ylabel('Max |delta residual| (%)');
title('Sensitivity summary');
grid on;

summary_headers = {'method_code','lambda_ref','min_delta_pct','max_delta_pct','max_abs_delta_pct'};
summary_data = [1, lambda_resid, resid_result.min_delta_pct, resid_result.max_delta_pct, resid_result.max_abs_delta_pct; ...
                2, lambda_lcurve, lcurve_result.min_delta_pct, lcurve_result.max_delta_pct, lcurve_result.max_abs_delta_pct];

resid_detail = [factors, resid_result.lambda_values, resid_result.mean_resid, resid_result.delta_pct];
lcurve_detail = [factors, lcurve_result.lambda_values, lcurve_result.mean_resid, lcurve_result.delta_pct];

xlswrite(out_file, summary_headers, 'summary', 'A1');
xlswrite(out_file, summary_data, 'summary', 'A2');

xlswrite(out_file, {'factor','lambda','mean_resid','delta_pct'}, 'residual_method', 'A1');
xlswrite(out_file, resid_detail, 'residual_method', 'A2');

xlswrite(out_file, {'factor','lambda','mean_resid','delta_pct'}, 'lcurve_method', 'A1');
xlswrite(out_file, lcurve_detail, 'lcurve_method', 'A2');

fprintf('Saved comparison file: %s\n', out_file);
end

function res = perturbation_result_local(freq_vec, Z_exp, lambda_ref, factors)
lambda_values = lambda_ref * factors;
n = numel(lambda_values);
mean_resid = zeros(n,1);

for i = 1:n
    [gamma_i, R_i] = tr_drt_local(freq_vec, Z_exp, lambda_values(i));
    Z_fit_i = calc_eis_local(freq_vec, gamma_i, R_i);
    mean_resid(i) = mean(abs(Z_exp - Z_fit_i));
end

ref_idx = find(abs(factors - 1) < 1e-12, 1);
if isempty(ref_idx)
    ref_idx = 1;
end
ref_resid = mean_resid(ref_idx);

delta_pct = 100 * (mean_resid - ref_resid) / max(ref_resid, eps);

res = struct('lambda_values', lambda_values, ...
             'mean_resid', mean_resid, ...
             'delta_pct', delta_pct, ...
             'min_delta_pct', min(delta_pct), ...
             'max_delta_pct', max(delta_pct), ...
             'max_abs_delta_pct', max(abs(delta_pct)));
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
