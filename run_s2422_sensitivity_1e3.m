function run_s2422_sensitivity_1e3()
% run_s2422_sensitivity_1e3
% Local lambda sensitivity around 1e-3 for S2422 dataset.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 's2422_sensitivity_lambda_1e3.xlsx');
out_png = fullfile(repo_root, 's2422_sensitivity_lambda_1e3.png');

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

lambda_ref = 1e-3;
factors = [0.3; 0.5; 0.7; 1.0; 1.5; 2.0; 3.0];
lambda_vals = lambda_ref * factors;

n = numel(lambda_vals);
mean_resid = zeros(n,1);
max_resid = zeros(n,1);

for i = 1:n
    [gamma_i, Rinf_i] = tr_drt_local(freq_vec, Z_exp, lambda_vals(i));
    Z_fit_i = calc_eis_local(freq_vec, gamma_i, Rinf_i);
    r = abs(Z_exp - Z_fit_i);
    mean_resid(i) = mean(r);
    max_resid(i) = max(r);
end

ref_idx = find(abs(factors - 1) < 1e-12, 1);
ref_mean = mean_resid(ref_idx);
delta_pct = 100 * (mean_resid - ref_mean) / max(ref_mean, eps);

[min_delta, min_i] = min(delta_pct);
[max_delta, max_i] = max(delta_pct);
max_abs_delta = max(abs(delta_pct));

fprintf('Sensitivity around lambda = %.3e\n', lambda_ref);
fprintf('Reference mean residual: %.6f\n', ref_mean);
fprintf('Min delta(%%): %.4f at factor %.2f (lambda %.3e)\n', min_delta, factors(min_i), lambda_vals(min_i));
fprintf('Max delta(%%): %.4f at factor %.2f (lambda %.3e)\n', max_delta, factors(max_i), lambda_vals(max_i));
fprintf('Max |delta|(%%): %.4f\n', max_abs_delta);

figure('Color','w','Name','S2422 sensitivity around 1e-3');
subplot(2,1,1);
semilogx(lambda_vals, mean_resid, 'bo-', 'LineWidth', 1.2); hold on;
plot(lambda_ref, ref_mean, 'rs', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('lambda');
ylabel('Mean absolute residual');
grid on;
legend('Perturbed lambdas', 'Reference 1e-3', 'Location', 'best');

subplot(2,1,2);
semilogx(lambda_vals, delta_pct, 'm.-', 'LineWidth', 1.2); hold on;
plot(lambda_ref, 0, 'ks', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('lambda');
ylabel('Residual change (%)');
grid on;
legend('Delta vs 1e-3', 'Reference', 'Location', 'best');

saveas(gcf, out_png);

xlswrite(out_xlsx, {'factor','lambda','mean_residual','max_residual','delta_pct'}, 'details', 'A1');
xlswrite(out_xlsx, [factors, lambda_vals, mean_resid, max_resid, delta_pct], 'details', 'A2');

xlswrite(out_xlsx, {'lambda_ref','reference_mean_residual','min_delta_pct','max_delta_pct','max_abs_delta_pct'}, 'summary', 'A1');
xlswrite(out_xlsx, [lambda_ref, ref_mean, min_delta, max_delta, max_abs_delta], 'summary', 'A2');

fprintf('Saved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_png);
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
