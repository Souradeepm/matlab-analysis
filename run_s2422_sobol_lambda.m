function run_s2422_sobol_lambda()
% run_s2422_sobol_lambda
% Sobol sensitivity analysis for lambda impact on DRT fitting residual (S2422).
% MATLAB 2011-compatible implementation with Saltelli estimators.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 's2422_sobol_lambda.xlsx');
out_png = fullfile(repo_root, 's2422_sobol_lambda.png');

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

% Analyze around practical lambda region centered near 1e-3.
lambda_min = 3e-4;
lambda_max = 3e-3;
log_min = log10(lambda_min);
log_max = log10(lambda_max);

d = 1;
N = 64;

fprintf('Sobol analysis on lambda in [%.3e, %.3e], N=%d\n', lambda_min, lambda_max, N);

% Build Sobol A and B matrices in [0,1].
if exist('sobolset', 'file') == 2
    p = sobolset(2 * d);
    U = net(p, 2 * N);
    used_sobol = 1;
else
    U = rand(2 * N, 2 * d);
    used_sobol = 0;
end

A = U(1:N, 1:d);
B = U(N+1:2*N, d+1:2*d);

% Map to lambda via log-space transform.
lamA = 10 .^ (log_min + (log_max - log_min) * A(:,1));
lamB = 10 .^ (log_min + (log_max - log_min) * B(:,1));

YA = zeros(N,1);
YB = zeros(N,1);
YAB = zeros(N,1);

for k = 1:N
    YA(k) = mean_residual_for_lambda(freq_vec, Z_exp, lamA(k));
    YB(k) = mean_residual_for_lambda(freq_vec, Z_exp, lamB(k));
end

% For d=1, AB is simply B replacing A's only column.
for k = 1:N
    YAB(k) = mean_residual_for_lambda(freq_vec, Z_exp, lamB(k));
end

Yall = [YA; YB];
VY = var(Yall, 1);
if VY <= eps
    S1 = 0;
    ST = 0;
else
    % Jansen estimators (numerically stable and bounded for practical use)
    S1 = 1 - mean((YB - YAB).^2) / (2 * VY);
    ST = mean((YA - YAB).^2) / (2 * VY);
end

fprintf('Sobol first-order S1: %.6f\n', S1);
fprintf('Sobol total-order ST: %.6f\n', ST);

if used_sobol == 1
    fprintf('Sampling: Sobol sequence\n');
else
    fprintf('Sampling: random fallback (sobolset unavailable)\n');
end

% Dense response curve for interpretation.
lambda_grid = logspace(log_min, log_max, 40).';
res_grid = zeros(size(lambda_grid));
for i = 1:numel(lambda_grid)
    res_grid(i) = mean_residual_for_lambda(freq_vec, Z_exp, lambda_grid(i));
end

figure('Color', 'w', 'Name', 'S2422 Sobol lambda sensitivity');
subplot(2,1,1);
semilogx(lambda_grid, res_grid, 'b-', 'LineWidth', 1.3); hold on;
scatter(lamA, YA, 20, 'r', 'filled');
scatter(lamB, YB, 20, 'g', 'filled');
xlabel('lambda');
ylabel('Mean absolute residual');
title('Response curve and Sobol sample evaluations');
legend('Response curve', 'A samples', 'B samples', 'Location', 'best');
grid on;

subplot(2,1,2);
bar([S1, ST]);
set(gca, 'XTickLabel', {'S1', 'ST'});
ylabel('Sobol index');
title('Sobol sensitivity indices for lambda');
grid on;

saveas(gcf, out_png);

details = [lamA, YA, lamB, YB, YAB];
xlswrite(out_xlsx, {'lambda_A','mean_resid_A','lambda_B','mean_resid_B','mean_resid_AB'}, 'samples', 'A1');
xlswrite(out_xlsx, details, 'samples', 'A2');

summary = [lambda_min, lambda_max, N, used_sobol, S1, ST, mean(Yall), std(Yall,1), VY];
xlswrite(out_xlsx, {'lambda_min','lambda_max','N','used_sobol','S1','ST','mean_output','std_output','var_output'}, 'summary', 'A1');
xlswrite(out_xlsx, summary, 'summary', 'A2');

xlswrite(out_xlsx, {'lambda','mean_residual'}, 'response_curve', 'A1');
xlswrite(out_xlsx, [lambda_grid, res_grid], 'response_curve', 'A2');

fprintf('Saved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_png);
end

function mr = mean_residual_for_lambda(freq_vec, Z_exp, lambda)
[gamma, Rinf] = tr_drt_local(freq_vec, Z_exp, lambda);
Z_fit = calc_eis_local(freq_vec, gamma, Rinf);
mr = mean(abs(Z_exp - Z_fit));
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
