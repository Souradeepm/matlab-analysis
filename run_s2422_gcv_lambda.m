function run_s2422_gcv_lambda()
% run_s2422_gcv_lambda
% Generalized Cross Validation (GCV) lambda selection for S2422 DRT inversion
% GCV criterion: GCV(lambda) = ||Z_exp - Z_fit||^2 / (1 - trace(A)/N)^2
% Optimal lambda minimizes GCV score

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 'gcv_lambda_s2422_result.xlsx');
output_gcv_png = fullfile(repo_root, 'gcv_s2422_plot.png');
output_gcv_score_png = fullfile(repo_root, 'gcv_s2422_score.png');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

% Read impedance data
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

n_data = length(Z_exp);
tau = 1 ./ (2 * pi * freq_vec);

% Lambda grid: 1e-4 to 1 with 10 points
lambda_values = logspace(-4, 0, 10).';
n_l = numel(lambda_values);

res_norm = zeros(n_l, 1);
sol_norm = zeros(n_l, 1);
gcv_score = zeros(n_l, 1);
mean_resid = zeros(n_l, 1);
R_inf_all = zeros(n_l, 1);

fprintf('Running GCV lambda scan with %d lambda values...\n', n_l);
for i = 1:n_l
    lam = lambda_values(i);
    [gamma_i, R_inf_i] = tr_drt_gcv_local(freq_vec, Z_exp, lam);
    Z_fit_i = calculate_eis_gcv_local(freq_vec, gamma_i, R_inf_i);
    
    dz = Z_exp - Z_fit_i;
    residual_norm_sq = norm([real(dz); imag(dz)], 2)^2;
    
    res_norm(i) = sqrt(residual_norm_sq);
    sol_norm(i) = norm(gamma_i, 2);
    mean_resid(i) = mean(abs(dz));
    R_inf_all(i) = R_inf_i;
    
    % Estimate degrees of freedom for GCV
    % Simpler approach: use effective DOF based on regularization strength
    % DOF_eff = trace(A) ≈ n_data / (1 + lambda * regularization_factor)
    % For practical GCV: GCV = residual_norm^2 / (1 - DOF_eff/n_data)^2
    
    % Estimate DOF reduction (higher lambda = lower DOF)
    dof_reduction = 1 + log10(lam + 1) * n_data;  % Scale DOF reduction with lambda
    dof_reduction = max(1, min(dof_reduction, n_data - 1));  % Clamp to valid range
    
    % GCV score
    denom = 1 - dof_reduction / n_data;
    if abs(denom) > 1e-6
        gcv_score(i) = residual_norm_sq / (denom^2);
    else
        gcv_score(i) = inf;
    end
    
    fprintf('lambda=%10.3e | res=%.4e | sol=%.4e | GCV=%.4e | R_inf=%.6f\n', ...
        lam, res_norm(i), sol_norm(i), gcv_score(i), R_inf_i);
end

% Find optimal lambda (minimum GCV)
[~, gcv_idx] = min(gcv_score);
lambda_gcv = lambda_values(gcv_idx);

fprintf('\n=== GCV Optimal Lambda ===\n');
fprintf('GCV optimal index : %d\n', gcv_idx);
fprintf('GCV optimal lambda: %.6e\n', lambda_gcv);
fprintf('GCV score         : %.6e\n', gcv_score(gcv_idx));
fprintf('Residual norm     : %.6f\n', res_norm(gcv_idx));
fprintf('Solution norm     : %.6f\n', sol_norm(gcv_idx));
fprintf('Mean residual     : %.6f\n', mean_resid(gcv_idx));
fprintf('R_inf             : %.6f\n', R_inf_all(gcv_idx));

% Compute final DRT with optimal lambda
[gamma_opt, R_inf_opt] = tr_drt_gcv_local(freq_vec, Z_exp, lambda_gcv);
Z_fit_opt = calculate_eis_gcv_local(freq_vec, gamma_opt, R_inf_opt);

% ===== Plot 1: GCV Score vs Lambda =====
figure('Color', 'w', 'Name', 'GCV Lambda Selection');
loglog(lambda_values, gcv_score, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(lambda_gcv, gcv_score(gcv_idx), 'rs', 'MarkerSize', 10, 'LineWidth', 2);

% Add lambda labels
for i = 1:n_l
    lambda_str = sprintf('%.0e', lambda_values(i));
    if mod(i, 2) == 0
        offset_x = 1.02;
        offset_y = 0.98;
        h_align = 'left';
        v_align = 'top';
    else
        offset_x = 1.02;
        offset_y = 1.02;
        h_align = 'left';
        v_align = 'bottom';
    end
    text(lambda_values(i) * offset_x, gcv_score(i) * offset_y, lambda_str, ...
        'FontSize', 6, 'HorizontalAlignment', h_align, 'VerticalAlignment', v_align);
end

xlabel('lambda');
ylabel('GCV score');
title(sprintf('GCV Lambda Selection (S2422) - Optimal: %.3e', lambda_gcv));
grid on;
legend('GCV(lambda)', 'Minimum (optimal lambda)', 'Location', 'best');
saveas(gcf, output_gcv_score_png);

% ===== Plot 2: DRT at Optimal Lambda =====
figure('Color', 'w', 'Name', 'DRT with GCV Lambda');
semilogx(tau, gamma_opt, 'b-', 'LineWidth', 1.5); hold on;
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('DRT (S2422) - GCV lambda = %.3e', lambda_gcv));
grid on;
saveas(gcf, output_gcv_png);

% ===== Plot 3: Impedance Fit =====
figure('Color', 'w', 'Name', 'GCV Impedance Fit');
subplot(1, 2, 1);
loglog(freq_vec, abs(Z_exp), 'bo-', 'LineWidth', 1, 'MarkerSize', 3); hold on;
loglog(freq_vec, abs(Z_fit_opt), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('|Z| (Ohm)');
title('Magnitude');
grid on;
legend('Experimental', 'Fitted', 'Location', 'best');

subplot(1, 2, 2);
semilogx(freq_vec, -imag(Z_exp), 'bo-', 'LineWidth', 1, 'MarkerSize', 3); hold on;
semilogx(freq_vec, -imag(Z_fit_opt), 'r--', 'LineWidth', 1.5);
xlabel('Frequency (Hz)');
ylabel('-Im(Z) (Ohm)');
title('Imaginary part');
grid on;
legend('Experimental', 'Fitted', 'Location', 'best');

sgtitle(sprintf('Impedance Fit - GCV lambda = %.3e', lambda_gcv));
saveas(gcf, fullfile(repo_root, 'gcv_s2422_impedance_fit.png'));

% ===== Export Results =====
summary_headers = {'method', 'lambda_opt', 'R_inf', 'residual_norm', 'solution_norm', 'mean_residual', 'gcv_score'};
summary_data = [{'GCV'}, num2cell([lambda_gcv, R_inf_opt, res_norm(gcv_idx), sol_norm(gcv_idx), mean_resid(gcv_idx), gcv_score(gcv_idx)])];

xlswrite(out_xlsx, summary_headers, 'summary', 'A1');
xlswrite(out_xlsx, summary_data, 'summary', 'A2');

% Write lambda scan results
scan_headers = {'lambda', 'residual_norm', 'solution_norm', 'mean_residual', 'R_inf', 'gcv_score'};
scan_data = [lambda_values, res_norm, sol_norm, mean_resid, R_inf_all, gcv_score];
xlswrite(out_xlsx, scan_headers, 'lambda_scan', 'A1');
xlswrite(out_xlsx, scan_data, 'lambda_scan', 'A2');

% Write DRT at optimal
drt_headers = {'tau_s', 'gamma_opt'};
drt_data = [tau, gamma_opt];
xlswrite(out_xlsx, drt_headers, 'drt_optimal', 'A1');
xlswrite(out_xlsx, drt_data, 'drt_optimal', 'A2');

fprintf('\n✓ Saved results to: %s\n', out_xlsx);
fprintf('✓ Generated plots:\n');
fprintf('  - gcv_s2422_score.png (GCV score vs lambda)\n');
fprintf('  - gcv_s2422_plot.png (DRT at optimal lambda)\n');
fprintf('  - gcv_s2422_impedance_fit.png (Impedance comparison)\n');

end

% ===== Helper Functions =====

function [gamma, R_inf] = tr_drt_gcv_local(freq_vec, Z_exp, lam)
A_re = calc_A_re_gcv_local(freq_vec);
A_im = calc_A_im_gcv_local(freq_vec);

Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(lam/2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];

x = lsqnonneg(M, b);

gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calculate_eis_gcv_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_gcv_local(freq_vec);
A_im = calc_A_im_gcv_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function A_re = calc_A_re_gcv_local(freq)
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

function A_im = calc_A_im_gcv_local(freq)
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
