function run_s2422_lcurve()
% run_s2422_lcurve - L-curve analysis for S2422 DRT inversion (MATLAB 2011 compatible)

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
output_xlsx = fullfile(repo_root, 'lcurve_s2422_result.xlsx');
output_lcurve_png = fullfile(repo_root, 'lcurve_s2422_plot.png');
output_curvature_png = fullfile(repo_root, 'lcurve_s2422_curvature.png');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

data = dlmread(input_path, '\t');
if size(data, 2) < 3
    error('S2422 input must contain frequency, magnitude, phase columns');
end

freq_vec = data(:,1);
mag = data(:,2);
phase_deg = data(:,3);

Z_exp = mag .* exp(1i * phase_deg * pi / 180);

valid = freq_vec > 0;
freq_vec = freq_vec(valid);
Z_exp = Z_exp(valid);

[freq_vec, idx] = sort(freq_vec);
Z_exp = Z_exp(idx);

lambda_values = logspace(-4, 0, 10).';
n_l = numel(lambda_values);

res_norm = zeros(n_l, 1);
sol_norm = zeros(n_l, 1);
mean_resid = zeros(n_l, 1);
R_inf_all = zeros(n_l, 1);

fprintf('Running L-curve scan with %d lambda values...\n', n_l);
for i = 1:n_l
    lam = lambda_values(i);
    [gamma_i, R_inf_i] = tr_drt_lcurve_local(freq_vec, Z_exp, lam);
    Z_fit_i = calculate_eis_local(freq_vec, gamma_i, R_inf_i);

    dz = Z_exp - Z_fit_i;
    res_norm(i) = norm([real(dz); imag(dz)], 2);
    sol_norm(i) = norm(gamma_i, 2);
    mean_resid(i) = mean(abs(dz));
    R_inf_all(i) = R_inf_i;

    fprintf('lambda=%9.3e | ||res||2=%10.4f | ||gamma||2=%10.4f | mean resid=%10.4f\n', ...
        lam, res_norm(i), sol_norm(i), mean_resid(i));
end

kappa = zeros(n_l, 1);
log_res = log10(max(res_norm, eps));
log_sol = log10(max(sol_norm, eps));
log_lam = log10(lambda_values);

for i = 2:n_l-1
    dt1 = log_lam(i) - log_lam(i-1);
    dt2 = log_lam(i+1) - log_lam(i);
    dt = log_lam(i+1) - log_lam(i-1);

    if dt1 <= 0 || dt2 <= 0 || dt <= 0
        continue;
    end

    dx = (log_res(i+1) - log_res(i-1)) / dt;
    dy = (log_sol(i+1) - log_sol(i-1)) / dt;

    d2x = 2 * ((log_res(i+1)-log_res(i))/dt2 - (log_res(i)-log_res(i-1))/dt1) / dt;
    d2y = 2 * ((log_sol(i+1)-log_sol(i))/dt2 - (log_sol(i)-log_sol(i-1))/dt1) / dt;

    denom = (dx*dx + dy*dy)^(3/2);
    if denom > 0
        kappa(i) = abs(dx * d2y - dy * d2x) / denom;
    end
end

[~, corner_idx] = max(kappa);
lambda_corner = lambda_values(corner_idx);

fprintf('\nL-curve corner index : %d\n', corner_idx);
fprintf('L-curve corner lambda: %.6e\n', lambda_corner);
fprintf('Corner ||res||2      : %.6f\n', res_norm(corner_idx));
fprintf('Corner ||gamma||2    : %.6f\n', sol_norm(corner_idx));

figure('Color', 'w', 'Name', 'L-curve (S2422)');
loglog(res_norm, sol_norm, 'bo-'); hold on;
plot(res_norm(corner_idx), sol_norm(corner_idx), 'rs', 'MarkerSize', 9, 'LineWidth', 1.5);

% Label each point with its lambda value
% Use adaptive positioning to keep labels close to data points
for i = 1:numel(lambda_values)
    lambda_str = sprintf('%.0e', lambda_values(i));
    
    % Determine label position (alternate above/below for minimal spacing)
    if mod(i, 2) == 0
        offset_x = 1.02;  % Slightly right
        offset_y = 0.98;  % Slightly below
        h_align = 'left';
        v_align = 'top';
    else
        offset_x = 1.02;  % Slightly right
        offset_y = 1.02;  % Slightly above
        h_align = 'left';
        v_align = 'bottom';
    end
    
    text(res_norm(i) * offset_x, sol_norm(i) * offset_y, lambda_str, ...
        'FontSize', 6, 'HorizontalAlignment', h_align, 'VerticalAlignment', v_align);
end

xlabel('||Z_{exp} - Z_{fit}||_2');
ylabel('||gamma||_2');
title('L-curve for S2422 DRT inversion (lambda values labeled)');
grid on;
legend('L-curve points', 'Corner (max curvature)', 'Location', 'best');
saveas(gcf, output_lcurve_png);

figure('Color', 'w', 'Name', 'L-curve curvature (S2422)');
semilogx(lambda_values, kappa, 'm.-'); hold on;
plot(lambda_corner, kappa(corner_idx), 'ks', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('lambda');
ylabel('Curvature');
title('L-curve curvature vs lambda (S2422)');
grid on;
legend('Curvature', 'Selected corner', 'Location', 'best');
saveas(gcf, output_curvature_png);

details = [lambda_values, res_norm, sol_norm, mean_resid, R_inf_all, kappa];
summary = [lambda_corner, res_norm(corner_idx), sol_norm(corner_idx), mean_resid(corner_idx), R_inf_all(corner_idx), kappa(corner_idx)];

xlswrite(output_xlsx, {'lambda','res_norm_l2','gamma_norm_l2','mean_abs_residual','R_inf','curvature'}, 'details', 'A1');
xlswrite(output_xlsx, details, 'details', 'A2');

xlswrite(output_xlsx, {'corner_lambda','corner_res_norm_l2','corner_gamma_norm_l2','corner_mean_abs_residual','corner_R_inf','corner_curvature'}, 'summary', 'A1');
xlswrite(output_xlsx, summary, 'summary', 'A2');

fprintf('L-curve outputs written to: %s\n', output_xlsx);
fprintf('Saved plot: %s\n', output_lcurve_png);
fprintf('Saved plot: %s\n', output_curvature_png);
end

function [gamma, R_inf] = tr_drt_lcurve_local(freq_vec, Z_exp, lam)
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

function Z_cal = calculate_eis_local(freq_vec, gamma, R_inf)
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
