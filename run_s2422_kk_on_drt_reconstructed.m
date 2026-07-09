function run_s2422_kk_on_drt_reconstructed()
% run_s2422_kk_on_drt_reconstructed
% Perform K-K test on both measured impedance and DRT-reconstructed impedance.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 's2422_kk_on_drt.xlsx');
out_png = fullfile(repo_root, 's2422_kk_on_drt.png');

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

lambda_values = [1e-4, 1e-3, 1e-2, 1e-1, 1e0];
best_lambda = lambda_values(1);
best_mean_resid = inf;
best_gamma = [];
best_rinf = 0;
best_z_drt = [];

fprintf('Running S2422 K-K on DRT-reconstructed impedance\n');
fprintf('Scanning lambda values for DRT fit selection:\n');
for i = 1:numel(lambda_values)
    lam = lambda_values(i);
    [gamma_i, rinf_i] = tr_drt_local(freq_vec, Z_exp, lam);
    z_drt_i = calc_eis_local(freq_vec, gamma_i, rinf_i);
    mean_resid_i = mean(abs(Z_exp - z_drt_i));
    fprintf('  lambda=%8.1e | mean residual=%10.6f\n', lam, mean_resid_i);
    if mean_resid_i < best_mean_resid
        best_mean_resid = mean_resid_i;
        best_lambda = lam;
        best_gamma = gamma_i;
        best_rinf = rinf_i;
        best_z_drt = z_drt_i;
    end
end

kk_exp = kk_test_local(freq_vec, Z_exp);
kk_drt = kk_test_local(freq_vec, best_z_drt);

fprintf('\nSelected lambda: %.3e\n', best_lambda);
fprintf('Selected R_inf : %.6f\n', best_rinf);
fprintf('Mean fit resid : %.6f\n', best_mean_resid);

fprintf('\nK-K on measured impedance:\n');
fprintf('  Real residual  : %.4f %%\n', kk_exp.real_residual_pct);
fprintf('  Imag residual  : %.4f %%\n', kk_exp.imag_residual_pct);
fprintf('  Total residual : %.4f %%\n', kk_exp.total_residual_pct);
fprintf('  Status         : %s\n', kk_exp.status);

fprintf('\nK-K on DRT-reconstructed impedance:\n');
fprintf('  Real residual  : %.4f %%\n', kk_drt.real_residual_pct);
fprintf('  Imag residual  : %.4f %%\n', kk_drt.imag_residual_pct);
fprintf('  Total residual : %.4f %%\n', kk_drt.total_residual_pct);
fprintf('  Status         : %s\n', kk_drt.status);

summary_header = {'dataset','selected_lambda','selected_R_inf','drt_fit_mean_residual', ...
    'kk_meas_real_pct','kk_meas_imag_pct','kk_meas_total_pct','kk_meas_max_pct','kk_meas_status_code', ...
    'kk_drt_real_pct','kk_drt_imag_pct','kk_drt_total_pct','kk_drt_max_pct','kk_drt_status_code'};
summary_row = {'S2422', best_lambda, best_rinf, best_mean_resid, ...
    kk_exp.real_residual_pct, kk_exp.imag_residual_pct, kk_exp.total_residual_pct, kk_exp.max_residual_pct, kk_exp.status_code, ...
    kk_drt.real_residual_pct, kk_drt.imag_residual_pct, kk_drt.total_residual_pct, kk_drt.max_residual_pct, kk_drt.status_code};

xlswrite(out_xlsx, summary_header, 'summary', 'A1');
xlswrite(out_xlsx, summary_row, 'summary', 'A2');

detail_header = {'freq_Hz','Zre_input','Zim_input','Zre_kk_fit','Zim_kk_fit','abs_residual','rel_residual_pct'};
exp_detail = [freq_vec, real(Z_exp), imag(Z_exp), real(kk_exp.Z_fit), imag(kk_exp.Z_fit), abs(kk_exp.residual), kk_exp.relative_residual_pct];
drt_detail = [freq_vec, real(best_z_drt), imag(best_z_drt), real(kk_drt.Z_fit), imag(kk_drt.Z_fit), abs(kk_drt.residual), kk_drt.relative_residual_pct];

xlswrite(out_xlsx, detail_header, 'measured_detail', 'A1');
xlswrite(out_xlsx, exp_detail, 'measured_detail', 'A2');
xlswrite(out_xlsx, detail_header, 'drt_reconstructed_detail', 'A1');
xlswrite(out_xlsx, drt_detail, 'drt_reconstructed_detail', 'A2');

figure('Color','w','Name','S2422 K-K measured vs DRT-reconstructed');
subplot(2,1,1);
semilogx(freq_vec, kk_exp.relative_residual_pct, 'k.-'); hold on;
semilogx(freq_vec, kk_drt.relative_residual_pct, 'b.-');
xlabel('Frequency (Hz)');
ylabel('Pointwise residual (%)');
title('K-K relative residual by frequency');
grid on;
legend('Measured input', 'DRT-reconstructed input', 'Location', 'best');

subplot(2,1,2);
plot(real(Z_exp), -imag(Z_exp), 'ko-'); hold on;
plot(real(best_z_drt), -imag(best_z_drt), 'r.-');
xlabel('Re(Z)');
ylabel('-Im(Z)');
title('Nyquist: measured vs DRT reconstructed');
grid on;
legend('Measured', 'DRT reconstructed', 'Location', 'best');

saveas(gcf, out_png);

fprintf('\nSaved: %s\n', out_xlsx);
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

function kk_result = kk_test_local(freq_vec, Z_in)
freq_vec = freq_vec(:);
Z_in = Z_in(:);
omega = 2 * pi * freq_vec;
n_freq = numel(freq_vec);

n_basis = min(max(12, ceil(n_freq / 2)), 30);
tau_min = 1 / (2 * pi * max(freq_vec));
tau_max = 1 / (2 * pi * min(freq_vec));
tau_basis = logspace(log10(tau_min / 5), log10(tau_max * 5), n_basis).';

B_re = zeros(n_freq, n_basis);
B_im = zeros(n_freq, n_basis);
for p = 1:n_freq
    for q = 1:n_basis
        omega_tau = omega(p) * tau_basis(q);
        den = 1 + omega_tau^2;
        B_re(p, q) = 1 / den;
        B_im(p, q) = -omega_tau / den;
    end
end

reg_weight = 1e-4;
M = [B_re, ones(n_freq,1), zeros(n_freq,1); ...
    B_im, zeros(n_freq,1), omega; ...
    sqrt(reg_weight) * eye(n_basis), zeros(n_basis, 2)];
b = [real(Z_in); imag(Z_in); zeros(n_basis,1)];

x = lsqnonneg(M, b);
R_branch = x(1:n_basis);
R_inf = x(n_basis + 1);
L_series = x(n_basis + 2);

Z_fit = R_inf + B_re * R_branch + 1i * (B_im * R_branch + omega * L_series);
residual = Z_in - Z_fit;

real_norm = max(norm(real(Z_in)), eps);
imag_norm = max(norm(imag(Z_in)), eps);
total_norm = max(norm([real(Z_in); imag(Z_in)]), eps);
point_scale = max(abs(Z_in), eps);

real_pct = 100 * norm(real(residual)) / real_norm;
imag_pct = 100 * norm(imag(residual)) / imag_norm;
total_pct = 100 * norm([real(residual); imag(residual)]) / total_norm;
rel_point_pct = 100 * abs(residual) ./ point_scale;
max_pct = max(rel_point_pct);

if total_pct <= 5
    status = 'pass';
    status_code = 1;
elseif total_pct <= 10
    status = 'warning';
    status_code = 0;
else
    status = 'fail';
    status_code = -1;
end

kk_result = struct('tau_basis', tau_basis, 'R_branch', R_branch, 'R_inf', R_inf, ...
    'L_series', L_series, 'Z_fit', Z_fit, 'residual', residual, ...
    'real_residual_pct', real_pct, 'imag_residual_pct', imag_pct, ...
    'total_residual_pct', total_pct, 'max_residual_pct', max_pct, ...
    'relative_residual_pct', rel_point_pct, 'status', status, ...
    'status_code', status_code);
end
