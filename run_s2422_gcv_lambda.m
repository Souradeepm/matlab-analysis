function run_s2422_gcv_lambda()
% run_s2422_gcv_lambda
% Generalized Cross Validation (GCV) lambda selection for S2422 DRT inversion.
%
% Mathematical model:
%   Re(Z) = R_inf + A_re * gamma
%   Im(Z) = A_im * gamma
%
% Strategy used here:
%   1) Invert gamma from Im(Z) only with nonnegative Tikhonov regularization.
%   2) Estimate R_inf from the real-channel mean residual.
%   3) Score each lambda using normalized GCV on BOTH Re and Im residuals.
%
% nGCV(lambda) = (||r(lambda)||_2^2 / m) / (1 - trace(H_lambda)/n)^2
%   where r stacks [Re residual; Im residual], m = 2*n, and
%   H_lambda is the ridge hat matrix of the imaginary inversion operator.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
out_xlsx = fullfile(repo_root, 'gcv_lambda_s2422_result.xlsx');
output_gcv_png = fullfile(repo_root, 'gcv_s2422_plot.png');
output_gcv_score_png = fullfile(repo_root, 'gcv_s2422_score.png');
output_gcv_deconv_png = fullfile(repo_root, 'gcv_s2422_deconv_overlay.png');
output_gcv_slope_png = fullfile(repo_root, 'gcv_s2422_slope.png');

% Gaussian-RBF peak fitting controls.
rbf_cfg.sigma_frac = 0.03;  % Gaussian RBF width as a fraction of the uniform grid.
rbf_cfg.alpha = 1e-3;       % Smoothness regularization for the coefficient fit.

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

% Remove high-frequency inductive tail (Imag(Z) > 0) before inversion.
[freq_vec, Z_exp, n_removed_ind] = remove_inductive_effects_local(freq_vec, Z_exp);
if n_removed_ind > 0
    fprintf('Removed %d high-frequency inductive point(s) before fitting.\n', n_removed_ind);
end

n_data = length(Z_exp);
tau = 1 ./ (2 * pi * freq_vec);
m_data = 2 * n_data;

% Lambda grid (log-uniform): 1e-6 to 1 with 30 points.
lambda_values = logspace(-6, 0, 30).';
n_l = numel(lambda_values);

res_norm = zeros(n_l, 1);
sol_norm = zeros(n_l, 1);
gcv_score = zeros(n_l, 1);
mean_resid = zeros(n_l, 1);
R_inf_all = zeros(n_l, 1);
hat_trace_all = zeros(n_l, 1);

% Real-only inversion pathway (for GCV comparison overlay).
res_norm_real = zeros(n_l, 1);
sol_norm_real = zeros(n_l, 1);
gcv_score_real = zeros(n_l, 1);
mean_resid_real = zeros(n_l, 1);
R_inf_all_real = zeros(n_l, 1);
hat_trace_all_real = zeros(n_l, 1);

fprintf('Running GCV lambda scan with %d lambda values...\n', n_l);
fprintf('GCV Strategy: Invert with IMAGINARY only, score with equal-weight REAL+IMAGINARY\n\n');
fprintf('Additional pathway: Invert with REAL only, score with equal-weight REAL+IMAGINARY\n\n');

% Precompute discretized DRT kernels once and reuse for all lambda values.
A_re = calc_A_re_gcv_local(freq_vec);
A_im = calc_A_im_gcv_local(freq_vec);

for i = 1:n_l
    lam = lambda_values(i);
    
    % Inversion uses IMAGINARY channel only:
    %   gamma = argmin_{gamma>=0} ||A_im*gamma - Im(Z)||_2^2 + (lambda/2)||gamma||_2^2
    [gamma_i, R_inf_i] = tr_drt_gcv_impart_local(Z_exp, A_re, A_im, lam);
    
    % Build fitted real/imaginary responses from the recovered gamma and R_inf.
    Z_re_fit = R_inf_i + A_re * gamma_i;
    Z_im_fit = A_im * gamma_i;
    
    % Residual vector used by GCV:
    %   r = [Re(Z_exp)-Re(Z_fit); Im(Z_exp)-Im(Z_fit)]
    dz_re = real(Z_exp) - Z_re_fit;
    dz_im = imag(Z_exp) - Z_im_fit;
    residual_norm_sq = norm([dz_re; dz_im], 2)^2;
    
    res_norm(i) = sqrt(residual_norm_sq);
    sol_norm(i) = norm(gamma_i, 2);
    mean_resid(i) = mean(abs([dz_re; dz_im]));
    R_inf_all(i) = R_inf_i;
    
    % Hat matrix complexity term for imaginary ridge inversion:
    %   H = A_im * inv(A_im'*A_im + (lambda/2)I) * A_im'
    tr_h = hat_trace_impart_local(A_im, lam);
    hat_trace_all(i) = tr_h;
    
    % Normalized GCV using combined residual and hat-matrix complexity penalty
    denom = 1 - tr_h / n_data;
    if abs(denom) > 1e-6
        gcv_score(i) = (residual_norm_sq / m_data) / (denom^2);
    else
        gcv_score(i) = inf;
    end
    
    fprintf('lambda=%10.3e | res_cplx=%.4e | sol=%.4e | tr(H)=%.3f | nGCV=%.4e | R_inf=%.6f\n', ...
        lam, res_norm(i), sol_norm(i), tr_h, gcv_score(i), R_inf_i);

    % Real-only inversion pathway (for GCV overlay comparison).
    [gamma_r, R_inf_r] = tr_drt_gcv_repart_local(Z_exp, A_re, lam);
    Z_re_fit_r = R_inf_r + A_re * gamma_r;
    Z_im_fit_r = A_im * gamma_r;
    dz_re_r = real(Z_exp) - Z_re_fit_r;
    dz_im_r = imag(Z_exp) - Z_im_fit_r;
    residual_norm_sq_r = norm([dz_re_r; dz_im_r], 2)^2;

    res_norm_real(i) = sqrt(residual_norm_sq_r);
    sol_norm_real(i) = norm(gamma_r, 2);
    mean_resid_real(i) = mean(abs([dz_re_r; dz_im_r]));
    R_inf_all_real(i) = R_inf_r;

    tr_h_r = hat_trace_repart_local(A_re, lam);
    hat_trace_all_real(i) = tr_h_r;

    denom_r = 1 - tr_h_r / n_data;
    if abs(denom_r) > 1e-6
        gcv_score_real(i) = (residual_norm_sq_r / m_data) / (denom_r^2);
    else
        gcv_score_real(i) = inf;
    end

    fprintf('                      REAL-only | res_cplx=%.4e | sol=%.4e | tr(H)=%.3f | nGCV=%.4e | R_inf=%.6f\n', ...
        res_norm_real(i), sol_norm_real(i), tr_h_r, gcv_score_real(i), R_inf_r);
end

% Select lambda from GCV slope flattening (max lambda in flat region).
log_lambda = log(lambda_values);
slope_seg = diff(log(max(gcv_score, realmin))) ./ diff(log_lambda);  % dlnG/dlnlambda
flat_slope_thresh = 0.002;
flat_seg = abs(slope_seg) <= flat_slope_thresh;

if any(flat_seg)
    flat_idx = find(flat_seg) + 1;
    gcv_idx = flat_idx(end);
    selection_mode = 'max_lambda_on_flat_slope';
else
    [~, gcv_idx] = min(gcv_score);
    selection_mode = 'fallback_min_gcv';
end
lambda_gcv = lambda_values(gcv_idx);

% Keep raw minima for reporting.
[~, gcv_min_idx] = min(gcv_score);
[~, gcv_min_idx_real] = min(gcv_score_real);

fprintf('\n=== GCV Lambda Selection (Imag Inversion, Equal Re+Im Validation) ===\n');
fprintf('selection mode    : %s\n', selection_mode);
fprintf('flat slope thresh : %.6f\n', flat_slope_thresh);
if any(flat_seg)
    fprintf('flat region start : %.6e\n', lambda_values(find(flat_seg, 1, 'first') + 1));
    fprintf('flat region end   : %.6e\n', lambda_values(find(flat_seg, 1, 'last') + 1));
end
fprintf('GCV optimal index : %d\n', gcv_idx);
fprintf('GCV optimal lambda: %.6e\n', lambda_gcv);
fprintf('GCV score         : %.6e\n', gcv_score(gcv_idx));
fprintf('min-GCV index     : %d\n', gcv_min_idx);
fprintf('min-GCV lambda    : %.6e\n', lambda_values(gcv_min_idx));
fprintf('min-GCV score     : %.6e\n', gcv_score(gcv_min_idx));
fprintf('real-only min idx : %d\n', gcv_min_idx_real);
fprintf('real-only min lam : %.6e\n', lambda_values(gcv_min_idx_real));
fprintf('real-only min nGCV: %.6e\n', gcv_score_real(gcv_min_idx_real));
fprintf('trace(H)          : %.6f\n', hat_trace_all(gcv_idx));
fprintf('Complex residual  : %.6f\n', res_norm(gcv_idx));
fprintf('Solution norm     : %.6f\n', sol_norm(gcv_idx));
fprintf('Mean residual     : %.6f\n', mean_resid(gcv_idx));
fprintf('R_inf             : %.6f\n', R_inf_all(gcv_idx));

% Compute final DRT with optimal lambda.
[gamma_opt, R_inf_opt] = tr_drt_gcv_impart_local(Z_exp, A_re, A_im, lambda_gcv);

% For impedance fit visualization, use full impedance
Z_fit_opt = calculate_eis_gcv_local(freq_vec, gamma_opt, R_inf_opt);

% ===== Plot 1: GCV Score vs Lambda =====
fig_score = figure('Color', 'w', 'Name', 'GCV Lambda Selection');
loglog(lambda_values, gcv_score, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
loglog(lambda_values, gcv_score_real, 'm^-', 'LineWidth', 1.3, 'MarkerSize', 5);
plot(lambda_gcv, gcv_score(gcv_idx), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
plot(lambda_values(gcv_min_idx_real), gcv_score_real(gcv_min_idx_real), 'kd', 'MarkerSize', 8, 'LineWidth', 1.5);

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
title(sprintf('Hat-matrix nGCV Overlay (S2422) - Imag selected: %.3e', lambda_gcv));
grid on;
legend('Imag-only nGCV', 'Real-only nGCV', 'Imag selected lambda', 'Real-only min nGCV', 'Location', 'best');
safe_saveas_local(fig_score, output_gcv_score_png);

% ===== Plot 1b: GCV slope dlnG/dlnlambda =====
fig_slope = figure('Color', 'w', 'Name', 'GCV Slope dlnG_dlnlambda');
lambda_seg = lambda_values(2:end);
semilogx(lambda_seg, slope_seg, 'k-o', 'LineWidth', 1.4, 'MarkerSize', 5); hold on;
yline(0, 'Color', [0.45 0.45 0.45], 'LineStyle', '--', 'LineWidth', 1.0);
seg_mask = false(size(lambda_seg));
seg_mask(flat_seg) = true;
if any(seg_mask)
    semilogx(lambda_seg(seg_mask), slope_seg(seg_mask), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
end
xlabel('lambda');
ylabel('dlnG / dlnlambda');
title('GCV slope vs lambda (dlnG/dlnlambda, S2422)');
grid on;
if any(seg_mask)
    legend('dlnG/dlnlambda', 'zero slope', 'flat-segment points', 'Location', 'best');
else
    legend('dlnG/dlnlambda', 'zero slope', 'Location', 'best');
end
safe_saveas_local(fig_slope, output_gcv_slope_png);

% ===== Plot 2: DRT at Optimal Lambda with Gaussian-RBF peak markers =====
fig_drt = figure('Color', 'w', 'Name', 'DRT with GCV Lambda');
semilogx(tau, gamma_opt, 'b-', 'LineWidth', 1.5); hold on;

% Fit a Gaussian-RBF peak model on the log-time axis and extract peak centers.
[pk_tau, pk_g, gamma_peak_scaled] = detect_rbf_peaks_local(tau, gamma_opt, rbf_cfg);
if ~isempty(pk_tau)
    semilogx(tau, gamma_peak_scaled, 'm--', 'LineWidth', 1.2);
    plot(pk_tau, pk_g, 'rv', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    for k = 1:numel(pk_tau)
        text(pk_tau(k), pk_g(k) * 1.05, sprintf('\\tau=%.2e', pk_tau(k)), ...
            'FontSize', 7, 'HorizontalAlignment', 'center');
    end
    n_sig_peaks = numel(pk_tau);
else
    semilogx(tau, gamma_peak_scaled, 'm--', 'LineWidth', 1.2);
    n_sig_peaks = 0;
end

set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('DRT (S2422) - Gaussian-RBF peaks at GCV \\lambda=%.3e - %d peak(s)', lambda_gcv, n_sig_peaks));
grid on;
if ~isempty(pk_tau)
    legend({'DRT', 'Scaled Gaussian-RBF fit', 'RBF peak centers'}, 'Location', 'best');
else
    legend({'DRT', 'Scaled Gaussian-RBF fit'}, 'Location', 'best');
end
safe_saveas_local(fig_drt, output_gcv_png);

% ===== Plot 3: DRT overlaid with Gaussian-RBF fit =====
fig_overlay = figure('Color', 'w', 'Name', 'DRT and Deconvolved Overlay');
semilogx(tau, gamma_opt, 'b-', 'LineWidth', 1.6); hold on;
semilogx(tau, gamma_peak_scaled, 'm--', 'LineWidth', 1.4);
if ~isempty(pk_tau)
    plot(pk_tau, pk_g, 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'y');
end
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('DRT vs Gaussian-RBF Fit Overlay (\\lambda=%.3e)', lambda_gcv));
grid on;
if ~isempty(pk_tau)
    legend({'DRT', 'Scaled Gaussian-RBF fit', 'Detected RBF peaks'}, 'Location', 'best');
else
    legend({'DRT', 'Scaled Gaussian-RBF fit'}, 'Location', 'best');
end
safe_saveas_local(fig_overlay, output_gcv_deconv_png);

% ===== Plot 4: Nyquist Overlay =====
fig_nyquist = figure('Color', 'w', 'Name', 'GCV Nyquist Overlay');
plot(real(Z_exp), -imag(Z_exp), 'bo-', 'LineWidth', 1, 'MarkerSize', 3); hold on;
plot(real(Z_fit_opt), -imag(Z_fit_opt), 'r--', 'LineWidth', 1.5);
axis equal;
xlabel('Re(Z) (Ohm)');
ylabel('-Im(Z) (Ohm)');
title(sprintf('Nyquist Overlay - GCV lambda = %.3e', lambda_gcv));
grid on;
legend('Experimental', 'DRT fit', 'Location', 'best');
safe_saveas_local(fig_nyquist, fullfile(repo_root, 'gcv_s2422_impedance_fit.png'));

% ===== Export Results =====
summary_headers = {'method', 'lambda_opt', 'R_inf', 'residual_norm', 'solution_norm', 'mean_residual', 'gcv_score'};
summary_data = [{'GCV_IMAG_INV_REIM_SCORE'}, num2cell([lambda_gcv, R_inf_opt, res_norm(gcv_idx), sol_norm(gcv_idx), mean_resid(gcv_idx), gcv_score(gcv_idx)])];

xlswrite(out_xlsx, summary_headers, 'summary', 'A1');
xlswrite(out_xlsx, summary_data, 'summary', 'A2');

% Write lambda scan results
scan_headers = {'lambda', 'residual_norm_im', 'solution_norm_im', 'mean_residual_im', 'R_inf_im', 'trace_H_im', 'gcv_score_im', 'residual_norm_re', 'solution_norm_re', 'mean_residual_re', 'R_inf_re', 'trace_H_re', 'gcv_score_re'};
scan_data = [lambda_values, res_norm, sol_norm, mean_resid, R_inf_all, hat_trace_all, gcv_score, res_norm_real, sol_norm_real, mean_resid_real, R_inf_all_real, hat_trace_all_real, gcv_score_real];
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
fprintf('  - gcv_s2422_slope.png (dlnG/dlnlambda vs lambda)\n');
fprintf('  - gcv_s2422_plot.png (DRT with Gaussian-RBF peaks)\n');
fprintf('  - gcv_s2422_deconv_overlay.png (DRT and Gaussian-RBF fit overlay)\n');
fprintf('  - gcv_s2422_impedance_fit.png (Nyquist overlay)\n');

end

% ===== Helper Functions =====

function [gamma, R_inf] = tr_drt_gcv_impart_local(Z_exp, A_re, A_im, lam)
% Solve gamma from the imaginary channel only.
%
% Optimization problem:
%   min_{gamma>=0} ||A_im*gamma - Im(Z_exp)||_2^2 + (lambda/2)||gamma||_2^2
%
% This is implemented via augmented least squares:
%   [A_im           ] gamma ~= [Im(Z_exp)]
%   [sqrt(lambda/2)I]         [0        ]

Z_im = imag(Z_exp(:));
n = numel(Z_im);

% Solve gamma from imaginary channel with Tikhonov regularization.
M = [A_im; sqrt(lam/2) * eye(n)];
b = [Z_im; zeros(n,1)];

x = lsqnonneg(M, b);

gamma = x;

% After gamma is fixed, estimate the constant real offset R_inf by least squares:
%   R_inf = mean(Re(Z_exp) - A_re*gamma)
R_inf = mean(real(Z_exp(:)) - A_re * gamma);
end

function tr_h = hat_trace_impart_local(A_im, lam)
% Compute trace(H_lambda) for imaginary ridge inversion.
%   H_lambda = A_im * inv(A_im'*A_im + (lambda/2)I) * A_im'
alpha = lam / 2;
S = (A_im' * A_im) + alpha * eye(size(A_im, 2));

% Equivalent and numerically stable trace identity:
%   trace(H_lambda) = trace((A_im'*A_im + alpha I)^(-1) * (A_im'*A_im))
tr_h = trace(S \ (A_im' * A_im));

% Keep trace in valid numeric range for denominator stability.
tr_h = max(0, min(tr_h, size(A_im, 1) - 1e-9));
end

function [gamma, R_inf] = tr_drt_gcv_repart_local(Z_exp, A_re, lam)
% Solve gamma from the real channel only.
Z_re = real(Z_exp(:));
n = numel(Z_re);

M = [A_re; sqrt(lam/2) * eye(n)];
b = [Z_re; zeros(n,1)];

x = lsqnonneg(M, b);
gamma = x;
R_inf = mean(real(Z_exp(:)) - A_re * gamma);
end

function tr_h = hat_trace_repart_local(A_re, lam)
% Compute trace(H_lambda) for real-channel ridge inversion.
alpha = lam / 2;
S = (A_re' * A_re) + alpha * eye(size(A_re, 2));
tr_h = trace(S \ (A_re' * A_re));
tr_h = max(0, min(tr_h, size(A_re, 1) - 1e-9));
end

function [freq_out, Z_out, n_removed] = remove_inductive_effects_local(freq_in, Z_in)
% Remove contiguous high-frequency inductive tail where Imag(Z) > 0.
freq_out = freq_in(:);
Z_out = Z_in(:);
n_removed = 0;

if isempty(freq_out)
    return;
end

tail_start = numel(freq_out) + 1;
for k = numel(freq_out):-1:1
    if imag(Z_out(k)) > 0
        tail_start = k;
    else
        break;
    end
end

if tail_start <= numel(freq_out)
    n_removed = numel(freq_out) - tail_start + 1;
    keep_idx = 1:(tail_start - 1);
    freq_out = freq_out(keep_idx);
    Z_out = Z_out(keep_idx);
end

if isempty(freq_out)
    error('All points were removed by inductive-tail filtering. Check data quality.');
end
end

function safe_saveas_local(fig_handle, out_path)
% Save figure robustly across MATLAB graphics-handle modes.
try
    if ishghandle(fig_handle)
        saveas(fig_handle, out_path);
    else
        saveas(gcf, out_path);
    end
catch
    try
        saveas(gcf, out_path);
    catch ME2
        warning('Failed to save figure %s: %s', out_path, ME2.message);
    end
end
end

function [pk_tau, pk_gamma, gamma_peak_scaled] = detect_rbf_peaks_local(tau, gamma, rbf_cfg)
% Gaussian-RBF peak fitting on a uniform log10(tau) grid.
%
% Observation model on log-time x = log10(tau):
%   gamma_obs(x) = (h * s)(x) + noise,
% where h is a Gaussian blur kernel and s is a sparse/nonnegative peak profile.
%
% Gaussian-RBF model:
%   s(x) = sum_j a_j * exp(-(x-c_j)^2 / (2*sigma^2)),   a_j >= 0
%   min_{a>=0} ||Phi*a - gamma_obs||_2^2 + alpha*||D*a||_2^2

% The fitted Gaussian-RBF profile is used only as a smooth surrogate for
% peak identification and overlay visualization.

tau = tau(:);
gamma = gamma(:);

% Sort by log-time ascending for stable interpolation and deconvolution.
x = log10(tau);
[x_sorted, ord] = sort(x, 'ascend');
g_sorted = gamma(ord);
g_sorted = max(g_sorted, 0);

% interp1 with 'pchip' requires strictly increasing sample points.
% Merge duplicate log10(tau) nodes by averaging gamma at repeated points.
[x_unique, ~, grp] = unique(x_sorted);
g_unique = accumarray(grp, g_sorted, [], @mean);

if numel(x_unique) < 2
    % Not enough distinct support points for interpolation/deconvolution.
    pk_tau = [];
    pk_gamma = [];
    gamma_peak_scaled = zeros(size(gamma));
    return;
end

n_uniform = max(300, 4 * numel(g_unique));
xq = linspace(min(x_unique), max(x_unique), n_uniform).';
gq = interp1(x_unique, g_unique, xq, 'pchip');
gq = max(gq, 0);

% Gaussian RBF basis centered on the uniform log-time grid.
% Sigma is measured in log10(tau) units, not sample-count units.
log_span = max(xq) - min(xq);
sigma_rbf = max(3 * median(diff(xq)), rbf_cfg.sigma_frac * log_span);
C = xq - xq.';
Phi = exp(-0.5 * (C / sigma_rbf).^2);
D = diff(eye(n_uniform), 2);
alpha = rbf_cfg.alpha;
M = [Phi; sqrt(alpha) * D];
b = [gq; zeros(size(D, 1), 1)];
a = lsqnonneg(M, b);

% Reconstructed Gaussian-RBF peak profile.
s = Phi * a;

% Local-max detector on the Gaussian-RBF fitted signal (no findpeaks).
if numel(s) >= 3
    cand = find(s(2:end-1) > s(1:end-2) & s(2:end-1) >= s(3:end)) + 1;
else
    cand = [];
end

if isempty(cand) || max(s) <= 0
    pk_tau = [];
    pk_gamma = [];
else
    amp_thresh = 0.05 * max(s);
    cand = cand(s(cand) >= amp_thresh);

    % Enforce a minimum separation in log-time to avoid duplicate micro-peaks.
    min_sep = max(3, round(0.03 * n_uniform));
    keep = select_spaced_peaks_local(cand, s, min_sep);

    x_pk = xq(keep);
    pk_tau = 10 .^ x_pk;
    pk_gamma = interp1(xq, s, x_pk, 'pchip');
end

% Scale fitted profile to DRT amplitude for visual overlay.
s_on_sorted = interp1(xq, s, x_sorted, 'pchip', 0);
if max(s_on_sorted) > 0
    scale = max(g_unique) / max(s_on_sorted);
else
    scale = 1;
end
g_model_sorted = scale * s_on_sorted;
gamma_peak_scaled = zeros(size(gamma));
gamma_peak_scaled(ord) = g_model_sorted;
end

function keep = select_spaced_peaks_local(cand, s, min_sep)
% Greedy non-maximum suppression in index space.
if isempty(cand)
    keep = cand;
    return;
end

[~, order] = sort(s(cand), 'descend');
cand_sorted = cand(order);
keep = zeros(0,1);
for i = 1:numel(cand_sorted)
    idx = cand_sorted(i);
    if isempty(keep) || all(abs(idx - keep) >= min_sep)
        keep(end+1,1) = idx; %#ok<AGROW>
    end
end
keep = sort(keep, 'ascend');
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
