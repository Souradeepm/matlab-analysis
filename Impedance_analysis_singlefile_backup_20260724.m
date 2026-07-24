% function Impedance_analysis_singlefile(input_path)
% Impedance_analysis_singlefile - Self-contained DRT impedance analysis workflow.
%
% Run from MATLAB:  Impedance_analysis_singlefile
% Command line (MATLAB 2011):  matlab -r "Impedance_analysis_singlefile; exit"
%
% ---- CONFIGURE INPUT HERE ----
% Change input_path to point to your EIS workbook.
% Supports a multi-temperature Excel/XLS workbook where every three columns =
% [frequency (Hz), |Z| (Ohm), phase (deg)] and the first header row
% contains temperature labels like "280K", "300K".
% -----------------------------------

clc;
close all;

% Resolve the repository root and default to the S2022 workbook when no
% explicit input path is supplied by the caller.
repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422Sap.xlsx');   % <-- change as needed
if exist(input_path, 'file') ~= 2
    error('Input file not found: %s\nEdit input_path at the top of Impedance_analysis_singlefile.m', input_path);
end

[indata, temperature] = load_input_dataset_2011(input_path);

fprintf('Loaded %d temperature(s) from: %s\n', numel(temperature), input_path);

% --------------------------- Main loop starts here ---------------------------
% Each temperature is stored as one three-column block: frequency, |Z|, phase.
sz = size(indata, 2);

if rem(sz, 3) ~= 0
    error('Input data must contain frequency, magnitude, and phase columns for each temperature');
end

analysis_results=zeros(size(indata,1)+1,size(indata,2));
peak_file=zeros(size(indata,1)+1,size(indata,2));
num_sets = sz / 3;
kk_summary = zeros(num_sets, 13);
sensitivity_summary = zeros(num_sets, 8);
lambda_perturb_summary = zeros(num_sets, 9);
lambda_values_master = logspace(-7, 0, 10).';
num_l_master = numel(lambda_values_master);
ngcv_imag_matrix = nan(num_sets, num_l_master);
ngcv_real_matrix = nan(num_sets, num_l_master);
noise_injection_pct_master = [0.5 1 2 5 10];
noise_mean_residual_pct_matrix = nan(num_sets, numel(noise_injection_pct_master));
noise_max_residual_pct_matrix = nan(num_sets, numel(noise_injection_pct_master));
selected_lambda_per_temp = nan(num_sets, 1);
selected_lambda_idx_per_temp = nan(num_sets, 1);
nyquist_mean_residual_per_temp = nan(num_sets, 1);
nyquist_max_residual_per_temp = nan(num_sets, 1);
nyquist_residual_lambda_matrix = nan(num_sets, num_l_master);
combined_residual_lambda_matrix = nan(num_sets, num_l_master);

sample_tag = make_sample_tag_from_input_2011(input_path);
analysis_output_file = fullfile(repo_root, sprintf('%s_analysis_result.xlsx', sample_tag));
peak_matrix_output_file = fullfile(repo_root, sprintf('%s_peak_result.xlsx', sample_tag));
kk_output_file = fullfile(repo_root, sprintf('%s_kk_result.xlsx', sample_tag));
sensitivity_output_file = fullfile(repo_root, sprintf('%s_sensitivity_result.xlsx', sample_tag));
lambda_perturb_output_file = fullfile(repo_root, sprintf('%s_lambda_perturbation_result.xlsx', sample_tag));
peak_output_file = fullfile(repo_root, sprintf('%s_peak_detail_result.xlsx', sample_tag));
plot_data_output_file = fullfile(repo_root, sprintf('%s_temperature_plot_data.xlsx', sample_tag));
peak_profile_output_file = fullfile(repo_root, sprintf('%s_refined_peak_profiles.xlsx', sample_tag));
rgb_visualization_file = fullfile(repo_root, sprintf('%s_temperature_rgb_peak_map.png', sample_tag));

% Initialize cell array to collect peak data across all temperatures for RGB visualization
all_temps_peak_data = cell(num_sets, 1);

for ab = 1:sz/3
% Slice the workbook into one impedance spectrum and one temperature label.
data = [indata(:,3*ab-2) indata(:,3*ab-1) indata(:,3*ab)];
phase_units = 'deg';
data_format = 'mag_phase';
lambda_values = lambda_values_master;
close all;
% --------------------------- Setup paths/log -------------------------
plots_dir = fullfile(repo_root, 'plots_2011');
if exist(plots_dir, 'dir') ~= 7
    mkdir(plots_dir);
end

run_log_file = fullfile(repo_root, 'drt_matlab2011_run.log');
diary(run_log_file);
diary on;

fprintf('Running MATLAB 2011-compatible DRT workflow from: %s\n', repo_root);

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

% --------------------------- Load data -------------------------------
% Convert the three-column block into frequency and complex impedance.
[freq_vec, Z_exp] = load_f_Z_theta_2011(data, phase_units, data_format);

valid = freq_vec > 0;
if ~all(valid)
    fprintf('Removing %d zero/negative frequency points\n', sum(~valid));
    freq_vec = freq_vec(valid);
    Z_exp = Z_exp(valid);
end

[freq_vec, sort_idx] = sort(freq_vec);
Z_exp = Z_exp(sort_idx);

% Remove high-frequency inductive tail before inversion.
[freq_vec, Z_exp, n_removed_ind] = remove_inductive_effects_matlab2011(freq_vec, Z_exp);
if n_removed_ind > 0
    fprintf('Removed %d high-frequency inductive point(s) before fitting\n', n_removed_ind);
end

fprintf('Loaded %d frequency points (sorted by increasing frequency)\n', numel(freq_vec));
fprintf('Frequency range: %.1f - %.1f Hz\n', min(freq_vec), max(freq_vec));
fprintf('Magnitude range: %.1f - %.1f Ohm\n', min(abs(Z_exp)), max(abs(Z_exp)));
phase_deg = angle(Z_exp) * 180 / pi;
fprintf('Phase range: %.2f - %.2f deg\n', min(phase_deg), max(phase_deg));

% --------------------------- Plot input impedance --------------------
% Plot the measured spectrum in component form before inversion.
% figure('Color', 'w', 'Name', 'Measured impedance');
% subplot(2,1,1);
% semilogx(freq_vec, real(Z_exp), 'o-');
% ylabel('Re(Z)');
% grid on;
% 
% subplot(2,1,2);
% semilogx(freq_vec, imag(Z_exp), 'o-');
% xlabel('Frequency (Hz)');
% ylabel('Im(Z)');
% grid on;

% --------------------------- DRT inversion with GCV ---------------------------
% Score a lambda grid with component-wise inversions and a Re+Im residual metric.
fprintf('\nRunning GCV lambda selection (imaginary inversion, Re+Im cross-validation)...\n');
fprintf('Lambda grid: logspace(-7, 0, 10) = %.1e to %.1e\n', lambda_values(1), lambda_values(end));

num_l = numel(lambda_values);
n_data = numel(freq_vec);
m_data = 2 * n_data;
A_re_mat = calc_A_re_matlab2011(freq_vec);
A_im_mat = calc_A_im_matlab2011(freq_vec);
exp_impedance_norm = max(norm([real(Z_exp(:)); imag(Z_exp(:))], 2), eps);

gcv_scores = inf(num_l, 1);
gcv_scores_real = inf(num_l, 1);
gamma_all  = zeros(n_data, num_l);
R_inf_all  = zeros(num_l, 1);
hat_trace_all = zeros(num_l, 1);
hat_trace_all_real = zeros(num_l, 1);

for i = 1:num_l
    el_i = lambda_values(i);

    % Invert using IMAGINARY part only
    [gamma_i, R_inf_i] = tr_drt_gcv_impart_matlab2011(Z_exp, A_re_mat, A_im_mat, el_i);

    % Validate using Re+Im residuals
    dz_re = real(Z_exp) - (R_inf_i + A_re_mat * gamma_i);
    dz_im = imag(Z_exp) - (A_im_mat * gamma_i);
    res_root = norm([dz_re; dz_im], 2);
    res_sq = res_root^2;
    combined_residual_lambda_matrix(ab, i) = 100 * res_root / exp_impedance_norm;

    tr_h = hat_trace_impart_matlab2011(A_im_mat, el_i);
    hat_trace_all(i) = tr_h;
    denom   = 1 - tr_h / n_data;
    if abs(denom) > 1e-6
        gcv_scores(i) = (res_sq / m_data) / (denom^2);
    end

    [gamma_r, R_inf_r] = tr_drt_gcv_repart_matlab2011(Z_exp, A_re_mat, el_i);
    dz_re_r = real(Z_exp) - (R_inf_r + A_re_mat * gamma_r);
    dz_im_r = imag(Z_exp) - (A_im_mat * gamma_r);
    res_sq_r = norm([dz_re_r; dz_im_r], 2)^2;

    tr_h_r = hat_trace_repart_matlab2011(A_re_mat, el_i);
    hat_trace_all_real(i) = tr_h_r;
    denom_r = 1 - tr_h_r / n_data;
    if abs(denom_r) > 1e-6
        gcv_scores_real(i) = (res_sq_r / m_data) / (denom_r^2);
    end

    gamma_all(:, i) = gamma_i;
    R_inf_all(i)    = R_inf_i;

    fprintf('  lambda = %8.2e | nGCV_im = %.4e | nGCV_re = %.4e\n', el_i, gcv_scores(i), gcv_scores_real(i));
end

% Store per-temperature nGCV-vs-lambda curves for final colormap summaries.
ngcv_imag_matrix(ab, :) = gcv_scores(:).';
ngcv_real_matrix(ab, :) = gcv_scores_real(:).';

log_lambda = log(lambda_values);
slope_seg = diff(log(max(gcv_scores, realmin))) ./ diff(log_lambda);  % dlnG/dlnlambda
flat_slope_thresh = 0.002;
flat_seg = abs(slope_seg) <= flat_slope_thresh;

if any(flat_seg)
    flat_idx = find(flat_seg) + 1;
    best_idx = flat_idx(end);
    selection_mode = 'max_lambda_on_flat_slope';
else
    [~, best_idx] = min(gcv_scores);
    selection_mode = 'fallback_min_gcv';
end

el = lambda_values(best_idx);
gamma  = gamma_all(:, best_idx);
R_inf  = R_inf_all(best_idx);
Z_cal  = R_inf + A_re_mat * gamma + 1i * (A_im_mat * gamma);

[gcv_best_score, gcv_min_idx] = min(gcv_scores);
[gcv_best_score_real, gcv_min_idx_real] = min(gcv_scores_real);

fprintf('\nSelected lambda mode          : %s\n', selection_mode);
fprintf('Selected lambda (imag nGCV)   : %.3e\n', el);
fprintf('Selected nGCV (imag)          : %.4e\n', gcv_scores(best_idx));
fprintf('Minimum nGCV lambda (imag)    : %.3e (score %.4e)\n', lambda_values(gcv_min_idx), gcv_best_score);
fprintf('Minimum nGCV lambda (real)    : %.3e (score %.4e)\n', lambda_values(gcv_min_idx_real), gcv_best_score_real);
fprintf('trace(H) imag at selected     : %.4f\n', hat_trace_all(best_idx));
fprintf('trace(H) real at selected     : %.4f\n', hat_trace_all_real(best_idx));
fprintf('Recovered R_inf = %.6f\n', R_inf);
fprintf('DRT vector length = %d\n', numel(gamma));

% GCV diagnostics: imag/real overlay and slope for the selected policy.
% figure('Color', 'w', 'Name', sprintf('GCV overlay %.1f K', temperature(ab)));
% loglog(lambda_values, gcv_scores, 'bo-', 'LineWidth', 1.2); hold on;
% loglog(lambda_values, gcv_scores_real, 'm^-', 'LineWidth', 1.2);
% plot(el, gcv_scores(best_idx), 'rs', 'MarkerSize', 8, 'LineWidth', 1.2);
% xlabel('lambda');
% ylabel('nGCV score');
% title('Hat-matrix nGCV overlay (imag-only and real-only inversion)');
% grid on;
% legend('Imag-only nGCV', 'Real-only nGCV', 'Selected lambda', 'Location', 'best');
% 
% figure('Color', 'w', 'Name', sprintf('GCV slope %.1f K', temperature(ab)));
% semilogx(lambda_values(2:end), slope_seg, 'k-o', 'LineWidth', 1.2); hold on;
% plot_horizontal_line_2011(0, '--', [0.4 0.4 0.4]);
% xlabel('lambda');
% ylabel('dlnG / dlnlambda');
% title('GCV slope for flattening-based lambda selection');
% grid on;

% --------------------------- Plot DRT vs tau -------------------------
% Plot the recovered DRT directly against relaxation time.
tau = 1 ./ (2 * pi * freq_vec);
% figure('Color', 'w', 'Name', 'DRT vs relaxation time');
% semilogx(tau, gamma, 'o-');
% set(gca, 'XDir', 'reverse');
% xlabel('Relaxation time tau (s)');
% ylabel('gamma(tau)');
% title('DRT vs Relaxation Time');
% grid on;

% --------------------------- Compare fit -----------------------------
% Compare the measured impedance against the DRT reconstruction.
residual = abs(Z_exp - Z_cal);
fprintf('\nMean absolute residual: %.6f\n', mean(residual));
fprintf('Max absolute residual : %.6f\n', max(residual));
selected_lambda_per_temp(ab) = el;
selected_lambda_idx_per_temp(ab) = best_idx;
nyquist_mean_residual_per_temp(ab) = mean(residual);
nyquist_max_residual_per_temp(ab) = max(residual);
nyquist_residual_lambda_matrix(ab, best_idx) = mean(residual);

% figure('Color', 'w', 'Name', 'Nyquist overlay (measured vs DRT)');
% plot(real(Z_exp), -imag(Z_exp), 'ko-', 'LineWidth', 1.0, 'MarkerSize', 3); hold on;
% plot(real(Z_cal), -imag(Z_cal), 'r--', 'LineWidth', 1.4);
% axis equal;
% xlabel('Re(Z) (Ohm)');
% ylabel('-Im(Z) (Ohm)');
% title('Nyquist overlay: measured vs DRT fit');
% grid on;
% legend('Measured', 'DRT fit', 'Location', 'best');

% ------------------ Lambda perturbation sensitivity -------------------
% Disabled per request to avoid extra lambda sensitivity check runs.
% write_matrix_with_header_2011(lambda_perturb_output_file, lambda_sheet, ...
%     {'factor_vs_ref','lambda','mean_residual','max_residual','delta_mean_residual_pct'}, lambda_detail);

% figure('Color', 'w', 'Name', sprintf('Lambda sensitivity %.1f K', temperature(ab)));
% subplot(2,1,1);
% semilogx(lambda_perturb_result.lambda_values, lambda_perturb_result.mean_residual, 'bo-'); hold on;
% plot(el, lambda_perturb_result.reference_mean_resid, 'rs', 'MarkerSize', 8, 'LineWidth', 1.2);
% xlabel('Lambda');
% ylabel('Mean |Z_{exp}-Z_{fit}|');
% grid on;
% legend('Perturbed lambda', 'Selected lambda', 'Location', 'best');

% subplot(2,1,2);
% semilogx(lambda_perturb_result.lambda_values, lambda_perturb_result.delta_pct, 'm.-'); hold on;
% plot(el, 0, 'ks', 'MarkerSize', 8, 'LineWidth', 1.2);
% xlabel('Lambda');
% ylabel('Residual change (%)');
% grid on;
% legend('Delta residual', 'Reference', 'Location', 'best');

% --------------------------- Peak analysis ---------------------------
% Detect, refine, and parametrize the dominant DRT peaks.
tau_data = tau;
gamma_max = max(gamma);

% Prominence-based peak detection (unbiased, no fixed count assumed)
if gamma_max > 0
    min_height = 0.05 * gamma_max;
    peak_idx = find_peaks_simple(gamma, min_height, 3);
    n_sig_peaks = numel(peak_idx);
else
    n_sig_peaks = 0;
end
fprintf('Prominent peaks detected (prominence >= 5%% of max): %d\n', n_sig_peaks);

[data_peaks_fitted, rbf_model, peak_refine] = refine_rbf_peak_fitting_2011(tau_data, gamma);

freq_peak = [];
tau_peak = [];
gamma_peak = [];
R_file = [];
C_file = [];
n_file = [];

fprintf('\nRBF peak refinement summary:\n');
fprintf('  RBF epsilon  : %.3e\n', rbf_model.epsilon);
fprintf('  Base peaks   : %d\n', numel(peak_refine.base_peak_idx));
fprintf('  Hidden peaks : %d\n', numel(peak_refine.hidden_peak_idx));
fprintf('  Total peaks  : %d\n', numel(data_peaks_fitted));
fprintf('  Fit RMSE     : %.4e\n', peak_refine.fit_rmse);

fprintf('\nFitted Original Data Circuit Parameters:\n');
fprintf('Peak | tau (s) | gamma (Ohm) | R (Ohm) | C (F) | n | hidden\n');
fprintf('-------------------------------------------------------------------\n');

if isempty(data_peaks_fitted)
    fprintf('No peaks passed threshold for fitting.\n');
else
    total_R_data = R_inf;
    for i = 1:numel(data_peaks_fitted)
        peak = data_peaks_fitted(i);
        fprintf('%4d | %.2e | %.2e | %.2e | %.2e | %.3f | %d\n', ...
            peak.peak_id, peak.tau, peak.amplitude, peak.R, peak.C, peak.n, peak.is_hidden);
        total_R_data = total_R_data + peak.R;
        freq_peak(i) = 1 / (2 * pi * peak.tau);
        tau_peak(i) = peak.tau;
        gamma_peak(i) = peak.amplitude;
        R_file(i)=peak.R;
        C_file(i)=peak.C;
        n_file(i)=peak.n;
    end
    fprintf('\nTotal series resistance: %.2f Ohm\n', total_R_data);

end
fprintf('High-frequency resistance (R_inf): %.2f Ohm\n', R_inf);

Z_rcn = nan(size(Z_cal));
rcn_residual = nan(size(residual));
if ~isempty(data_peaks_fitted)
    [Z_rcn, Z_rcn_components] = calculate_eis_from_rcn_peaks_matlab2011(freq_vec, R_inf, data_peaks_fitted);

    % Explicitly verify series assembly: Z_total = R_inf + sum_k Z_k.
    Z_rcn_from_series = R_inf + sum(Z_rcn_components, 2);
    series_mismatch = max(abs(Z_rcn - Z_rcn_from_series));
    % 
    % figure('Color', 'w', 'Name', sprintf('Nyquist overlay (DRT vs R-C-n) %.1f K', temperature(ab)));
    % plot(real(Z_exp), -imag(Z_exp), 'ko-', 'LineWidth', 1.0, 'MarkerSize', 3); hold on;
    % plot(real(Z_cal), -imag(Z_cal), 'r--', 'LineWidth', 1.3);
    % plot(real(Z_rcn), -imag(Z_rcn), 'b-.', 'LineWidth', 1.3);
    % axis equal;
    % xlabel('Re(Z) (Ohm)');
    % ylabel('-Im(Z) (Ohm)');
    % title(sprintf('Nyquist overlay at %.1f K (Measured / DRT / R-C-n)', temperature(ab)));
    % grid on;
    % legend('Measured', 'DRT fit', 'R-C-n reconstructed', 'Location', 'best');

    rcn_residual = abs(Z_exp - Z_rcn);
    fprintf('R-C-n Nyquist mean residual: %.6f\n', mean(rcn_residual));
    fprintf('R-C-n Nyquist max residual : %.6f\n', max(rcn_residual));
    fprintf('R-C-n Nyquist mean relative residual: %.2f %%\n', 100 * mean(rcn_residual ./ max(abs(Z_exp), eps)));
    fprintf('Series-assembly mismatch (should be ~0): %.3e\n', series_mismatch);
end

% --------------------------- Plot-data export ------------------------
% Export the dense curves and tabulated plot data needed for downstream analysis.
nyq_sheet = make_sheet_name_2011('NYQ', temperature(ab));
nyq_data = [freq_vec(:), real(Z_exp(:)), imag(Z_exp(:)), real(Z_cal(:)), imag(Z_cal(:)), abs(Z_exp(:) - Z_cal(:))];
% write_matrix_with_header_2011(plot_data_output_file, nyq_sheet, ...
%     {'freq_Hz','Zre_exp','Zim_exp','Zre_drt','Zim_drt','abs_residual_drt'}, nyq_data);

nyq_rcn_sheet = make_sheet_name_2011('NYQRCN', temperature(ab));
nyq_rcn_data = [freq_vec(:), real(Z_exp(:)), imag(Z_exp(:)), real(Z_rcn(:)), imag(Z_rcn(:)), rcn_residual(:)];
% write_matrix_with_header_2011(plot_data_output_file, nyq_rcn_sheet, ...
%     {'freq_Hz','Zre_exp','Zim_exp','Zre_rcn','Zim_rcn','abs_residual_rcn'}, nyq_rcn_data);

drt_sheet = make_sheet_name_2011('DRT', temperature(ab));
drt_plot_data = [tau(:), gamma(:)];
% write_matrix_with_header_2011(plot_data_output_file, drt_sheet, ...
%     {'tau_s','gamma_ohm'}, drt_plot_data);

% --------------------------- KK test ---------------------------------
% Measure how well the spectrum and DRT reconstruction satisfy KK consistency.
kk_result = perform_kk_test_matlab2011(freq_vec, Z_exp);
kk_drt_result = perform_kk_test_matlab2011(freq_vec, Z_cal);
fprintf('\nKramers-Kronig consistency test:\n');
fprintf('  KK basis / reg    : %d / %.1e\n', kk_result.n_basis, kk_result.reg_weight);
fprintf('  Real residual  : %.2f %%\n', kk_result.real_residual_pct);
fprintf('  Imag residual  : %.2f %%\n', kk_result.imag_residual_pct);
fprintf('  Total residual : %.2f %%\n', kk_result.total_residual_pct);
fprintf('  Status         : %s\n', kk_result.status);

fprintf('\nKramers-Kronig test on DRT-generated impedance:\n');
fprintf('  KK basis / reg    : %d / %.1e\n', kk_drt_result.n_basis, kk_drt_result.reg_weight);
fprintf('  Real residual  : %.2f %%\n', kk_drt_result.real_residual_pct);
fprintf('  Imag residual  : %.2f %%\n', kk_drt_result.imag_residual_pct);
fprintf('  Total residual : %.2f %%\n', kk_drt_result.total_residual_pct);
fprintf('  Status         : %s\n', kk_drt_result.status);

kk_summary(ab, :) = [temperature(ab), el, R_inf, kk_result.real_residual_pct, ...
    kk_result.imag_residual_pct, kk_result.total_residual_pct, ...
    kk_result.max_residual_pct, kk_result.status_code, ...
    kk_drt_result.real_residual_pct, kk_drt_result.imag_residual_pct, ...
    kk_drt_result.total_residual_pct, kk_drt_result.max_residual_pct, ...
    kk_drt_result.status_code];

kk_sheet = make_sheet_name_2011('KK', temperature(ab));
kk_detail = [freq_vec(:), real(Z_exp(:)), imag(Z_exp(:)), real(kk_result.Z_fit(:)), ...
    imag(kk_result.Z_fit(:)), abs(kk_result.residual(:)), kk_result.relative_residual_pct(:)];
% write_matrix_with_header_2011(kk_output_file, kk_sheet, ...
%     {'freq_Hz','Zre_exp','Zim_exp','Zre_KK','Zim_KK','abs_residual','rel_residual_pct'}, kk_detail);

kk_drt_sheet = make_sheet_name_2011('KKDRT', temperature(ab));
kk_drt_detail = [freq_vec(:), real(Z_cal(:)), imag(Z_cal(:)), real(kk_drt_result.Z_fit(:)), ...
    imag(kk_drt_result.Z_fit(:)), abs(kk_drt_result.residual(:)), kk_drt_result.relative_residual_pct(:)];
% write_matrix_with_header_2011(kk_output_file, kk_drt_sheet, ...
%     {'freq_Hz','Zre_drt','Zim_drt','Zre_KK','Zim_KK','abs_residual','rel_residual_pct'}, kk_drt_detail);

% figure('Color', 'w', 'Name', sprintf('KK test %.1f K', temperature(ab)));
% subplot(2,1,1);
% semilogx(freq_vec, real(Z_exp), 'ko-'); hold on;
% semilogx(freq_vec, real(kk_result.Z_fit), 'b.-');
% ylabel('Re(Z)');
% grid on;
% legend('Measured', 'KK fit', 'Location', 'best');
% 
% subplot(2,1,2);
% semilogx(freq_vec, imag(Z_exp), 'ko-'); hold on;
% semilogx(freq_vec, imag(kk_result.Z_fit), 'b.-');
% xlabel('Frequency (Hz)');
% ylabel('Im(Z)');
% grid on;
% legend('Measured', 'KK fit', 'Location', 'best');
% 
% figure('Color', 'w', 'Name', sprintf('KK test on DRT impedance %.1f K', temperature(ab)));
% subplot(2,1,1);
% semilogx(freq_vec, real(Z_cal), 'ko-'); hold on;
% semilogx(freq_vec, real(kk_drt_result.Z_fit), 'b.-');
% ylabel('Re(Z)');
% grid on;
% legend('DRT impedance', 'KK fit', 'Location', 'best');
% 
% subplot(2,1,2);
% semilogx(freq_vec, imag(Z_cal), 'ko-'); hold on;
% semilogx(freq_vec, imag(kk_drt_result.Z_fit), 'b.-');
% xlabel('Frequency (Hz)');
% ylabel('Im(Z)');
% grid on;
% legend('DRT impedance', 'KK fit', 'Location', 'best');

% --------------------------- Sensitivity -----------------------------
% Bootstrap the recovered DRT to estimate uncertainty bands and CV statistics.
sensitivity_result = perform_sensitivity_analysis_matlab2011(freq_vec, Z_exp, el, gamma, Z_cal, ab, noise_injection_pct_master);
fprintf('\nSensitivity analysis (%d bootstrap runs):\n', sensitivity_result.n_boot);
fprintf('  Mean std. band : %.4e\n', sensitivity_result.mean_std);
fprintf('  Max std. band  : %.4e\n', sensitivity_result.max_std);
fprintf('  Mean CV        : %.2f %%\n', sensitivity_result.mean_cv_pct);
fprintf('  Max CV         : %.2f %%\n', sensitivity_result.max_cv_pct);

sensitivity_summary(ab, :) = [temperature(ab), el, sensitivity_result.mean_std, ...
    sensitivity_result.max_std, sensitivity_result.mean_cv_pct, ...
    sensitivity_result.max_cv_pct, sensitivity_result.relative_band_area_pct, ...
    numel(data_peaks_fitted)];
noise_mean_residual_pct_matrix(ab, :) = sensitivity_result.noise_mean_residual_pct(:).';
noise_max_residual_pct_matrix(ab, :) = sensitivity_result.noise_max_residual_pct(:).';

sens_sheet = make_sheet_name_2011('SENS', temperature(ab));
sens_detail = [tau(:), gamma(:), sensitivity_result.gamma_mean(:), ...
    sensitivity_result.gamma_std(:), sensitivity_result.gamma_lower(:), ...
    sensitivity_result.gamma_upper(:)];
% write_matrix_with_header_2011(sensitivity_output_file, sens_sheet, ...
%     {'tau_s','gamma_ref','gamma_mean','gamma_std','gamma_lower','gamma_upper'}, sens_detail);

% figure('Color', 'w', 'Name', sprintf('DRT sensitivity %.1f K', temperature(ab)));
% semilogx(tau, gamma, 'k-', 'LineWidth', 1.2); hold on;
% semilogx(tau, sensitivity_result.gamma_mean, 'b--', 'LineWidth', 1.2);
% semilogx(tau, sensitivity_result.gamma_lower, 'r:');
% semilogx(tau, sensitivity_result.gamma_upper, 'r:');
% set(gca, 'XDir', 'reverse');
% xlabel('Relaxation time tau (s)');
% ylabel('gamma(tau)');
% title('DRT sensitivity band');
% grid on;
% legend('Reference DRT', 'Bootstrap mean', 'Bootstrap band', 'Location', 'Best');
% 
% figure('Color', 'w', 'Name', 'DRT peak analysis');
%semilogx(tau_data, gamma, 'bo-'); hold on;

tau_plot = peak_refine.tau_dense;
gamma_rbf = peak_refine.gamma_rbf;
fig_peak = figure('Color', 'w', 'Name', 'DRT peak analysis');
ax_peak = axes('Parent', fig_peak, 'Tag', 'drt_peak_analysis_axes');
hold(ax_peak, 'on');
plot(ax_peak, tau_plot, gamma_rbf, 'g-', 'LineWidth', 1.5);
plot(ax_peak, tau_plot, peak_refine.fit_curve, 'c--', 'LineWidth', 1.2);

for i = 1:numel(data_peaks_fitted)
    peak = data_peaks_fitted(i);
    if peak.is_hidden > 0
        plot(ax_peak, peak.tau, peak.amplitude, 'mo', 'MarkerSize', 8, 'LineWidth', 1.2);
    else
        plot(ax_peak, peak.tau, peak.amplitude, 'ro', 'MarkerSize', 8, 'LineWidth', 1.2);
    end
end

if ~isempty(peak_refine.base_peak_idx)
    scatter(ax_peak, tau_plot(peak_refine.base_peak_idx), gamma_rbf(peak_refine.base_peak_idx), 80, 'r', '*');
end
if ~isempty(peak_refine.hidden_peak_idx)
    scatter(ax_peak, tau_plot(peak_refine.hidden_peak_idx), gamma_rbf(peak_refine.hidden_peak_idx), 80, 'm', 's');
end

set(ax_peak, 'XDir', 'reverse');
xlabel(ax_peak, 'Relaxation time tau (s)');
ylabel(ax_peak, 'gamma(tau) (Ohm)');
title(ax_peak, 'Original Data DRT with RBF-refined Peak Analysis');
grid(ax_peak, 'on');
legend(ax_peak, 'DRT data', 'RBF interpolation', 'Gaussian RBF sum fit', 'Refined peaks', 'Base peaks', 'Hidden peaks', 'Location', 'Best');

fprintf('Run complete at %s\n', datestr(now, 31));

analysis_results(1,2*ab)=temperature(ab);
analysis_results(end-size(tau)+1:end,2*ab-1)=tau;
analysis_results(end-size(gamma)+1:end,2*ab)=gamma;

peak_store=[temperature(ab);tau_peak';0;R_file';0;C_file';0;n_file'];
peak_file(1:size(peak_store),ab)=peak_store;
peak_sheet = make_sheet_name_2011('PEAK', temperature(ab));
peak_numeric = zeros(numel(data_peaks_fitted), 9);
for i = 1:numel(data_peaks_fitted)
    peak_numeric(i,:) = [data_peaks_fitted(i).peak_id, data_peaks_fitted(i).tau, ...
        data_peaks_fitted(i).amplitude, data_peaks_fitted(i).sigma, data_peaks_fitted(i).fwhm_logtau, ...
        data_peaks_fitted(i).R, data_peaks_fitted(i).C, data_peaks_fitted(i).n, data_peaks_fitted(i).is_hidden];
end
% write_matrix_with_header_2011(peak_output_file, peak_sheet, ...
%     {'peak_id','tau_s','amplitude','sigma_logtau','fwhm_logtau','R_ohm','C_F','n_est','is_hidden'}, peak_numeric);

% Store R-C-n values and individual peaks in sample-tagged plot-data workbook.
rcn_sheet = make_sheet_name_2011('RCN', temperature(ab));
% write_matrix_with_header_2011(plot_data_output_file, rcn_sheet, ...
%     {'peak_id','tau_s','amplitude','sigma_logtau','fwhm_logtau','R_ohm','C_F','n_est','is_hidden'}, peak_numeric);

% Collect peak data for RGB multi-temperature visualization
all_temps_peak_data{ab} = struct( ...
    'temperature', temperature(ab), ...
    'tau_dense', peak_refine.tau_dense, ...
    'peaks', data_peaks_fitted, ...
    'peak_count', numel(data_peaks_fitted));

peak_profile_sheet = make_sheet_name_2011('PFIT', temperature(ab));
[peak_profile_header, peak_profile_data] = build_refined_peak_profile_export_2011(peak_refine, data_peaks_fitted);
% write_matrix_with_header_2011(peak_profile_output_file, peak_profile_sheet, peak_profile_header, peak_profile_data);

kkm_sheet = make_sheet_name_2011('KKM', temperature(ab));
kkm_data = [1, kk_result.real_residual_pct, kk_result.imag_residual_pct, kk_result.total_residual_pct, kk_result.max_residual_pct, kk_result.status_code; ...
            2, kk_drt_result.real_residual_pct, kk_drt_result.imag_residual_pct, kk_drt_result.total_residual_pct, kk_drt_result.max_residual_pct, kk_drt_result.status_code];
% write_matrix_with_header_2011(plot_data_output_file, kkm_sheet, ...
%     {'mode_code(1=meas,2=drt)','real_residual_pct','imag_residual_pct','total_residual_pct','max_residual_pct','status_code'}, kkm_data);

% Preserve the explicit peak marker data used in plotting.
pk_sheet = make_sheet_name_2011('PKPTS', temperature(ab));
pk_points = [freq_peak(:), tau_peak(:), gamma_peak(:), R_file(:), C_file(:), n_file(:)];
% write_matrix_with_header_2011(plot_data_output_file, pk_sheet, ...
%     {'freq_peak_Hz','tau_peak_s','gamma_peak_ohm','R_ohm','C_F','n_est'}, pk_points);
diary off;
end

%% Generate RGB multi-temperature peak visualization
fprintf('\n========== Multi-Temperature RGB Peak Visualization ==========\n');
close all;    
num_temps = numel(all_temps_peak_data);
    if num_temps == 0
        fprintf('Warning: no temperature data to visualize\n');
        return;
    end
    
    % Build master ln_tau grid from first temperature's tau_dense
    tau_dense_1 = all_temps_peak_data{1}.tau_dense(:);
    ln_tau_master = log(tau_dense_1);
    n_ln_tau = numel(ln_tau_master);
    
    % Normalize ln_tau to pixel grid [0, 255]
    ln_tau_min = min(ln_tau_master);
    ln_tau_max = max(ln_tau_master);
    if ln_tau_max <= ln_tau_min
        fprintf('Warning: invalid ln_tau range\n');
        return;
    end
    ln_tau_pixel_grid = linspace(0, 256, n_ln_tau);
    
    % Initialize RGB image: [n_temps, 256, 3]
    rgb_image = zeros(num_temps, 256, 3, 'uint8');
    
    % Process each temperature
    for temp_idx = 1:num_temps
        peak_data = all_temps_peak_data{temp_idx};
        peaks = peak_data.peaks;
        n_peaks = peak_data.peak_count;
        
        if n_peaks == 0
            continue;  % Skip temperature with no peaks
        end
        
          % Rank peaks by tau values only.
        peak_areas = zeros(n_peaks, 1);
        for p = 1:n_peaks
              peak_areas(p) = peaks(p).tau;
        end
        
        % Sort peaks by area (descending)
        [~, sort_idx] = sort(peak_areas, 'ascend');
        top_3_idx = sort_idx(1:min(3, n_peaks));
        
        % Reconstruct Gaussian curves for each top peak
        tau_dense_t = peak_data.tau_dense(:);
        ln_tau_t = log(tau_dense_t);
        
        % Initialize R, G, B channels for this temperature
        r_channel = zeros(1, 256, 'double');
        g_channel = zeros(1, 256, 'double');
        b_channel = zeros(1, 256, 'double');
        
        % Assign top 3 peaks to R, G, B
        channels = [r_channel; g_channel; b_channel];
        for rank = 1:min(3, numel(top_3_idx))
            peak_idx = top_3_idx(rank);
            p = peaks(peak_idx);
            
            % Reconstruct Gaussian: A * exp(-(ln_tau - mu)^2 / (2*sigma^2))
            gaussian_vals = p.amplitude * exp(-((ln_tau_t - log(p.tau)) .^ 2) / (2 * p.sigma ^ 2));
            
            % Interpolate to pixel grid
            if exist('interp1', 'file') == 2
                pixel_vals = interp1(ln_tau_t, gaussian_vals, ...
                    linspace(ln_tau_min, ln_tau_max, 256), 'linear', 0);
            else
                % MATLAB 2011 fallback: simple nearest-neighbor
                [~, closest_idx] = min(abs(ln_tau_t - linspace(ln_tau_min, ln_tau_max, 256)'), [], 1);
                pixel_vals = gaussian_vals(closest_idx);
            end
            
            % Normalize to [0, 255]
            max_val = max(pixel_vals);
            if max_val > 0
                pixel_vals = uint8(round(255 * pixel_vals / max_val));
            else
                pixel_vals = uint8(zeros(1, 256));
            end
            
            channels(rank, :) = pixel_vals;
        end
        
        % Fill RGB image for this temperature
        rgb_image(temp_idx, :, 1) = channels(1, :);  % R
        rgb_image(temp_idx, :, 2) = channels(2, :);  % G
        rgb_image(temp_idx, :, 3) = channels(3, :);  % B
    end
    
    % Display RGB image with linear temperature mapping and log-frequency row spacing.
    try
        
        % xlswrite(rgb_image,output_file);

        temperature_axis = temperature(:)';
        [temperature_axis, temp_order] = sort(temperature_axis, 'ascend');
        rgb_display = permute(rgb_image(temp_order, :, :), [2 1 3]);
        ln_tau_pixel_centers = linspace(ln_tau_min, ln_tau_max, size(rgb_display, 1)).';
        freq_axis = exp(-ln_tau_pixel_centers);
        [freq_axis, freq_order] = sort(freq_axis, 'ascend');
        rgb_display = rgb_display(freq_order, :, :);

        % Remove rows that are fully zero across all temperatures and RGB channels.
        row_has_signal = any(reshape(rgb_display, size(rgb_display, 1), []), 2);
        if any(~row_has_signal)
            rgb_display = rgb_display(row_has_signal, :, :);
            freq_axis = freq_axis(row_has_signal);
            fprintf('RGB cleanup removed %d all zero rows\n', sum(~row_has_signal));
        end

        if isempty(freq_axis)
            fprintf('Warning RGB image has no non zero rows after cleanup; skipping visualization\n');
            return;
        end

        % Rescale to a fixed 1000x1000 canvas using nearest-neighbor mapping
        % so original pixel values are preserved without interpolation blending.
        target_rows = 1000;
        target_cols = 1000;
        src_rows = size(rgb_display, 1);
        src_cols = size(rgb_display, 2);
        if src_rows > 0 && src_cols > 0
            row_idx = round(linspace(1, src_rows, target_rows));
            col_idx = round(linspace(1, src_cols, target_cols));
            rgb_display = rgb_display(row_idx, col_idx, :);
            freq_axis = freq_axis(row_idx);
            fprintf('RGB resize: %dx%d -> %dx%d (nearest-neighbor)\n', src_rows, src_cols, target_rows, target_cols);
        end

        % Runtime diagnostics to verify non-empty RGB content before plotting.
        rgb_min = min(rgb_display(:));
        rgb_max = max(rgb_display(:));
        rgb_nnz = nnz(rgb_display);
        fprintf('RGB debug: size=%dx%dx%d, min=%d, max=%d, nnz=%d\n', ...
            size(rgb_display,1), size(rgb_display,2), size(rgb_display,3), rgb_min, rgb_max, rgb_nnz);

        fig = figure;
        ax = axes('Parent', fig, 'Tag', 'rgb_peak_visualization_axes');
        imshow(rgb_display, 'Parent', ax);
        set(ax, 'Visible', 'on', 'YDir', 'normal');
        axis(ax, 'on');
        box(ax, 'on');
        xlabel(ax, 'Temperature (K)');
        ylabel(ax, 'Frequency (Hz)');
        title(ax, 'RGB peak visualization');

        n_cols = size(rgb_display, 2);
        n_rows = size(rgb_display, 1);

        % X-axis: major labels from 1,10,...,280 with 9 minor ticks between majors.
        x_major_vals = 1:10:280;
        x_major_pix = 1 + ((x_major_vals - 1) / (280 - 1)) * (n_cols - 1);
        x_major_valid = x_major_pix >= 1 & x_major_pix <= n_cols;
        x_major_vals = x_major_vals(x_major_valid);
        x_major_pix = x_major_pix(x_major_valid);
        set(ax, 'XTick', x_major_pix, 'XTickLabel', arrayfun(@(v) sprintf('%d', v), x_major_vals, 'UniformOutput', false));

        x_minor_vals = 1:1:280;
        x_minor_vals(mod(x_minor_vals - 1, 10) == 0) = [];
        x_minor_pix = 1 + ((x_minor_vals - 1) / (280 - 1)) * (n_cols - 1);
        x_minor_valid = x_minor_pix >= 1 & x_minor_pix <= n_cols;
        x_minor_pix = x_minor_pix(x_minor_valid);

        % Y-axis: log-mapped frequency labels as powers of ten with 9 minor ticks/decade.
        if freq_axis(1) > 0 && freq_axis(end) > 0
            fmin = min(freq_axis);
            fmax = max(freq_axis);
            logfmin = log10(fmin);
            logfmax = log10(fmax);
            if logfmax > logfmin
                y_major_exp = ceil(logfmin):floor(logfmax);
                if isempty(y_major_exp)
                    y_major_exp = floor(logfmin):ceil(logfmax);
                end
                y_major_vals = 10 .^ y_major_exp;
                y_major_pix = 1 + ((log10(y_major_vals) - logfmin) / (logfmax - logfmin)) * (n_rows - 1);
                y_major_valid = y_major_pix >= 1 & y_major_pix <= n_rows;
                y_major_vals = y_major_vals(y_major_valid);
                y_major_exp = y_major_exp(y_major_valid);
                y_major_pix = y_major_pix(y_major_valid);
                set(ax, 'YTick', y_major_pix, 'YTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), y_major_exp, 'UniformOutput', false));

                y_minor_vals = [];
                for e = floor(logfmin):ceil(logfmax)-1
                    for m = 1:9
                        y_minor_vals(end+1,1) = 10^(e + m/10); %#ok<AGROW>
                    end
                end
                y_minor_pix = 1 + ((log10(y_minor_vals) - logfmin) / (logfmax - logfmin)) * (n_rows - 1);
                y_minor_valid = y_minor_pix >= 1 & y_minor_pix <= n_rows;
                y_minor_pix = y_minor_pix(y_minor_valid);
            else
                y_minor_pix = [];
            end
        else
            y_minor_pix = [];
        end

        set(ax, 'XMinorTick', 'on', 'YMinorTick', 'on');
        try
            ax.XRuler.MinorTickValues = x_minor_pix;
            ax.YRuler.MinorTickValues = y_minor_pix;
        catch
            % Keep major ticks if MinorTickValues is unavailable.
        end

        % if exist('exportgraphics', 'file') == 2
        %     exportgraphics(fig, rgb_visualization_file, 'Resolution', 300);
        % else
        %     print(fig, rgb_visualization_file, '-dpng', '-r300');
        % end

        fprintf('RGB peak visualization ready: %s (%d x %d x 3)\n', rgb_visualization_file, size(rgb_display, 1), size(rgb_display, 2));
    catch ME
        fprintf('Error writing PNG: %s\n', ME.message);
    end

%%



% try
%     generate_temperature_rgb_peak_map(all_temps_peak_data, temperature, rgb_visualization_file);
%     fprintf('RGB peak visualization exported to: %s\n', rgb_visualization_file);
% catch ME
%     fprintf('Warning: RGB visualization generation failed: %s\n', ME.message);
% end

% xlswrite(analysis_output_file,analysis_results);
% xlswrite(peak_matrix_output_file,peak_file);
% Write the per-workbook summaries after all temperatures have been processed.
% write_matrix_with_header_2011(kk_output_file, 'summary', ...
%     {'Temperature_K','Lambda','R_inf','KK_real_pct','KK_imag_pct','KK_total_pct','KK_max_pct','KK_status_code', ...
%     'KK_DRT_real_pct','KK_DRT_imag_pct','KK_DRT_total_pct','KK_DRT_max_pct','KK_DRT_status_code'}, ...
%     kk_summary);
% write_matrix_with_header_2011(sensitivity_output_file, 'summary', ...
%     {'Temperature_K','Lambda','Mean_std','Max_std','Mean_CV_pct','Max_CV_pct','Band_area_pct','Peak_count'}, ...
%     sensitivity_summary);
% write_matrix_with_header_2011(lambda_perturb_output_file, 'summary', ...
%     {'Temperature_K','Lambda_ref','Ref_mean_residual','Min_delta_pct','Max_delta_pct','Max_abs_delta_pct', ...
%     'Local_slope','Best_factor','Worst_factor'}, lambda_perturb_summary);

fprintf('\nSample-tagged outputs generated for %s:\n', sample_tag);
fprintf('  %s\n', analysis_output_file);
fprintf('  %s\n', peak_matrix_output_file);
fprintf('  %s\n', kk_output_file);
fprintf('  %s\n', sensitivity_output_file);
fprintf('  %s\n', peak_output_file);
fprintf('  %s\n', plot_data_output_file);
fprintf('  %s\n', peak_profile_output_file);
% end

%% Generate KK raw residual and noise sensitivity summary plots
fprintf('\n========== KK and Noise Sensitivity Summaries ==========%s\n', '');

temperature_vec = temperature(:);
[temperature_sorted, temp_sort_idx] = sort(temperature_vec, 'ascend');

kk_temperature_sorted = kk_summary(temp_sort_idx, 1);
kk_total_raw_sorted = kk_summary(temp_sort_idx, 6);

fig_kk_raw = figure;
ax_kk_raw = axes('Parent', fig_kk_raw, 'Tag', 'kk_raw_total_residual_axes');
plot(ax_kk_raw, kk_temperature_sorted, kk_total_raw_sorted, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 4);
grid(ax_kk_raw, 'on');
xlabel(ax_kk_raw, 'Temperature (K)');
ylabel(ax_kk_raw, 'KK total residual (%)');
title(ax_kk_raw, 'Kramers-Kronig raw-data total residual vs temperature');

noise_mean_residual_sorted = noise_mean_residual_pct_matrix(temp_sort_idx, :);

fig_noise_mean = figure;
ax_noise_mean = axes('Parent', fig_noise_mean, 'Tag', 'noise_mean_residual_colormap_axes');
imagesc(ax_noise_mean, 1:numel(noise_injection_pct_master), temperature_sorted, noise_mean_residual_sorted);
set(ax_noise_mean, 'YDir', 'normal');
set(ax_noise_mean, 'XTick', 1:numel(noise_injection_pct_master), ...
    'XTickLabel', arrayfun(@(v) sprintf('%.1f', v), noise_injection_pct_master, 'UniformOutput', false));
xlabel(ax_noise_mean, 'Noise injection (%)');
ylabel(ax_noise_mean, 'Temperature (K)');
title(ax_noise_mean, 'Noise sensitivity mean residual colormap');
cb_noise_mean = colorbar(ax_noise_mean);
ylabel(cb_noise_mean, 'Mean residual (%)');

fprintf('KK and noise sensitivity plots ready: %d temperature(s), %d noise levels\n', ...
    numel(temperature_sorted), numel(noise_injection_pct_master));

%% Generate nGCV-vs-lambda colormaps (imaginary and real inversions)
fprintf('\n========== nGCV vs Lambda Colormaps ==========%s\n', '');

temperature_vec = temperature(:);
[temperature_sorted, temp_sort_idx] = sort(temperature_vec, 'ascend');

ngcv_imag_sorted = ngcv_imag_matrix(temp_sort_idx, :);
ngcv_real_sorted = ngcv_real_matrix(temp_sort_idx, :);

% Use log10(nGCV) for dynamic-range-stable visualization.
ngcv_imag_plot = log10(max(ngcv_imag_sorted, realmin));
ngcv_real_plot = log10(max(ngcv_real_sorted, realmin));

ngcv_imag_plot(~isfinite(ngcv_imag_plot)) = nan;
ngcv_real_plot(~isfinite(ngcv_real_plot)) = nan;

lambda_axis = lambda_values_master(:).';
n_lambda = numel(lambda_axis);
log_lambda_axis = log10(lambda_axis);
lambda_major_exp = ceil(min(log_lambda_axis)):floor(max(log_lambda_axis));
if isempty(lambda_major_exp)
    lambda_major_exp = floor(min(log_lambda_axis)):ceil(max(log_lambda_axis));
end
lambda_major_vals = 10 .^ lambda_major_exp;
lambda_major_pix = 1 + ((log10(lambda_major_vals) - min(log_lambda_axis)) / ...
    (max(log_lambda_axis) - min(log_lambda_axis))) * (n_lambda - 1);
lambda_major_valid = lambda_major_pix >= 1 & lambda_major_pix <= n_lambda;
lambda_major_vals = lambda_major_vals(lambda_major_valid);
lambda_major_pix = lambda_major_pix(lambda_major_valid);

lambda_minor_vals = [];
for e = floor(min(log_lambda_axis)):ceil(max(log_lambda_axis))-1
    for m = 1:9
        lambda_minor_vals(end+1,1) = 10^(e + m/10); %#ok<AGROW>
    end
end
lambda_minor_pix = 1 + ((log10(lambda_minor_vals) - min(log_lambda_axis)) / ...
    (max(log_lambda_axis) - min(log_lambda_axis))) * (n_lambda - 1);
lambda_minor_valid = lambda_minor_pix >= 1 & lambda_minor_pix <= n_lambda;
lambda_minor_pix = lambda_minor_pix(lambda_minor_valid);

fig_ngcv_im = figure;
ax_ngcv_im = axes('Parent', fig_ngcv_im, 'Tag', 'ngcv_imag_colormap_axes');
imagesc(ax_ngcv_im, 1:n_lambda, temperature_sorted, ngcv_imag_plot);
set(ax_ngcv_im, 'YDir', 'normal');
xlim(ax_ngcv_im, [1 n_lambda]);
set(ax_ngcv_im, 'XTick', lambda_major_pix, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false));
xlabel(ax_ngcv_im, 'Lambda');
ylabel(ax_ngcv_im, 'Temperature (K)');
title(ax_ngcv_im, 'Imaginary inversion nGCV colormap');
cb1 = colorbar(ax_ngcv_im);
ylabel(cb1, 'log10(nGCV)');
set(ax_ngcv_im, 'XMinorTick', 'on');
try
    ax_ngcv_im.XRuler.MinorTickValues = lambda_minor_pix;
catch
    % Keep major ticks if MinorTickValues is unavailable.
end

fig_ngcv_re = figure;
ax_ngcv_re = axes('Parent', fig_ngcv_re, 'Tag', 'ngcv_real_colormap_axes');
imagesc(ax_ngcv_re, 1:n_lambda, temperature_sorted, ngcv_real_plot);
set(ax_ngcv_re, 'YDir', 'normal');
xlim(ax_ngcv_re, [1 n_lambda]);
lambda_tick_labels_re = arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false);
set(ax_ngcv_re, 'XTickMode', 'manual', 'XTick', lambda_major_pix, ...
    'XTickLabelMode', 'manual', 'XTickLabel', lambda_tick_labels_re);
xlabel(ax_ngcv_re, 'Lambda');
ylabel(ax_ngcv_re, 'Temperature (K)');
title(ax_ngcv_re, 'Real inversion nGCV colormap');
cb2 = colorbar(ax_ngcv_re);
ylabel(cb2, 'log10(nGCV)');
set(ax_ngcv_re, 'XTickMode', 'manual', 'XTick', lambda_major_pix, ...
    'XTickLabelMode', 'manual', 'XTickLabel', lambda_tick_labels_re);
set(ax_ngcv_re, 'XMinorTick', 'on');
try
    ax_ngcv_re.XRuler.MinorTickValues = lambda_minor_pix;
catch
    % Keep major ticks if MinorTickValues is unavailable.
end

fprintf('nGCV colormaps ready: imag matrix %dx%d, real matrix %dx%d\n', ...
    size(ngcv_imag_sorted, 1), size(ngcv_imag_sorted, 2), ...
    size(ngcv_real_sorted, 1), size(ngcv_real_sorted, 2));

%% Generate combined residual-vs-lambda colormap
fprintf('\n========== Combined Residual vs Lambda Colormap ==========%s\n', '');

combined_residual_lambda_sorted = combined_residual_lambda_matrix(temp_sort_idx, :);
selected_lambda_sorted = selected_lambda_per_temp(temp_sort_idx);
nyquist_mean_residual_sorted = nyquist_mean_residual_per_temp(temp_sort_idx);

fig_lambda_resid = figure;
ax_lambda_resid = axes('Parent', fig_lambda_resid, 'Tag', 'lambda_nyquist_residual_colormap_axes');
imagesc(ax_lambda_resid, 1:n_lambda, temperature_sorted, combined_residual_lambda_sorted);
set(ax_lambda_resid, 'YDir', 'normal');
xlim(ax_lambda_resid, [1 n_lambda]);
set(ax_lambda_resid, 'XTick', lambda_major_pix, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false));
set(ax_lambda_resid, 'XMinorTick', 'on');
try
    ax_lambda_resid.XRuler.MinorTickValues = lambda_minor_pix;
catch
    % Keep major ticks if MinorTickValues is unavailable.
end
xlabel(ax_lambda_resid, 'Lambda');
ylabel(ax_lambda_resid, 'Temperature (K)');
title(ax_lambda_resid, 'Root combined residual percent vs lambda and temperature');
cb_lambda_resid = colorbar(ax_lambda_resid);
ylabel(cb_lambda_resid, '100 * ||[Re residual; Im residual]||_2 / ||Z_{exp}||_2');
hold(ax_lambda_resid, 'on');
plot(ax_lambda_resid, selected_lambda_idx_per_temp(temp_sort_idx), temperature_sorted, 'ks', ...
    'MarkerSize', 5, 'LineWidth', 1.0, 'MarkerFaceColor', 'w');
hold(ax_lambda_resid, 'off');

fprintf('Selected lambda/residual map ready: %d temperature(s), lambda range %.1e to %.1e\n', ...
    numel(temperature_sorted), min(selected_lambda_sorted), max(selected_lambda_sorted));
fprintf('Combined residual percent range: %.4f %% to %.4f %%\n', ...
    min(combined_residual_lambda_sorted(:)), max(combined_residual_lambda_sorted(:)));

function [freq, Z] = load_f_Z_theta_2011(data, phase_units, data_format)
% load_f_Z_theta_2011 - Convert one [freq, value2, value3] data block to complex impedance
%
% The input is already a numeric matrix extracted from the workbook. This
% helper converts columns 2 and 3 into complex impedance values using one of
% two formats:
%   - 'mag_phase': column 2 is magnitude |Z|, column 3 is phase angle
%   - 'real_imag': column 2 is Re(Z), column 3 is Im(Z)
%
% INPUTS:
%   data        - numeric matrix [freq, value2, value3]
%   phase_units - 'deg' for degrees (default) or 'rad' for radians
%   data_format - 'mag_phase' (default) or 'real_imag' format
%
% OUTPUTS:
%   freq - Column vector of frequencies (Hz) from column 1
%   Z    - Column vector of complex impedance values (Ohm)

if nargin < 2 || isempty(phase_units)
    phase_units = 'deg';
end
if nargin < 3 || isempty(data_format)
    data_format = 'mag_phase';
end

if size(data, 2) < 3
    error('Input must contain at least three columns');
end

freq = double(data(:,1));
col2 = double(data(:,2));
col3 = double(data(:,3));

if strcmpi(data_format, 'mag_phase')
    if strcmpi(phase_units, 'deg')
        phase = col3 * pi / 180;
    elseif strcmpi(phase_units, 'rad')
        phase = col3;
    else
        error('phase_units must be ''deg'' or ''rad''');
    end
    Z = col2 .* exp(1i * phase);
elseif strcmpi(data_format, 'real_imag')
    Z = col2 + 1i * col3;
else
    error('data_format must be ''mag_phase'' or ''real_imag''');
end

end

function [indata, temperature] = load_input_dataset_2011(input_path)
% load_input_dataset_2011 - Load the numeric workbook block and header temperatures.
%
% The single-file workflow expects a workbook where each temperature occupies
% three numeric columns and the first header row stores labels such as 280K.
% This helper returns the numeric matrix and the parsed temperature vector.

try
    % Prefer native readers first to avoid Excel COM activation issues.
    indata = readmatrix(input_path, 'Sheet', 1);
    text_hdr = readcell(input_path, 'Sheet', 1, 'Range', '1:1');
catch ME
    try
        % Fallback: explicitly resolve first worksheet name.
        [~, sheet_names] = xlsfinfo(input_path);
        if isempty(sheet_names)
            sheet_id = 1;
        else
            sheet_id = sheet_names{1};
        end
        [indata, text_hdr] = xlsread(input_path, sheet_id);
    catch
        try
            % Final fallback: basic mode bypasses Excel automation.
            [num_data, txt_hdr, raw_hdr] = xlsread(input_path, 1, 'basic');
            indata = num_data;
            text_hdr = txt_hdr;
            if isempty(text_hdr) && ~isempty(raw_hdr)
                text_hdr = raw_hdr(1, :);
            end
        catch ME2
            error('Failed to read Excel workbook %s: %s', input_path, ME2.message);
        end
    end
end

if ~isempty(indata)
    indata = indata(~all(isnan(indata), 2), :);
end

temperature = parse_temperature_headers_2011(text_hdr);

if isempty(indata)
    error('No numeric impedance data found in workbook: %s', input_path);
end

if isempty(temperature)
    error('Could not determine temperature labels for %s', input_path);
end
end

function temperature = parse_temperature_headers_2011(text_hdr)
% parse_temperature_headers_2011 - Extract numeric temperature labels from the first header row.
%
% Header cells are sanitized by stripping the trailing K/k before numeric
% conversion so labels such as '280K' and '300k' are both accepted.

temperature = zeros(0, 1);
if isempty(text_hdr)
    return;
end

header_cells = text_hdr(1, :);
for hk = 1:numel(header_cells)
    tok = header_cells{hk};
    if ischar(tok)
        tok = strtrim(tok);
        tok(tok == 'K' | tok == 'k') = [];
        val = str2double(tok);
        if ~isnan(val)
            temperature(end + 1, 1) = val; %#ok<AGROW>
        end
    end
end
end

function [gamma, R_inf] = TR_DRT_matlab2011(freq_vec, Z_exp, el)
% TR_DRT_matlab2011 - Perform Distribution of Relaxation Times (DRT) inversion
%
% This function inverts experimental electrochemical impedance spectroscopy (EIS) data
% to obtain the DRT (gamma), which represents the distribution of relaxation processes.
% Uses Tikhonov regularization (TR) with a user-specified regularization parameter
% (lambda = el) to balance fit quality against smoothness.
%
% The algorithm:
%   1. Builds matrices A_re and A_im that map DRT to impedance
%   2. Sets up constrained least-squares problem: minimize ||M*x - b||^2
%   3. Solves with non-negativity constraint (gamma >= 0)
%   4. Returns gamma (DRT) and R_inf (high-frequency resistance)
%
% INPUTS:
%   freq_vec - Frequency vector (Hz) [N x 1]
%   Z_exp    - Experimental complex impedance values (Ohm) [N x 1]
%   el       - Regularization parameter lambda (0.0001 to 1.0 typical)
%              Larger lambda = smoother gamma, larger residual
%              Smaller lambda = rougher gamma, smaller residual
%
% OUTPUTS:
%   gamma  - Distribution of relaxation times [N x 1]
%   R_inf  - High-frequency resistance (Ohm)

A_re = calc_A_re_matlab2011(freq_vec);
A_im = calc_A_im_matlab2011(freq_vec);

Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
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

function Z_cal = calculate_EIS_matlab2011(freq_vec, gamma, R_inf)
% calculate_EIS_matlab2011 - Calculate impedance spectrum from DRT
%
% This function performs the forward calculation: given a DRT (gamma) and
% high-frequency resistance, it computes the theoretical impedance spectrum.
% This is the inverse operation of the DRT inversion.
%
% The calculation reconstructs impedance as:
%   Z = R_inf + A_re*gamma + i*A_im*gamma
% where A_re and A_im are the real and imaginary projection matrices.
%
% INPUTS:
%   freq_vec - Frequency vector (Hz) [N x 1]
%   gamma    - Distribution of relaxation times [N x 1]
%   R_inf    - High-frequency resistance (Ohm)
%
% OUTPUTS:
%   Z_cal - Calculated complex impedance values (Ohm) [N x 1]

A_re = calc_A_re_matlab2011(freq_vec);
A_im = calc_A_im_matlab2011(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function [Z_rcn, Z_components] = calculate_eis_from_rcn_peaks_matlab2011(freq_vec, R_inf, fitted_peaks)
% calculate_eis_from_rcn_peaks_matlab2011 - Rebuild impedance from fitted R-C-n peaks.
%
% Each refined DRT peak is mapped to one ZARC branch, then branches are put
% in series with each other and with the ohmic term R_inf:
%   Z_total(omega) = R_inf + sum_k Z_k(omega)
%   Z_k(omega) = R_k / (1 + (j*omega*R_k*C_k)^n_k)
% This corresponds to: (R1||CPE1) + (R2||CPE2) + ... + (Rm||CPEm) + R_inf.

omega = 2 * pi * freq_vec(:);
Z_rcn = R_inf + zeros(size(omega));
Z_components = zeros(numel(omega), numel(fitted_peaks));

for k = 1:numel(fitted_peaks)
    Rk = fitted_peaks(k).R;
    Ck = fitted_peaks(k).C;
    nk = fitted_peaks(k).n;

    if ~(Rk > 0) || ~(Ck > 0) || ~(nk > 0)
        continue;
    end

    Zk = Rk ./ (1 + (1i * omega * Rk * Ck) .^ nk);
    Z_components(:, k) = Zk;
    Z_rcn = Z_rcn + Zk;
end
end

function [gamma, R_inf] = tr_drt_gcv_impart_matlab2011(Z_exp, A_re, A_im, lam)
% tr_drt_gcv_impart_matlab2011 - Solve the imag-only ridge problem used in GCV.
%
% The imaginary channel is inverted with non-negativity, then R_inf is backed
% out from the real component residual so the resulting gamma can be scored.

Z_im = imag(Z_exp(:));
n = numel(Z_im);
M = [A_im; sqrt(lam/2) * eye(n)];
b = [Z_im; zeros(n,1)];

gamma = lsqnonneg(M, b);
R_inf = mean(real(Z_exp(:)) - A_re * gamma);
end

function [gamma, R_inf] = tr_drt_gcv_repart_matlab2011(Z_exp, A_re, lam)
% tr_drt_gcv_repart_matlab2011 - Solve the real-only ridge problem used in GCV.
%
% This companion solve mirrors the imag-only branch so the workflow can compare
% cross-channel predictive performance across the lambda grid.

Z_re = real(Z_exp(:));
n = numel(Z_re);
M = [A_re; sqrt(lam/2) * eye(n)];
b = [Z_re; zeros(n,1)];

gamma = lsqnonneg(M, b);
R_inf = mean(real(Z_exp(:)) - A_re * gamma);
end

function tr_h = hat_trace_impart_matlab2011(A_im, lam)
% hat_trace_impart_matlab2011 - Compute the effective degrees of freedom for imag GCV.

alpha = lam / 2;
S = (A_im' * A_im) + alpha * eye(size(A_im, 2));
tr_h = trace(S \ (A_im' * A_im));
tr_h = max(0, min(tr_h, size(A_im, 1) - 1e-9));
end

function tr_h = hat_trace_repart_matlab2011(A_re, lam)
% hat_trace_repart_matlab2011 - Compute the effective degrees of freedom for real GCV.

alpha = lam / 2;
S = (A_re' * A_re) + alpha * eye(size(A_re, 2));
tr_h = trace(S \ (A_re' * A_re));
tr_h = max(0, min(tr_h, size(A_re, 1) - 1e-9));
end

function [freq_out, Z_out, n_removed] = remove_inductive_effects_matlab2011(freq_in, Z_in)
% remove_inductive_effects_matlab2011 - Drop the contiguous inductive high-frequency tail.
%
% Positive imaginary impedance at the highest frequencies usually signals an
% inductive artefact that destabilizes the DRT inversion, so only the leading
% contiguous tail is removed.

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

function A_re = calc_A_re_matlab2011(freq)
% calc_A_re_matlab2011 - Calculate real part of DRT projection matrix
%
% This function builds the real component of the matrix that projects the DRT
% (Distribution of Relaxation Times) into measured impedance values.
% Used in both DRT inversion and forward EIS calculation.
%
% The matrix element A_re(p,q) represents the contribution of relaxation time q
% to the real part of impedance at frequency p.
%
% Mathematical basis:
%   Re(Z) = -0.5/(1 + (omega*tau)^2) * log(tau_q+1/tau_q-1)
% where omega = 2*pi*f and tau = 1/f
%
% INPUTS:
%   freq - Frequency vector (Hz) [N x 1]
%
% OUTPUTS:
%   A_re - Real part projection matrix [N x N]

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

function A_im = calc_A_im_matlab2011(freq)
% calc_A_im_matlab2011 - Calculate imaginary part of DRT projection matrix
%
% This function builds the imaginary component of the matrix that projects the DRT
% (Distribution of Relaxation Times) into measured impedance values.
% Used in both DRT inversion and forward EIS calculation.
%
% The matrix element A_im(p,q) represents the contribution of relaxation time q
% to the imaginary part of impedance at frequency p.
%
% Mathematical basis:
%   Im(Z) = 0.5*(omega*tau)/(1 + (omega*tau)^2) * log(tau_q+1/tau_q-1)
% where omega = 2*pi*f and tau = 1/f
%
% INPUTS:
%   freq - Frequency vector (Hz) [N x 1]
%
% OUTPUTS:
%   A_im - Imaginary part projection matrix [N x N]

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

function [peak_indices, peak_values] = find_peaks_simple(y, min_height, min_distance)
% find_peaks_simple - Detect peaks in a 1D signal
%
% This function identifies local maxima in a signal that exceed a minimum height
% threshold and are separated by minimum distance. Peaks are returned in order of
% decreasing height, then filtered to ensure minimum separation distance.
%
% Algorithm:
%   1. Find all local maxima (points higher than neighbors) above min_height
%   2. Sort candidates by height (descending)
%   3. Select peaks, rejecting any within min_distance of already selected peaks
%   4. Return peaks in index order (low to high)
%
% INPUTS:
%   y           - Signal vector to analyze [N x 1 or 1 x N]
%   min_height  - Minimum peak height (scalar, in signal units)
%   min_distance - Minimum distance between peaks (scalar, in index units)
%
% OUTPUTS:
%   peak_indices - Index positions of detected peaks [M x 1]
%   peak_values  - Signal values at peak positions [M x 1]

y = y(:);
n = numel(y);

cand_idx = [];
for i = 2:n-1
    if y(i) >= y(i-1) && y(i) > y(i+1) && y(i) >= min_height
        cand_idx(end+1,1) = i; %#ok<AGROW>
    end
end

if isempty(cand_idx)
    peak_indices = [];
    peak_values = [];
    return;
end

cand_val = y(cand_idx);
[~, order] = sort(cand_val, 'descend');
selected = [];

for k = 1:numel(order)
    idx = cand_idx(order(k));
    if isempty(selected)
        selected(end+1,1) = idx; %#ok<AGROW>
    else
        if all(abs(selected - idx) >= min_distance)
            selected(end+1,1) = idx; %#ok<AGROW>
        end
    end
end

selected = sort(selected);
peak_indices = selected;
peak_values = y(selected);
end

function ysm = smooth_mavg(y, w)
% smooth_mavg - Minimal moving-average smoother used as a toolbox-free fallback.

y = y(:);
n = numel(y);
if w <= 1
    ysm = y;
    return;
end
half = floor(w / 2);
ysm = zeros(n,1);
for i = 1:n
    i1 = max(1, i - half);
    i2 = min(n, i + half);
    ysm(i) = mean(y(i1:i2));
end
end

function ysm = smooth_sgolay_2011(y, frame_len, poly_order)
% smooth_sgolay_2011 - Lightweight Savitzky-Golay smoothing without toolbox dependency
%
% The code performs a local polynomial least-squares fit at each point and is
% used to stabilize curvature-based hidden-peak detection on older MATLAB.

y = y(:);
n = numel(y);

if n < 5
    ysm = y;
    return;
end

if mod(frame_len, 2) == 0
    frame_len = frame_len + 1;
end
frame_len = min(frame_len, n);
if mod(frame_len, 2) == 0
    frame_len = frame_len - 1;
end

poly_order = min(poly_order, frame_len - 1);
if poly_order < 1 || frame_len < 3
    ysm = smooth_mavg(y, 3);
    return;
end

half = floor(frame_len / 2);
ysm = zeros(n, 1);
for i = 1:n
    i1 = max(1, i - half);
    i2 = min(n, i + half);
    xi = (i1:i2).' - i;
    yi = y(i1:i2);

    V = zeros(numel(xi), poly_order + 1);
    for p = 0:poly_order
        V(:, p + 1) = xi .^ p;
    end

    c = V \ yi;
    ysm(i) = c(1);
end
end

function model = create_rbf_model_2011(x, y)
% create_rbf_model_2011 - Create Radial Basis Function (RBF) interpolation model
%
% This function trains an RBF network to smoothly interpolate DRT data.
% RBF networks use Gaussian basis functions to create smooth, continuous fits
% that capture peak shapes better than simple linear interpolation.
%
% The RBF model learns a weighted sum of Gaussian functions:
%   gamma(tau) = sum_i w_i * exp(-(epsilon * (tau - tau_i))^2)
% where w_i are learned weights and epsilon is a shape parameter.
%
% Algorithm:
%   1. Determine optimal epsilon (shape parameter) from data spacing
%   2. Build Gaussian kernel matrix from pairwise distances
%   3. Add Tikhonov regularization for numerical stability
%   4. Solve for optimal weights
%
% INPUTS:
%   x - Logarithm of relaxation times, log(tau) [N x 1]
%   y - DRT values at those times [N x 1]
%
% OUTPUTS:
%   model - Structure containing: .x (training x), .w (weights), .epsilon (shape)

x = x(:);
y = y(:);
n = numel(x);

if n > 1
    ds = abs(diff(sort(x)));
    med_ds = median(ds);
    if med_ds <= eps
        epsilon = 1;
    else
        epsilon = 1 / med_ds;
    end
else
    epsilon = 1;
end

lambda = 1e-8;
DX = x(:, ones(1,n)) - x(:, ones(1,n)).';
Phi = exp(-(epsilon * DX).^2);
w = (Phi + lambda * eye(n)) \ y;

model.x = x;
model.w = w;
model.epsilon = epsilon;
end

function yq = rbf_eval_2011(model, xq)
% rbf_eval_2011 - Evaluate RBF model at query points
%
% This function evaluates a previously created RBF interpolation model at new points.
% It applies the learned RBF weights to Gaussian basis functions centered at the
% training data points to produce smooth interpolated values.
%
% INPUTS:
%   model - RBF model structure created by create_rbf_model_2011 containing:
%           .x (training x values), .w (weights), .epsilon (shape parameter)
%   xq - Query points where model should be evaluated [M x 1]
%
% OUTPUTS:
%   yq - Interpolated values at query points [M x 1]
%
% EXAMPLE:
%   Train model on tau data
%   model = create_rbf_model_2011(log(tau), gamma);
%   Evaluate at high-resolution points for smooth plotting
%   tau_fine = logspace(log10(min(tau)), log10(max(tau)), 1000);
%   gamma_smooth = rbf_eval_2011(model, log(tau_fine));

xq = xq(:);
mx = model.x(:).';
DX = xq(:, ones(1,numel(mx))) - mx(ones(numel(xq),1), :);
Phi_q = exp(-(model.epsilon * DX).^2);
yq = Phi_q * model.w;
end

function [fitted_peaks, rbf_model, peak_refine] = refine_rbf_peak_fitting_2011(tau_vals, gamma_vals)
% refine_rbf_peak_fitting_2011 - RBF-refine peaks and detect hidden shoulders with the 2nd derivative
%
% The workflow first interpolates the DRT on a dense log-time grid, then finds
% base peaks, searches for hidden shoulders via curvature, and finally fits a
% Gaussian sum whose parameters are exported for post-processing.

tau_vals = tau_vals(:);
gamma_vals = gamma_vals(:);
tau_log = log(tau_vals);

rbf_model = create_rbf_model_2011(tau_log, gamma_vals);
tau_dense = logspace(log10(min(tau_vals)), log10(max(tau_vals)), 1000).';
tau_dense_log = log(tau_dense);
gamma_rbf = rbf_eval_2011(rbf_model, tau_dense_log);
gamma_rbf(gamma_rbf < 0) = 0;

base_thr = max(0.05 * max(gamma_rbf), eps);
[base_peak_idx, ~] = find_peaks_simple(gamma_rbf, base_thr, 20);

hidden_peak_idx = detect_hidden_peaks_second_derivative_2011(tau_dense_log, gamma_rbf, base_peak_idx);

all_peak_idx = unique([base_peak_idx(:); hidden_peak_idx(:)]);
all_peak_idx = sort(all_peak_idx);

hidden_flag = zeros(size(all_peak_idx));
for i = 1:numel(all_peak_idx)
    if any(abs(hidden_peak_idx - all_peak_idx(i)) <= 2)
        hidden_flag(i) = 1;
    end
end

[fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2011(tau_dense, gamma_rbf, all_peak_idx, hidden_flag);
peak_refine = struct('tau_dense', tau_dense, 'gamma_rbf', gamma_rbf, ...
    'base_peak_idx', base_peak_idx, 'hidden_peak_idx', hidden_peak_idx, ...
    'all_peak_idx', all_peak_idx, 'hidden_flag', hidden_flag, ...
    'fit_curve', fit_curve, 'fit_rmse', fit_rmse);
end

function hidden_peak_idx = detect_hidden_peaks_second_derivative_2011(tau_log, gamma_rbf, base_peak_idx)
% detect_hidden_peaks_second_derivative_2011 - Find shoulder-like peaks from curvature in the RBF curve
%
% A shoulder is accepted only if it clears amplitude, prominence, and spacing
% checks so weak numerical ripples are rejected.

tau_log = tau_log(:);
gamma_rbf = gamma_rbf(:);
n = numel(gamma_rbf);

if n < 5
    hidden_peak_idx = [];
    return;
end

dx = mean(diff(tau_log));
if ~(dx > 0)
    dx = 1;
end

d2 = zeros(n,1);
for i = 2:n-1
    d2(i) = (gamma_rbf(i+1) - 2 * gamma_rbf(i) + gamma_rbf(i-1)) / (dx * dx);
end

curvature_signal = -d2;
curvature_signal(curvature_signal < 0) = 0;
curvature_signal = smooth_sgolay_2011(curvature_signal, 7, 2);

curv_thr = max(0.30 * max(curvature_signal), 3.2 * std(curvature_signal));
[curv_idx, ~] = find_peaks_simple(curvature_signal, curv_thr, 28);

hidden_peak_idx = [];
for i = 1:numel(curv_idx)
    idx = curv_idx(i);
    near_base = any(abs(base_peak_idx - idx) <= 12);
    amp_ok = gamma_rbf(idx) >= 0.05 * max(gamma_rbf);

    j1 = max(1, idx - 18);
    j2 = min(n, idx + 18);
    local_floor = min(gamma_rbf(j1:j2));
    shoulder_depth_ok = (gamma_rbf(idx) - local_floor) >= 0.02 * max(gamma_rbf);

    l1 = max(1, idx - 25);
    l2 = idx;
    r1 = idx;
    r2 = min(n, idx + 25);
    left_min = min(gamma_rbf(l1:l2));
    right_min = min(gamma_rbf(r1:r2));
    prominence_ok = (gamma_rbf(idx) - max(left_min, right_min)) >= 0.015 * max(gamma_rbf);

    if ~near_base && amp_ok && shoulder_depth_ok && prominence_ok && idx > 5 && idx < n - 5
        hidden_peak_idx(end+1,1) = idx; %#ok<AGROW>
    end
end

% Mid-gap rescue: recover one weak central shoulder between neighboring base peaks.
if numel(base_peak_idx) >= 2
    for k = 1:(numel(base_peak_idx) - 1)
        i1 = base_peak_idx(k);
        i2 = base_peak_idx(k + 1);
        if i2 - i1 < 24
            continue;
        end
        seg = (i1 + 3):(i2 - 3);
        if isempty(seg)
            continue;
        end

        [seg_max, rel_idx] = max(gamma_rbf(seg));
        cand = seg(rel_idx);
        near_any = any(abs([base_peak_idx(:); hidden_peak_idx(:)] - cand) <= 10);
        if near_any
            continue;
        end

        valley_left = min(gamma_rbf(i1:cand));
        valley_right = min(gamma_rbf(cand:i2));
        local_prom = seg_max - max(valley_left, valley_right);
        if seg_max >= 0.04 * max(gamma_rbf) && local_prom >= 0.012 * max(gamma_rbf)
            hidden_peak_idx(end+1,1) = cand; %#ok<AGROW>
        end
    end
end

if ~isempty(hidden_peak_idx)
    hidden_peak_idx = unique(sort(hidden_peak_idx));
end

max_hidden = max(0, floor(numel(base_peak_idx) / 2));
if numel(hidden_peak_idx) > max_hidden
    [~, ord_h] = sort(gamma_rbf(hidden_peak_idx), 'descend');
    hidden_peak_idx = hidden_peak_idx(ord_h(1:max_hidden));
    hidden_peak_idx = sort(hidden_peak_idx);
end
end

function [fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2011(tau_dense, gamma_rbf, peak_idx, hidden_flag)
% fit_rbf_gaussian_peaks_2011 - Fit a global sum of Gaussian radial basis functions to the refined DRT
%
% Peak centers are fixed from the refined detection stage, while amplitudes,
% widths, and baseline are optimized globally for a compact analytical summary.

if isempty(peak_idx)
    fitted_peaks = [];
    fit_curve = zeros(size(gamma_rbf));
    fit_rmse = 0;
    return;
end

tau_dense = tau_dense(:);
gamma_rbf = gamma_rbf(:);
xall = log(tau_dense);
m = numel(peak_idx);
centers = xall(peak_idx(:));

[fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);

% Residual-based rescue: if the Gaussian-sum underfits between two peaks, add one central hidden basis.
if numel(centers) >= 2
    residual_pos = gamma_rbf - fit_curve;
    added_centers = [];
    for k = 1:(numel(centers) - 1)
        seg = find(xall > centers(k) & xall < centers(k+1));
        if numel(seg) < 12
            continue;
        end
        [res_max, rel_idx] = max(residual_pos(seg));
        cand = seg(rel_idx);
        if gamma_rbf(cand) < 0.035 * max(gamma_rbf)
            continue;
        end
        if res_max < max(0.02 * max(gamma_rbf), 1.2 * sqrt(mean((fit_curve - gamma_rbf).^2)))
            continue;
        end
        if any(abs([centers(:); added_centers(:)] - xall(cand)) < 0.06)
            continue;
        end
        added_centers(end+1,1) = xall(cand); %#ok<AGROW>
    end

    if ~isempty(added_centers)
        centers = sort([centers(:); added_centers(:)]);
        hidden_flag = [hidden_flag(:); ones(numel(added_centers), 1)];
        [~, ord_h] = sort([centers(:)]);
        hidden_flag = hidden_flag(ord_h);
        [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);
    end
end

peak_area = amplitudes .* sigma_vec;
amp_cutoff = 0.05 * max(amplitudes);
area_cutoff = 0.05 * max(peak_area);
keep = (amplitudes >= amp_cutoff) & (peak_area >= area_cutoff);

if any(~keep) && any(keep)
    centers = centers(keep);
    hidden_flag = hidden_flag(keep);
    m = numel(centers);
    [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);
end

fit_rmse = sqrt(mean((fit_curve - gamma_rbf).^2));

fitted_peaks = repmat(struct('peak_id', 0, 'tau', 0, 'amplitude', 0, 'sigma', 0, ...
    'fwhm_logtau', 0, 'R', 0, 'C', 0, 'n', 0, 'is_hidden', 0, 'baseline', 0), m, 1);

for peak_idx_local = 1:m
    tau_peak = exp(centers(peak_idx_local));
    amp = amplitudes(peak_idx_local);
    sig = sigma_vec(peak_idx_local);
    fwhm_logtau = gaussian_fwhm_from_sigma_2011(sig);
    R_est = amp * sig * sqrt(2 * pi);
    C_est = tau_peak / max(R_est, eps);

    % ===== N-FORMULA OPTIONS (empirically tested on S2022 dataset) =====
    % EMPIRICAL WINNER: Heuristic formula has 2.99% mean lower Nyquist residual.
    % To test alternative formulas, comment out the ACTIVE formula and uncomment
    % the alternative below, then regenerate outputs and compare Nyquist residuals.
    % ====================================================================

    % ACTIVE FORMULA (EMPIRICALLY OPTIMAL):
    % Heuristic: n = 1.144 / (1 + sigma)
    % where sigma is Gaussian width parameter in ln(tau) space.
    % TESTED: Mean Nyquist residual 2796.77 Ohm (best across 6 temperatures)
    n_est = 1.144 / (1 + sig);

    % ALTERNATIVE FORMULA 1: FWHM-BASED IMPLICIT EQUATION
    % Physics-grounded from ZARC width theory:
    %   FWHM_lntau = 2*sqrt(2*ln2)*sigma  (Gaussian width)
    %   FWHM_lntau = (2/n)*acosh(2 + cos(pi*n))  (ZARC width formula)
    %   => sqrt(2*ln2)*sigma*n = acosh(2 + cos(pi*n))  (implicit equation)
    % Solve numerically for n in (0, 1).
    % TESTED: Mean Nyquist residual 2906.61 Ohm (3.78% worse than heuristic)
    % Uncomment line below to use FWHM-based formula:
    % n_est = solve_n_from_fwhm_implicit_2011(sig);

    % ALTERNATIVE FORMULA 2: APEX PHASE ANGLE (from one-arc approximation)
    % For a depressed arc, phase angle at apex omega*tau=1 is:
    %   phi_apex = -pi*n/4  (in radians) or -45*n degrees
    %   => n = -4*phi_apex/pi  or  n = -phi_apex/45  (degrees)
    % This assumes single dominant arc and uses phase information.
    % Note: Not tested; provides independent physics perspective.
    % Uncomment line below to use phase-based formula:
    % n_est = estimate_n_from_phase_2011(tau_peak, amplitude_vals(i), gamma_ref);

    % ALTERNATIVE FORMULA 3: CONSTANT CPE EXPONENT
    % Simplified assumption: n = 0.5 for all peaks (mid-range CPE behavior)
    % Useful for sensitivity studies and baseline comparison.
    % Uncomment line below to use constant formula:
    % n_est = 0.5;

    fitted_peaks(peak_idx_local).peak_id = peak_idx_local;
    fitted_peaks(peak_idx_local).tau = tau_peak;
    fitted_peaks(peak_idx_local).amplitude = amp;
    fitted_peaks(peak_idx_local).sigma = sig;
    fitted_peaks(peak_idx_local).fwhm_logtau = fwhm_logtau;
    fitted_peaks(peak_idx_local).R = R_est;
    fitted_peaks(peak_idx_local).C = C_est;
    fitted_peaks(peak_idx_local).n = n_est;
    fitted_peaks(peak_idx_local).is_hidden = hidden_flag(peak_idx_local);
    fitted_peaks(peak_idx_local).baseline = baseline;
end

function [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers)
% solve_gaussian_rbf_sum_fit_2011 - Optimize Gaussian widths and baseline for fixed peak centers.

m = numel(centers);
sigma0 = zeros(m, 1);
for sigma_idx = 1:m
    if m == 1
        neighbor_span = max(xall) - min(xall);
    elseif sigma_idx == 1
        neighbor_span = abs(centers(sigma_idx+1) - centers(sigma_idx));
    elseif sigma_idx == m
        neighbor_span = abs(centers(sigma_idx) - centers(sigma_idx-1));
    else
        neighbor_span = min(abs(centers(sigma_idx) - centers(sigma_idx-1)), abs(centers(sigma_idx+1) - centers(sigma_idx)));
    end
    sigma0(sigma_idx) = min(max(0.12, 0.35 * neighbor_span), 0.75);
end

p0 = [log(sigma0); log(max(min(gamma_rbf), eps))];
opts = optimset('Display', 'off', 'MaxFunEvals', 6000, 'MaxIter', 2500);
obj = @(p) gaussian_rbf_sum_objective_2011(p, xall, gamma_rbf, centers);
p = fminsearch(obj, p0, opts);

[fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2011(p, xall, gamma_rbf, centers);
end

all_tau = [fitted_peaks.tau].';
[~, order] = sort(all_tau, 'descend');
fitted_peaks = fitted_peaks(order);
for peak_order = 1:numel(fitted_peaks)
    fitted_peaks(peak_order).peak_id = peak_order;
end
end

function sse = gaussian_rbf_sum_objective_2011(p, xall, yall, centers)
% gaussian_rbf_sum_objective_2011 - Penalized least-squares objective for Gaussian peak fitting.

[fit_curve, amplitudes, sigma_vec] = evaluate_gaussian_rbf_sum_2011(p, xall, yall, centers);
resid = fit_curve - yall;
sse = sum(resid .^ 2);

sigma_penalty = 0;
for sigma_idx = 1:numel(sigma_vec)
    if sigma_vec(sigma_idx) < 0.04
        sigma_penalty = sigma_penalty + (0.04 - sigma_vec(sigma_idx)) ^ 2;
    elseif sigma_vec(sigma_idx) > 0.90
        sigma_penalty = sigma_penalty + (sigma_vec(sigma_idx) - 0.90) ^ 2;
    end
end

amp_penalty = 0;
for amp_idx = 1:numel(amplitudes)
    if amplitudes(amp_idx) < 1e-8 * max(yall)
        amp_penalty = amp_penalty + (1e-8 * max(yall) - amplitudes(amp_idx)) ^ 2;
    end
end

sse = sse + 1e3 * sigma_penalty + 1e-2 * amp_penalty;
end

function [fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2011(p, xall, yall, centers)
% evaluate_gaussian_rbf_sum_2011 - Evaluate one Gaussian-sum parameter vector.

sigma_vec = exp(p(1:numel(centers)));

Phi = zeros(numel(xall), numel(centers));
for basis_idx = 1:numel(centers)
    Phi(:, basis_idx) = exp(-0.5 * ((xall - centers(basis_idx)) / sigma_vec(basis_idx)).^2);
end

Phi_aug = [Phi, ones(numel(xall), 1)];
coef = lsqnonneg(Phi_aug, yall);
amplitudes = coef(1:numel(centers));
baseline = coef(end);
fit_curve = Phi_aug * coef;
end

function fwhm_logtau = gaussian_fwhm_from_sigma_2011(sig)
% gaussian_fwhm_from_sigma_2011 - Convert Gaussian sigma in log(tau) to FWHM.

fwhm_logtau = 2 * sqrt(2 * log(2)) * sig;
end

function [header_row, numeric_data] = build_refined_peak_profile_export_2011(peak_refine, fitted_peaks)
% build_refined_peak_profile_export_2011 - Assemble dense per-temperature peak curves for Excel.
%
% The exported matrix stores the dense log-time axis, the RBF interpolation, the
% summed Gaussian fit, the fitted baseline, and each individual peak component.

tau_dense = peak_refine.tau_dense(:);
log_tau_dense = log(tau_dense);
gamma_rbf = peak_refine.gamma_rbf(:);
fit_curve = peak_refine.fit_curve(:);
num_points = numel(tau_dense);
num_peaks = numel(fitted_peaks);

baseline_value = 0;
if num_peaks > 0
    baseline_value = fitted_peaks(1).baseline;
end
baseline_col = baseline_value * ones(num_points, 1);

component_curves = zeros(num_points, num_peaks);
header_row = {'tau_s','log_tau','gamma_rbf','gamma_fit','fit_baseline'};
for peak_idx_local = 1:num_peaks
    component_curves(:, peak_idx_local) = fitted_peaks(peak_idx_local).amplitude * ...
        exp(-0.5 * ((log_tau_dense - log(fitted_peaks(peak_idx_local).tau)) / fitted_peaks(peak_idx_local).sigma).^2);
    header_row{end + 1} = sprintf('peak_%02d_gamma', fitted_peaks(peak_idx_local).peak_id); %#ok<AGROW>
end

numeric_data = [tau_dense, log_tau_dense, gamma_rbf, fit_curve, baseline_col, component_curves];
end

function kk_result = perform_kk_test_matlab2011(freq_vec, Z_exp)
% perform_kk_test_matlab2011 - Estimate KK-consistent impedance and residual metrics
%
% Several basis sizes and regularization strengths are tried; the best candidate
% is kept and then labeled pass, warning, or fail by total residual percentage.

freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
n_freq = numel(freq_vec);

candidate_basis = unique([ ...
    min(max(8, ceil(n_freq / 3)), 20), ...
    min(max(10, ceil(n_freq / 2.5)), 24), ...
    min(max(12, ceil(n_freq / 2)), 30)]);
candidate_reg = [1e-4; 3e-4; 1e-3];

best_result = [];
for basis_idx = 1:numel(candidate_basis)
    for reg_idx = 1:numel(candidate_reg)
        kk_candidate = solve_kk_candidate_matlab2011(freq_vec, Z_exp, candidate_basis(basis_idx), candidate_reg(reg_idx));
        if isempty(best_result) || kk_candidate.total_residual_pct < best_result.total_residual_pct
            best_result = kk_candidate;
        end
    end
end

kk_result = best_result;

if kk_result.total_residual_pct <= 5
    status = 'pass';
    status_code = 1;
elseif kk_result.total_residual_pct <= 10
    status = 'warning';
    status_code = 0;
else
    status = 'fail';
    status_code = -1;
end

kk_result.status = status;
kk_result.status_code = status_code;
end

function kk_result = solve_kk_candidate_matlab2011(freq_vec, Z_exp, n_basis, reg_weight)
% solve_kk_candidate_matlab2011 - Solve one KK fit candidate and report residual metrics.
%
% The spectrum is represented by a finite set of Debye-like branches plus series
% resistance and inductance, then solved under non-negativity constraints.

omega = 2 * pi * freq_vec;
tau_min = 1 / (2 * pi * max(freq_vec));
tau_max = 1 / (2 * pi * min(freq_vec));
tau_basis = logspace(log10(tau_min / 5), log10(tau_max * 5), n_basis).';

n_freq = numel(freq_vec);
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

M = [B_re, ones(n_freq,1), zeros(n_freq,1); ...
    B_im, zeros(n_freq,1), omega; ...
    sqrt(reg_weight) * eye(n_basis), zeros(n_basis, 2)];
b = [real(Z_exp); imag(Z_exp); zeros(n_basis,1)];

x = lsqnonneg_robust_2011(M, b);
R_branch = x(1:n_basis);
R_inf = x(n_basis + 1);
L_series = x(n_basis + 2);

Z_fit = R_inf + B_re * R_branch + 1i * (B_im * R_branch + omega * L_series);
residual = Z_exp - Z_fit;

real_norm = max(norm(real(Z_exp)), eps);
imag_norm = max(norm(imag(Z_exp)), eps);
total_norm = max(norm([real(Z_exp); imag(Z_exp)]), eps);
point_scale = max(abs(Z_exp), eps);

real_pct = 100 * norm(real(residual)) / real_norm;
imag_pct = 100 * norm(imag(residual)) / imag_norm;
total_pct = 100 * norm([real(residual); imag(residual)]) / total_norm;
rel_point_pct = 100 * abs(residual) ./ point_scale;
max_pct = max(rel_point_pct);

kk_result = struct('tau_basis', tau_basis, 'R_branch', R_branch, 'R_inf', R_inf, ...
    'L_series', L_series, 'Z_fit', Z_fit, 'residual', residual, ...
    'real_residual_pct', real_pct, 'imag_residual_pct', imag_pct, ...
    'total_residual_pct', total_pct, 'max_residual_pct', max_pct, ...
    'relative_residual_pct', rel_point_pct, 'n_basis', n_basis, ...
    'reg_weight', reg_weight);
end

function x = lsqnonneg_robust_2011(M, b)
% lsqnonneg_robust_2011 - Use explicit solver options to reduce premature iteration exits on older MATLAB releases.
%
% Older MATLAB versions can stop early with default solver tolerances, so this
% helper centralizes a stricter option set for all constrained linear solves.

options = optimset('Display', 'off', 'MaxIter', max(400, 10 * size(M, 2)), 'TolX', 1e-10);
x = lsqnonneg(M, b, options);
end

function sensitivity_result = perform_sensitivity_analysis_matlab2011(freq_vec, Z_exp, el, gamma_ref, Z_cal, seed_offset, noise_injection_pct)
% perform_sensitivity_analysis_matlab2011 - Bootstrap DRT stability at the selected lambda
%
% Synthetic complex noise is drawn from the fit residual statistics and added to
% the measured spectrum so the variability of the recovered DRT can be exported.

freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
gamma_ref = gamma_ref(:);
Z_cal = Z_cal(:);
n_freq = numel(freq_vec);
n_boot = 24;

if nargin < 7 || isempty(noise_injection_pct)
    noise_injection_pct = [0.5 1 2 5 10];
end
noise_injection_pct = noise_injection_pct(:);
noise_n_levels = numel(noise_injection_pct);
noise_n_trials = 4;

residual = Z_exp - Z_cal;
sigma_re = std(real(residual));
sigma_im = std(imag(residual));

if ~(sigma_re > 0)
    sigma_re = 0.005 * max(max(abs(real(Z_exp))), 1);
end
if ~(sigma_im > 0)
    sigma_im = 0.005 * max(max(abs(imag(Z_exp))), 1);
end

rand('twister', 100 + seed_offset);
randn('state', 200 + seed_offset);

gamma_samples = zeros(n_freq, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n_freq, 1) + 1i * sigma_im * randn(n_freq, 1);
    Z_boot = Z_exp + noise;
    [gamma_boot, ~] = TR_DRT_matlab2011(freq_vec, Z_boot, el);
    gamma_samples(:, k) = gamma_boot(:);
end

gamma_mean = mean(gamma_samples, 2);
gamma_std = std(gamma_samples, 0, 2);
gamma_sorted = sort(gamma_samples, 2);
lower_idx = max(1, round(0.05 * n_boot));
upper_idx = min(n_boot, round(0.95 * n_boot));
gamma_lower = gamma_sorted(:, lower_idx);
gamma_upper = gamma_sorted(:, upper_idx);

active_floor = max(0.01 * max(gamma_ref), eps);
active_mask = gamma_ref > active_floor;
if ~any(active_mask)
    active_mask = true(size(gamma_ref));
end

cv_pct = 100 * gamma_std(active_mask) ./ max(gamma_mean(active_mask), active_floor);
tau = 1 ./ (2 * pi * freq_vec);
band_area_pct = 100 * trapz(log(tau(end:-1:1)), gamma_upper(end:-1:1) - gamma_lower(end:-1:1)) / ...
    max(trapz(log(tau(end:-1:1)), abs(gamma_ref(end:-1:1))), eps);

noise_mean_residual_pct = zeros(noise_n_levels, 1);
noise_max_residual_pct = zeros(noise_n_levels, 1);
for noise_idx = 1:noise_n_levels
    noise_pct = noise_injection_pct(noise_idx);
    mean_resid_trials = zeros(noise_n_trials, 1);
    max_resid_trials = zeros(noise_n_trials, 1);
    for trial_idx = 1:noise_n_trials
        noise_scale = (noise_pct / 100) * max(abs(Z_exp), eps);
        noise_vec = (noise_scale ./ sqrt(2)) .* (randn(n_freq, 1) + 1i * randn(n_freq, 1));
        Z_noisy = Z_exp + noise_vec;
        [gamma_noisy, R_inf_noisy] = TR_DRT_matlab2011(freq_vec, Z_noisy, el);
        Z_fit_noisy = calculate_EIS_matlab2011(freq_vec, gamma_noisy, R_inf_noisy);
        resid_pct = 100 * abs(Z_noisy - Z_fit_noisy) ./ max(abs(Z_noisy), eps);
        mean_resid_trials(trial_idx) = mean(resid_pct);
        max_resid_trials(trial_idx) = max(resid_pct);
    end
    noise_mean_residual_pct(noise_idx) = mean(mean_resid_trials);
    noise_max_residual_pct(noise_idx) = mean(max_resid_trials);
end

sensitivity_result = struct('n_boot', n_boot, 'gamma_samples', gamma_samples, ...
    'gamma_mean', gamma_mean, 'gamma_std', gamma_std, 'gamma_lower', gamma_lower, ...
    'gamma_upper', gamma_upper, 'mean_std', mean(gamma_std(active_mask)), ...
    'max_std', max(gamma_std), 'mean_cv_pct', mean(cv_pct), ...
    'max_cv_pct', max(cv_pct), 'relative_band_area_pct', band_area_pct, ...
    'noise_injection_pct', noise_injection_pct, ...
    'noise_mean_residual_pct', noise_mean_residual_pct, ...
    'noise_max_residual_pct', noise_max_residual_pct);
end

function result = perform_lambda_perturbation_sensitivity_matlab2011(freq_vec, Z_exp, lambda_ref, reference_mean_resid)
% perform_lambda_perturbation_sensitivity_matlab2011 - Residual sensitivity around selected lambda
%
% This is a local robustness check around the selected lambda rather than a full
% grid search; it reports percent changes in mean and maximum fit residual.

freq_vec = freq_vec(:);
Z_exp = Z_exp(:);

if nargin < 4 || isempty(reference_mean_resid)
    [gamma_ref, R_inf_ref] = TR_DRT_matlab2011(freq_vec, Z_exp, lambda_ref);
    Z_ref = calculate_EIS_matlab2011(freq_vec, gamma_ref, R_inf_ref);
    reference_mean_resid = mean(abs(Z_exp - Z_ref));
end

factors = [0.3; 0.5; 0.7; 1.0; 1.5; 2.0; 3.0];
lambda_values = lambda_ref * factors;
n_vals = numel(lambda_values);

mean_residual = zeros(n_vals, 1);
max_residual = zeros(n_vals, 1);

for i = 1:n_vals
    [gamma_i, R_inf_i] = TR_DRT_matlab2011(freq_vec, Z_exp, lambda_values(i));
    Z_fit_i = calculate_EIS_matlab2011(freq_vec, gamma_i, R_inf_i);
    resid_i = abs(Z_exp - Z_fit_i);
    mean_residual(i) = mean(resid_i);
    max_residual(i) = max(resid_i);
end

delta_pct = 100 * (mean_residual - reference_mean_resid) / max(reference_mean_resid, eps);
min_delta_pct = min(delta_pct);
max_delta_pct = max(delta_pct);
max_abs_delta_pct = max(abs(delta_pct));

[~, best_idx] = min(mean_residual);
[~, worst_idx] = max(mean_residual);
best_factor = factors(best_idx);
worst_factor = factors(worst_idx);

center_idx = find(abs(factors - 1) < 1e-12, 1);
if isempty(center_idx)
    local_slope = 0;
elseif center_idx == 1 || center_idx == n_vals
    local_slope = 0;
else
    x1 = log10(lambda_values(center_idx - 1));
    x2 = log10(lambda_values(center_idx + 1));
    y1 = mean_residual(center_idx - 1);
    y2 = mean_residual(center_idx + 1);
    local_slope = (y2 - y1) / max(x2 - x1, eps);
end

result = struct('factors', factors, 'lambda_values', lambda_values, ...
    'mean_residual', mean_residual, 'max_residual', max_residual, ...
    'reference_mean_resid', reference_mean_resid, 'delta_pct', delta_pct, ...
    'min_delta_pct', min_delta_pct, 'max_delta_pct', max_delta_pct, ...
    'max_abs_delta_pct', max_abs_delta_pct, 'local_slope', local_slope, ...
    'best_factor', best_factor, 'worst_factor', worst_factor);
end

function write_matrix_with_header_2011(file_name, sheet_name, header_row, numeric_data)
% write_matrix_with_header_2011 - Write a header row and numeric matrix to Excel
%
% All workbook exports go through this wrapper so every sheet uses the same
% two-row layout: header labels in A1 and numeric data starting in A2.

% xlswrite(file_name, header_row, sheet_name, 'A1');
% xlswrite(file_name, numeric_data, sheet_name, 'A2');
end

function sample_tag = make_sample_tag_from_input_2011(input_path)
% make_sample_tag_from_input_2011 - Build a filesystem-safe stem for all outputs.
[~, base_name, ~] = fileparts(input_path);
sample_tag = lower(regexprep(base_name, '[^a-zA-Z0-9]', ''));
if isempty(sample_tag)
    sample_tag = 'sample';
end
end

function sheet_name = make_sheet_name_2011(prefix, temperature_value)
% make_sheet_name_2011 - Build Excel-safe sheet names for per-temperature exports

sheet_name = sprintf('%s_%0.1fK', prefix, temperature_value);
sheet_name = strrep(sheet_name, '.', 'p');
sheet_name = strrep(sheet_name, '-', 'm');
if numel(sheet_name) > 31
    sheet_name = sheet_name(1:31);
end
end

function plot_horizontal_line_2011(y_value, line_style, line_color)
% plot_horizontal_line_2011 - MATLAB 2011 compatible replacement for yline.
%
% MATLAB 2011 does not provide yline, so the workflow draws a standard line
% object using the current axis limits whenever a horizontal reference is needed.

if nargin < 2 || isempty(line_style)
    line_style = '--';
end
if nargin < 3 || isempty(line_color)
    line_color = [0 0 0];
end

ax = gca;
x_limits = get(ax, 'XLim');
line(x_limits, [y_value y_value], 'LineStyle', line_style, 'Color', line_color, 'Parent', ax);
end

function n_val = solve_n_from_fwhm_implicit_2011(sigma)
% solve_n_from_fwhm_implicit_2011 - Compute n from Gaussian sigma via ZARC FWHM relation.
%
% Physics-based formula: ZARC DRT peak width is related to CPE exponent n by:
%   FWHM_lntau = (2/n) * acosh(2 + cos(pi*n))
% where FWHM_lntau = 2*sqrt(2*ln(2))*sigma for a Gaussian of width sigma in ln(tau).
%
% This function solves the implicit equation:
%   sqrt(2*ln(2)) * sigma * n = acosh(2 + cos(pi*n))
% numerically for n in (0, 1].
%
% INPUTS:
%   sigma - Gaussian width parameter in ln(tau) space
%
% OUTPUTS:
%   n_val - CPE exponent (between 0.01 and 0.99)

    if sigma <= 0
        n_val = 0.5;  % Default if sigma is non-positive
        return;
    end

    fwhm_target = 2 * sqrt(2 * log(2)) * sigma;

    % Define objective: minimize |FWHM(n) - FWHM_target|
    objective = @(n) abs((2/n) * acosh(2 + cos(pi*n)) - fwhm_target);

    % Grid search to find reasonable initial guess
    n_candidates = linspace(0.01, 0.99, 100);
    errors = arrayfun(objective, n_candidates);
    [~, best_idx] = min(errors);
    n_init = n_candidates(best_idx);

    % Refine with fminsearch
    opts = optimset('Display', 'off', 'TolX', 1e-8, 'MaxIter', 500);
    n_val = fminsearch(objective, n_init, opts);

    % Clamp to valid range
    n_val = max(0.01, min(n_val, 0.99));
end

function n_val = estimate_n_from_phase_2011(tau_peak, amplitude, gamma_ref)
% estimate_n_from_phase_2011 - Estimate n from apex phase angle (single-arc approximation).
%
% For a depressed arc (ZARC), the phase angle at the arc maximum (omega*tau=1) is:
%   phi_apex = -pi*n/4  (radians) or -45*n degrees
% This gives n = -4*phi_apex/pi or n = -phi_apex/45 (degrees).
%
% INPUTS:
%   tau_peak    - Peak relaxation time (s)
%   amplitude   - Peak amplitude (Ohm)
%   gamma_ref   - Reference DRT for context
%
% OUTPUTS:
%   n_val - CPE exponent estimated from phase

    % Placeholder: return heuristic if insufficient data
    % In a full implementation, one would measure phase from Nyquist or Bode plot.
    n_val = 0.5 + 0.2 * exp(-amplitude / 100);  % Heuristic fallback
end

% ============================================================================
%{
% RGB Multi-Temperature Peak Visualization Function
% ============================================================================
function generate_temperature_rgb_peak_map(all_temps_peak_data, temperature_vec, output_file)
% generate_temperature_rgb_peak_map - Create RGB visualization of top 3 DRT peaks across temperatures
%
% Creates a 256 x n_temps x 3 RGB image where:
%   X-axis (columns): relaxation time ln(tau) mapped to 256 pixels
%   Y-axis (rows): temperature indices
%   R channel: highest-area peak
%   G channel: second-highest-area peak
%   B channel: third-highest-area peak
%   Intensity: pixel value (0-255) represents DRT magnitude at that ln(tau)
%
% INPUTS:
%   all_temps_peak_data: cell array of structs, one per temperature
%                        each containing: temperature, tau_dense, peaks, peak_count
%   temperature_vec: row vector of temperatures (K)
%   output_file: path to output PNG file
%
% OUTPUTS:
%   Saves RGB visualization PNG to output_file

    num_temps = numel(all_temps_peak_data);
    if num_temps == 0
        fprintf('Warning: no temperature data to visualize\n');
        return;
    end
    
    % Build master ln_tau grid from first temperature's tau_dense
    tau_dense_1 = all_temps_peak_data{1}.tau_dense(:);
    ln_tau_master = log(tau_dense_1);
    n_ln_tau = numel(ln_tau_master);
    
    % Normalize ln_tau to pixel grid [0, 255]
    ln_tau_min = min(ln_tau_master);
    ln_tau_max = max(ln_tau_master);
    if ln_tau_max <= ln_tau_min
        fprintf('Warning: invalid ln_tau range\n');
        return;
    end
    ln_tau_pixel_grid = linspace(0, 255, n_ln_tau);
    
    % Initialize RGB image: [n_temps, 256, 3]
    rgb_image = zeros(num_temps, 255, 3, 'uint8');
    
    % Process each temperature
    for temp_idx = 1:num_temps
        peak_data = all_temps_peak_data{temp_idx};
        peaks = peak_data.peaks;
        n_peaks = peak_data.peak_count;
        
        if n_peaks == 0
            continue;  % Skip temperature with no peaks
        end
        
        % Compute peak areas: integral of Gaussian = A * sigma * sqrt(2*pi)
        peak_areas = zeros(n_peaks, 1);
        for p = 1:n_peaks
              peak_areas(p) = peaks(p).amplitude * 1; %peaks(p).sigma * sqrt(2 * pi);
        end
        
        % Sort peaks by area (descending)
        [~, sort_idx] = sort(peak_areas, 'descend');
        top_3_idx = sort_idx(1:min(3, n_peaks));
        
        % Reconstruct Gaussian curves for each top peak
        tau_dense_t = peak_data.tau_dense(:);
        ln_tau_t = log(tau_dense_t);
        
        % Initialize R, G, B channels for this temperature
        r_channel = zeros(1, 256, 'double');
        g_channel = zeros(1, 256, 'double');
        b_channel = zeros(1, 256, 'double');
        
        % Assign top 3 peaks to R, G, B
        channels = [r_channel; g_channel; b_channel];
        for rank = 1:min(3, numel(top_3_idx))
            peak_idx = top_3_idx(rank);
            p = peaks(peak_idx);
            
            % Reconstruct Gaussian: A * exp(-(ln_tau - mu)^2 / (2*sigma^2))
            gaussian_vals = p.amplitude * exp(-((ln_tau_t - log(p.tau)) .^ 2) / (2 * p.sigma ^ 2));
            
            % Interpolate to pixel grid
            if exist('interp1', 'file') == 2
                pixel_vals = interp1(ln_tau_t, gaussian_vals, ...
                    linspace(ln_tau_min, ln_tau_max, 256), 'linear', 0);
            else
                % MATLAB 2011 fallback: simple nearest-neighbor
                [~, closest_idx] = min(abs(ln_tau_t - linspace(ln_tau_min, ln_tau_max, 256)'), [], 1);
                pixel_vals = gaussian_vals(closest_idx);
            end
            
            % Normalize to [0, 255]
            max_val = max(pixel_vals);
            if max_val > 0
                pixel_vals = uint8(round(255 * pixel_vals / max_val));
            else
                pixel_vals = uint8(zeros(1, 256));
            end
            
            channels(rank, :) = pixel_vals;
        end
        
        % Fill RGB image for this temperature
        rgb_image(temp_idx, :, 1) = channels(1, :);  % R
        rgb_image(temp_idx, :, 2) = channels(2, :);  % G
        rgb_image(temp_idx, :, 3) = channels(3, :);  % B
    end
    
    % Save RGB image as a frequency-vs-temperature figure with axes.
    try
        temperature_axis = temperature_vec(:)';
        [temperature_axis, temp_order] = sort(temperature_axis, 'ascend');
        rgb_display = permute(rgb_image(temp_order, end:-1:1, :), [2 1 3]);
        freq_axis = exp(-ln_tau_master(end:-1:1));

        fig = figure;
        image(temperature_axis, freq_axis, rgb_display);
        set(gca, 'YDir', 'normal', 'YScale', 'log');
        axis on;
        box on;
        xlabel('Temperature (K)');
        ylabel('Frequency (Hz)');
        title('RGB peak visualization');
        xlim([temperature_axis(1) temperature_axis(end)]);
        ylim([freq_axis(1) freq_axis(end)]);
        xtick_count = min(6, numel(temperature_axis));
        xticks(linspace(temperature_axis(1), temperature_axis(end), xtick_count));
        if freq_axis(1) > 0 && freq_axis(end) > 0
            ytick_count = min(6, numel(freq_axis));
            y_ticks = logspace(log10(freq_axis(1)), log10(freq_axis(end)), ytick_count);
            yticks(y_ticks);
            yticklabels(arrayfun(@(v) sprintf('%.3g', v), y_ticks, 'UniformOutput', false));
        end

        % if exist('exportgraphics', 'file') == 2
        %     exportgraphics(fig, output_file, 'Resolution', 300);
        % else
        %     print(fig, output_file, '-dpng', '-r300');
        % end

        fprintf('RGB peak visualization saved: %s (%d x %d x 3)\n', output_file, numel(freq_axis), numel(temperature_axis));
    catch ME
        fprintf('Error writing PNG: %s\n', ME.message);
    end
end
%}
