function drt_input_analysis_matlab2011(repo_root,input_path,temperature,indata)
% drt_input_analysis_matlab2011 - Complete DRT (Distribution of Relaxation Times) Analysis Workflow
%
% MAIN FUNCTION - Performs end-to-end DRT analysis on electrochemical impedance spectroscopy (EIS) data
%
% WORKFLOW OVERVIEW:
%   1. LOAD DATA: Read impedance data (frequency, magnitude, phase) from Excel file
%   2. VALIDATE: Remove zero/negative frequencies, sort data
%   3. INVERT: Test multiple lambda values and select lambda using
%      ChemElectroChem-style criteria (RID + Re/Im CV + resampling variance)
%   4. CALCULATE: Compute fitted impedance and residual metrics from recovered DRT
%   5. VALIDATE: Run Kramers-Kronig (KK) consistency tests on measured and DRT-fitted impedance
%   6. ANALYZE PEAKS: Detect peaks in DRT, characterize them with RBF interpolation
%   7. EXTRACT PARAMETERS: Calculate equivalent circuit R and C for each peak
%   8. VISUALIZE: Generate plots (console only, no PNG files for easy script copying)
%
% KEY VARIABLES:
%   freq_vec - Measurement frequencies (Hz), log-spaced preferred
%   Z_exp    - Experimental complex impedance (Ohm)
%   gamma    - Distribution of Relaxation Times [N x 1], main output
%   R_inf    - High-frequency limit resistance (Ohm)
%   Z_cal    - Fitted impedance spectrum from DRT
%   tau      - Relaxation times corresponding to DRT (seconds)
%   lambda_values - Regularization parameters to test (typical: 1e-4 to 1.0)
%
% OUTPUT FILES:
%   - drt_matlab2011_run.log : Text log of all console output
%   - plots_2011/ directory : Figures displayed on screen (no PNG files)
%
% CUSTOMIZATION:
%   Edit 'User settings' section below to change:
%   - input_path: Path to Excel file with EIS data
%   - lambda_values: Regularization parameters to test
%   - phase_units: 'deg' (degrees) or 'rad' (radians)
%   - data_format: 'mag_phase' or 'real_imag'
%
% EXPECTED DATA FORMAT:
%   Column 1 = frequency (Hz)
%   Column 2 = |Z| (Ohm) or Re(Z)
%   Column 3 = theta (degrees/radians) or Im(Z)

clc;


% --------------------------- User settings ---------------------------
sz=size(indata,2);

if rem(sz, 3) ~= 0
    error('Input data must contain frequency, magnitude, and phase columns for each temperature');
end

analysis_results=zeros(size(indata,1)+1,size(indata,2));
peak_file=zeros(size(indata,1)+1,size(indata,2));
num_sets = sz / 3;
kk_summary = zeros(num_sets, 13);
sensitivity_summary = zeros(num_sets, 8);
lambda_perturb_summary = zeros(num_sets, 9);
kk_output_file = fullfile(repo_root, 'kk_result.xlsx');
sensitivity_output_file = fullfile(repo_root, 'sensitivity_result.xlsx');
lambda_perturb_output_file = fullfile(repo_root, 'lambda_perturbation_result.xlsx');
peak_output_file = fullfile(repo_root, 'peak_result.xlsx');

for ab=1:sz/3
data=[indata(:,3*ab-2) indata(:,3*ab-1) indata(:,3*ab)];     
phase_units = 'deg';
data_format = 'mag_phase';
lambda_values = [1e-4, 1e-3, 1e-2, 1e-1, 1e0];
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
[freq_vec, Z_exp] = load_f_Z_theta_2011(data, phase_units, data_format);

valid = freq_vec > 0;
if ~all(valid)
    fprintf('Removing %d zero/negative frequency points\n', sum(~valid));
    freq_vec = freq_vec(valid);
    Z_exp = Z_exp(valid);
end

[freq_vec, sort_idx] = sort(freq_vec);
Z_exp = Z_exp(sort_idx);

fprintf('Loaded %d frequency points (sorted by increasing frequency)\n', numel(freq_vec));
fprintf('Frequency range: %.1f - %.1f Hz\n', min(freq_vec), max(freq_vec));
fprintf('Magnitude range: %.1f - %.1f Ohm\n', min(abs(Z_exp)), max(abs(Z_exp)));
phase_deg = angle(Z_exp) * 180 / pi;
fprintf('Phase range: %.2f - %.2f deg\n', min(phase_deg), max(phase_deg));

% --------------------------- Plot input impedance --------------------
figure('Color', 'w', 'Name', 'Measured impedance');
subplot(2,1,1);
semilogx(freq_vec, real(Z_exp), 'o-');
ylabel('Re(Z)');
grid on;

subplot(2,1,2);
semilogx(freq_vec, imag(Z_exp), 'o-');
xlabel('Frequency (Hz)');
ylabel('Im(Z)');
grid on;

% --------------------------- DRT inversion ---------------------------
fprintf('\nTesting different lambda values...\n');
num_l = numel(lambda_values);
results(num_l) = struct('lambda', 0, 'gamma', [], 'R_inf', 0, 'mean_resid', 0, 'Z_cal', [], ...
    'rid', 0, 'cv_real', 0, 'cv_imag', 0, 'resample_var', 0, ...
    'gamma_real_only', [], 'gamma_imag_only', []);

for i = 1:num_l
    el = lambda_values(i);
    [gamma_test, R_inf_test] = TR_DRT_matlab2011(freq_vec, Z_exp, el);
    Z_cal_test = calculate_EIS_matlab2011(freq_vec, gamma_test, R_inf_test);

    residual_test = abs(Z_exp - Z_cal_test);
    mean_resid = mean(residual_test);

    % ChemElectroChem (Schlueter et al., 2019) style lambda metrics:
    % RID compares DRT from Re-only and Im-only inversions.
    [gamma_re_only, R_re_only] = TR_DRT_component_matlab2011(freq_vec, Z_exp, el, 'real');
    [gamma_im_only, ~] = TR_DRT_component_matlab2011(freq_vec, Z_exp, el, 'imag');
    rid_val = mean((gamma_re_only - gamma_im_only).^2);

    % Cross-validation terms: predict Im from Re-only DRT and Re from Im-only DRT.
    Z_from_re_only = calculate_EIS_matlab2011(freq_vec, gamma_re_only, R_re_only);

    R_from_im = estimate_Rinf_for_gamma_matlab2011(freq_vec, Z_exp, gamma_im_only);
    Z_from_im_only = calculate_EIS_matlab2011(freq_vec, gamma_im_only, R_from_im);

    cv_imag_val = mean((imag(Z_exp) - imag(Z_from_re_only)).^2);
    cv_real_val = mean((real(Z_exp) - real(Z_from_im_only)).^2);

    % Resampling variance proxy for stability (paper's repetitive/resampling idea).
    resample_var_val = estimate_resample_variance_matlab2011(freq_vec, Z_exp, el, 8);

    fprintf('lambda = %6.1e: R_inf = %8.3f, Mean residual = %8.3f, RID = %.3e, CVre = %.3e, CVim = %.3e\n', ...
        el, R_inf_test, mean_resid, rid_val, cv_real_val, cv_imag_val);

    results(i).lambda = el;
    results(i).gamma = gamma_test;
    results(i).R_inf = R_inf_test;
    results(i).mean_resid = mean_resid;
    results(i).Z_cal = Z_cal_test;
    results(i).rid = rid_val;
    results(i).cv_real = cv_real_val;
    results(i).cv_imag = cv_imag_val;
    results(i).resample_var = resample_var_val;
    results(i).gamma_real_only = gamma_re_only;
    results(i).gamma_imag_only = gamma_im_only;
end

all_resid = zeros(num_l,1);
all_rid = zeros(num_l,1);
all_cv = zeros(num_l,1);
all_var = zeros(num_l,1);
for i = 1:num_l
    all_resid(i) = results(i).mean_resid;
    all_rid(i) = results(i).rid;
    all_cv(i) = results(i).cv_real + results(i).cv_imag;
    all_var(i) = results(i).resample_var;
end

% ChemElectroChem-inspired selection: lower/upper lambda boundary from CV and variance,
% then select minimum RID inside the admissible interval.
[~, idx_cv_min] = min(all_cv);
[~, idx_var_min] = min(all_var);
idx_low = min(idx_cv_min, idx_var_min);
idx_high = max(idx_cv_min, idx_var_min);
if idx_low == idx_high
    idx_candidates = (1:num_l).';
else
    idx_candidates = (idx_low:idx_high).';
end

rid_cand = all_rid(idx_candidates);
[~, rid_rel_idx] = min(rid_cand);
best_idx = idx_candidates(rid_rel_idx);

score = normalize_metric_2011(all_rid) + normalize_metric_2011(all_cv) + normalize_metric_2011(all_var);
[~, score_best_idx] = min(score);
if score_best_idx ~= best_idx
    best_idx = score_best_idx;
end

el = results(best_idx).lambda;
gamma = results(best_idx).gamma;
R_inf = results(best_idx).R_inf;
Z_cal = results(best_idx).Z_cal;

fprintf('\nSelected lambda (ChemElectroChem optimization) = %.1e\n', el);
fprintf('  Boundary index from CV/variance: [%d, %d]\n', idx_low, idx_high);
fprintf('  Selected RID = %.3e | CV(total) = %.3e | Var = %.3e\n', ...
    all_rid(best_idx), all_cv(best_idx), all_var(best_idx));
fprintf('Recovered R_inf = %.6f\n', R_inf);
fprintf('DRT vector length = %d\n', numel(gamma));

% --------------------------- Plot DRT vs tau -------------------------
tau = 1 ./ (2 * pi * freq_vec);
figure('Color', 'w', 'Name', 'DRT vs relaxation time');
semilogx(tau, gamma, 'o-');
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title('DRT vs Relaxation Time');
grid on;

% --------------------------- Compare fit -----------------------------
residual = abs(Z_exp - Z_cal);
fprintf('\nMean absolute residual: %.6f\n', mean(residual));
fprintf('Max absolute residual : %.6f\n', max(residual));

figure('Color', 'w', 'Name', 'Impedance comparison');
subplot(2,1,1);
semilogx(freq_vec, real(Z_exp), 'ko-'); hold on;
semilogx(freq_vec, real(Z_cal), 'r.-');
ylabel('Re(Z)');
grid on;
legend('Re(Z) exp', 'Re(Z) fit', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vec, imag(Z_exp), 'ko-'); hold on;
semilogx(freq_vec, imag(Z_cal), 'r.-');
xlabel('Frequency (Hz)');
ylabel('Im(Z)');
grid on;
legend('Im(Z) exp', 'Im(Z) fit', 'Location', 'best');

% ------------------ Lambda perturbation sensitivity -------------------
lambda_perturb_result = perform_lambda_perturbation_sensitivity_matlab2011( ...
    freq_vec, Z_exp, el, mean(residual));

fprintf('\nLambda perturbation sensitivity around selected lambda:\n');
fprintf('  Reference lambda      : %.3e\n', el);
fprintf('  Reference mean resid  : %.6f\n', lambda_perturb_result.reference_mean_resid);
fprintf('  Min residual change   : %.2f %%\n', lambda_perturb_result.min_delta_pct);
fprintf('  Max residual change   : %.2f %%\n', lambda_perturb_result.max_delta_pct);
fprintf('  Max abs change        : %.2f %%\n', lambda_perturb_result.max_abs_delta_pct);

lambda_perturb_summary(ab, :) = [temperature(ab), el, ...
    lambda_perturb_result.reference_mean_resid, lambda_perturb_result.min_delta_pct, ...
    lambda_perturb_result.max_delta_pct, lambda_perturb_result.max_abs_delta_pct, ...
    lambda_perturb_result.local_slope, lambda_perturb_result.best_factor, ...
    lambda_perturb_result.worst_factor];

lambda_sheet = make_sheet_name_2011('LAMBDA', temperature(ab));
lambda_detail = [lambda_perturb_result.factors(:), lambda_perturb_result.lambda_values(:), ...
    lambda_perturb_result.mean_residual(:), lambda_perturb_result.max_residual(:), ...
    lambda_perturb_result.delta_pct(:)];
write_matrix_with_header_2011(lambda_perturb_output_file, lambda_sheet, ...
    {'factor_vs_ref','lambda','mean_residual','max_residual','delta_mean_residual_pct'}, lambda_detail);

figure('Color', 'w', 'Name', sprintf('Lambda sensitivity %.1f K', temperature(ab)));
subplot(2,1,1);
semilogx(lambda_perturb_result.lambda_values, lambda_perturb_result.mean_residual, 'bo-'); hold on;
plot(el, lambda_perturb_result.reference_mean_resid, 'rs', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('Lambda');
ylabel('Mean |Z_{exp}-Z_{fit}|');
grid on;
legend('Perturbed lambda', 'Selected lambda', 'Location', 'best');

subplot(2,1,2);
semilogx(lambda_perturb_result.lambda_values, lambda_perturb_result.delta_pct, 'm.-'); hold on;
plot(el, 0, 'ks', 'MarkerSize', 8, 'LineWidth', 1.2);
xlabel('Lambda');
ylabel('Residual change (%)');
grid on;
legend('Delta residual', 'Reference', 'Location', 'best');

% --------------------------- Peak analysis ---------------------------
tau_data = tau;
gamma_max = max(gamma);
thresholds = [0.01, 0.05, 0.10, 0.15];

for t = 1:numel(thresholds)
    thr = thresholds(t);
    min_height = thr * gamma_max;
    [locs_t, pks_t] = find_peaks_simple(gamma, min_height, 10);
    fprintf('Threshold %.2f: Found %d peaks\n', thr, numel(locs_t));
    for k = 1:numel(locs_t)
        fprintf('  Peak %d: gamma = %.1f\n', k, pks_t(k));
    end
end

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
fprintf('Peak | tau (s) | gamma (1/Ohm) | R (Ohm) | C (F) | n | hidden\n');
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

% --------------------------- KK test ---------------------------------
kk_result = perform_kk_test_matlab2011(freq_vec, Z_exp);
kk_drt_result = perform_kk_test_matlab2011(freq_vec, Z_cal);
fprintf('\nKramers-Kronig consistency test:\n');
fprintf('  Real residual  : %.2f %%\n', kk_result.real_residual_pct);
fprintf('  Imag residual  : %.2f %%\n', kk_result.imag_residual_pct);
fprintf('  Total residual : %.2f %%\n', kk_result.total_residual_pct);
fprintf('  Status         : %s\n', kk_result.status);

fprintf('\nKramers-Kronig test on DRT-generated impedance:\n');
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
write_matrix_with_header_2011(kk_output_file, kk_sheet, ...
    {'freq_Hz','Zre_exp','Zim_exp','Zre_KK','Zim_KK','abs_residual','rel_residual_pct'}, kk_detail);

kk_drt_sheet = make_sheet_name_2011('KKDRT', temperature(ab));
kk_drt_detail = [freq_vec(:), real(Z_cal(:)), imag(Z_cal(:)), real(kk_drt_result.Z_fit(:)), ...
    imag(kk_drt_result.Z_fit(:)), abs(kk_drt_result.residual(:)), kk_drt_result.relative_residual_pct(:)];
write_matrix_with_header_2011(kk_output_file, kk_drt_sheet, ...
    {'freq_Hz','Zre_drt','Zim_drt','Zre_KK','Zim_KK','abs_residual','rel_residual_pct'}, kk_drt_detail);

figure('Color', 'w', 'Name', sprintf('KK test %.1f K', temperature(ab)));
subplot(2,1,1);
semilogx(freq_vec, real(Z_exp), 'ko-'); hold on;
semilogx(freq_vec, real(kk_result.Z_fit), 'b.-');
ylabel('Re(Z)');
grid on;
legend('Measured', 'KK fit', 'Location', 'best');

figure('Color', 'w', 'Name', sprintf('KK test on DRT impedance %.1f K', temperature(ab)));
subplot(2,1,1);
semilogx(freq_vec, real(Z_cal), 'ko-'); hold on;
semilogx(freq_vec, real(kk_drt_result.Z_fit), 'b.-');
ylabel('Re(Z)');
grid on;
legend('DRT impedance', 'KK fit', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vec, imag(Z_cal), 'ko-'); hold on;
semilogx(freq_vec, imag(kk_drt_result.Z_fit), 'b.-');
xlabel('Frequency (Hz)');
ylabel('Im(Z)');
grid on;
legend('DRT impedance', 'KK fit', 'Location', 'best');

subplot(2,1,2);
semilogx(freq_vec, imag(Z_exp), 'ko-'); hold on;
semilogx(freq_vec, imag(kk_result.Z_fit), 'b.-');
xlabel('Frequency (Hz)');
ylabel('Im(Z)');
grid on;
legend('Measured', 'KK fit', 'Location', 'best');

% --------------------------- Sensitivity -----------------------------
sensitivity_result = perform_sensitivity_analysis_matlab2011(freq_vec, Z_exp, el, gamma, Z_cal, ab);
fprintf('\nSensitivity analysis (%d bootstrap runs):\n', sensitivity_result.n_boot);
fprintf('  Mean std. band : %.4e\n', sensitivity_result.mean_std);
fprintf('  Max std. band  : %.4e\n', sensitivity_result.max_std);
fprintf('  Mean CV        : %.2f %%\n', sensitivity_result.mean_cv_pct);
fprintf('  Max CV         : %.2f %%\n', sensitivity_result.max_cv_pct);

sensitivity_summary(ab, :) = [temperature(ab), el, sensitivity_result.mean_std, ...
    sensitivity_result.max_std, sensitivity_result.mean_cv_pct, ...
    sensitivity_result.max_cv_pct, sensitivity_result.relative_band_area_pct, ...
    numel(data_peaks_fitted)];

sens_sheet = make_sheet_name_2011('SENS', temperature(ab));
sens_detail = [tau(:), gamma(:), sensitivity_result.gamma_mean(:), ...
    sensitivity_result.gamma_std(:), sensitivity_result.gamma_lower(:), ...
    sensitivity_result.gamma_upper(:)];
write_matrix_with_header_2011(sensitivity_output_file, sens_sheet, ...
    {'tau_s','gamma_ref','gamma_mean','gamma_std','gamma_lower','gamma_upper'}, sens_detail);

figure('Color', 'w', 'Name', sprintf('DRT sensitivity %.1f K', temperature(ab)));
semilogx(tau, gamma, 'k-', 'LineWidth', 1.2); hold on;
semilogx(tau, sensitivity_result.gamma_mean, 'b--', 'LineWidth', 1.2);
semilogx(tau, sensitivity_result.gamma_lower, 'r:');
semilogx(tau, sensitivity_result.gamma_upper, 'r:');
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title('DRT sensitivity band');
grid on;
legend('Reference DRT', 'Bootstrap mean', 'Bootstrap band', 'Location', 'Best');

figure('Color', 'w', 'Name', 'DRT peak analysis');
semilogx(tau_data, gamma, 'bo-'); hold on;

tau_plot = peak_refine.tau_dense;
gamma_rbf = peak_refine.gamma_rbf;
plot(tau_plot, gamma_rbf, 'g-', 'LineWidth', 1.5);
plot(tau_plot, peak_refine.fit_curve, 'c--', 'LineWidth', 1.2);

for i = 1:numel(data_peaks_fitted)
    peak = data_peaks_fitted(i);
    if peak.is_hidden > 0
        plot(peak.tau, peak.amplitude, 'mo', 'MarkerSize', 8, 'LineWidth', 1.2);
    else
        plot(peak.tau, peak.amplitude, 'ro', 'MarkerSize', 8, 'LineWidth', 1.2);
    end
end

if ~isempty(peak_refine.base_peak_idx)
    scatter(tau_plot(peak_refine.base_peak_idx), gamma_rbf(peak_refine.base_peak_idx), 80, 'r', '*');
end
if ~isempty(peak_refine.hidden_peak_idx)
    scatter(tau_plot(peak_refine.hidden_peak_idx), gamma_rbf(peak_refine.hidden_peak_idx), 80, 'm', 's');
end

set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau) (1/Ohm)');
title('Original Data DRT with RBF-refined Peak Analysis');
grid on;
legend('DRT data', 'RBF interpolation', 'Gaussian RBF sum fit', 'Refined peaks', 'Base peaks', 'Hidden peaks', 'Location', 'Best');

fprintf('Run complete at %s\n', datestr(now, 31));

analysis_results(1,2*ab)=temperature(ab);
analysis_results(end-size(tau)+1:end,2*ab-1)=tau;
analysis_results(end-size(gamma)+1:end,2*ab)=gamma;

peak_store=[temperature(ab);tau_peak';0;R_file';0;C_file';0;n_file'];
peak_file(1:size(peak_store),ab)=peak_store;
peak_sheet = make_sheet_name_2011('PEAK', temperature(ab));
peak_numeric = zeros(numel(data_peaks_fitted), 8);
for i = 1:numel(data_peaks_fitted)
    peak_numeric(i,:) = [data_peaks_fitted(i).peak_id, data_peaks_fitted(i).tau, ...
        data_peaks_fitted(i).amplitude, data_peaks_fitted(i).sigma, data_peaks_fitted(i).R, ...
        data_peaks_fitted(i).C, data_peaks_fitted(i).n, data_peaks_fitted(i).is_hidden];
end
write_matrix_with_header_2011(peak_output_file, peak_sheet, ...
    {'peak_id','tau_s','amplitude','sigma_logtau','R_ohm','C_F','n_est','is_hidden'}, peak_numeric);
diary off;
end

xlswrite('analysis_result.xlsx',analysis_results);
xlswrite('peak_result.xlsx',peak_file);
write_matrix_with_header_2011(kk_output_file, 'summary', ...
    {'Temperature_K','Lambda','R_inf','KK_real_pct','KK_imag_pct','KK_total_pct','KK_max_pct','KK_status_code', ...
    'KK_DRT_real_pct','KK_DRT_imag_pct','KK_DRT_total_pct','KK_DRT_max_pct','KK_DRT_status_code'}, ...
    kk_summary);
write_matrix_with_header_2011(sensitivity_output_file, 'summary', ...
    {'Temperature_K','Lambda','Mean_std','Max_std','Mean_CV_pct','Max_CV_pct','Band_area_pct','Peak_count'}, ...
    sensitivity_summary);
write_matrix_with_header_2011(lambda_perturb_output_file, 'summary', ...
    {'Temperature_K','Lambda_ref','Ref_mean_residual','Min_delta_pct','Max_delta_pct','Max_abs_delta_pct', ...
    'Local_slope','Best_factor','Worst_factor'}, lambda_perturb_summary);
end

function [freq, Z] = load_f_Z_theta_2011(data, phase_units, data_format)
% load_f_Z_theta_2011 - Load impedance data from file and convert to complex impedance
%
% This function reads impedance data from a text file (handles .xls, .xlsx, .txt, .csv)
% and converts it to complex impedance format. It handles two input formats:
%   - 'mag_phase': column 2 is magnitude |Z|, column 3 is phase angle
%   - 'real_imag': column 2 is Re(Z), column 3 is Im(Z)
%
% INPUTS:
%   data        - data wrt temperature in f |Z| phase in degrees for each
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

function [gamma, R_inf] = TR_DRT_component_matlab2011(freq_vec, Z_exp, el, component)
% TR_DRT_component_matlab2011 - DRT inversion using only real or imaginary data.

A_re = calc_A_re_matlab2011(freq_vec);
A_im = calc_A_im_matlab2011(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

switch lower(component)
    case 'real'
        M = [A_re, ones(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
        b = [Z_re; zeros(n,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
        b = [Z_im; zeros(n,1)];
    otherwise
        error('Unknown component mode: %s', component);
end

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

function R_inf = estimate_Rinf_for_gamma_matlab2011(freq_vec, Z_exp, gamma)
% estimate_Rinf_for_gamma_matlab2011 - Best nonnegative R_inf for fixed gamma from Re(Z).

A_re = calc_A_re_matlab2011(freq_vec);
re_model_wo_r = A_re * gamma(:);
R_inf = mean(real(Z_exp(:)) - re_model_wo_r);
if ~(R_inf > 0)
    R_inf = 0;
end
end

function mean_var = estimate_resample_variance_matlab2011(freq_vec, Z_exp, el, n_boot)
% estimate_resample_variance_matlab2011 - Bootstrap noise-resampling variance of gamma.

n = numel(freq_vec);
if nargin < 4 || isempty(n_boot)
    n_boot = 8;
end

sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);
sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);

gamma_samples = zeros(n, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n,1) + 1i * sigma_im * randn(n,1);
    Z_boot = Z_exp(:) + noise;
    [gk, ~] = TR_DRT_matlab2011(freq_vec, Z_boot, el);
    gamma_samples(:, k) = gk(:);
end

v = var(gamma_samples, 0, 2);
mean_var = mean(v);
end

function nvec = normalize_metric_2011(vec)
% normalize_metric_2011 - Min-max normalization with robust fallback.

vec = vec(:);
vmin = min(vec);
vmax = max(vec);
if vmax > vmin
    nvec = (vec - vmin) / (vmax - vmin);
else
    nvec = zeros(size(vec));
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

function [fitted_peaks, rbf_model] = rbf_peak_fitting_2011(tau_vals, gamma_vals, peaks, n_peaks)
% rbf_peak_fitting_2011 - Fit RBF model to detected peaks and extract circuit parameters
%
% This function performs detailed analysis of detected peaks in the DRT curve:
%   1. Creates a radial basis function (RBF) interpolant of the entire DRT
%   2. For each detected peak, refines its position and extracts electrical parameters
%   3. Calculates equivalent circuit resistance R and capacitance C from peak properties
%
% The RBF model provides smooth interpolation between measurement points, allowing
% better peak characterization and fine-tuning of peak positions/widths.
%
% INPUTS:
%   tau_vals - Relaxation times (seconds) [N x 1]
%   gamma_vals - DRT values corresponding to tau_vals [N x 1]
%   peaks - Indices of detected peaks in gamma_vals [M x 1]
%   n_peaks - Number of top peaks to fit (default: all peaks in 'peaks')
%
% OUTPUTS:
%   fitted_peaks - Struct array with fields: peak_id, amplitude, tau, width, R, C
%   rbf_model - RBF model structure for future evaluation
%
% PEAK OUTPUT STRUCTURE:
%   .peak_id - Peak identifier (1, 2, 3, ...)
%   .amplitude - Peak height in DRT (Ohm)
%   .tau - Relaxation time at refined peak position (seconds)
%   .width - Peak width estimate (log-space standard deviation)
%   .R - Estimated resistance of this process (Ohm)
%   .C - Estimated capacitance of this process (Farad)

if nargin < 4 || isempty(n_peaks)
    n_peaks = numel(peaks);
end

tau_vals = tau_vals(:);
gamma_vals = gamma_vals(:);
tau_log = log(tau_vals);

rbf_model = create_rbf_model_2011(tau_log, gamma_vals);

n_fit = min(n_peaks, numel(peaks));
if n_fit == 0
    fitted_peaks = [];
    return;
end

fitted_peaks = repmat(struct('peak_id', 0, 'amplitude', 0, 'tau', 0, 'width', 0.5, 'R', 0, 'C', 0), n_fit, 1);

for i = 1:n_fit
    peak_idx = peaks(i);
    tau_peak = tau_vals(peak_idx);
    gamma_peak = gamma_vals(peak_idx);

    tau_local = logspace(log10(tau_peak/5), log10(tau_peak*5), 50).';
    tau_local_log = log(tau_local);

    gamma_interp = rbf_eval_2011(rbf_model, tau_local_log);
    [gamma_max, max_idx] = max(gamma_interp);
    tau_max = tau_local(max_idx);

    peak_region = gamma_interp > (0.1 * gamma_max);
    if sum(peak_region) > 3
        width = std(log(tau_local(peak_region)));
    else
        width = 0.5;
    end

    amp = max(gamma_max, gamma_peak);
    if amp > 0
        R_val = rbf_model.w(peak_idx)*sqrt(pi)./abs(rbf_model.epsilon);
        C_val = tau_peak/R_val;
        n_val = 1.144./(1+width);
    else
        R_val = 1e6;
        C_val = 0;
        n_val = 1;
    end

    fitted_peaks(i).peak_id = i;
    fitted_peaks(i).amplitude = amp;
    fitted_peaks(i).tau = tau_max;
    fitted_peaks(i).width = n_val;
    fitted_peaks(i).R = R_val;
    fitted_peaks(i).C = C_val;
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
    'R', 0, 'C', 0, 'n', 0, 'is_hidden', 0, 'baseline', 0), m, 1);

for peak_idx_local = 1:m
    tau_peak = exp(centers(peak_idx_local));
    amp = amplitudes(peak_idx_local);
    sig = sigma_vec(peak_idx_local);
    R_est = amp * sig * sqrt(2 * pi);
    C_est = tau_peak / max(R_est, eps);
    n_est = 1.144 / (1 + sig);

    fitted_peaks(peak_idx_local).peak_id = peak_idx_local;
    fitted_peaks(peak_idx_local).tau = tau_peak;
    fitted_peaks(peak_idx_local).amplitude = amp;
    fitted_peaks(peak_idx_local).sigma = sig;
    fitted_peaks(peak_idx_local).R = R_est;
    fitted_peaks(peak_idx_local).C = C_est;
    fitted_peaks(peak_idx_local).n = n_est;
    fitted_peaks(peak_idx_local).is_hidden = hidden_flag(peak_idx_local);
    fitted_peaks(peak_idx_local).baseline = baseline;
end

function [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers)
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

function kk_result = perform_kk_test_matlab2011(freq_vec, Z_exp)
% perform_kk_test_matlab2011 - Estimate KK-consistent impedance and residual metrics

freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
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
b = [real(Z_exp); imag(Z_exp); zeros(n_basis,1)];

x = lsqnonneg(M, b);
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

function sensitivity_result = perform_sensitivity_analysis_matlab2011(freq_vec, Z_exp, el, gamma_ref, Z_cal, seed_offset)
% perform_sensitivity_analysis_matlab2011 - Bootstrap DRT stability at the selected lambda

freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
gamma_ref = gamma_ref(:);
Z_cal = Z_cal(:);
n_freq = numel(freq_vec);
n_boot = 24;

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

sensitivity_result = struct('n_boot', n_boot, 'gamma_samples', gamma_samples, ...
    'gamma_mean', gamma_mean, 'gamma_std', gamma_std, 'gamma_lower', gamma_lower, ...
    'gamma_upper', gamma_upper, 'mean_std', mean(gamma_std(active_mask)), ...
    'max_std', max(gamma_std), 'mean_cv_pct', mean(cv_pct), ...
    'max_cv_pct', max(cv_pct), 'relative_band_area_pct', band_area_pct);
end

function result = perform_lambda_perturbation_sensitivity_matlab2011(freq_vec, Z_exp, lambda_ref, reference_mean_resid)
% perform_lambda_perturbation_sensitivity_matlab2011 - Residual sensitivity around selected lambda

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

xlswrite(file_name, header_row, sheet_name, 'A1');
xlswrite(file_name, numeric_data, sheet_name, 'A2');
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




