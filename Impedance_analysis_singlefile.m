% Impedance_analysis_singlefile - streamlined MATLAB 2011 workflow
% Single-file DRT analysis with temperature-tagged post-processing exports.

clc;
close all;

%% Configuration
% Purpose:
% - Resolve paths, derive a clean sample tag from the input file name, and
%   create output folders/workbook paths used by the rest of the workflow.
% Outputs created here are naming-only; no analysis runs in this cell.
% Inputs:
% - Hardcoded workbook name in this script (input_path).
% - Current script location (mfilename) for repo_root.
% Outputs:
% - repo_root, sample_tag, output_dir, fig_dir, excel_output_file.
% - Verified existence check for input_path.
repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422Sap.xlsx');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[~, sample_base, ~] = fileparts(input_path);
sample_tag = lower(regexprep(sample_base, '[^a-zA-Z0-9]', ''));
if isempty(sample_tag)
    sample_tag = 'sample';
end

output_dir = fullfile(repo_root, [sample_tag '_postprocess']);
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

fig_dir = fullfile(output_dir, 'figures_eps');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

excel_output_file = fullfile(output_dir, [sample_tag '_postprocess_data.xlsx']);

fprintf('Input workbook : %s\n', input_path);
fprintf('Output folder  : %s\n', output_dir);
fprintf('Output workbook: %s\n', excel_output_file);
fprintf('Running MATLAB 2011-compatible DRT workflow from: %s\n', repo_root);

%% Load workbook and parse temperatures
% Purpose:
% - Read numeric EIS blocks and parse temperature labels from the header row.
% - Validate that columns are grouped as [freq, |Z|, phase] per temperature.
% This cell establishes the dataset count and temperature vector used globally.
% Inputs:
% - input_path from the Configuration cell.
% - Workbook format assumption: repeating [freq, |Z|, phase] columns.
% Outputs:
% - indata numeric matrix, temperature vector, and num_sets consistency checks.
try
    % Prefer MATLAB-native readers first (no Excel COM dependency).
    indata = readmatrix(input_path, 'Sheet', 1);
    text_hdr = readcell(input_path, 'Sheet', 1, 'Range', '1:1');
catch ME
    try
        % Fallback: resolve first sheet name and read via xlsread.
        [~, sheet_names] = xlsfinfo(input_path);
        if isempty(sheet_names)
            sheet_id = 1;
        else
            sheet_id = sheet_names{1};
        end
        [indata, text_hdr] = xlsread(input_path, sheet_id);
    catch
        try
            % Final fallback: basic mode (avoids Excel automation issues).
            [num_data, txt_hdr, raw_hdr] = xlsread(input_path, 1, 'basic');
            indata = num_data;
            text_hdr = txt_hdr;
            if isempty(text_hdr) && ~isempty(raw_hdr)
                text_hdr = raw_hdr(1, :);
            end
        catch ME2
            error('Failed to read workbook: %s', ME2.message);
        end
    end
end

% Drop rows that are completely non-numeric after import.
if ~isempty(indata)
    indata = indata(~all(isnan(indata), 2), :);
end

if isempty(indata)
    error('No numeric data found in workbook.');
end

temperature = zeros(0, 1);
if ~isempty(text_hdr)
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

sz = size(indata, 2);
if rem(sz, 3) ~= 0
    error('Data must be organized as [freq, |Z|, phase] per temperature block.');
end

num_sets = sz / 3;
if numel(temperature) < num_sets
    error('Temperature headers are incomplete (%d found, %d expected).', numel(temperature), num_sets);
end
temperature = temperature(1:num_sets);
fprintf('Loaded %d temperature(s) from: %s\n', numel(temperature), input_path);

%% Master arrays for active summary outputs
% Purpose:
% - Preallocate all matrices that collect per-temperature results.
% - Define the fixed lambda grid and noise levels used by summary figures.
% Keeping these arrays preallocated avoids repeated memory growth in loops.
% Inputs:
% - num_sets and temperature from the previous cell.
% - Chosen analysis settings (lambda grid and noise levels) defined here.
% Outputs:
% - Preallocated summary/result matrices and all_temps_peak_data container.
% - Shared configuration vectors: lambda_values_master, noise_injection_pct_master.
lambda_values_master = logspace(-7, 0, 10).';
num_l_master = numel(lambda_values_master);
noise_injection_pct_master = [0.5 1 2 5 10];

ngcv_imag_matrix = nan(num_sets, num_l_master);
ngcv_real_matrix = nan(num_sets, num_l_master);
combined_residual_lambda_matrix = nan(num_sets, num_l_master);
selected_lambda_per_temp = nan(num_sets, 1);
selected_lambda_idx_per_temp = nan(num_sets, 1);
nyquist_mean_residual_per_temp = nan(num_sets, 1);
nyquist_max_residual_per_temp = nan(num_sets, 1);
kk_summary = nan(num_sets, 13);
noise_mean_residual_pct_matrix = nan(num_sets, numel(noise_injection_pct_master));

% Kept because active RGB plot uses per-temperature peak results.
all_temps_peak_data = cell(num_sets, 1);

% Deferred storage buffers: collect all per-temperature outputs first, then write once at end.
raw_all = [];
drt_all = [];
peak_all = [];
rcn_all = [];

raw_header = {'Temperature_K','freq_Hz','Zmag_ohm','phase_deg','Zre_ohm','Zim_ohm'};
drt_header = {'Temperature_K','tau_s','gamma_ohm','freq_Hz_equiv'};
peak_header = {'Temperature_K','peak_id','tau_s','amplitude','sigma_logtau','fwhm_logtau','R_ohm','C_F','n_est','is_hidden'};
rcn_header = {'Temperature_K','freq_Hz','Zre_exp','Zim_exp','Zre_drt','Zim_drt','Zre_rcn','Zim_rcn','abs_residual_drt','abs_residual_rcn'};

figure_export_jobs = cell(0, 2);
rgb_export_ready = false;
rgb_export_freq_axis = [];
rgb_export_temperature_axis = [];
rgb_export_r = [];
rgb_export_g = [];
rgb_export_b = [];

%% Main per-temperature workflow
% Purpose:
% - Process each temperature block end-to-end: cleaning, DRT inversion,
%   lambda selection, peak extraction, R-C-n reconstruction, KK metrics,
%   noise sensitivity, per-temperature figure export, and tagged Excel writes.
% This is the core compute cell and fills all summary matrices.
% Inputs:
% - indata, temperature, num_sets from workbook parsing.
% - lambda_values_master and noise_injection_pct_master from preallocation cell.
% - output paths: fig_dir and excel_output_file from Configuration.
% Outputs:
% - Populated per-temperature metrics in ngcv/KK/residual/noise summary arrays.
% - Temperature-tagged sheet exports (RAW/DRT/PEAK/RCN).
% - all_temps_peak_data entries used by the global RGB visualization cell.
for ab = 1:num_sets
    fprintf('\n================ Dataset %d/%d ================\n', ab, num_sets);
    fprintf('Temperature: %.1f K\n', temperature(ab));

    data = [indata(:,3*ab-2), indata(:,3*ab-1), indata(:,3*ab)];

    % Convert [freq, |Z|, phase(deg)] to complex impedance and clean rows.
    freq_vec = double(data(:,1));
    Z_mag = double(data(:,2));
    phase_deg = double(data(:,3));
    Z_exp = Z_mag .* exp(1i * (phase_deg * pi / 180));

    valid = freq_vec > 0;
    if ~all(valid)
        fprintf('Removing %d zero/negative frequency points\n', sum(~valid));
    end
    freq_vec = freq_vec(valid);
    Z_exp = Z_exp(valid);

    [freq_vec, sort_idx] = sort(freq_vec);
    Z_exp = Z_exp(sort_idx);

    % Remove contiguous inductive high-frequency tail.
    tail_start = numel(freq_vec) + 1;
    for k = numel(freq_vec):-1:1
        if imag(Z_exp(k)) > 0
            tail_start = k;
        else
            break;
        end
    end
    if tail_start <= numel(freq_vec)
        fprintf('Removed %d high-frequency inductive point(s) before fitting\n', numel(freq_vec) - tail_start + 1);
        freq_vec = freq_vec(1:tail_start-1);
        Z_exp = Z_exp(1:tail_start-1);
    end

    if isempty(freq_vec)
        warning('All points removed after inductive-tail filtering at %.1f K.', temperature(ab));
        continue;
    end

    n_data = numel(freq_vec);
    m_data = 2 * n_data;
    lsq_opts = optimset('Display', 'off', 'MaxIter', max(400, 10 * n_data), 'TolX', 1e-10);
    fprintf('Loaded %d frequency points (sorted by increasing frequency)\n', n_data);
    fprintf('Frequency range: %.1f - %.1f Hz\n', min(freq_vec), max(freq_vec));
    fprintf('Magnitude range: %.1f - %.1f Ohm\n', min(abs(Z_exp)), max(abs(Z_exp)));
    phase_deg_now = angle(Z_exp) * 180 / pi;
    fprintf('Phase range: %.2f - %.2f deg\n', min(phase_deg_now), max(phase_deg_now));

    % Precompute DRT projection matrices once per temperature.
    omega = 2 * pi * freq_vec(:);
    tau_base = 1 ./ freq_vec(:);
    A_re_mat = zeros(n_data, n_data);
    A_im_mat = zeros(n_data, n_data);
    for p = 1:n_data
        for q = 1:n_data
            if q == 1
                log_term = log(tau_base(q+1) / tau_base(q));
            elseif q == n_data
                log_term = log(tau_base(q) / tau_base(q-1));
            else
                log_term = log(tau_base(q+1) / tau_base(q-1));
            end
            den = 1 + (omega(p) * tau_base(q))^2;
            A_re_mat(p, q) = -0.5 * log_term / den;
            A_im_mat(p, q) = 0.5 * (omega(p) * tau_base(q)) * log_term / den;
        end
    end

    % Lambda sweep with imag-only and real-only nGCV scoring.
    fprintf('\nRunning GCV lambda selection (imaginary inversion, Re+Im cross-validation)...\n');
    fprintf('Lambda grid: logspace(-7, 0, 10) = %.1e to %.1e\n', lambda_values_master(1), lambda_values_master(end));
    num_l = numel(lambda_values_master);
    gcv_scores = inf(num_l, 1);
    gcv_scores_real = inf(num_l, 1);
    gamma_all = zeros(n_data, num_l);
    R_inf_all = zeros(num_l, 1);
    hat_trace_all = zeros(num_l, 1);
    hat_trace_all_real = zeros(num_l, 1);
    exp_impedance_norm = max(norm([real(Z_exp(:)); imag(Z_exp(:))], 2), eps);
    eye_n = eye(n_data);

    for i = 1:num_l
        lam = lambda_values_master(i);

        M_im = [A_im_mat; sqrt(lam/2) * eye_n];
        b_im = [imag(Z_exp(:)); zeros(n_data, 1)];
        gamma_i = lsqnonneg(M_im, b_im);
        R_inf_i = mean(real(Z_exp(:)) - A_re_mat * gamma_i);

        dz_re = real(Z_exp) - (R_inf_i + A_re_mat * gamma_i);
        dz_im = imag(Z_exp) - (A_im_mat * gamma_i);
        res_root = norm([dz_re; dz_im], 2);
        res_sq = res_root^2;
        combined_residual_lambda_matrix(ab, i) = 100 * res_root / exp_impedance_norm;

        alpha = lam / 2;
        S_im = (A_im_mat' * A_im_mat) + alpha * eye_n;
        tr_h = trace(S_im \ (A_im_mat' * A_im_mat));
        tr_h = max(0, min(tr_h, n_data - 1e-9));
        hat_trace_all(i) = tr_h;
        denom = 1 - tr_h / n_data;
        if abs(denom) > 1e-6
            gcv_scores(i) = (res_sq / m_data) / (denom^2);
        end

        M_re = [A_re_mat; sqrt(lam/2) * eye_n];
        b_re = [real(Z_exp(:)); zeros(n_data, 1)];
        gamma_r = lsqnonneg(M_re, b_re);
        R_inf_r = mean(real(Z_exp(:)) - A_re_mat * gamma_r);
        dz_re_r = real(Z_exp) - (R_inf_r + A_re_mat * gamma_r);
        dz_im_r = imag(Z_exp) - (A_im_mat * gamma_r);
        res_sq_r = norm([dz_re_r; dz_im_r], 2)^2;

        S_re = (A_re_mat' * A_re_mat) + alpha * eye_n;
        tr_h_r = trace(S_re \ (A_re_mat' * A_re_mat));
        tr_h_r = max(0, min(tr_h_r, n_data - 1e-9));
        hat_trace_all_real(i) = tr_h_r;
        denom_r = 1 - tr_h_r / n_data;
        if abs(denom_r) > 1e-6
            gcv_scores_real(i) = (res_sq_r / m_data) / (denom_r^2);
        end

        gamma_all(:, i) = gamma_i;
        R_inf_all(i) = R_inf_i;
        fprintf('  lambda = %8.2e | nGCV_im = %.4e | nGCV_re = %.4e\n', lam, gcv_scores(i), gcv_scores_real(i));
    end

    ngcv_imag_matrix(ab, :) = gcv_scores(:).';
    ngcv_real_matrix(ab, :) = gcv_scores_real(:).';

    % Flat-slope lambda selection (fallback to min nGCV).
    log_lambda = log(lambda_values_master);
    slope_seg = diff(log(max(gcv_scores, realmin))) ./ diff(log_lambda);
    flat_slope_thresh = 0.002;
    flat_seg = abs(slope_seg) <= flat_slope_thresh;
    if any(flat_seg)
        flat_idx = find(flat_seg) + 1;
        best_idx = flat_idx(end);
    else
        [~, best_idx] = min(gcv_scores);
    end

    el = lambda_values_master(best_idx);
    gamma = gamma_all(:, best_idx);
    R_inf = R_inf_all(best_idx);
    Z_cal = R_inf + A_re_mat * gamma + 1i * (A_im_mat * gamma);

    [gcv_best_score, gcv_min_idx] = min(gcv_scores);
    [gcv_best_score_real, gcv_min_idx_real] = min(gcv_scores_real);
    if any(flat_seg)
        selection_mode = 'max_lambda_on_flat_slope';
    else
        selection_mode = 'fallback_min_gcv';
    end

    residual = abs(Z_exp - Z_cal);
    selected_lambda_per_temp(ab) = el;
    selected_lambda_idx_per_temp(ab) = best_idx;
    nyquist_mean_residual_per_temp(ab) = mean(residual);
    nyquist_max_residual_per_temp(ab) = max(residual);

    fprintf('\nSelected lambda mode          : %s\n', selection_mode);
    fprintf('Selected lambda (imag nGCV)   : %.3e\n', el);
    fprintf('Selected nGCV (imag)          : %.4e\n', gcv_scores(best_idx));
    fprintf('Minimum nGCV lambda (imag)    : %.3e (score %.4e)\n', lambda_values_master(gcv_min_idx), gcv_best_score);
    fprintf('Minimum nGCV lambda (real)    : %.3e (score %.4e)\n', lambda_values_master(gcv_min_idx_real), gcv_best_score_real);
    fprintf('trace(H) imag at selected     : %.4f\n', hat_trace_all(best_idx));
    fprintf('trace(H) real at selected     : %.4f\n', hat_trace_all_real(best_idx));
    fprintf('Recovered R_inf = %.6f\n', R_inf);
    fprintf('DRT vector length = %d\n', numel(gamma));
    fprintf('\nMean absolute residual: %.6f\n', nyquist_mean_residual_per_temp(ab));
    fprintf('Max absolute residual : %.6f\n', nyquist_max_residual_per_temp(ab));

    % --------------------------- Peak analysis ---------------------------
    % Detect, refine, and parametrize the dominant DRT peaks.
    tau_data = 1 ./ (2 * pi * freq_vec(:));
    gamma_max = max(gamma);

    % Prominence-based peak detection (unbiased, no fixed count assumed)
    if gamma_max > 0
        min_height = 0.05 * gamma_max;
        peak_idx = find_peaks_simple_2023(gamma, min_height, 3);
        n_sig_peaks = numel(peak_idx);
    else
        n_sig_peaks = 0;
    end
    fprintf('Prominent peaks detected (prominence >= 5%% of max): %d\n', n_sig_peaks);

    [data_peaks_fitted, rbf_model, peak_refine] = refine_rbf_peak_fitting_2023(tau_data, gamma);

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
            freq_peak(i) = 1 / (2 * pi * peak.tau); %#ok<AGROW>
            tau_peak(i) = peak.tau; %#ok<AGROW>
            gamma_peak(i) = peak.amplitude; %#ok<AGROW>
            R_file(i) = peak.R; %#ok<AGROW>
            C_file(i) = peak.C; %#ok<AGROW>
            n_file(i) = peak.n; %#ok<AGROW>
        end
        fprintf('\nTotal series resistance: %.2f Ohm\n', total_R_data);

    end
    fprintf('High-frequency resistance (R_inf): %.2f Ohm\n', R_inf);

    n_peaks = numel(data_peaks_fitted);
    peak_numeric = zeros(n_peaks, 9);
    for i = 1:n_peaks
        peak_numeric(i,:) = [data_peaks_fitted(i).peak_id, data_peaks_fitted(i).tau, ...
            data_peaks_fitted(i).amplitude, data_peaks_fitted(i).sigma, data_peaks_fitted(i).fwhm_logtau, ...
            data_peaks_fitted(i).R, data_peaks_fitted(i).C, data_peaks_fitted(i).n, data_peaks_fitted(i).is_hidden];
    end

    tau_dense = peak_refine.tau_dense(:);
    gamma_rbf = peak_refine.gamma_rbf(:);

    Z_rcn = nan(size(Z_cal));
    rcn_residual = nan(size(residual));
    if ~isempty(data_peaks_fitted)
        [Z_rcn, Z_rcn_components] = calculate_eis_from_rcn_peaks_2023(freq_vec, R_inf, data_peaks_fitted);

        % Explicitly verify series assembly: Z_total = R_inf + sum_k Z_k.
        Z_rcn_from_series = R_inf + sum(Z_rcn_components, 2);
        series_mismatch = max(abs(Z_rcn - Z_rcn_from_series));

        rcn_residual = abs(Z_exp - Z_rcn);
        fprintf('R-C-n Nyquist mean residual: %.6f\n', mean(rcn_residual));
        fprintf('R-C-n Nyquist max residual : %.6f\n', max(rcn_residual));
        fprintf('R-C-n Nyquist mean relative residual: %.2f %%\n', 100 * mean(rcn_residual ./ max(abs(Z_exp), eps)));
        fprintf('Series-assembly mismatch (should be ~0): %.3e\n', series_mismatch);
    end

    % KK candidate sweep (compact but stable).
    n_freq = n_data;
    candidate_basis = unique([min(max(8, ceil(n_freq / 3)), 20), ...
                              min(max(10, ceil(n_freq / 2.5)), 24), ...
                              min(max(12, ceil(n_freq / 2)), 30)]);
    candidate_reg = [1e-4; 3e-4; 1e-3];
    best_total = inf;
    kk_best = struct();

    for ib = 1:numel(candidate_basis)
        n_basis = candidate_basis(ib);
        tau_min = 1 / (2 * pi * max(freq_vec));
        tau_max = 1 / (2 * pi * min(freq_vec));
        tau_basis = logspace(log10(tau_min / 5), log10(tau_max * 5), n_basis).';

        B_re = zeros(n_freq, n_basis);
        B_im = zeros(n_freq, n_basis);
        for p = 1:n_freq
            for q = 1:n_basis
                ot = omega(p) * tau_basis(q);
                den = 1 + ot^2;
                B_re(p, q) = 1 / den;
                B_im(p, q) = -ot / den;
            end
        end

        for ir = 1:numel(candidate_reg)
            regw = candidate_reg(ir);
            M = [B_re, ones(n_freq,1), zeros(n_freq,1); ...
                 B_im, zeros(n_freq,1), omega; ...
                 sqrt(regw) * eye(n_basis), zeros(n_basis,2)];
            b = [real(Z_exp); imag(Z_exp); zeros(n_basis,1)];
            x = lsqnonneg(M, b, lsq_opts);

            R_branch = x(1:n_basis);
            R_inf_kk = x(n_basis + 1);
            L_series = x(n_basis + 2);
            Z_fit = R_inf_kk + B_re * R_branch + 1i * (B_im * R_branch + omega * L_series);
            res_kk = Z_exp - Z_fit;

            total_norm = max(norm([real(Z_exp); imag(Z_exp)]), eps);
            total_pct = 100 * norm([real(res_kk); imag(res_kk)]) / total_norm;

            if total_pct < best_total
                best_total = total_pct;
                kk_best.Z_fit = Z_fit;
                kk_best.residual = res_kk;
                kk_best.reg_weight = regw;
                kk_best.n_basis = n_basis;
            end
        end
    end

    real_norm = max(norm(real(Z_exp)), eps);
    imag_norm = max(norm(imag(Z_exp)), eps);
    total_norm = max(norm([real(Z_exp); imag(Z_exp)]), eps);
    rel_point_pct = 100 * abs(kk_best.residual) ./ max(abs(Z_exp), eps);

    kk_real_pct = 100 * norm(real(kk_best.residual)) / real_norm;
    kk_imag_pct = 100 * norm(imag(kk_best.residual)) / imag_norm;
    kk_total_pct = 100 * norm([real(kk_best.residual); imag(kk_best.residual)]) / total_norm;
    kk_max_pct = max(rel_point_pct);

    kk_drt_res = Z_cal - kk_best.Z_fit;
    kk_drt_real_pct = 100 * norm(real(kk_drt_res)) / max(norm(real(Z_cal)), eps);
    kk_drt_imag_pct = 100 * norm(imag(kk_drt_res)) / max(norm(imag(Z_cal)), eps);
    kk_drt_total_pct = 100 * norm([real(kk_drt_res); imag(kk_drt_res)]) / max(norm([real(Z_cal); imag(Z_cal)]), eps);
    kk_drt_max_pct = max(100 * abs(kk_drt_res) ./ max(abs(Z_cal), eps));

    kk_status_code = -1;
    if kk_total_pct <= 5
        kk_status_code = 1;
    elseif kk_total_pct <= 10
        kk_status_code = 0;
    end

    kk_drt_status_code = -1;
    if kk_drt_total_pct <= 5
        kk_drt_status_code = 1;
    elseif kk_drt_total_pct <= 10
        kk_drt_status_code = 0;
    end

    kk_summary(ab, :) = [temperature(ab), el, R_inf, kk_real_pct, kk_imag_pct, kk_total_pct, kk_max_pct, kk_status_code, ...
                         kk_drt_real_pct, kk_drt_imag_pct, kk_drt_total_pct, kk_drt_max_pct, kk_drt_status_code];
    kk_status = 'fail';
    if kk_status_code == 1
        kk_status = 'pass';
    elseif kk_status_code == 0
        kk_status = 'warning';
    end

    kk_drt_status = 'fail';
    if kk_drt_status_code == 1
        kk_drt_status = 'pass';
    elseif kk_drt_status_code == 0
        kk_drt_status = 'warning';
    end

    fprintf('\nKramers-Kronig consistency test:\n');
    fprintf('  KK basis / reg    : %d / %.1e\n', kk_best.n_basis, kk_best.reg_weight);
    fprintf('  Real residual  : %.2f %%\n', kk_real_pct);
    fprintf('  Imag residual  : %.2f %%\n', kk_imag_pct);
    fprintf('  Total residual : %.2f %%\n', kk_total_pct);
    fprintf('  Status         : %s\n', kk_status);

    fprintf('\nKramers-Kronig test on DRT-generated impedance:\n');
    fprintf('  KK basis / reg    : %d / %.1e\n', kk_best.n_basis, kk_best.reg_weight);
    fprintf('  Real residual  : %.2f %%\n', kk_drt_real_pct);
    fprintf('  Imag residual  : %.2f %%\n', kk_drt_imag_pct);
    fprintf('  Total residual : %.2f %%\n', kk_drt_total_pct);
    fprintf('  Status         : %s\n', kk_drt_status);

    % Lightweight noise sensitivity around selected fit (faster than bootstrap).
    n_trials = 4;
    n_levels = numel(noise_injection_pct_master);
    noise_mean_residual_pct = zeros(n_levels, 1);
    rand('twister', 100 + ab);
    randn('state', 200 + ab);
    for nl = 1:n_levels
        noise_pct = noise_injection_pct_master(nl);
        mean_trials = zeros(n_trials, 1);
        for tr = 1:n_trials
            noise_scale = (noise_pct / 100) * max(abs(Z_exp), eps);
            noise_vec = (noise_scale ./ sqrt(2)) .* (randn(n_data, 1) + 1i * randn(n_data, 1));
            Z_noisy = Z_exp + noise_vec;
            % Reuse selected lambda and prebuilt matrices for speed.
            M_im = [A_im_mat; sqrt(el/2) * eye_n];
            b_im = [imag(Z_noisy(:)); zeros(n_data, 1)];
            gamma_noisy = lsqnonneg(M_im, b_im);
            R_inf_noisy = mean(real(Z_noisy(:)) - A_re_mat * gamma_noisy);
            Z_fit_noisy = R_inf_noisy + A_re_mat * gamma_noisy + 1i * (A_im_mat * gamma_noisy);
            resid_pct = 100 * abs(Z_noisy - Z_fit_noisy) ./ max(abs(Z_noisy), eps);
            mean_trials(tr) = mean(resid_pct);
        end
        noise_mean_residual_pct(nl) = mean(mean_trials);
    end
    noise_mean_residual_pct_matrix(ab, :) = noise_mean_residual_pct(:).';
    fprintf('\nSensitivity analysis (%d runs):\n', n_trials * n_levels);
    fprintf('  Mean residual at %.1f%% noise : %.4f %%\n', noise_injection_pct_master(1), noise_mean_residual_pct(1));
    fprintf('  Mean residual at %.1f%% noise : %.4f %%\n', noise_injection_pct_master(end), noise_mean_residual_pct(end));
    fprintf('Run complete at %s\n', datestr(now, 31));

    %% Store data for RGB visualization
    % Purpose:
    % - Convert per-temperature peak outputs into a compact structure used
    %   later to build the cross-temperature RGB map.
    % A dense tau grid is stored so RGB reconstruction stays smooth.
    % Inputs:
    % - tau, gamma, peak_numeric, n_peaks from current loop pass.
    % - temperature(ab) to tag this temperature entry.
    % Outputs:
    % - all_temps_peak_data{ab} populated with dense tau/gamma and peak structs.
    % Use dense tau grid for smoother per-temperature color profiles.
    gamma_dense = gamma_rbf;

    peaks_struct = repmat(struct('tau', 0, 'amplitude', 0, 'sigma', 0), n_peaks, 1);
    for i = 1:n_peaks
        peaks_struct(i).tau = peak_numeric(i,2);
        peaks_struct(i).amplitude = peak_numeric(i,3);
        peaks_struct(i).sigma = peak_numeric(i,4);
    end

    all_temps_peak_data{ab} = struct('temperature', temperature(ab), ...
                                     'tau_dense', tau_dense, ...
                                     'gamma_dense', gamma_dense, ...
                                     'peaks', peaks_struct, ...
                                     'peak_count', n_peaks);

    %% Temperature-tagged Excel export buffers for post processing
    % Purpose:
    % - Accumulate raw impedance, DRT values, fitted peaks, and R-C-n
    %   reconstruction rows into consolidated arrays.
    % - Final workbook write is deferred to the end; no temperature-wise sheets.
    % Inputs:
    % - Current-loop vectors/matrices: freq_vec, Z_exp, tau, gamma,
    %   peak_numeric, Z_cal, Z_rcn, residual, rcn_residual, temperature(ab).
    % Outputs:
    % - Consolidated arrays raw_all, drt_all, peak_all, rcn_all.
    tcol_raw = repmat(temperature(ab), numel(freq_vec), 1);
    raw_all = [raw_all; tcol_raw, freq_vec(:), abs(Z_exp(:)), angle(Z_exp(:))*180/pi, real(Z_exp(:)), imag(Z_exp(:))]; %#ok<AGROW>

    tcol_drt = repmat(temperature(ab), numel(tau_data), 1);
    drt_all = [drt_all; tcol_drt, tau_data(:), gamma(:), 1 ./ (2*pi*tau_data(:))]; %#ok<AGROW>

    if ~isempty(peak_numeric)
        tcol_peak = repmat(temperature(ab), size(peak_numeric, 1), 1);
        peak_all = [peak_all; tcol_peak, peak_numeric]; %#ok<AGROW>
    end

    tcol_rcn = repmat(temperature(ab), numel(freq_vec), 1);
    rcn_all = [rcn_all; tcol_rcn, freq_vec(:), real(Z_exp(:)), imag(Z_exp(:)), real(Z_cal(:)), imag(Z_cal(:)), real(Z_rcn(:)), imag(Z_rcn(:)), residual(:), rcn_residual(:)]; %#ok<AGROW>
end

%% Active figure 2: RGB multi-temperature peak visualization
% Purpose:
% - Build one RGB image over all temperatures where color channels encode
%   top-ranked peaks; map axes to temperature and log-frequency view.
% - Export vectorized EPS and store RGB channel matrices for post-processing.
% Inputs:
% - all_temps_peak_data collected in the main loop.
% - temperature vector for sorting/axis labels.
% - fig_dir, sample_tag, excel_output_file for exporting image and RGB arrays.
% Outputs:
% - RGB EPS figure and RGB_map_data sheet containing axes and R/G/B matrices.
fprintf('\n========== Multi-Temperature RGB Peak Visualization ==========\n');

num_temps = numel(all_temps_peak_data);
if num_temps == 0 || isempty(all_temps_peak_data{1})
    fprintf('Warning: no temperature data to visualize\n');
else
    tau_dense_1 = all_temps_peak_data{1}.tau_dense(:);
    ln_tau_master = log(tau_dense_1);
    ln_tau_min = min(ln_tau_master);
    ln_tau_max = max(ln_tau_master);

    rgb_image = zeros(num_temps, 256, 3, 'uint8');
    ln_tau_pixels = linspace(ln_tau_min, ln_tau_max, 256);

    for temp_idx = 1:num_temps
        peak_data = all_temps_peak_data{temp_idx};
        if isempty(peak_data)
            continue;
        end
        peaks = peak_data.peaks;
        n_peaks = peak_data.peak_count;
        if n_peaks == 0
            continue;
        end

        peak_rank_metric = zeros(n_peaks, 1);
        for p = 1:n_peaks
            peak_rank_metric(p) = peaks(p).tau;
        end
        [~, ord] = sort(peak_rank_metric, 'ascend');
        top_3_idx = ord(1:min(3, n_peaks));

        ln_tau_t = log(peak_data.tau_dense(:));
        channels = zeros(3, 256);
        for rank = 1:min(3, numel(top_3_idx))
            pk = peaks(top_3_idx(rank));
            gvals = pk.amplitude * exp(-((ln_tau_t - log(pk.tau)).^2) / (2 * pk.sigma^2));
            px = interp1(ln_tau_t, gvals, ln_tau_pixels, 'linear', 0);
            pxmax = max(px);
            if pxmax > 0
                channels(rank, :) = 255 * px / pxmax;
            end
        end

        rgb_image(temp_idx, :, 1) = uint8(round(channels(1, :)));
        rgb_image(temp_idx, :, 2) = uint8(round(channels(2, :)));
        rgb_image(temp_idx, :, 3) = uint8(round(channels(3, :)));
    end

    temperature_axis = temperature(:)';
    [temperature_axis, temp_order] = sort(temperature_axis, 'ascend');
    rgb_display = permute(rgb_image(temp_order, :, :), [2 1 3]);

    ln_tau_pixel_centers = linspace(ln_tau_min, ln_tau_max, size(rgb_display, 1)).';
    freq_axis = exp(-ln_tau_pixel_centers);
    [freq_axis, freq_order] = sort(freq_axis, 'ascend');
    rgb_display = rgb_display(freq_order, :, :);

    row_has_signal = any(reshape(rgb_display, size(rgb_display, 1), []), 2);
    removed_rows = sum(~row_has_signal);
    rgb_display = rgb_display(row_has_signal, :, :);
    freq_axis = freq_axis(row_has_signal);
    if removed_rows > 0
        fprintf('RGB cleanup removed %d all zero rows\n', removed_rows);
    end

    if ~isempty(freq_axis)
        target_rows = 1000;
        target_cols = 1000;
        src_rows = size(rgb_display, 1);
        src_cols = size(rgb_display, 2);
        row_idx = round(linspace(1, src_rows, target_rows));
        col_idx = round(linspace(1, src_cols, target_cols));
        rgb_display = rgb_display(row_idx, col_idx, :);
        freq_axis = freq_axis(row_idx);
        fprintf('RGB resize: %dx%d -> %dx%d (nearest-neighbor)\n', src_rows, src_cols, target_rows, target_cols);
        fprintf('RGB debug: size=%dx%dx%d, min=%d, max=%d, nnz=%d\n', ...
            size(rgb_display,1), size(rgb_display,2), size(rgb_display,3), min(rgb_display(:)), max(rgb_display(:)), nnz(rgb_display));

        fig_rgb = figure('Color', 'w', 'Name', 'RGB peak visualization');
        ax_rgb = axes('Parent', fig_rgb);
        imshow(rgb_display, 'Parent', ax_rgb);
        set(ax_rgb, 'Visible', 'on', 'YDir', 'normal');
        axis(ax_rgb, 'on');
        box(ax_rgb, 'on');
        xlabel(ax_rgb, 'Temperature (K)');
        ylabel(ax_rgb, 'Frequency (Hz)');
        title(ax_rgb, 'RGB peak visualization');

        n_cols = size(rgb_display, 2);
        n_rows = size(rgb_display, 1);

        x_major_vals = 1:10:280;
        x_major_pix = 1 + ((x_major_vals - 1) / (280 - 1)) * (n_cols - 1);
        x_major_valid = x_major_pix >= 1 & x_major_pix <= n_cols;
        x_major_vals = x_major_vals(x_major_valid);
        x_major_pix = x_major_pix(x_major_valid);
        set(ax_rgb, 'XTick', x_major_pix, 'XTickLabel', arrayfun(@(v) sprintf('%d', v), x_major_vals, 'UniformOutput', false));

        x_minor_vals = 1:1:280;
        x_minor_vals(mod(x_minor_vals - 1, 10) == 0) = [];
        x_minor_pix = 1 + ((x_minor_vals - 1) / (280 - 1)) * (n_cols - 1);
        x_minor_valid = x_minor_pix >= 1 & x_minor_pix <= n_cols;
        x_minor_pix = x_minor_pix(x_minor_valid);

        fmin = min(freq_axis);
        fmax = max(freq_axis);
        y_minor_pix = [];
        if fmin > 0 && fmax > fmin
            logfmin = log10(fmin);
            logfmax = log10(fmax);
            y_major_exp = ceil(logfmin):floor(logfmax);
            if isempty(y_major_exp)
                y_major_exp = floor(logfmin):ceil(logfmax);
            end
            y_major_vals = 10 .^ y_major_exp;
            y_major_pix = 1 + ((log10(y_major_vals) - logfmin) / (logfmax - logfmin)) * (n_rows - 1);
            y_major_valid = y_major_pix >= 1 & y_major_pix <= n_rows;
            y_major_exp = y_major_exp(y_major_valid);
            y_major_pix = y_major_pix(y_major_valid);
            set(ax_rgb, 'YTick', y_major_pix, 'YTickLabel', arrayfun(@(e) sprintf('10^{%d}', e), y_major_exp, 'UniformOutput', false));

            y_minor_vals = [];
            for e = floor(logfmin):ceil(logfmax)-1
                for m = 1:9
                    y_minor_vals(end+1,1) = 10^(e + m/10); %#ok<AGROW>
                end
            end
            y_minor_pix = 1 + ((log10(y_minor_vals) - logfmin) / (logfmax - logfmin)) * (n_rows - 1);
            y_minor_valid = y_minor_pix >= 1 & y_minor_pix <= n_rows;
            y_minor_pix = y_minor_pix(y_minor_valid);
        end

        set(ax_rgb, 'XMinorTick', 'on', 'YMinorTick', 'on');
        try
            ax_rgb.XRuler.MinorTickValues = x_minor_pix;
            ax_rgb.YRuler.MinorTickValues = y_minor_pix;
        catch
        end

        rgb_eps_file = fullfile(fig_dir, [sample_tag '_rgb_peak_map.eps']);
        figure_export_jobs(end+1, :) = {fig_rgb, rgb_eps_file}; %#ok<AGROW>
        fprintf('RGB peak visualization ready: %s (%d x %d x 3)\n', rgb_eps_file, size(rgb_display, 1), size(rgb_display, 2));

        % Defer RGB matrix write until final export stage.
        rgb_export_ready = true;
        rgb_export_freq_axis = freq_axis(:);
        rgb_export_temperature_axis = temperature_axis(:);
        rgb_export_r = double(rgb_display(:,:,1));
        rgb_export_g = double(rgb_display(:,:,2));
        rgb_export_b = double(rgb_display(:,:,3));
    end
end

%% Active figure 3: KK raw residual vs temperature
% Purpose:
% - Plot how KK total residual on raw data changes with temperature.
% This gives a quick quality trend check across the full sweep.
% Inputs:
% - kk_summary matrix and temperature vector from the main loop outputs.
% - fig_dir and sample_tag for vectorized figure export.
% Outputs:
% - KK residual trend figure exported as EPS.
fprintf('\n========== KK and Noise Summaries ==========\n');
[temperature_sorted, temp_sort_idx] = sort(temperature(:), 'ascend');
kk_temperature_sorted = kk_summary(temp_sort_idx, 1);
kk_total_raw_sorted = kk_summary(temp_sort_idx, 6);

fig_kk_raw = figure('Color', 'w', 'Name', 'KK raw residual vs temperature');
ax_kk_raw = axes('Parent', fig_kk_raw);
plot(ax_kk_raw, kk_temperature_sorted, kk_total_raw_sorted, 'ko-', 'LineWidth', 1.2, 'MarkerSize', 4);
grid(ax_kk_raw, 'on');
xlabel(ax_kk_raw, 'Temperature (K)');
ylabel(ax_kk_raw, 'KK total residual (%)');
title(ax_kk_raw, 'Kramers-Kronig raw-data total residual vs temperature');
figure_export_jobs(end+1, :) = {fig_kk_raw, fullfile(fig_dir, [sample_tag '_kk_raw_residual_vs_temp.eps'])}; %#ok<AGROW>

%% Active figure 4: noise mean residual colormap
% Purpose:
% - Visualize mean residual sensitivity versus injected noise level and
%   temperature using a compact 2D colormap.
% This summarizes robustness of selected-lambda fits.
% Inputs:
% - noise_mean_residual_pct_matrix from per-temperature calculations.
% - noise_injection_pct_master and temperature_sorted axis vectors.
% - fig_dir and sample_tag for export.
% Outputs:
% - Noise mean residual colormap exported as EPS.
noise_mean_residual_sorted = noise_mean_residual_pct_matrix(temp_sort_idx, :);
fig_noise_mean = figure('Color', 'w', 'Name', 'Noise mean residual colormap');
ax_noise_mean = axes('Parent', fig_noise_mean);
imagesc(ax_noise_mean, 1:numel(noise_injection_pct_master), temperature_sorted, noise_mean_residual_sorted);
set(ax_noise_mean, 'YDir', 'normal');
set(ax_noise_mean, 'XTick', 1:numel(noise_injection_pct_master), ...
    'XTickLabel', arrayfun(@(v) sprintf('%.1f', v), noise_injection_pct_master, 'UniformOutput', false));
xlabel(ax_noise_mean, 'Noise injection (%)');
ylabel(ax_noise_mean, 'Temperature (K)');
title(ax_noise_mean, 'Noise sensitivity mean residual colormap');
cb_noise = colorbar(ax_noise_mean);
ylabel(cb_noise, 'Mean residual (%)');
figure_export_jobs(end+1, :) = {fig_noise_mean, fullfile(fig_dir, [sample_tag '_noise_mean_residual_colormap.eps'])}; %#ok<AGROW>
fprintf('KK and noise sensitivity plots ready: %d temperature(s), %d noise levels\n', ...
    numel(temperature_sorted), numel(noise_injection_pct_master));

%% Active figure 5 and 6: nGCV imag/real colormaps
% Purpose:
% - Render temperature-vs-lambda nGCV maps for imaginary-only and real-only
%   inversion branches on the same lambda axis convention.
% These two figures are used to compare regularization behavior by method.
% Inputs:
% - ngcv_imag_matrix, ngcv_real_matrix from the main loop.
% - lambda_values_master and temperature_sorted for axes/ticks.
% - fig_dir and sample_tag for EPS exports.
% Outputs:
% - Two EPS colormap figures: imag nGCV and real nGCV.
fprintf('\n========== nGCV vs Lambda Colormaps ==========\n');
ngcv_imag_sorted = ngcv_imag_matrix(temp_sort_idx, :);
ngcv_real_sorted = ngcv_real_matrix(temp_sort_idx, :);
ngcv_imag_plot = log10(max(ngcv_imag_sorted, realmin));
ngcv_real_plot = log10(max(ngcv_real_sorted, realmin));

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

fig_ngcv_im = figure('Color', 'w', 'Name', 'Imag nGCV colormap');
ax_ngcv_im = axes('Parent', fig_ngcv_im);
imagesc(ax_ngcv_im, 1:n_lambda, temperature_sorted, ngcv_imag_plot);
set(ax_ngcv_im, 'YDir', 'normal');
xlim(ax_ngcv_im, [1 n_lambda]);
set(ax_ngcv_im, 'XTick', lambda_major_pix, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false));
set(ax_ngcv_im, 'XMinorTick', 'on');
try
    ax_ngcv_im.XRuler.MinorTickValues = lambda_minor_pix;
catch
end
xlabel(ax_ngcv_im, 'Lambda');
ylabel(ax_ngcv_im, 'Temperature (K)');
title(ax_ngcv_im, 'Imaginary inversion nGCV colormap');
cb1 = colorbar(ax_ngcv_im);
ylabel(cb1, 'log10(nGCV)');
figure_export_jobs(end+1, :) = {fig_ngcv_im, fullfile(fig_dir, [sample_tag '_ngcv_imag_colormap.eps'])}; %#ok<AGROW>

fig_ngcv_re = figure('Color', 'w', 'Name', 'Real nGCV colormap');
ax_ngcv_re = axes('Parent', fig_ngcv_re);
imagesc(ax_ngcv_re, 1:n_lambda, temperature_sorted, ngcv_real_plot);
set(ax_ngcv_re, 'YDir', 'normal');
xlim(ax_ngcv_re, [1 n_lambda]);
set(ax_ngcv_re, 'XTick', lambda_major_pix, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false));
set(ax_ngcv_re, 'XMinorTick', 'on');
try
    ax_ngcv_re.XRuler.MinorTickValues = lambda_minor_pix;
catch
end
xlabel(ax_ngcv_re, 'Lambda');
ylabel(ax_ngcv_re, 'Temperature (K)');
title(ax_ngcv_re, 'Real inversion nGCV colormap');
cb2 = colorbar(ax_ngcv_re);
ylabel(cb2, 'log10(nGCV)');
figure_export_jobs(end+1, :) = {fig_ngcv_re, fullfile(fig_dir, [sample_tag '_ngcv_real_colormap.eps'])}; %#ok<AGROW>
fprintf('nGCV colormaps ready: imag matrix %dx%d, real matrix %dx%d\n', ...
    size(ngcv_imag_sorted, 1), size(ngcv_imag_sorted, 2), ...
    size(ngcv_real_sorted, 1), size(ngcv_real_sorted, 2));

%% Active figure 7: combined residual-vs-lambda colormap
% Purpose:
% - Show percent-normalized combined residual across lambda and temperature,
%   with selected lambda overlaid as markers.
% This is the main map for final lambda-vs-error interpretation.
% Inputs:
% - combined_residual_lambda_matrix and selected_lambda_idx_per_temp.
% - lambda axis/ticks built in the previous colormap cell.
% - fig_dir and sample_tag for EPS export.
% Outputs:
% - Combined residual map with selected-lambda overlay exported as EPS.
fprintf('\n========== Combined Residual vs Lambda Colormap ==========\n');
combined_residual_lambda_sorted = combined_residual_lambda_matrix(temp_sort_idx, :);
fig_lambda_resid = figure('Color', 'w', 'Name', 'Combined residual vs lambda');
ax_lambda_resid = axes('Parent', fig_lambda_resid);
imagesc(ax_lambda_resid, 1:n_lambda, temperature_sorted, combined_residual_lambda_sorted);
set(ax_lambda_resid, 'YDir', 'normal');
xlim(ax_lambda_resid, [1 n_lambda]);
set(ax_lambda_resid, 'XTick', lambda_major_pix, ...
    'XTickLabel', arrayfun(@(v) sprintf('%.0e', v), lambda_major_vals, 'UniformOutput', false));
set(ax_lambda_resid, 'XMinorTick', 'on');
try
    ax_lambda_resid.XRuler.MinorTickValues = lambda_minor_pix;
catch
end
xlabel(ax_lambda_resid, 'Lambda');
ylabel(ax_lambda_resid, 'Temperature (K)');
title(ax_lambda_resid, 'Root combined residual percent vs lambda and temperature');
cb_lambda = colorbar(ax_lambda_resid);
ylabel(cb_lambda, '100 * ||[Re residual; Im residual]||_2 / ||Z_{exp}||_2');
hold(ax_lambda_resid, 'on');
plot(ax_lambda_resid, selected_lambda_idx_per_temp(temp_sort_idx), temperature_sorted, 'ks', ...
    'MarkerSize', 5, 'LineWidth', 1.0, 'MarkerFaceColor', 'w');
hold(ax_lambda_resid, 'off');
figure_export_jobs(end+1, :) = {fig_lambda_resid, fullfile(fig_dir, [sample_tag '_combined_residual_lambda_temp.eps'])}; %#ok<AGROW>
fprintf('Selected lambda/residual map ready: %d temperature(s), lambda range %.1e to %.1e\n', ...
    numel(temperature_sorted), min(selected_lambda_per_temp), max(selected_lambda_per_temp));
fprintf('Combined residual percent range: %.4f %% to %.4f %%\n', ...
    min(combined_residual_lambda_sorted(:)), max(combined_residual_lambda_sorted(:)));

%% Export active summary matrices to Excel
% Purpose:
% - Persist cross-temperature summary tables used by downstream analysis.
% - Includes selected lambda trends, nGCV matrices, residual map values,
%   noise-response summaries, and KK diagnostics.
% This cell is the final aggregation/export stage.
% Inputs:
% - Final populated summary matrices: selected_lambda_per_temp,
%   nyquist residual vectors, ngcv matrices, combined residual map,
%   noise_mean_residual_pct_matrix, kk_summary.
% - excel_output_file destination path.
% Outputs:
% - SUMMARY, NGCV_IMAG, NGCV_REAL, RESID_LAMBDA, NOISE_MEAN, and
%   KK_SUMMARY sheets written to excel_output_file.

% Prepare output workbook path; if locked, switch to a timestamped fallback.
[excel_output_file, switched_output_file] = prepare_excel_output_path_2023(excel_output_file);
if switched_output_file
    fprintf('Output workbook was locked. Writing to fallback file instead:\n  %s\n', excel_output_file);
end

% Write consolidated tagged sheets (no temperature-wise sheet splitting).
write_excel_block_2023(excel_output_file, 'RAW_ALL', 'A1', raw_header);
if ~isempty(raw_all)
    write_excel_block_2023(excel_output_file, 'RAW_ALL', 'A2', raw_all);
end

write_excel_block_2023(excel_output_file, 'DRT_ALL', 'A1', drt_header);
if ~isempty(drt_all)
    write_excel_block_2023(excel_output_file, 'DRT_ALL', 'A2', drt_all);
end

write_excel_block_2023(excel_output_file, 'PEAK_ALL', 'A1', peak_header);
if ~isempty(peak_all)
    write_excel_block_2023(excel_output_file, 'PEAK_ALL', 'A2', peak_all);
end

write_excel_block_2023(excel_output_file, 'RCN_ALL', 'A1', rcn_header);
if ~isempty(rcn_all)
    write_excel_block_2023(excel_output_file, 'RCN_ALL', 'A2', rcn_all);
end

% Write deferred RGB matrix/axis export after all calculations.
if rgb_export_ready
    rgb_sheet = 'RGB_map_data';
    write_excel_block_2023(excel_output_file, rgb_sheet, 'A1', {'freq_Hz_axis'});
    write_excel_block_2023(excel_output_file, rgb_sheet, 'A2', rgb_export_freq_axis);

    write_excel_block_2023(excel_output_file, rgb_sheet, 'C1', {'temp_K_axis'});
    write_excel_block_2023(excel_output_file, rgb_sheet, 'C2', rgb_export_temperature_axis);

    write_excel_block_2023(excel_output_file, rgb_sheet, 'E1', {'RGB_R'});
    write_excel_block_2023(excel_output_file, rgb_sheet, 'E2', rgb_export_r);
    write_excel_block_2023(excel_output_file, rgb_sheet, 'AN1', {'RGB_G'});
    write_excel_block_2023(excel_output_file, rgb_sheet, 'AN2', rgb_export_g);
    write_excel_block_2023(excel_output_file, rgb_sheet, 'BW1', {'RGB_B'});
    write_excel_block_2023(excel_output_file, rgb_sheet, 'BW2', rgb_export_b);
end

write_excel_block_2023(excel_output_file, 'SUMMARY', 'A1', {'temperature_K','lambda_selected','nyq_mean_residual','nyq_max_residual'});
write_excel_block_2023(excel_output_file, 'SUMMARY', 'A2', [temperature(:), selected_lambda_per_temp(:), nyquist_mean_residual_per_temp(:), nyquist_max_residual_per_temp(:)]);

write_excel_block_2023(excel_output_file, 'NGCV_IMAG', 'A1', [{'temperature_K'}, arrayfun(@(v) sprintf('lam_%.0e', v), lambda_values_master(:).', 'UniformOutput', false)]);
write_excel_block_2023(excel_output_file, 'NGCV_IMAG', 'A2', [temperature(:), ngcv_imag_matrix]);

write_excel_block_2023(excel_output_file, 'NGCV_REAL', 'A1', [{'temperature_K'}, arrayfun(@(v) sprintf('lam_%.0e', v), lambda_values_master(:).', 'UniformOutput', false)]);
write_excel_block_2023(excel_output_file, 'NGCV_REAL', 'A2', [temperature(:), ngcv_real_matrix]);

write_excel_block_2023(excel_output_file, 'RESID_LAMBDA', 'A1', [{'temperature_K'}, arrayfun(@(v) sprintf('lam_%.0e', v), lambda_values_master(:).', 'UniformOutput', false)]);
write_excel_block_2023(excel_output_file, 'RESID_LAMBDA', 'A2', [temperature(:), combined_residual_lambda_matrix]);

write_excel_block_2023(excel_output_file, 'NOISE_MEAN', 'A1', [{'temperature_K'}, arrayfun(@(v) sprintf('noise_%.1fpct', v), noise_injection_pct_master, 'UniformOutput', false)]);
write_excel_block_2023(excel_output_file, 'NOISE_MEAN', 'A2', [temperature(:), noise_mean_residual_pct_matrix]);

write_excel_block_2023(excel_output_file, 'KK_SUMMARY', 'A1', {'Temperature_K','Lambda','R_inf','KK_real_pct','KK_imag_pct','KK_total_pct','KK_max_pct','KK_status_code', ...
    'KK_DRT_real_pct','KK_DRT_imag_pct','KK_DRT_total_pct','KK_DRT_max_pct','KK_DRT_status_code'});
write_excel_block_2023(excel_output_file, 'KK_SUMMARY', 'A2', kk_summary);

% Save all deferred figures at the very end.
for k = 1:size(figure_export_jobs, 1)
    fig_h = figure_export_jobs{k, 1};
    fig_file = figure_export_jobs{k, 2};
    if isgraphics(fig_h)
        print(fig_h, fig_file, '-depsc');
    end
end

fprintf('\nSample-tagged outputs generated for %s:\n', sample_tag);
fprintf('  %s\n', excel_output_file);
fprintf('  %s\n', fig_dir);

function write_excel_block_2023(file_name, sheet_name, range_ref, payload)
try
    if iscell(payload)
        writecell(payload, file_name, 'Sheet', sheet_name, 'Range', range_ref, 'UseExcel', false);
    else
        writematrix(payload, file_name, 'Sheet', sheet_name, 'Range', range_ref, 'UseExcel', false);
    end
catch ME
    error('Failed to write workbook %s sheet %s (%s): %s', file_name, sheet_name, range_ref, ME.message);
end
end

function [resolved_file, switched] = prepare_excel_output_path_2023(target_file)
resolved_file = target_file;
switched = false;

if exist(target_file, 'file') == 2
    try
        delete(target_file);
    catch
        [out_dir, base_name, ext] = fileparts(target_file);
        if isempty(ext)
            ext = '.xlsx';
        end
        stamp = datestr(now, 'yyyymmdd_HHMMSS');
        resolved_file = fullfile(out_dir, [base_name '_autosave_' stamp ext]);
        switched = true;
    end
end
end

function [Z_rcn, Z_components] = calculate_eis_from_rcn_peaks_2023(freq_vec, R_inf, fitted_peaks)
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

function [peak_indices, peak_values] = find_peaks_simple_2023(y, min_height, min_distance)
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
    if isempty(selected) || all(abs(selected - idx) >= min_distance)
        selected(end+1,1) = idx; %#ok<AGROW>
    end
end

selected = sort(selected);
peak_indices = selected;
peak_values = y(selected);
end

function ysm = smooth_mavg_2023(y, w)
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

function ysm = smooth_sgolay_2023(y, frame_len, poly_order)
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
    ysm = smooth_mavg_2023(y, 3);
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

function model = create_rbf_model_2023(x, y)
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

function yq = rbf_eval_2023(model, xq)
xq = xq(:);
mx = model.x(:).';
DX = xq(:, ones(1,numel(mx))) - mx(ones(numel(xq),1), :);
Phi_q = exp(-(model.epsilon * DX).^2);
yq = Phi_q * model.w;
end

function [fitted_peaks, rbf_model, peak_refine] = refine_rbf_peak_fitting_2023(tau_vals, gamma_vals)
tau_vals = tau_vals(:);
gamma_vals = gamma_vals(:);
tau_log = log(tau_vals);

rbf_model = create_rbf_model_2023(tau_log, gamma_vals);
tau_dense = logspace(log10(min(tau_vals)), log10(max(tau_vals)), 1000).';
tau_dense_log = log(tau_dense);
gamma_rbf = rbf_eval_2023(rbf_model, tau_dense_log);
gamma_rbf(gamma_rbf < 0) = 0;

base_thr = max(0.05 * max(gamma_rbf), eps);
[base_peak_idx, ~] = find_peaks_simple_2023(gamma_rbf, base_thr, 20);

hidden_peak_idx = detect_hidden_peaks_second_derivative_2023(tau_dense_log, gamma_rbf, base_peak_idx);

all_peak_idx = unique([base_peak_idx(:); hidden_peak_idx(:)]);
all_peak_idx = sort(all_peak_idx);

hidden_flag = zeros(size(all_peak_idx));
for i = 1:numel(all_peak_idx)
    if any(abs(hidden_peak_idx - all_peak_idx(i)) <= 2)
        hidden_flag(i) = 1;
    end
end

[fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2023(tau_dense, gamma_rbf, all_peak_idx, hidden_flag);
peak_refine = struct('tau_dense', tau_dense, 'gamma_rbf', gamma_rbf, ...
    'base_peak_idx', base_peak_idx, 'hidden_peak_idx', hidden_peak_idx, ...
    'all_peak_idx', all_peak_idx, 'hidden_flag', hidden_flag, ...
    'fit_curve', fit_curve, 'fit_rmse', fit_rmse);
end

function hidden_peak_idx = detect_hidden_peaks_second_derivative_2023(tau_log, gamma_rbf, base_peak_idx)
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
curvature_signal = smooth_sgolay_2023(curvature_signal, 7, 2);

curv_thr = max(0.30 * max(curvature_signal), 3.2 * std(curvature_signal));
[curv_idx, ~] = find_peaks_simple_2023(curvature_signal, curv_thr, 28);

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

function [fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2023(tau_dense, gamma_rbf, peak_idx, hidden_flag)
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

[fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2023(xall, gamma_rbf, centers);

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
        [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2023(xall, gamma_rbf, centers);
        m = numel(centers);
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
    [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2023(xall, gamma_rbf, centers);
end

fit_rmse = sqrt(mean((fit_curve - gamma_rbf).^2));
fitted_peaks = repmat(struct('peak_id', 0, 'tau', 0, 'amplitude', 0, 'sigma', 0, ...
    'fwhm_logtau', 0, 'R', 0, 'C', 0, 'n', 0, 'is_hidden', 0, 'baseline', 0), m, 1);

for peak_idx_local = 1:m
    tau_peak = exp(centers(peak_idx_local));
    amp = amplitudes(peak_idx_local);
    sig = sigma_vec(peak_idx_local);
    fwhm_logtau = gaussian_fwhm_from_sigma_2023(sig);
    R_est = amp * sig * sqrt(2 * pi);
    C_est = tau_peak / max(R_est, eps);
    n_est = 1.144 / (1 + sig);

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

all_tau = [fitted_peaks.tau].';
[~, order] = sort(all_tau, 'descend');
fitted_peaks = fitted_peaks(order);
for peak_order = 1:numel(fitted_peaks)
    fitted_peaks(peak_order).peak_id = peak_order;
end
end

function [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2023(xall, gamma_rbf, centers)
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
obj = @(p) gaussian_rbf_sum_objective_2023(p, xall, gamma_rbf, centers);
p = fminsearch(obj, p0, opts);

[fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2023(p, xall, gamma_rbf, centers);
end

function sse = gaussian_rbf_sum_objective_2023(p, xall, yall, centers)
[fit_curve, amplitudes, sigma_vec] = evaluate_gaussian_rbf_sum_2023(p, xall, yall, centers);
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

function [fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2023(p, xall, yall, centers)
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

function fwhm_logtau = gaussian_fwhm_from_sigma_2023(sig)
fwhm_logtau = 2 * sqrt(2 * log(2)) * sig;
end
