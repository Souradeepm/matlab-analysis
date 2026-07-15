function result = run_bayes_drt_workflow_matlab2011(dataset_file)
% run_bayes_drt_workflow_matlab2011
% Bayes-DRT-style MATLAB 2011 workflow.
%
% This separate workflow follows the Bayes-DRT ridge Re/Im cross-validation
% idea: choose lambda by minimizing the sum of real-part and imaginary-part
% cross-validation errors over a lambda grid, then fit the full spectrum with
% the selected lambda.
%
% Default scope: first 10 temperatures of a dataset workbook.

clc;

if nargin < 1 || isempty(dataset_file)
    dataset_file = getenv('BAYES_DRT_DATASET');
end
if isempty(dataset_file)
    dataset_file = 'S2422Al.xlsx';
end

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, dataset_file);
if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

max_temps = 10;
tmp = str2double(getenv('BAYES_DRT_MAX_TEMPS'));
if ~isnan(tmp) && tmp >= 1
    max_temps = round(tmp);
end

start_idx = 1;
tmp = str2double(getenv('BAYES_DRT_START_IDX'));
if ~isnan(tmp) && tmp >= 1
    start_idx = round(tmp);
end

out_tag = strtrim(getenv('BAYES_DRT_OUT_TAG'));

lambda_min_exp = -4;   % 1e-4 minimum (prevent underfitting)
lambda_max_exp = -1;   % 1e-1 maximum (prevent overfitting)
lambda_count = 31;     % Same number of points for resolution
tmp = str2double(getenv('BAYES_DRT_LAMBDA_MIN_EXP'));
if ~isnan(tmp)
    lambda_min_exp = tmp;
end
tmp = str2double(getenv('BAYES_DRT_LAMBDA_MAX_EXP'));
if ~isnan(tmp)
    lambda_max_exp = tmp;
end
tmp = str2double(getenv('BAYES_DRT_LAMBDA_COUNT'));
if ~isnan(tmp) && tmp >= 3
    lambda_count = round(tmp);
end
lambda_values = logspace(lambda_min_exp, lambda_max_exp, lambda_count).';

[~, base_name, ~] = fileparts(dataset_file);
if isempty(out_tag)
    out_txt = fullfile(repo_root, [lower(base_name) '_bayes_drt_matlab2011_10temp.txt']);
else
    out_txt = fullfile(repo_root, [lower(base_name) '_bayes_drt_matlab2011_' out_tag '.txt']);
end

[data, text, raw] = xlsread(input_path); %#ok<ASGLU>
n_sets = detect_n_sets_local(text, raw, size(data,2));
if n_sets < 1
    n_sets = floor(size(data,2) / 3);
end
max_cols = min(size(data,2), 3 * n_sets);
data = data(:, 1:max_cols);

temperature = parse_temperature_labels_local(text, raw, n_sets);
if isempty(temperature)
    temperature = (1:n_sets).';
end
temperature = temperature(1:min(numel(temperature), n_sets));
if isempty(temperature)
    error('No temperatures parsed in %s', dataset_file);
end

total_available = numel(temperature);
if start_idx > total_available
    fprintf('No temperatures left for %s at start index %d (available=%d)\n', dataset_file, start_idx, total_available);
    result = struct();
    result.dataset = dataset_file;
    result.output_file = '';
    result.total_temperatures = 0;
    result.total_available_temperatures = total_available;
    return;
end
end_idx = min(total_available, start_idx + max_temps - 1);
if size(data,2) >= 3 * end_idx
    data = data(:, (3*start_idx-2):(3*end_idx));
else
    data = data(:, (3*start_idx-2):size(data,2));
end
n_temp = numel(temperature(start_idx:end_idx));
selected_temperature = temperature(start_idx:end_idx);

selected_lambda = zeros(n_temp,1);
re_cv = zeros(n_temp,1);
im_cv = zeros(n_temp,1);
tot_cv = zeros(n_temp,1);
mean_resid = zeros(n_temp,1);
peak_count = zeros(n_temp,1);
residual_peak_ratio = zeros(n_temp,1);

for t = 1:n_temp
    cols = (3*t-2):(3*t);
    block = data(:, cols);
    freq_vec = block(:,1);
    mag = block(:,2);
    phase_deg = block(:,3);

    valid = isfinite(freq_vec) & isfinite(mag) & isfinite(phase_deg) & (freq_vec > 0);
    freq_vec = freq_vec(valid);
    mag = mag(valid);
    phase_deg = phase_deg(valid);

    Z_exp = mag .* exp(1i * phase_deg * pi / 180);
    [freq_vec, ord] = sort(freq_vec);
    Z_exp = Z_exp(ord);

    [idx_best, lam_best, cv_real_best, cv_imag_best, cv_tot_best] = select_lambda_reimcv_local(freq_vec, Z_exp, lambda_values);

    [gamma_full, R_inf_full] = tr_drt_full_local(freq_vec, Z_exp, lam_best);
    Z_fit_full = calculate_eis_local(freq_vec, gamma_full, R_inf_full);

    res_full = abs(Z_exp - Z_fit_full);
    mean_resid(t) = mean(res_full);
    selected_lambda(t) = lam_best;
    re_cv(t) = cv_real_best;
    im_cv(t) = cv_imag_best;
    tot_cv(t) = cv_tot_best;

    thr = 0.05 * max(gamma_full);
    if ~(thr > 0)
        thr = 0;
    end
    min_dist = max(3, round(numel(freq_vec) / 50));
    [pk_idx, ~] = find_peaks_simple_local(gamma_full, thr, min_dist);
    peak_count(t) = numel(pk_idx);
    residual_peak_ratio(t) = mean_resid(t) / max(peak_count(t), 1);

    fprintf('Dataset %s temp %d/%d: lambda=%.3e, total CV=%.6e, mean resid=%.6e, peaks=%d\n', ...
        dataset_file, t, n_temp, lam_best, cv_tot_best, mean_resid(t), peak_count(t));
    fprintf('  lambda index %d of %d\n', idx_best, numel(lambda_values));
end

fid = fopen(out_txt, 'w');
if fid < 0
    error('Failed to open output file: %s', out_txt);
end

fprintf(fid, 'Bayes-DRT MATLAB 2011 workflow (Re/Im CV lambda selection)\n');
fprintf(fid, 'Dataset: %s\n', dataset_file);
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Lambda grid: %.8e to %.8e (%d values)\n', lambda_values(1), lambda_values(end), numel(lambda_values));
fprintf(fid, 'Batch start index: %d\n', start_idx);
fprintf(fid, 'Batch end index  : %d\n', end_idx);
fprintf(fid, 'Total available temperatures: %d\n', total_available);
fprintf(fid, 'Total processed temperatures: %d\n\n', n_temp);

fprintf(fid, 'Aggregate summary\n');
fprintf(fid, 'Mean selected lambda: %.8e\n', mean(selected_lambda));
fprintf(fid, 'Median selected lambda: %.8e\n', median(selected_lambda));
fprintf(fid, 'Mean real-part CV: %.8e\n', mean(re_cv));
fprintf(fid, 'Mean imag-part CV: %.8e\n', mean(im_cv));
fprintf(fid, 'Mean total CV: %.8e\n', mean(tot_cv));
fprintf(fid, 'Mean absolute residual: %.8e\n', mean(mean_resid));
fprintf(fid, 'Mean peak count: %.4f\n', mean(peak_count));
fprintf(fid, 'Max peak count: %d\n\n', max(peak_count));

fprintf(fid, 'Per-temperature detail\n');
fprintf(fid, 'TemperatureK,SelectedLambda,RealCV,ImagCV,TotalCV,MeanAbsResidual,PeakCount,ResidualPerPeak\n');
for t = 1:n_temp
    fprintf(fid, '%.6f,%.8e,%.8e,%.8e,%.8e,%.8e,%d,%.8e\n', ...
        selected_temperature(t), selected_lambda(t), re_cv(t), im_cv(t), tot_cv(t), mean_resid(t), peak_count(t), residual_peak_ratio(t));
end
fclose(fid);

fprintf('Wrote Bayes-DRT workflow report: %s\n', out_txt);

result = struct();
result.dataset = dataset_file;
result.output_file = out_txt;
result.total_available_temperatures = total_available;
result.batch_start_index = start_idx;
result.batch_end_index = end_idx;
result.total_temperatures = n_temp;
result.lambda_grid = lambda_values;
result.selected_lambda = selected_lambda;
result.real_cv = re_cv;
result.imag_cv = im_cv;
result.total_cv = tot_cv;
result.mean_residual = mean_resid;
result.peak_count = peak_count;
end

function [idx_best, lam_best, cv_real_best, cv_imag_best, cv_total_best] = select_lambda_reimcv_local(freq_vec, Z_exp, lambda_values)
num_l = numel(lambda_values);
cv_real = zeros(num_l,1);
cv_imag = zeros(num_l,1);
cv_total = zeros(num_l,1);

for i = 1:num_l
    lam = lambda_values(i);
    [gamma_re, R_re] = tr_drt_component_local(freq_vec, Z_exp, lam, 'real');
    [gamma_im, ~] = tr_drt_component_local(freq_vec, Z_exp, lam, 'imag');

    Z_from_re = calculate_eis_local(freq_vec, gamma_re, R_re);
    R_from_im = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, gamma_im);
    Z_from_im = calculate_eis_local(freq_vec, gamma_im, R_from_im);

    cv_imag(i) = mean((imag(Z_exp) - imag(Z_from_re)).^2);
    cv_real(i) = mean((real(Z_exp) - real(Z_from_im)).^2);
    cv_total(i) = cv_real(i) + cv_imag(i);
end

[cv_total_best, idx_best] = min(cv_total);
lam_best = lambda_values(idx_best);
cv_real_best = cv_real(idx_best);
cv_imag_best = cv_imag(idx_best);
end

function [gamma, R_inf] = tr_drt_full_local(freq_vec, Z_exp, el)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);
M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];
[gamma, R_inf] = solve_nonnegative_ls_local(M, b, n + 1, Z_exp);
end

function [gamma, R_inf] = tr_drt_component_local(freq_vec, Z_exp, el, component)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
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
[gamma, R_inf] = solve_nonnegative_ls_local(M, b, n + 1, Z_exp);
end

function [gamma, R_inf] = solve_nonnegative_ls_local(M, b, n_gamma, Z_exp)
ub = max(abs(Z_exp)) * ones(n_gamma, 1);
if exist('lsqlin', 'file') == 2
    lb = zeros(n_gamma, 1);
    options = optimset('Display', 'off');
    x = lsqlin(M, b, [], [], [], [], lb, ub, [], options);
else
    x = lsqnonneg(M, b);
    x = min(x, ub);
end
gamma = x(1:end-1);
R_inf = x(end);
if ~(R_inf > 0)
    R_inf = 0;
end
end

function R_inf = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, gamma)
A_re = calc_A_re_local(freq_vec);
re_model_wo_r = A_re * gamma(:);
R_inf = mean(real(Z_exp(:)) - re_model_wo_r);
if ~(R_inf > 0)
    R_inf = 0;
end
end

function Z_cal = calculate_eis_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma(:) + 1i * (A_im * gamma(:));
end

function [peak_indices, peak_values] = find_peaks_simple_local(y, min_height, min_distance)
y = y(:);
n = numel(y);
peak_indices = [];
peak_values = [];
if n < 3
    return;
end
if nargin < 2 || isempty(min_height)
    min_height = -inf;
end
if nargin < 3 || isempty(min_distance)
    min_distance = 1;
end
cand = [];
for i = 2:(n-1)
    if y(i) >= y(i-1) && y(i) > y(i+1) && y(i) >= min_height
        cand(end+1,1) = i; %#ok<AGROW>
    end
end
if isempty(cand)
    return;
end
[~, ord] = sort(y(cand), 'descend');
selected = [];
for k = 1:numel(ord)
    idx = cand(ord(k));
    if isempty(selected)
        selected(end+1,1) = idx; %#ok<AGROW>
    elseif min(abs(selected - idx)) >= min_distance
        selected(end+1,1) = idx; %#ok<AGROW>
    end
end
selected = sort(selected);
peak_indices = selected;
peak_values = y(selected);
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

function n_sets = detect_n_sets_local(text, raw, n_data_cols)
n_sets = 0;
if ~isempty(text) && iscell(text)
    n_sets = count_header_triplets_local(text(1,:));
end
if n_sets < 1 && ~isempty(raw)
    n_sets = count_header_triplets_local(raw(1,:));
end
if n_sets < 1
    n_sets = floor(n_data_cols / 3);
end
end

function n_triplets = count_header_triplets_local(row_cells)
n_triplets = 0;
if isempty(row_cells), return; end
n_cols = numel(row_cells);
max_start = n_cols - 2;
for c = 1:3:max_start
    a = row_cells{c}; b = row_cells{c+1}; d = row_cells{c+2};
    if ischar(a) && ischar(b) && ischar(d)
        if strcmpi(strtrim(a), 'freq') && strcmpi(strtrim(b), 'z') && strcmpi(strtrim(d), 'theta')
            n_triplets = n_triplets + 1;
        end
    end
end
end

function temperature = parse_temperature_labels_local(text, raw, n_sets)
temperature = [];
if ~isempty(text)
    if iscell(text)
        header = text(1, :);
    else
        header = cellstr(char(text(1, :)));
    end
    for k = 1:numel(header)
        tok = header{k};
        if ischar(tok)
            tok(tok == 'K' | tok == 'k') = [];
            val = str2double(strtrim(tok));
            if ~isnan(val)
                temperature(end+1,1) = val; %#ok<AGROW>
            end
        end
    end
end
if numel(temperature) >= n_sets
    temperature = temperature(1:n_sets);
    return;
end
if isempty(raw)
    temperature = [];
    return;
end
scan_rows = min(size(raw,1), 12);
best_row = 0; best_count = 0; best_offset = 0;
for r = 1:scan_rows
    for off = 0:2
        count = 0;
        for c = (1+off):3:(3*n_sets)
            v = raw{r,c};
            if ischar(v) && (~isempty(strfind(v, 'K')) || ~isempty(strfind(v, 'k')))
                count = count + 1;
            end
        end
        if count > best_count
            best_count = count;
            best_row = r;
            best_offset = off;
        end
    end
end
if best_row == 0
    temperature = [];
    return;
end
temperature = nan(n_sets,1);
for i = 1:n_sets
    cols3 = (3*i-2):(3*i);
    found = NaN;
    c0 = 3*(i-1) + (1 + best_offset);
    if c0 >= 1 && c0 <= size(raw,2)
        v0 = raw{best_row,c0};
        if ischar(v0)
            tok0 = strtrim(v0);
            tok0(tok0 == 'K' | tok0 == 'k') = [];
            val0 = str2double(tok0);
            if ~isnan(val0)
                found = val0;
            end
        end
    end
    for c = cols3
        if ~isnan(found)
            break;
        end
        v = raw{best_row,c};
        if ischar(v)
            tok = strtrim(v);
            tok(tok == 'K' | tok == 'k') = [];
            val = str2double(tok);
            if ~isnan(val)
                found = val;
                break;
            end
        end
    end
    temperature(i) = found;
end
temperature = temperature(~isnan(temperature));
temperature = temperature(1:min(numel(temperature), n_sets));
end
