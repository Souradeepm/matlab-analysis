function result = run_paper_vs_residual_peak_cv_compare(dataset_file)
% Compare paper-style lambda selection vs minimum-residual selection.
% Sensitivity metric: residual change after random 15% data removal.

clc;

if nargin < 1 || isempty(dataset_file)
    dataset_file = getenv('COMPARE_DATASET');
end
if isempty(dataset_file)
    dataset_file = 'S2422Al.xlsx';
end

result = struct();
repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, dataset_file);
if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[~, base_name, ~] = fileparts(dataset_file);

lambda_values = [1e-4; 1e-3; 1e-2; 1e-1];  % Range [1e-4, 1e-1] constrained
n_boot = 3;
max_temps = 10;
start_idx = 1;
drop_ratio = 0.15;
n_sens_repeat = 5;

tmp = str2double(getenv('COMPARE_NBOOT')); if ~isnan(tmp) && tmp >= 1, n_boot = round(tmp); end
tmp = str2double(getenv('COMPARE_MAX_TEMPS')); if ~isnan(tmp) && tmp >= 1, max_temps = round(tmp); end
tmp = str2double(getenv('COMPARE_START_IDX')); if ~isnan(tmp) && tmp >= 1, start_idx = round(tmp); end
tmp = str2double(getenv('COMPARE_DROP_RATIO')); if ~isnan(tmp) && tmp > 0 && tmp < 1, drop_ratio = tmp; end
tmp = str2double(getenv('COMPARE_NSENS')); if ~isnan(tmp) && tmp >= 1, n_sens_repeat = round(tmp); end
out_tag = strtrim(getenv('COMPARE_OUT_TAG'));

if isempty(out_tag)
    out_txt = fullfile(repo_root, [lower(base_name) '_paper_vs_residual_peak_cv_comparison.txt']);
else
    out_txt = fullfile(repo_root, [lower(base_name) '_paper_vs_residual_peak_cv_comparison_' out_tag '.txt']);
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
total_available = numel(temperature);
if start_idx > total_available
    fprintf('No temperatures left for %s at start index %d (available=%d)\n', dataset_file, start_idx, total_available);
    result.dataset = dataset_file;
    result.output_file = '';
    result.total_temperatures = 0;
    result.total_available_temperatures = total_available;
    return;
end

end_idx = min(total_available, start_idx + max_temps - 1);
temperature = temperature(start_idx:end_idx);
if size(data,2) >= 3 * end_idx
    data = data(:, (3*start_idx-2):(3*end_idx));
else
    data = data(:, (3*start_idx-2):size(data,2));
end
n_temp = numel(temperature);
if n_temp < 1
    error('No temperature columns parsed in %s', dataset_file);
end

res_lambda = zeros(n_temp,1);
paper_lambda = zeros(n_temp,1);
res_peak_count = zeros(n_temp,1);
paper_peak_count = zeros(n_temp,1);
res_sens = zeros(n_temp,1);
paper_sens = zeros(n_temp,1);

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

    n_l = numel(lambda_values);
    mean_abs_resid = zeros(n_l, 1);
    rid_metric = zeros(n_l, 1);
    sens_metric = zeros(n_l, 1);
    var_metric = zeros(n_l, 1);
    gamma_all = zeros(numel(freq_vec), n_l);

    for i = 1:n_l
        lam = lambda_values(i);
        [gamma_i, R_inf_i] = tr_drt_local(freq_vec, Z_exp, lam);
        Z_fit_i = calculate_eis_local(freq_vec, gamma_i, R_inf_i);
        mean_abs_resid(i) = mean(abs(Z_exp - Z_fit_i));
        gamma_all(:, i) = gamma_i;

        [g_re, ~] = tr_drt_component_local(freq_vec, Z_exp, lam, 'real');
        [g_im, ~] = tr_drt_component_local(freq_vec, Z_exp, lam, 'imag');
        rid_metric(i) = mean((g_re - g_im).^2);

        sens_metric(i) = estimate_residual_sensitivity_local(freq_vec, Z_exp, lam, ...
            mean_abs_resid(i), drop_ratio, n_sens_repeat);
        var_metric(i) = estimate_resample_variance_local(freq_vec, Z_exp, lam, n_boot);
    end

    [~, idx_res] = min(mean_abs_resid);
    idx_paper = select_lambda_idx_local(rid_metric, sens_metric, var_metric);

    res_lambda(t) = lambda_values(idx_res);
    paper_lambda(t) = lambda_values(idx_paper);
    res_sens(t) = sens_metric(idx_res);
    paper_sens(t) = sens_metric(idx_paper);

    thr_res = 0.05 * max(gamma_all(:, idx_res)); if ~(thr_res > 0), thr_res = 0; end
    thr_paper = 0.05 * max(gamma_all(:, idx_paper)); if ~(thr_paper > 0), thr_paper = 0; end
    min_dist = max(3, round(numel(freq_vec) / 50));

    [pk_res, ~] = find_peaks_simple_local(gamma_all(:, idx_res), thr_res, min_dist);
    [pk_paper, ~] = find_peaks_simple_local(gamma_all(:, idx_paper), thr_paper, min_dist);
    res_peak_count(t) = numel(pk_res);
    paper_peak_count(t) = numel(pk_paper);
end

delta_peak_signed = paper_peak_count - res_peak_count;
delta_sens_signed = paper_sens - res_sens;
delta_peak = abs(delta_peak_signed);
delta_sens = abs(delta_sens_signed);

fid = fopen(out_txt, 'w');
if fid < 0
    error('Failed to open output file: %s', out_txt);
end

fprintf(fid, 'Paper vs minimum-residual DRT comparison\n');
fprintf(fid, 'Dataset: %s\n', dataset_file);
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Batch start index: %d\n', start_idx);
fprintf(fid, 'Batch end index  : %d\n', end_idx);
fprintf(fid, 'Total available temperatures: %d\n', total_available);
fprintf(fid, 'Total temperatures: %d\n\n', n_temp);
fprintf(fid, 'Bootstrap repeats for variance metric: %d\n\n', n_boot);
fprintf(fid, 'Drop ratio for sensitivity metric: %.2f\n', drop_ratio);
fprintf(fid, 'Sensitivity repeats per lambda: %d\n\n', n_sens_repeat);

fprintf(fid, 'Aggregate comparison\n');
fprintf(fid, 'Mean peaks (residual method): %.4f\n', mean(res_peak_count));
fprintf(fid, 'Mean peaks (paper method)   : %.4f\n', mean(paper_peak_count));
fprintf(fid, 'Mean peak delta abs |paper-res| : %.4f\n', mean(delta_peak));
fprintf(fid, 'Median peak delta abs |paper-res|: %.4f\n', median(delta_peak));
fprintf(fid, 'Temps paper has more peaks      : %d\n', sum(delta_peak_signed > 0));
fprintf(fid, 'Temps equal peaks               : %d\n', sum(delta_peak_signed == 0));
fprintf(fid, 'Temps paper has fewer peaks     : %d\n\n', sum(delta_peak_signed < 0));

fprintf(fid, 'Mean residual-change sensitivity %% (residual-selected lambda): %.8e\n', mean(res_sens));
fprintf(fid, 'Mean residual-change sensitivity %% (paper-selected lambda)   : %.8e\n', mean(paper_sens));
fprintf(fid, 'Mean sensitivity delta abs %% |paper-res|                    : %.8e\n', mean(delta_sens));
fprintf(fid, 'Median sensitivity delta abs %% |paper-res|                  : %.8e\n', median(delta_sens));
fprintf(fid, 'Temps paper lower sensitivity                             : %d\n', sum(delta_sens_signed < 0));
fprintf(fid, 'Temps equal sensitivity                                   : %d\n', sum(delta_sens_signed == 0));
fprintf(fid, 'Temps paper higher sensitivity                            : %d\n\n', sum(delta_sens_signed > 0));

fprintf(fid, 'Per-temperature detail\n');
fprintf(fid, 'TemperatureK,ResidualLambda,PaperLambda,ResidualPeaks,PaperPeaks,PeakDeltaAbs_PaperMinusRes,ResidualSensitivity,PaperSensitivity,SensitivityDeltaAbs_PaperMinusRes\n');
for t = 1:n_temp
    fprintf(fid, '%.6f,%.8e,%.8e,%d,%d,%d,%.8e,%.8e,%.8e\n', ...
        temperature(t), res_lambda(t), paper_lambda(t), ...
        res_peak_count(t), paper_peak_count(t), delta_peak(t), ...
        res_sens(t), paper_sens(t), delta_sens(t));
end
fclose(fid);

fprintf('Wrote comparison report: %s\n', out_txt);

result.dataset = dataset_file;
result.output_file = out_txt;
result.total_temperatures = n_temp;
result.batch_start_index = start_idx;
result.batch_end_index = end_idx;
result.total_available_temperatures = total_available;
result.mean_residual_peaks = mean(res_peak_count);
result.mean_paper_peaks = mean(paper_peak_count);
result.mean_peak_delta = mean(delta_peak);
result.mean_residual_sensitivity = mean(res_sens);
result.mean_paper_sensitivity = mean(paper_sens);
result.mean_sensitivity_delta = mean(delta_sens);
result.paper_lower_sensitivity_count = sum(delta_sens_signed < 0);
result.paper_higher_sensitivity_count = sum(delta_sens_signed > 0);
end

function idx_best = select_lambda_idx_local(rid_vec, sens_vec, var_vec)
% Paper Method v2: Score-based selection (FIXED Issue #3)
% Combines three metrics equally (normalized):
% - RID: Re/Im divergence
% - Sensitivity: Perturbation robustness
% - Var: Bootstrap resampling variance
%
% Lambda selected: argmin(norm(RID) + norm(Sensitivity) + norm(Var))

norm_rid = normalize_metric_local(rid_vec);
norm_sens = normalize_metric_local(sens_vec);
norm_var = normalize_metric_local(var_vec);

score = norm_rid + norm_sens + norm_var;
[~, idx_best] = min(score);
end

function nvec = normalize_metric_local(vec)
vec = vec(:); vmin = min(vec); vmax = max(vec);
if vmax > vmin
    nvec = (vec - vmin) ./ (vmax - vmin);
else
    nvec = zeros(size(vec));
end
end

function [gamma, R_inf] = tr_drt_local(freq_vec, Z_exp, el)
freq_vec = freq_vec(:); Z_exp = Z_exp(:); n = numel(freq_vec);
A_re = calc_A_re_local(freq_vec); A_im = calc_A_im_local(freq_vec);
M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
b = [real(Z_exp); imag(Z_exp); zeros(n,1)];
opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n); R_inf = x(n+1);
end

function [gamma, R_inf] = tr_drt_component_local(freq_vec, Z_exp, el, component)
freq_vec = freq_vec(:); Z_exp = Z_exp(:); n = numel(freq_vec);
A_re = calc_A_re_local(freq_vec); A_im = calc_A_im_local(freq_vec);
switch lower(component)
    case 'real'
        M = [A_re, ones(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
        b = [real(Z_exp); zeros(n,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(el/2) * eye(n), zeros(n,1)];
        b = [imag(Z_exp); zeros(n,1)];
    otherwise
        error('Unknown component mode: %s', component);
end
opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n); R_inf = x(n+1);
end

function Z_cal = calculate_eis_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec); A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, el, n_boot)
n = numel(freq_vec); if nargin < 4 || isempty(n_boot), n_boot = 8; end
sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);
sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);
gamma_samples = zeros(n, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n,1) + 1i * sigma_im * randn(n,1);
    Z_boot = Z_exp(:) + noise;
    [gk, ~] = tr_drt_local(freq_vec, Z_boot, el);
    gamma_samples(:, k) = gk(:);
end
mean_var = mean(var(gamma_samples, 0, 2));
end

function sens_val = estimate_residual_sensitivity_local(freq_vec, Z_exp, el, base_residual, drop_ratio, n_repeat)
n = numel(freq_vec);
if n < 10
    sens_val = 0;
    return;
end
n_drop = max(1, round(drop_ratio * n));
delta = zeros(n_repeat,1);
for k = 1:n_repeat
    keep_mask = true(n,1);
    drop_idx = randperm(n, n_drop);
    keep_mask(drop_idx) = false;
    freq_sub = freq_vec(keep_mask);
    Z_sub = Z_exp(keep_mask);
    [gamma_sub, R_sub] = tr_drt_local(freq_sub, Z_sub, el);
    Z_fit_sub = calculate_eis_local(freq_sub, gamma_sub, R_sub);
    resid_sub = mean(abs(Z_sub - Z_fit_sub));
    base_safe = max(base_residual, eps);
    delta(k) = 100 * abs(resid_sub - base_residual) / base_safe;
end
sens_val = mean(delta);
end

function A_re = calc_A_re_local(freq)
omega = 2 * pi * freq(:); tau = 1 ./ freq(:); n = numel(freq); A_re = zeros(n, n);
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
omega = 2 * pi * freq(:); tau = 1 ./ freq(:); n = numel(freq); A_im = zeros(n, n);
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

function [peak_indices, peak_values] = find_peaks_simple_local(y, min_height, min_distance)
y = y(:); n = numel(y);
if n < 3
    peak_indices = []; peak_values = []; return;
end
if nargin < 2 || isempty(min_height), min_height = -inf; end
if nargin < 3 || isempty(min_distance), min_distance = 1; end
candidates = [];
for i = 2:(n-1)
    if y(i) >= y(i-1) && y(i) >= y(i+1) && y(i) >= min_height
        candidates(end+1,1) = i; %#ok<AGROW>
    end
end
if isempty(candidates)
    peak_indices = []; peak_values = []; return;
end
[~, ord] = sort(y(candidates), 'descend');
selected = [];
for k = 1:numel(ord)
    idx = candidates(ord(k));
    if isempty(selected)
        selected = idx;
    else
        if min(abs(selected - idx)) >= min_distance
            selected(end+1,1) = idx; %#ok<AGROW>
        end
    end
end
selected = sort(selected);
peak_indices = selected;
peak_values = y(selected);
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
n_cols = numel(row_cells); max_start = n_cols - 2;
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
    if iscell(text), header = text(1, :); else, header = cellstr(char(text(1, :))); end
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
temperature = [];
if isempty(raw), return; end
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
            best_count = count; best_row = r; best_offset = off;
        end
    end
end
if best_row == 0, return; end
temperature = nan(n_sets,1);
for i = 1:n_sets
    cols3 = (3*i-2):(3*i); found = NaN;
    c0 = 3*(i-1) + (1 + best_offset);
    if c0 >= 1 && c0 <= size(raw,2)
        v0 = raw{best_row,c0};
        if ischar(v0)
            tok0 = strtrim(v0);
            tok0(tok0 == 'K' | tok0 == 'k') = [];
            val0 = str2double(tok0);
            if ~isnan(val0), found = val0; end
        end
    end
    for c = cols3
        if ~isnan(found), break; end
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
