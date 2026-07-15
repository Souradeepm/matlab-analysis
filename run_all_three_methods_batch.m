function run_all_three_methods_batch(dataset_file, start_temp_idx, end_temp_idx)
% run_all_three_methods_batch: Process a batch of temperatures using all three lambda selection methods
%
% Syntax:
%   run_all_three_methods_batch(dataset_file, start_temp_idx, end_temp_idx)
%
% Inputs:
%   dataset_file   = Excel file path (e.g., 'S2022Sap.xlsx')
%   start_temp_idx = Starting temperature index (1-based)
%   end_temp_idx   = Ending temperature index (inclusive)
%
% Outputs:
%   Creates batch-specific output files:
%   - *_bayes_drt_matlab2011_b{batch}.txt
%   - *_gcv_method_b{batch}.txt
%   - *_residual_method_b{batch}.txt
%
% Environment Variables:
%   MATLAB_ANALYSIS_REPO_ROOT  = Repo base directory
%   BAYES_DRT_LAMBDA_MIN_EXP   = Min lambda exponent (default -4)
%   BAYES_DRT_LAMBDA_MAX_EXP   = Max lambda exponent (default -1)
%   BAYES_DRT_LAMBDA_COUNT     = Lambda grid points (default 31)

repo_root = getenv('MATLAB_ANALYSIS_REPO_ROOT');
if isempty(repo_root)
    repo_root = pwd;
end

% Ensure dataset file exists
if ~isfile(dataset_file) && isfile(fullfile(repo_root, dataset_file))
    dataset_file = fullfile(repo_root, dataset_file);
end
if ~isfile(dataset_file)
    error('Dataset file not found: %s', dataset_file);
end

fprintf('[All Methods] Processing %s, temps %d to %d\n', dataset_file, start_temp_idx, end_temp_idx);

% ---- BAYES-DRT METHOD ----
fprintf('[Bayes-DRT] Running batch...\n');
try
    result_bayes = run_bayes_drt_batch_local(dataset_file, start_temp_idx, end_temp_idx);
    fprintf('[Bayes-DRT] ✓ Batch complete. %d temps processed.\n', result_bayes.n_temps);
catch ME
    fprintf('[Bayes-DRT] ✗ Error: %s\n', ME.message);
    result_bayes.n_temps = 0;
    result_bayes.success = false;
end

% ---- GCV METHOD ----
fprintf('[GCV Method] Running batch...\n');
try
    result_gcv = run_gcv_method_batch_local(dataset_file, start_temp_idx, end_temp_idx);
    fprintf('[GCV Method] ✓ Batch complete. %d temps processed.\n', result_gcv.n_temps);
catch ME
    fprintf('[GCV Method] ✗ Error: %s\n', ME.message);
    result_gcv.n_temps = 0;
    result_gcv.success = false;
end

% ---- RESIDUAL METHOD ----
fprintf('[Residual Method] Running batch...\n');
try
    result_residual = run_residual_method_batch_local(dataset_file, start_temp_idx, end_temp_idx);
    fprintf('[Residual Method] ✓ Batch complete. %d temps processed.\n', result_residual.n_temps);
catch ME
    fprintf('[Residual Method] ✗ Error: %s\n', ME.message);
    result_residual.n_temps = 0;
    result_residual.success = false;
end

fprintf('[All Methods] Batch complete. Total output files: 3\n');
% NOTE: no exit here - caller (run_all_datasets_comprehensive.m) controls termination
end

%% ---- BAYES-DRT BATCH HANDLER ----
function result = run_bayes_drt_batch_local(dataset_file, start_idx, end_idx)
[~, base_name, ~] = fileparts(dataset_file);
base_name = lower(base_name);
repo_root = getenv('MATLAB_ANALYSIS_REPO_ROOT');
if isempty(repo_root), repo_root = pwd; end

lambda_min_exp = -4;
lambda_max_exp = -1;
lambda_count = 31;

tmp = str2double(getenv('BAYES_DRT_LAMBDA_MIN_EXP'));
if ~isnan(tmp), lambda_min_exp = tmp; end
tmp = str2double(getenv('BAYES_DRT_LAMBDA_MAX_EXP'));
if ~isnan(tmp), lambda_max_exp = tmp; end
tmp = str2double(getenv('BAYES_DRT_LAMBDA_COUNT'));
if ~isnan(tmp) && tmp >= 3, lambda_count = round(tmp); end

lambda_values = logspace(lambda_min_exp, lambda_max_exp, lambda_count).';

% Read dataset using xlsread
try
    [data, ~, raw] = xlsread(dataset_file); %#ok<XLSRD>
    [data, temperature] = load_dataset_local(data, raw);
    if isempty(temperature)
        result.n_temps = 0; result.success = false; return;
    end
catch ME
    fprintf('Bayes read error: %s\n', ME.message);
    result.n_temps = 0; result.success = false; return;
end

end_idx = min(end_idx, numel(temperature));
if start_idx > end_idx
    result.n_temps = 0;
    result.success = true;
    return;
end

% Process batch
n_points = 100;
tau_max = 1e0;
tau_min = 1e-9;
tau_vec = logspace(log10(tau_min), log10(tau_max), n_points).';

output_file = fullfile(repo_root, sprintf('%s_bayes_drt_matlab2011_b%d.txt', base_name, ceil(start_idx/10)));
fid = fopen(output_file, 'w');
if fid < 0
    result.n_temps = 0;
    result.success = false;
    return;
end

fprintf(fid, 'Temperature,SelectedLambda,RealCV,ImagCV,TotalCV,MeanAbsResidual,PeakCount\n');

count = 0;
for t = start_idx:end_idx
    try
        block = data(:, (t-1)*3 + 1 : (t-1)*3 + 3);
        valid = all(isfinite(block), 2) & block(:,1) > 0;
        block = block(valid, :);
        if size(block, 1) < 5, continue; end
        freq_vec = block(:, 1);
        Z_mag    = block(:, 2);
        Z_phase  = block(:, 3);
        Z_exp = Z_mag .* exp(1i * Z_phase * pi / 180);
        
        [lam_best, cv_real, cv_imag, cv_total, mean_residual, n_peaks] = ...
            select_lambda_reimcv_local(freq_vec, Z_exp, lambda_values);
        
        fprintf(fid, '%.6f,%.8e,%.6e,%.6e,%.6e,%.6e,%d\n', ...
            temperature(t), lam_best, cv_real, cv_imag, cv_total, mean_residual, n_peaks);
        
        count = count + 1;
    catch ME
        fprintf(fid, '%% t=%d ERROR: %s\n', t, ME.message);
        continue;
    end
end

fclose(fid);
result.n_temps = count;
result.success = true;
result.output_file = output_file;
end

%% ---- GCV METHOD BATCH HANDLER ----
function result = run_gcv_method_batch_local(dataset_file, start_idx, end_idx)
[~, base_name, ~] = fileparts(dataset_file);
base_name = lower(base_name);
repo_root = getenv('MATLAB_ANALYSIS_REPO_ROOT');
if isempty(repo_root), repo_root = pwd; end

% GCV lambda grid: imaginary inversion, Re+Im cross-validation
lambda_values = logspace(-4, 0, 10).';

% Read dataset using xlsread
try
    [data, ~, raw] = xlsread(dataset_file); %#ok<XLSRD>
    [data, temperature] = load_dataset_local(data, raw);
    if isempty(temperature)
        result.n_temps = 0; result.success = false; return;
    end
catch ME
    fprintf('GCV read error: %s\n', ME.message);
    result.n_temps = 0; result.success = false; return;
end

end_idx = min(end_idx, numel(temperature));
if start_idx > end_idx
    result.n_temps = 0;
    result.success = true;
    return;
end

output_file = fullfile(repo_root, sprintf('%s_gcv_method_b%d.txt', base_name, ceil(start_idx/10)));
fid = fopen(output_file, 'w');
if fid < 0
    result.n_temps = 0;
    result.success = false;
    return;
end

fprintf(fid, 'Temperature,SelectedLambda,GCVScore,MeanResidual,PeakCount\n');

count = 0;
for t = start_idx:end_idx
    try
        block = data(:, (t-1)*3 + 1 : (t-1)*3 + 3);
        valid = all(isfinite(block), 2) & block(:,1) > 0;
        block = block(valid, :);
        if size(block, 1) < 5, continue; end
        freq_vec = block(:, 1);
        Z_mag    = block(:, 2);
        Z_phase  = block(:, 3);
        Z_exp = Z_mag .* exp(1i * Z_phase * pi / 180);

        [lam_best, gcv_best, mean_resid, n_peaks] = ...
            select_lambda_gcv_local(freq_vec, Z_exp, lambda_values);

        fprintf(fid, '%.6f,%.8e,%.6e,%.6e,%d\n', ...
            temperature(t), lam_best, gcv_best, mean_resid, n_peaks);

        count = count + 1;
    catch
        continue;
    end
end

fclose(fid);
result.n_temps = count;
result.success = true;
result.output_file = output_file;
end

%% ---- RESIDUAL METHOD BATCH HANDLER ----
function result = run_residual_method_batch_local(dataset_file, start_idx, end_idx)
[~, base_name, ~] = fileparts(dataset_file);
base_name = lower(base_name);
repo_root = getenv('MATLAB_ANALYSIS_REPO_ROOT');
if isempty(repo_root), repo_root = pwd; end

lambda_values = logspace(-4, -1, 11).';

% Read dataset using xlsread
try
    [data, ~, raw] = xlsread(dataset_file); %#ok<XLSRD>
    [data, temperature] = load_dataset_local(data, raw);
    if isempty(temperature)
        result.n_temps = 0; result.success = false; return;
    end
catch ME
    fprintf('Residual read error: %s\n', ME.message);
    result.n_temps = 0; result.success = false; return;
end

end_idx = min(end_idx, numel(temperature));
if start_idx > end_idx
    result.n_temps = 0;
    result.success = true;
    return;
end

output_file = fullfile(repo_root, sprintf('%s_residual_method_b%d.txt', base_name, ceil(start_idx/10)));
fid = fopen(output_file, 'w');
if fid < 0
    result.n_temps = 0;
    result.success = false;
    return;
end

fprintf(fid, 'Temperature,SelectedLambda,MeanResidual,PeakCount\n');

count = 0;
for t = start_idx:end_idx
    try
        block = data(:, (t-1)*3 + 1 : (t-1)*3 + 3);
        valid = all(isfinite(block), 2) & block(:,1) > 0;
        block = block(valid, :);
        if size(block, 1) < 5, continue; end
        freq_vec = block(:, 1);
        Z_mag    = block(:, 2);
        Z_phase  = block(:, 3);
        Z_exp = Z_mag .* exp(1i * Z_phase * pi / 180);
        
        residual_best = inf;
        idx_best = 1;
        
        for i = 1:numel(lambda_values)
            [gamma, R_inf] = tr_drt_local(freq_vec, Z_exp, lambda_values(i));
            A_re = calc_A_re_local(freq_vec);
            A_im = calc_A_im_local(freq_vec);
            Z_fit = R_inf + A_re * gamma + 1i * (A_im * gamma);
            residual = mean(abs(Z_exp - Z_fit));
            
            if residual < residual_best
                residual_best = residual;
                idx_best = i;
            end
        end
        
        [gamma, ~] = tr_drt_local(freq_vec, Z_exp, lambda_values(idx_best));
        n_peaks = find_peaks_simple_local(gamma);
        
        fprintf(fid, '%.6f,%.8e,%.6e,%d\n', ...
            temperature(t), lambda_values(idx_best), residual_best, n_peaks);
        
        count = count + 1;
    catch
        continue;
    end
end

fclose(fid);
result.n_temps = count;
result.success = true;
result.output_file = output_file;
end

%% ---- HELPER: Process loaded xlsread data into clean matrix + temperature vector ----
function [data_out, temperature] = load_dataset_local(data_in, raw)
% data_in: numeric matrix from xlsread
% raw: cell array from xlsread
data_out = data_in;
n_cols = size(data_in, 2);
n_sets = floor(n_cols / 3);

temperature = [];

% Try row 1 first (simple parsing)
if ~isempty(raw) && size(raw, 1) >= 1
    header = raw(1, :);
    for k = 1:numel(header)
        tok = header{k};
        if isnumeric(tok) && isscalar(tok) && isfinite(tok) && tok > 0 && tok < 1000
            temperature(end + 1, 1) = tok; %#ok<AGROW>
        elseif ischar(tok)
            tok(tok == 'K' | tok == 'k') = [];
            val = str2double(strtrim(tok));
            if ~isnan(val)
                temperature(end + 1, 1) = val; %#ok<AGROW>
            end
        end
    end
end

if numel(temperature) >= n_sets
    temperature = temperature(1:n_sets);
    return;
end

% If row 1 didn't work, scan rows 1-12 with offsets 0-2 (for datasets like S2422Al)
if isempty(raw)
    return;
end

best_row = 0;
best_count = 0;
best_offset = 0;
scan_rows = min(size(raw, 1), 12);

for r = 1:scan_rows
    for off = 0:2
        count = 0;
        for c = (1+off):3:(3*n_sets)
            if c <= size(raw, 2)
                v = raw{r, c};
                if ischar(v) && (~isempty(strfind(v, 'K')) || ~isempty(strfind(v, 'k')))
                    count = count + 1;
                end
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
    return;
end

% Extract temperatures from the best row/offset
temperature = nan(n_sets, 1);
for i = 1:n_sets
    c0 = 3 * (i - 1) + (1 + best_offset);
    if c0 >= 1 && c0 <= size(raw, 2)
        v0 = raw{best_row, c0};
        if ischar(v0)
            tok0 = strtrim(v0);
            tok0(tok0 == 'K' | tok0 == 'k') = [];
            val = str2double(strtrim(tok0));
            if ~isnan(val) && isfinite(val) && val > 0 && val < 1000
                temperature(i) = val;
            end
        elseif isnumeric(v0) && isscalar(v0) && isfinite(v0) && v0 > 0 && v0 < 1000
            temperature(i) = v0;
        end
    end
end

% Remove NaN entries
temperature = temperature(~isnan(temperature));
end

%% ---- HELPER: GCV lambda selection (Im inversion, Re+Im cross-validation) ----
function [lambda_best, gcv_best, mean_resid, n_peaks] = ...
    select_lambda_gcv_local(freq_vec, Z_exp, lambda_values)

n = numel(freq_vec);
n_lambda = numel(lambda_values);
n_data = n;

A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);

gcv_scores = inf(n_lambda, 1);

for i = 1:n_lambda
    lam = lambda_values(i);

    % Invert using imaginary part only
    M = [A_im, zeros(n,1); sqrt(lam/2) * eye(n), zeros(n,1)];
    b = [imag(Z_exp); zeros(n,1)];
    opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
    x = lsqnonneg(M, b, [], opts);
    gamma_i = x(1:n);
    R_inf_i = x(n+1);

    % Validate using Re+Im residuals
    dz_re = real(Z_exp) - (R_inf_i + A_re * gamma_i);
    dz_im = imag(Z_exp) - (A_im * gamma_i);
    res_sq = norm([dz_re; dz_im], 2)^2;

    dof_red = 1 + log10(lam + 1) * n_data;
    dof_red = max(1, min(dof_red, n_data - 1));
    denom = 1 - dof_red / n_data;
    if abs(denom) > 1e-6
        gcv_scores(i) = res_sq / (denom^2);
    end
end

[gcv_best, idx_best] = min(gcv_scores);
lambda_best = lambda_values(idx_best);

% Final DRT at optimal lambda (imaginary inversion)
M = [A_im, zeros(n,1); sqrt(lambda_best/2) * eye(n), zeros(n,1)];
b = [imag(Z_exp); zeros(n,1)];
opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
x = lsqnonneg(M, b, [], opts);
gamma_opt = x(1:n);
R_inf_opt = x(n+1);

Z_fit = R_inf_opt + A_re * gamma_opt + 1i * (A_im * gamma_opt);
mean_resid = mean(abs(Z_exp - Z_fit));
n_peaks = find_peaks_simple_local(gamma_opt);
end

%% ---- HELPER: Bayes-DRT lambda selection (Re/Im CV) ----
function [lambda_best, cv_real, cv_imag, cv_total, mean_residual, n_peaks] = ...
    select_lambda_reimcv_local(freq_vec, Z_exp, lambda_values)

cv_real_all = [];
cv_imag_all = [];
n_lambda = numel(lambda_values);

for i = 1:n_lambda
    lam = lambda_values(i);
    [gamma_re, ~] = tr_drt_component_local(freq_vec, Z_exp, lam, 'real');
    [gamma_im, ~] = tr_drt_component_local(freq_vec, Z_exp, lam, 'imag');
    
    A_re = calc_A_re_local(freq_vec);
    A_im = calc_A_im_local(freq_vec);
    
    Z_from_re = A_re * gamma_re + 1i * (A_im * gamma_re);
    Z_from_im = A_re * gamma_im + 1i * (A_im * gamma_im);
    
    cv_imag_all(i) = mean(abs(imag(Z_exp) - imag(Z_from_re)).^2);
    cv_real_all(i) = mean(abs(real(Z_exp) - real(Z_from_im)).^2);
end

cv_total_all = cv_real_all + cv_imag_all;
[~, idx_best] = min(cv_total_all);
lambda_best = lambda_values(idx_best);
cv_real = cv_real_all(idx_best);
cv_imag = cv_imag_all(idx_best);
cv_total = cv_total_all(idx_best);

[gamma, ~] = tr_drt_local(freq_vec, Z_exp, lambda_best);
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_fit = A_re * gamma + 1i * (A_im * gamma);
mean_residual = mean(abs(Z_exp - Z_fit));
n_peaks = find_peaks_simple_local(gamma);
end

%% ---- HELPER: Ridge regression DRT solver ----
function [gamma, R_inf] = tr_drt_local(freq_vec, Z_exp, lambda)
freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
n = numel(freq_vec);

A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(lambda/2) * eye(n), zeros(n,1)];
b = [real(Z_exp); imag(Z_exp); zeros(n,1)];

opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n);
R_inf = x(n+1);
end

%% ---- HELPER: DRT component solver ----
function [gamma, R_inf] = tr_drt_component_local(freq_vec, Z_exp, lambda, component)
freq_vec = freq_vec(:);
Z_exp = Z_exp(:);
n = numel(freq_vec);

A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);

switch lower(component)
    case 'real'
        M = [A_re, ones(n,1); sqrt(lambda/2) * eye(n), zeros(n,1)];
        b = [real(Z_exp); zeros(n,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(lambda/2) * eye(n), zeros(n,1)];
        b = [imag(Z_exp); zeros(n,1)];
    otherwise
        error('Unknown component mode: %s', component);
end

opts = optimset('Display', 'off', 'MaxIter', 10000, 'TolX', 1e-10);
x = lsqnonneg(M, b, [], opts);
gamma = x(1:n);
R_inf = x(n+1);
end

%% ---- HELPER: Calculate A_re matrix ----
function A_re = calc_A_re_local(freq_vec)
freq_vec = freq_vec(:);
omega = 2 * pi * freq_vec;
tau = 1 ./ freq_vec;
n = numel(freq_vec);
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

%% ---- HELPER: Calculate A_im matrix ----
function A_im = calc_A_im_local(freq_vec)
freq_vec = freq_vec(:);
omega = 2 * pi * freq_vec;
tau = 1 ./ freq_vec;
n = numel(freq_vec);
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

%% ---- HELPER: Find peaks ----
function n_peaks = find_peaks_simple_local(gamma)
n = numel(gamma);
if n < 3
    n_peaks = 0;
    return;
end

n_peaks = 0;
for i = 2:(n-1)
    if gamma(i) > gamma(i-1) && gamma(i) > gamma(i+1)
        n_peaks = n_peaks + 1;
    end
end
end

%% ---- HELPER: Variance estimation ----
function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, lambda, n_boot)
n = numel(freq_vec);
if nargin < 4 || isempty(n_boot), n_boot = 3; end

sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);
sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);

gamma_samples = zeros(n, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n,1) + 1i * sigma_im * randn(n,1);
    Z_boot = Z_exp(:) + noise;
    [gk, ~] = tr_drt_local(freq_vec, Z_boot, lambda);
    gamma_samples(:, k) = gk(:);
end

mean_var = mean(var(gamma_samples, 0, 2));
end
