% Rerun S2322Al and S2422Al with fixed numeric header parsing
setenv('MATLAB_ANALYSIS_REPO_ROOT', pwd);

datasets = {
    'S2322Al.xlsx', 's2322al';
    'S2422Al.xlsx', 's2422al'
};

fprintf('\n====== RERUN: S2322Al and S2422Al with Fixed Parsing ======\n');
fprintf('Started: %s\n\n', datetime('now'));

output_file = fullfile(pwd, 'rerun_s2322al_s2422al_status.txt');
fid_status = fopen(output_file, 'w');
fprintf(fid_status, 'Rerun S2322Al and S2422Al - Fixed Numeric Header Parsing\n');
fprintf(fid_status, 'Started: %s\n', datetime('now'));
fprintf(fid_status, '=====================================\n\n');

for d = 1:size(datasets, 1)
    dataset_file = datasets{d, 1};
    ds_tag = datasets{d, 2};
    
    fprintf('Processing %s...\n', ds_tag);
    fprintf(fid_status, 'Dataset: %s\n', dataset_file);
    
    % Count temperatures
    try
        [data, ~, raw] = xlsread(dataset_file); %#ok<XLSRD>
        n_cols = size(data, 2);
        n_sets = floor(n_cols / 3);
        
        % Parse temperatures with FIXED parsing (handles numeric headers)
        temperature = [];
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
        end
        
        if isempty(temperature)
            fprintf(fid_status, '  ERROR: Found 0 temperatures\n');
            fprintf('  ERROR: Found 0 temperatures\n');
            continue;
        end
        
        n_temps = numel(temperature);
        fprintf(fid_status, '  Found %d temperatures\n', n_temps);
        fprintf('  Found %d temperatures\n', n_temps);
        
        % Calculate batches (10 temps per batch)
        n_batches = ceil(n_temps / 10);
        fprintf(fid_status, '  Batches: %d\n', n_batches);
        
        success_batches = 0;
        failed_batches = 0;
        
        % Process batches
        for b = 1:n_batches
            start_idx = (b-1)*10 + 1;
            end_idx = min(b*10, n_temps);
            
            fprintf('    Batch %d [%d-%d]... ', b, start_idx, end_idx);
            fprintf(fid_status, '    Batch %d [%d-%d]... ', b, start_idx, end_idx);
            
            try
                fprintf('[All Methods] Processing %s, temps %d to %d\n', dataset_file, start_idx, end_idx);
                fprintf(fid_status, '[All Methods] Processing %s, temps %d to %d\n', dataset_file, start_idx, end_idx);
                
                % Run all three methods
                fprintf('[Bayes-DRT] Running batch...\n');
                result_bayes = run_bayes_drt_batch_local(dataset_file, start_idx, end_idx);
                fprintf('[Bayes-DRT]  Batch complete. %d temps processed.\n', result_bayes.n_temps);
                fprintf(fid_status, '[Bayes-DRT]  Batch complete. %d temps processed.\n', result_bayes.n_temps);
                
                fprintf('[Paper Method] Running batch...\n');
                result_paper = run_paper_method_batch_local(dataset_file, start_idx, end_idx);
                fprintf('[Paper Method]  Batch complete. %d temps processed.\n', result_paper.n_temps);
                fprintf(fid_status, '[Paper Method]  Batch complete. %d temps processed.\n', result_paper.n_temps);
                
                fprintf('[Residual Method] Running batch...\n');
                result_residual = run_residual_method_batch_local(dataset_file, start_idx, end_idx);
                fprintf('[Residual Method]  Batch complete. %d temps processed.\n', result_residual.n_temps);
                fprintf(fid_status, '[Residual Method]  Batch complete. %d temps processed.\n', result_residual.n_temps);
                
                fprintf(fid_status, '  Batch OK\n');
                fprintf('OK\n');
                success_batches = success_batches + 1;
                
            catch ME
                fprintf(fid_status, '  FAIL: %s\n', ME.message);
                fprintf('FAIL: %s\n', ME.message);
                failed_batches = failed_batches + 1;
            end
        end
        
        fprintf(fid_status, '  Result: %d/%d batches OK\n\n', success_batches, n_batches);
        fprintf('  Result: %d/%d batches OK\n', success_batches, n_batches);
        
    catch ME
        fprintf(fid_status, '  Error: %s\n\n', ME.message);
        fprintf('  Error: %s\n', ME.message);
    end
end

fprintf(fid_status, '=====================================\n');
fprintf(fid_status, 'Completed: %s\n', datetime('now'));
fclose(fid_status);

fprintf('\nStatus file: %s\n', output_file);
exit;

%% ---- BAYES-DRT BATCH HANDLER (copied from run_all_three_methods_batch) ----
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

%% ---- PAPER METHOD BATCH HANDLER ----
function result = run_paper_method_batch_local(dataset_file, start_idx, end_idx)
[~, base_name, ~] = fileparts(dataset_file);
base_name = lower(base_name);
repo_root = getenv('MATLAB_ANALYSIS_REPO_ROOT');
if isempty(repo_root), repo_root = pwd; end

lambda_values = [1e-4; 1e-3; 1e-2; 1e-1];

% Read dataset using xlsread
try
    [data, ~, raw] = xlsread(dataset_file); %#ok<XLSRD>
    [data, temperature] = load_dataset_local(data, raw);
    if isempty(temperature)
        result.n_temps = 0; result.success = false; return;
    end
catch ME
    fprintf('Paper read error: %s\n', ME.message);
    result.n_temps = 0; result.success = false; return;
end

end_idx = min(end_idx, numel(temperature));
if start_idx > end_idx
    result.n_temps = 0;
    result.success = true;
    return;
end

output_file = fullfile(repo_root, sprintf('%s_paper_method_b%d.txt', base_name, ceil(start_idx/10)));
fid = fopen(output_file, 'w');
if fid < 0
    result.n_temps = 0;
    result.success = false;
    return;
end

fprintf(fid, 'Temperature,SelectedLambda,RID,CV,Variance,Score,PeakCount\n');

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
        
        [rid_vec, cv_vec, var_vec] = compute_paper_metrics_local(freq_vec, Z_exp, lambda_values);
        idx_best = select_lambda_idx_local(rid_vec, cv_vec, var_vec);
        lam_best = lambda_values(idx_best);
        
        % Fit with best lambda
        [gamma, ~] = tr_drt_local(freq_vec, Z_exp, lam_best);
        n_peaks = find_peaks_simple_local(gamma);
        
        norm_rid = normalize_metric_local(rid_vec);
        norm_cv = normalize_metric_local(cv_vec);
        norm_var = normalize_metric_local(var_vec);
        score = norm_rid + norm_cv + norm_var;
        
        fprintf(fid, '%.6f,%.8e,%.6e,%.6e,%.6e,%.6e,%d\n', ...
            temperature(t), lam_best, rid_vec(idx_best), cv_vec(idx_best), var_vec(idx_best), score(idx_best), n_peaks);
        
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

%% ---- SHARED HELPERS ----
function [data_out, temperature] = load_dataset_local(data_in, raw)
data_out = data_in;
n_cols = size(data_in, 2);
n_sets = floor(n_cols / 3);

temperature = [];
if ~isempty(raw) && size(raw, 1) >= 1
    header = raw(1, :);
    for k = 1:numel(header)
        tok = header{k};
        % Handle both text headers and numeric headers (FIXED)
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
end
end

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
[cv_total, idx_best] = min(cv_total_all);
cv_real = cv_real_all(idx_best);
cv_imag = cv_imag_all(idx_best);
lambda_best = lambda_values(idx_best);

[gamma, ~] = tr_drt_local(freq_vec, Z_exp, lambda_best);
n_peaks = find_peaks_simple_local(gamma);
mean_residual = mean(abs(Z_exp - (A_re*gamma + 1i*A_im*gamma)));
end

function [rid, cv, var] = compute_paper_metrics_local(freq_vec, Z_exp, lambda_values)
rid = [];
cv = [];
var = [];

for i = 1:numel(lambda_values)
    [gamma, ~] = tr_drt_local(freq_vec, Z_exp, lambda_values(i));
    A_re = calc_A_re_local(freq_vec);
    A_im = calc_A_im_local(freq_vec);
    Z_fit = A_re * gamma + 1i * A_im * gamma;
    Z_residual = Z_exp - Z_fit;
    
    rid(i) = mean(abs(Z_residual));
    cv(i) = mean(abs(Z_residual).^2);
    var(i) = std(abs(Z_residual).^2);
end
end

function idx = select_lambda_idx_local(rid_vec, cv_vec, var_vec)
n_rid = normalize_metric_local(rid_vec);
n_cv = normalize_metric_local(cv_vec);
n_var = normalize_metric_local(var_vec);
score = n_rid + n_cv + n_var;
[~, idx] = min(score);
end

function norm_v = normalize_metric_local(v)
v = v(:);
v_min = min(v);
v_max = max(v);
if v_max > v_min
    norm_v = (v - v_min) / (v_max - v_min);
else
    norm_v = zeros(size(v));
end
end

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

function n_peaks = find_peaks_simple_local(gamma)
if numel(gamma) < 3
    n_peaks = 0;
    return;
end
d1 = diff(gamma);
sign_changes = sum(diff(sign(d1)) ~= 0);
n_peaks = max(0, sign_changes);
end
