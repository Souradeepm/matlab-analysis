%% Simple K-K validation on raw data - All 11 datasets
% Uses working data-loading logic from Bayes-DRT workflow

clear all; close all; clc;

repo_root = fileparts(mfilename('fullpath'));

% Datasets
datasets = {
    'S2022Sap.xlsx', 's2022sap'
    'S2022Al.xlsx', 's2022al'
    'S2222Sap.xlsx', 's2222sap'
    'S2222Al.xlsx', 's2222al'
    'S2302Sap.xlsx', 's2302sap'
    'S2302Al.xlsx', 's2302al'
    'S2322Sap.xlsx', 's2322sap'
    'S2322Al.xlsx', 's2322al'
    'S2332Sap.xlsx', 's2332sap'
    'S2422Sap.xlsx', 's2422sap'
    'S2422Al.xlsx', 's2422al'
};

out_file = fullfile(repo_root, 'raw_data_kk_check.txt');
fid = fopen(out_file, 'w');

fprintf(fid, 'Raw Data K-K Validation\n');
fprintf(fid, 'Generated: %s\n', datetime('now'));
fprintf(fid, '======================\n\n');

all_results = [];

for d = 1:size(datasets, 1)
    xlsx_file = datasets{d, 1};
    ds_tag = datasets{d, 2};
    input_path = fullfile(repo_root, xlsx_file);
    
    if ~isfile(input_path)
        fprintf(fid, '%s: FILE NOT FOUND\n\n', ds_tag);
        continue;
    end
    
    fprintf('Processing %s...\n', ds_tag);
    fprintf(fid, '%s:\n', ds_tag);
    
    try
        [num, ~, raw] = xlsread(input_path);
        
        % Use the tested parsing functions from Bayes-DRT
        n_sets = detect_n_sets_local(raw, num, size(num,2));
        if n_sets < 1
            n_sets = floor(size(num,2) / 3);
        end
        
        temps = parse_temperature_labels_local(raw, num, n_sets);
        
        fprintf(fid, '  Temperatures found: %d\n', numel(temps));
        
        kk_pass = 0;
        kk_warn = 0;
        kk_fail = 0;
        
        for t = 1:numel(temps)
            temp_k = temps(t);
            col_start = 3*t - 2;
            
            if col_start + 2 > size(num, 2)
                break;
            end
            
            % Extract data
            freq = num(:, col_start);
            mag = num(:, col_start + 1);
            phase_deg = num(:, col_start + 2);
            
            % Clean
            valid = ~isnan(freq) & ~isnan(mag) & ~isnan(phase_deg) & freq > 0 & mag > 0;
            freq = freq(valid);
            mag = mag(valid);
            phase_deg = phase_deg(valid);
            
            if numel(freq) < 5
                continue;
            end
            
            Z = mag .* exp(1i * phase_deg * pi / 180);
            
            % K-K test
            kk_res = kk_test_simple(freq, Z);
            
            if kk_res.total_pct <= 5
                kk_pass = kk_pass + 1;
                status = 'PASS';
            elseif kk_res.total_pct <= 10
                kk_warn = kk_warn + 1;
                status = 'WARN';
            else
                kk_fail = kk_fail + 1;
                status = 'FAIL';
            end
            
            fprintf(fid, '    %.1f K: %s (total residual=%.2f%%)\n', temp_k, status, kk_res.total_pct);
            
            all_results = [all_results; table({ds_tag}, temp_k, kk_res.total_pct, {status}, ...
                'VariableNames', {'Dataset', 'TempK', 'TotalResidualPct', 'Status'})];
        end
        
        fprintf(fid, '  Summary: %d pass, %d warn, %d fail\n\n', kk_pass, kk_warn, kk_fail);
        
    catch ME
        fprintf(fid, '  ERROR: %s\n\n', ME.message);
    end
end

fprintf(fid, '\nOVERALL SUMMARY\n');
fprintf(fid, '===============\n');
if ~isempty(all_results)
    n_total = height(all_results);
    n_pass = sum(strcmp(all_results.Status, 'PASS'));
    n_warn = sum(strcmp(all_results.Status, 'WARN'));
    n_fail = sum(strcmp(all_results.Status, 'FAIL'));
    
    fprintf(fid, 'Total measurements: %d\n', n_total);
    fprintf(fid, 'PASS: %d (%.1f%%)\n', n_pass, 100*n_pass/max(1,n_total));
    fprintf(fid, 'WARN: %d (%.1f%%)\n', n_warn, 100*n_warn/max(1,n_total));
    fprintf(fid, 'FAIL: %d (%.1f%%)\n', n_fail, 100*n_fail/max(1,n_total));
    fprintf(fid, '\nCriteria: PASS<=5%%, WARN 5-10%%, FAIL>10%%\n');
end

fclose(fid);

fprintf('\n✓ K-K check complete: %s\n', out_file);

%% Simple K-K test
function kk_res = kk_test_simple(freq, Z)
freq = freq(:); Z = Z(:);
omega = 2*pi*freq;
n_freq = numel(freq);

if n_freq < 5
    kk_res = struct('total_pct', nan);
    return;
end

n_basis = min(max(8, ceil(n_freq/3)), 20);
tau_min = 1/(2*pi*max(freq));
tau_max = 1/(2*pi*min(freq));
tau_basis = logspace(log10(tau_min/5), log10(tau_max*5), n_basis).';

B_re = zeros(n_freq, n_basis);
B_im = zeros(n_freq, n_basis);
for p = 1:n_freq
    for q = 1:n_basis
        wt = omega(p)*tau_basis(q);
        B_re(p,q) = 1/(1+wt^2);
        B_im(p,q) = -wt/(1+wt^2);
    end
end

M = [B_re, ones(n_freq,1), zeros(n_freq,1); ...
     B_im, zeros(n_freq,1), omega; ...
     sqrt(1e-4)*eye(n_basis), zeros(n_basis,2)];
b = [real(Z); imag(Z); zeros(n_basis,1)];

x = lsqnonneg(M, b);
Z_fit = x(n_basis+1) + B_re*x(1:n_basis) + 1i*(B_im*x(1:n_basis) + omega*x(n_basis+2));
residual = Z - Z_fit;

real_pct = 100*norm(real(residual))/max(norm(real(Z)), eps);
imag_pct = 100*norm(imag(residual))/max(norm(imag(Z)), eps);
total_pct = 100*norm([real(residual); imag(residual)])/max(norm([real(Z); imag(Z)]), eps);

kk_res = struct('real_pct', real_pct, 'imag_pct', imag_pct, 'total_pct', total_pct);
end

%% Helper functions from Bayes-DRT (copied for compatibility)
function n_sets = detect_n_sets_local(raw, data, n_col_data)
n_sets = 0;
if isempty(raw) || size(raw,1) < 1
    n_sets = floor(n_col_data / 3);
    return;
end
header_row = raw(1, :);
for h = 1:numel(header_row)
    val = extract_temp_val_local(header_row{h});
    if ~isnan(val)
        n_sets = n_sets + 1;
    end
end
if n_sets < 1
    n_sets = floor(n_col_data / 3);
end
end

function temps = parse_temperature_labels_local(raw, data, n_sets)
temps = [];
if isempty(raw)
    return;
end

% Try row 1 first with simple parsing
if size(raw,1) >= 1
    header_row = raw(1, :);
    for h = 1:min(numel(header_row), 3*n_sets)
        val = extract_temp_val_local(header_row{h});
        if ~isnan(val)
            temps = [temps; val]; %#ok<AGROW>
        end
    end
end

if numel(temps) >= n_sets
    temps = temps(1:n_sets);
    return;
end

% If row 1 didn't work, scan rows 1-12 with offsets 0-2 (for datasets like S2422Al)
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
temps = nan(n_sets, 1);
for i = 1:n_sets
    c0 = 3 * (i - 1) + (1 + best_offset);
    if c0 >= 1 && c0 <= size(raw, 2)
        v0 = raw{best_row, c0};
        if ischar(v0)
            tok0 = strtrim(v0);
            tok0(tok0 == 'K' | tok0 == 'k') = [];
            val = str2double(strtrim(tok0));
            if ~isnan(val) && isfinite(val) && val > 0 && val < 1000
                temps(i) = val;
            end
        elseif isnumeric(v0) && isscalar(v0) && isfinite(v0) && v0 > 0 && v0 < 1000
            temps(i) = v0;
        end
    end
end

% Remove NaN entries
temps = temps(~isnan(temps));
end

function val = extract_temp_val_local(cell_val)
% Extracts a temperature value from either a string or numeric cell.
% Returns NaN if not a valid temperature (0-1000 K).
val = NaN;
if isnumeric(cell_val) && isscalar(cell_val)
    % Numeric cell (e.g. S2322Al stores temps as numbers, not strings)
    if isfinite(cell_val) && cell_val > 0 && cell_val < 1000
        val = cell_val;
    end
elseif ischar(cell_val)
    val_str = regexprep(cell_val, '[^\d.]', '');
    if ~isempty(val_str)
        v = str2double(val_str);
        if isfinite(v) && v > 0 && v < 1000
            val = v;
        end
    end
end
end
