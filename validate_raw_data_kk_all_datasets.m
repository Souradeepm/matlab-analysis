%% Validate raw data using Kramers-Kronig criterion for all 11 datasets
% This script checks K-K consistency of raw measured impedance before DRT analysis
% and generates a comprehensive report with data quality metrics

clear all; close all; clc;

repo_root = fileparts(mfilename('fullpath'));

% All 11 datasets
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

% Output files
kk_summary_file = fullfile(repo_root, 'raw_data_kk_validation_summary.txt');
kk_excel_file = fullfile(repo_root, 'raw_data_kk_validation.xlsx');

fprintf('Raw Data K-K Validation for All 11 Datasets\n');
fprintf('===========================================\n\n');

% Storage
all_kk_data = table();
kk_summary_lines = {};

fid_summary = fopen(kk_summary_file, 'w');
fprintf(fid_summary, 'Raw Data K-K Validation Report\n');
fprintf(fid_summary, 'Generated: %s\n', datetime('now'));
fprintf(fid_summary, '==================================\n\n');

% Process each dataset
for d = 1:size(datasets, 1)
    xlsx_file = datasets{d, 1};
    ds_tag = datasets{d, 2};
    input_path = fullfile(repo_root, xlsx_file);
    
    if ~isfile(input_path)
        fprintf('Skipping %s: file not found\n', xlsx_file);
        continue;
    end
    
    fprintf('Processing %s (%s)...\n', ds_tag, xlsx_file);
    fprintf(fid_summary, '\nDataset: %s\n', ds_tag);
    fprintf(fid_summary, '-----------------------------------\n');
    
    % Read Excel file
    try
        [num, ~, raw] = xlsread(input_path);
    catch ME
        fprintf('Error reading %s: %s\n', xlsx_file, ME.message);
        continue;
    end
    
    if isempty(num)
        fprintf('No numeric data in %s\n', xlsx_file);
        continue;
    end
    
    % Parse header to get temperatures and column structure
    header_row = raw(1, :);
    temps = [];
    col_groups = [];  % Store [start_col, end_col] for each temperature
    
    for h = 1:3:size(raw, 2)  % Temperature data comes in groups of 3 columns
        if h + 2 <= size(raw, 2)
            header_str = '';
            if ischar(header_row{h})
                header_str = header_row{h};
            end
            
            % Extract temperature from header
            temp_val_str = regexprep(header_str, '[^\d.]', '');
            if ~isempty(temp_val_str)
                temp_val = str2double(temp_val_str);
                if isfinite(temp_val) && temp_val > 0 && temp_val < 1000
                    temps = [temps; temp_val];
                    col_groups = [col_groups; h, h+2];
                end
            end
        end
    end
    
    fprintf('Found %d temperature columns (groups of 3)\n', numel(temps));
    
    % K-K validation for each temperature
    kk_pass = 0;
    kk_warning = 0;
    kk_fail = 0;
    pass_list = [];
    warning_list = [];
    fail_list = [];
    
    for t = 1:numel(temps)
        temp_k = temps(t);
        col_start = col_groups(t, 1);
        
        % Extract frequency, magnitude, phase for this temperature (3 consecutive columns)
        if col_start + 2 > size(num, 2)
            continue;
        end
        
        freq = num(:, col_start);
        mag = num(:, col_start + 1);
        phase_deg = num(:, col_start + 2);
        
        % Clean data
        valid = ~isnan(freq) & ~isnan(mag) & ~isnan(phase_deg) & freq > 0 & mag > 0;
        freq = freq(valid);
        mag = mag(valid);
        phase_deg = phase_deg(valid);
        
        if numel(freq) < 10
            continue;
        end
        
        % Convert to impedance
        Z_complex = mag .* exp(1i * phase_deg * pi / 180);
        
        % Perform K-K test
        kk_result = perform_kk_test_all_datasets(freq, Z_complex);
        
        % Log results
        if strcmp(kk_result.status, 'pass')
            kk_pass = kk_pass + 1;
            pass_list = [pass_list; temp_k];
        elseif strcmp(kk_result.status, 'warning')
            kk_warning = kk_warning + 1;
            warning_list = [warning_list; temp_k];
        else
            kk_fail = kk_fail + 1;
            fail_list = [fail_list; temp_k];
        end
        
        % Store in table
        all_kk_data = [all_kk_data; table(...
            {ds_tag}, temp_k, numel(freq), ...
            kk_result.real_residual_pct, kk_result.imag_residual_pct, kk_result.total_residual_pct, ...
            kk_result.max_residual_pct, kk_result.status_code, {kk_result.status}, ...
            'VariableNames', {'Dataset', 'TemperatureK', 'DataPoints', 'Real_ResidualPct', ...
                             'Imag_ResidualPct', 'Total_ResidualPct', 'Max_ResidualPct', 'StatusCode', 'Status'})];
    end
    
    % Summary for this dataset
    n_temps = numel(temps);
    pass_pct = 100 * kk_pass / max(1, n_temps);
    warn_pct = 100 * kk_warning / max(1, n_temps);
    fail_pct = 100 * kk_fail / max(1, n_temps);
    
    fprintf(fid_summary, 'Total temperatures: %d\n', n_temps);
    fprintf(fid_summary, 'K-K PASS:    %3d (%.1f%%)\n', kk_pass, pass_pct);
    fprintf(fid_summary, 'K-K WARNING: %3d (%.1f%%)\n', kk_warning, warn_pct);
    fprintf(fid_summary, 'K-K FAIL:    %3d (%.1f%%)\n', kk_fail, fail_pct);
    
    if ~isempty(pass_list)
        fprintf(fid_summary, 'PASS temps:    '); fprintf(fid_summary, '%.1f ', pass_list); fprintf(fid_summary, '\n');
    end
    if ~isempty(warning_list)
        fprintf(fid_summary, 'WARNING temps: '); fprintf(fid_summary, '%.1f ', warning_list); fprintf(fid_summary, '\n');
    end
    if ~isempty(fail_list)
        fprintf(fid_summary, 'FAIL temps:    '); fprintf(fid_summary, '%.1f ', fail_list); fprintf(fid_summary, '\n');
    end
    
    fprintf('%s: %d pass, %d warning, %d fail\n\n', ds_tag, kk_pass, kk_warning, kk_fail);
end

% Write Excel file with detailed results
if ~isempty(all_kk_data)
    writetable(all_kk_data, kk_excel_file, 'Sheet', 'AllResults');
    
    % Summary sheet
    dataset_names = unique(all_kk_data.Dataset);
    summary_table = table();
    
    for d = 1:numel(dataset_names)
        ds_mask = strcmp(all_kk_data.Dataset, dataset_names{d});
        ds_data = all_kk_data(ds_mask, :);
        
        n_pass = sum(strcmp(ds_data.Status, 'pass'));
        n_warn = sum(strcmp(ds_data.Status, 'warning'));
        n_fail = sum(strcmp(ds_data.Status, 'fail'));
        n_total = height(ds_data);
        
        pass_pct = 100 * n_pass / max(1, n_total);
        warn_pct = 100 * n_warn / max(1, n_total);
        fail_pct = 100 * n_fail / max(1, n_total);
        
        mean_total_residual = mean(ds_data.Total_ResidualPct);
        max_total_residual = max(ds_data.Total_ResidualPct);
        
        summary_table = [summary_table; table(...
            dataset_names(d), n_total, n_pass, n_warn, n_fail, ...
            pass_pct, warn_pct, fail_pct, mean_total_residual, max_total_residual, ...
            'VariableNames', {'Dataset', 'TotalTemps', 'Pass', 'Warning', 'Fail', ...
                             'PassPct', 'WarningPct', 'FailPct', 'MeanTotalResidualPct', 'MaxTotalResidualPct'})];
    end
    
    writetable(summary_table, kk_excel_file, 'Sheet', 'Summary');
end

fprintf(fid_summary, '\n\nOverall K-K Validation Statistics\n');
fprintf(fid_summary, '==================================\n');
fprintf(fid_summary, 'Total measurements: %d\n', height(all_kk_data));

if height(all_kk_data) > 0
    n_pass = sum(strcmp(all_kk_data.Status, 'pass'));
    n_warn = sum(strcmp(all_kk_data.Status, 'warning'));
    n_fail = sum(strcmp(all_kk_data.Status, 'fail'));
    n_total = height(all_kk_data);
    
    fprintf(fid_summary, 'K-K PASS:    %d (%.1f%%)\n', n_pass, 100*n_pass/max(1,n_total));
    fprintf(fid_summary, 'K-K WARNING: %d (%.1f%%)\n', n_warn, 100*n_warn/max(1,n_total));
    fprintf(fid_summary, 'K-K FAIL:    %d (%.1f%%)\n', n_fail, 100*n_fail/max(1,n_total));
end
fprintf(fid_summary, '\nCriteria: PASS <= 5%%, WARNING 5-10%%, FAIL > 10%%\n');

fclose(fid_summary);

fprintf('\n✓ K-K validation complete\n');
fprintf('  Summary: %s\n', kk_summary_file);
fprintf('  Excel:   %s\n', kk_excel_file);

%% Helper function: K-K test
function kk_result = perform_kk_test_all_datasets(freq_vec, Z_exp)
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

kk_result = struct(...
    'real_residual_pct', real_pct, ...
    'imag_residual_pct', imag_pct, ...
    'total_residual_pct', total_pct, ...
    'max_residual_pct', max_pct, ...
    'status', status, ...
    'status_code', status_code);
end
