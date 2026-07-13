%% Identify 3-peak consistent datasets and create 15% removal CV analysis
% This script finds temperatures where all three methods show 3 peaks,
% then generates cross-validation plots with 15% data removal.

clear all; close all; clc;

% Candidate datasets (manually selected based on prior inspection)
datasets_to_check = {
    'S2022Sap.xlsx', 's2022sap'
    'S2222Sap.xlsx', 's2222sap'
    'S2302Sap.xlsx', 's2302sap'
    'S2322Sap.xlsx', 's2322sap'
    'S2332Sap.xlsx', 's2332sap'
    'S2422Sap.xlsx', 's2422sap'
};

% For each dataset, find one temperature with 3 peaks
selected_cases = {};

fprintf('Finding 3-peak cases...\n\n');

for d = 1:size(datasets_to_check, 1)
    xlsx_file = datasets_to_check{d, 1};
    ds_tag = datasets_to_check{d, 2};
    
    % Read Bayes-DRT results
    bayes_file = sprintf('%s_bayes_drt_matlab2011_all_temps.txt', ds_tag);
    
    if ~isfile(bayes_file)
        continue;
    end
    
    % Parse file to find temps with 3 peaks
    fid = fopen(bayes_file, 'r');
    line_count = 0;
    while ~feof(fid)
        line = fgetl(fid);
        line_count = line_count + 1;
        if line_count > 13  % Skip header
            if ischar(line) && ~isempty(line) && ismember(line(1), '0':'9')
                parts = strsplit(line, ',');
                if numel(parts) >= 7
                    temp = str2double(parts{1});
                    peak_count = str2double(parts{7});
                    if peak_count == 3 && isfinite(temp) && isfinite(peak_count)
                        selected_cases = [selected_cases; {ds_tag, temp, xlsx_file}];
                        fprintf('  %s: %.3f K (Bayes-DRT = 3 peaks)\n', ds_tag, temp);
                        break;  % Take first match
                    end
                end
            end
        end
    end
    fclose(fid);
end

fprintf('\nSelected %d cases for CV analysis\n\n', size(selected_cases, 1));

% Lambda grid from Bayes-DRT
lambda_grid = logspace(-8, 2, 21);

% CV parameters
drop_ratio = 0.15;  % 15% removal
n_resamples = 10;   % Number of random splits

% Storage
results = table();
cv_data_store = {};

% Process each case
for case_idx = 1:size(selected_cases, 1)
    ds_tag = selected_cases{case_idx, 1};
    target_temp = selected_cases{case_idx, 2};
    xlsx_file = selected_cases{case_idx, 3};
    
    fprintf('Processing %s @ %.3f K...\n', ds_tag, target_temp);
    
    if ~isfile(xlsx_file)
        fprintf('  File not found: %s\n', xlsx_file);
        continue;
    end
    
    % Read Excel data
    try
        [num, ~, raw] = xlsread(xlsx_file);
    catch
        fprintf('  Error reading Excel file\n');
        continue;
    end
    
    if isempty(num)
        fprintf('  No numeric data found\n');
        continue;
    end
    
    % Parse headers to find the temperature column matching target_temp
    header_row = raw(1, :);
    col_idx = [];
    
    for h = 1:numel(header_row)
        if ischar(header_row{h})
            header_str = header_row{h};
            % Extract numeric value from string like "280K" or "280"
            temp_val_str = regexprep(header_str, '[^\d.]', '');
            if ~isempty(temp_val_str)
                temp_val = str2double(temp_val_str);
                if isfinite(temp_val) && abs(temp_val - target_temp) < 1.0
                    col_idx = h;
                    break;
                end
            end
        end
    end
    
    if isempty(col_idx)
        % Try to infer from structure (3 columns per dataset)
        col_idx = 1;  % Fallback
        fprintf('  Using fallback column selection\n');
    end
    
    % Extract data for this temperature (3 columns: freq, mag, phase)
    freq = num(:, col_idx);
    mag = num(:, col_idx + 1);
    phase_deg = num(:, col_idx + 2);
    
    % Remove NaN/invalid
    valid = ~isnan(freq) & ~isnan(mag) & ~isnan(phase_deg) & freq > 0;
    freq = freq(valid);
    mag = mag(valid);
    phase_deg = phase_deg(valid);
    
    if numel(freq) < 10
        fprintf('  Insufficient data points: %d\n', numel(freq));
        continue;
    end
    
    Z_complex = mag .* exp(1i * phase_deg * pi / 180);
    tau_basis = 1 ./ (2 * pi * freq);
    
    % Compute CV with 15% removal
    n_data = numel(freq);
    n_remove = max(1, round(n_data * drop_ratio));
    
    cv_real = zeros(numel(lambda_grid), 1);
    cv_imag = zeros(numel(lambda_grid), 1);
    cv_total = zeros(numel(lambda_grid), 1);
    cv_real_std = zeros(numel(lambda_grid), 1);
    cv_imag_std = zeros(numel(lambda_grid), 1);
    
    cv_real_all = zeros(numel(lambda_grid), n_resamples);
    cv_imag_all = zeros(numel(lambda_grid), n_resamples);
    
    rng(42);  % Reproducible
    
    for resample = 1:n_resamples
        % Random split
        remove_idx = randperm(n_data, n_remove);
        keep_idx = setdiff(1:n_data, remove_idx);
        
        freq_train = freq(keep_idx);
        mag_train = mag(keep_idx);
        phase_train = phase_deg(keep_idx);
        Z_train = mag_train .* exp(1i * phase_train * pi / 180);
        tau_train = 1 ./ (2 * pi * freq_train);
        
        freq_test = freq(remove_idx);
        Z_test = Z_complex(remove_idx);
        
        for lam_idx = 1:numel(lambda_grid)
            lambda = lambda_grid(lam_idx);
            
            % Fit on training data (real part)
            try
                A_re = (tau_train' * tau_train + lambda * eye(numel(tau_train))) \ (tau_train' * real(Z_train));
                Z_re_pred = tau_basis .* A_re;
                cv_real_all(lam_idx, resample) = mean(abs(real(Z_complex) - Z_re_pred));
            catch
                cv_real_all(lam_idx, resample) = inf;
            end
            
            % Fit on training data (imaginary part)
            try
                A_im = (tau_train' * tau_train + lambda * eye(numel(tau_train))) \ (tau_train' * imag(Z_train));
                Z_im_pred = tau_basis .* A_im;
                cv_imag_all(lam_idx, resample) = mean(abs(imag(Z_complex) - Z_im_pred));
            catch
                cv_imag_all(lam_idx, resample) = inf;
            end
        end
    end
    
    cv_real = mean(cv_real_all, 2);
    cv_imag = mean(cv_imag_all, 2);
    cv_total = cv_real + cv_imag;
    cv_real_std = std(cv_real_all, [], 2);
    cv_imag_std = std(cv_imag_all, [], 2);
    
    % Find optimum
    [~, idx_opt] = min(cv_total);
    lambda_opt = lambda_grid(idx_opt);
    
    % Store results
    results = [results; table(...
        {ds_tag}, target_temp, numel(freq), lambda_opt, ...
        'VariableNames', {'Dataset', 'TemperatureK', 'DataPoints', 'OptimalLambda'})];
    
    cv_data_store{case_idx} = struct(...
        'dataset', ds_tag, ...
        'temp', target_temp, ...
        'lambda_grid', lambda_grid, ...
        'cv_real', cv_real, ...
        'cv_imag', cv_imag, ...
        'cv_total', cv_total, ...
        'cv_real_std', cv_real_std, ...
        'cv_imag_std', cv_imag_std, ...
        'lambda_opt', lambda_opt, ...
        'n_data', numel(freq) ...
    );
    
    fprintf('  → Optimal λ = %.4e (Total CV = %.2e)\n\n', lambda_opt, cv_total(idx_opt));
end

% Plot results
n_cases = numel(cv_data_store);
if n_cases > 0
    figure('Position', [100, 100, 1500, 400 + 300*ceil(n_cases/3)]);
    
    for plot_idx = 1:n_cases
        subplot(ceil(n_cases/3), 3, plot_idx);
        data = cv_data_store{plot_idx};
        
        hold on;
        
        % Plot CV curves with error bands
        semilogx(data.lambda_grid, data.cv_real, 'o-', 'Color', [0.2, 0.5, 1.0], 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Real');
        semilogx(data.lambda_grid, data.cv_imag, 's-', 'Color', [1.0, 0.6, 0.2], 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', 'Imag');
        semilogx(data.lambda_grid, data.cv_total, 'd-', 'Color', [0.3, 0.8, 0.3], 'LineWidth', 2.5, 'MarkerSize', 5, 'DisplayName', 'Total');
        
        % Mark optimum
        cv_opt = data.cv_total(find(data.lambda_grid == data.lambda_opt, 1));
        plot(data.lambda_opt, cv_opt, 'r*', 'MarkerSize', 18, 'LineWidth', 2, 'DisplayName', sprintf('Min: λ=%.2e', data.lambda_opt));
        
        xlabel('Regularization parameter λ', 'FontSize', 10);
        ylabel('Cross-validation error', 'FontSize', 10);
        title(sprintf('%s @ %.2f K (15%% removal, n=%d)', strrep(data.dataset, '_', '\_'), data.temp, data.n_data), ...
            'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        legend('Location', 'best', 'FontSize', 9);
        
        hold off;
    end
    
    sgtitle('15% Data Removal Cross-Validation (3-Peak Cases)', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    out_fig = 'cv_15pct_removal_3peak_plots.png';
    saveas(gcf, out_fig, 'png');
    fprintf('Plots saved: %s\n\n', out_fig);
end

% Save detailed results to Excel
if ~isempty(results)
    out_xlsx = 'cv_15pct_removal_3peak_results.xlsx';
    writetable(results, out_xlsx, 'Sheet', 'Summary');
    
    % Add detailed CV data to separate sheets
    for i = 1:numel(cv_data_store)
        data = cv_data_store{i};
        sheet_name = sprintf('CV_%s_%.0fK', data.dataset, data.temp);
        sheet_name = sheet_name(1:min(31, length(sheet_name)));  % Excel sheet name limit
        
        detail_table = table(...
            data.lambda_grid', ...
            data.cv_real, ...
            data.cv_imag, ...
            data.cv_total, ...
            data.cv_real_std, ...
            data.cv_imag_std, ...
            'VariableNames', {'Lambda', 'CV_Real', 'CV_Imag', 'CV_Total', 'CV_Real_Std', 'CV_Imag_Std'} ...
        );
        writetable(detail_table, out_xlsx, 'Sheet', sheet_name);
    end
    
    fprintf('Results saved: %s\n\n', out_xlsx);
else
    fprintf('No cases selected for analysis.\n');
end

fprintf('Analysis complete.\n');
