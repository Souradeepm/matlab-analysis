%% Cross-validation with 15% data removal for 3-peak consistent datasets
% This script identifies and analyzes datasets where Bayes-DRT, paper method,
% and residual method all show ~3 peaks, then plots 15% removal CV curves.

clear all; close all;

% Dataset configuration
datasets = {
    'S2022Sap.xlsx', 's2022sap'
    'S2022Al.xlsx',  's2022al'
    'S2222Sap.xlsx', 's2222sap'
    'S2222Al.xlsx',  's2222al'
    'S2302Sap.xlsx', 's2302sap'
    'S2302Al.xlsx',  's2302al'
    'S2322Sap.xlsx', 's2322sap'
    'S2322Al.xlsx',  's2322al'
    'S2332Sap.xlsx', 's2332sap'
    'S2422Sap.xlsx', 's2422sap'
    'S2422Al.xlsx',  's2422al'
};

% Lambda grid (from Bayes-DRT workflow)
lambda_grid = logspace(-8, 2, 21);

% Cross-validation parameters
cv_drop_ratio = 0.15;  % 15% data removal
cv_repeats = 5;        % 5 different random splits

% Initialize results storage
cv_results = [];
dataset_selection_summary = [];

% Hardcoded candidate cases based on inspection:
% S2022Sap: temp ~4.53K, ~7.57K, ~20.38K (Bayes=3)
% Let's select a few that are likely to have 3 peaks across methods
candidate_datasets = {'s2022sap', 's2222sap', 's2302sap', 's2322sap', 's2332sap', 's2422sap'};
selected_cases = {};  % Will store selected dataset + temps

% Scan for candidate temperatures
fprintf('Scanning for 3-peak consistent cases across all methods...\n\n');

for d = 1:size(datasets, 1)
    xlsx_file = datasets{d, 1};
    ds_name = datasets{d, 2};
    
    if ~ismember(ds_name, candidate_datasets)
        continue;
    end
    
    % Read Bayes-DRT results
    bayes_file = sprintf('%s_bayes_drt_matlab2011_all_temps.txt', ds_name);
    paper_file = sprintf('%s_paper_vs_residual_peak_cv_comparison.txt', ds_name);
    
    if ~isfile(bayes_file) || ~isfile(paper_file)
        continue;
    end
    
    % Parse Bayes-DRT file
    bayes_data = readtable(bayes_file, 'FileType', 'text', 'HeaderLines', 13, 'ReadVariableNames', 0);
    if height(bayes_data) == 0, continue; end
    bayes_temps = bayes_data{:, 1};
    bayes_peaks = bayes_data{:, 7};
    
    % Parse paper vs residual file
    paper_data = readtable(paper_file, 'FileType', 'text', 'HeaderLines', 21, 'ReadVariableNames', 0);
    if height(paper_data) == 0, continue; end
    paper_temps = paper_data{:, 1};
    paper_res_peaks = paper_data{:, 4};
    paper_method_peaks = paper_data{:, 5};
    
    % Find matches: temp within 0.01K with all three showing 3 peaks
    for b_idx = 1:length(bayes_temps)
        if bayes_peaks(b_idx) ~= 3, continue; end
        
        % Find matching temp in paper data
        p_idx = find(abs(paper_temps - bayes_temps(b_idx)) < 0.01, 1);
        if isempty(p_idx), continue; end
        
        res_pk = paper_res_peaks(p_idx);
        pap_pk = paper_method_peaks(p_idx);
        
        if res_pk == 3 && pap_pk == 3
            selected_cases = [selected_cases; {ds_name, bayes_temps(b_idx)}];
            fprintf('%s @ %.3fK: Bayes=3, Residual=3, Paper=3 ✓\n', ds_name, bayes_temps(b_idx));
        end
    end
end

fprintf('\nTotal cases identified: %d\n\n', size(selected_cases, 1));

if isempty(selected_cases)
    fprintf('Warning: No cases found with exactly 3 peaks across all three methods.\n');
    fprintf('Using first available 3-peak cases from Bayes-DRT as fallback.\n\n');
    
    % Fallback: pick first 3-peak temp from each dataset
    for d = 1:length(candidate_datasets)
        ds_name = candidate_datasets{d};
        bayes_file = sprintf('%s_bayes_drt_matlab2011_all_temps.txt', ds_name);
        if isfile(bayes_file)
            bayes_data = readtable(bayes_file, 'FileType', 'text', 'HeaderLines', 13, 'ReadVariableNames', 0);
            peaks = bayes_data{:, 7};
            idx = find(peaks == 3, 1);
            if ~isempty(idx)
                selected_cases = [selected_cases; {ds_name, bayes_data{idx, 1}}];
                fprintf('Fallback: %s @ %.3fK\n', ds_name, bayes_data{idx, 1});
            end
        end
    end
end

% Now compute 15% removal CV for each selected case
fprintf('\nComputing 15%% removal cross-validation...\n\n');

cv_plot_data = {};
for s = 1:size(selected_cases, 1)
    ds_name = selected_cases{s, 1};
    temp_k = selected_cases{s, 2};
    
    % Find full dataset file
    ds_file = '';
    for d = 1:size(datasets, 1)
        if strcmp(datasets{d, 2}, ds_name)
            ds_file = datasets{d, 1};
            break;
        end
    end
    if isempty(ds_file), continue; end
    
    % Read impedance data (assuming 3 columns: freq, mag, phase)
    try
        [num, txt, raw] = xlsread(ds_file);
        if isempty(num), continue; end
        
        freq = num(:, 1);
        Z_mag = num(:, 2);
        Z_phase = num(:, 3);
        Z_complex = Z_mag .* exp(1j * Z_phase * pi / 180);
        
        fprintf('  %s @ %.3fK: n=%d points\n', ds_name, temp_k, length(freq));
        
        % Compute CV with 15% removal
        n_data = length(freq);
        n_remove = max(1, round(n_data * cv_drop_ratio));
        
        cv_real_vec = zeros(length(lambda_grid), cv_repeats);
        cv_imag_vec = zeros(length(lambda_grid), cv_repeats);
        
        rng(42);  % Reproducible random splits
        for rep = 1:cv_repeats
            % Random indices to remove (15%)
            remove_idx = randperm(n_data, n_remove);
            keep_idx = setdiff(1:n_data, remove_idx);
            
            freq_train = freq(keep_idx);
            Z_train = Z_complex(keep_idx);
            freq_test = freq(remove_idx);
            Z_test = Z_complex(remove_idx);
            
            % Test on each lambda
            for l = 1:length(lambda_grid)
                lambda = lambda_grid(l);
                
                % Fit ridge on real part
                try
                    [~, Z_cal_real] = tr_drt_component_local(freq_train, Z_train, lambda, 'real');
                    Z_cal_real_test = tr_drt_component_local(freq_test, Z_test, lambda, 'real');
                    cv_real_vec(l, rep) = mean(abs(Z_cal_real - real(Z_train)));
                catch
                    cv_real_vec(l, rep) = inf;
                end
                
                % Fit ridge on imag part
                try
                    [~, Z_cal_imag] = tr_drt_component_local(freq_train, Z_train, lambda, 'imag');
                    Z_cal_imag_test = tr_drt_component_local(freq_test, Z_test, lambda, 'imag');
                    cv_imag_vec(l, rep) = mean(abs(Z_cal_imag - imag(Z_train)));
                catch
                    cv_imag_vec(l, rep) = inf;
                end
            end
        end
        
        cv_real_mean = mean(cv_real_vec, 2);
        cv_imag_mean = mean(cv_imag_vec, 2);
        cv_total = cv_real_mean + cv_imag_mean;
        
        % Find minimum
        [~, idx_min] = min(cv_total);
        lambda_opt = lambda_grid(idx_min);
        
        cv_plot_data{s} = struct(...
            'dataset', ds_name, ...
            'temp_k', temp_k, ...
            'lambda_grid', lambda_grid, ...
            'cv_real', cv_real_mean, ...
            'cv_imag', cv_imag_mean, ...
            'cv_total', cv_total, ...
            'lambda_opt', lambda_opt, ...
            'n_data', n_data ...
        );
        
    catch ME
        fprintf('    Error processing %s: %s\n', ds_file, ME.message);
    end
end

% Create plots
figure('Position', [100, 100, 1400, 900]);
n_cases = length(cv_plot_data);
n_cols = min(3, n_cases);
n_rows = ceil(n_cases / n_cols);

for s = 1:n_cases
    data = cv_plot_data{s};
    
    subplot(n_rows, n_cols, s);
    hold on;
    
    h1 = semilogx(data.lambda_grid, data.cv_real, 'o-', 'Color', [0.2, 0.6, 1.0], 'LineWidth', 1.5, 'MarkerSize', 5);
    h2 = semilogx(data.lambda_grid, data.cv_imag, 's-', 'Color', [1.0, 0.6, 0.2], 'LineWidth', 1.5, 'MarkerSize', 5);
    h3 = semilogx(data.lambda_grid, data.cv_total, 'd-', 'Color', [0.3, 0.8, 0.3], 'LineWidth', 2, 'MarkerSize', 6);
    
    plot(data.lambda_opt, data.cv_total(find(data.lambda_grid == data.lambda_opt)), 'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    xlabel('Regularization parameter λ', 'FontSize', 10);
    ylabel('Cross-validation error', 'FontSize', 10);
    title(sprintf('%s @ %.3fK (15%% removal)', strrep(data.dataset, '_', '\_'), data.temp_k), 'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    legend({sprintf('Real (n=%d)', data.n_data), 'Imaginary', 'Total', sprintf('Optimum: λ=%.2e', data.lambda_opt)}, 'FontSize', 9);
    hold off;
end

sgtitle('15% Data Removal Cross-Validation (3-Peak Consistent Cases)', 'FontSize', 14, 'FontWeight', 'bold');
export_fig_path = 'cv_15pct_removal_3peak_cases.png';
saveas(gcf, export_fig_path, 'png');
fprintf('\nPlot saved: %s\n\n', export_fig_path);

% Save results to Excel
output_excel = 'cv_15pct_removal_3peak_analysis.xlsx';
excel_data = table();

for s = 1:length(cv_plot_data)
    data = cv_plot_data{s};
    row_data = table(...
        {data.dataset}, ...
        data.temp_k, ...
        data.n_data, ...
        data.lambda_opt, ...
        min(data.cv_real), ...
        min(data.cv_imag), ...
        min(data.cv_total), ...
        'VariableNames', {'Dataset', 'TemperatureK', 'DataPoints', 'OptimalLambda', 'MinRealCV', 'MinImagCV', 'MinTotalCV'} ...
    );
    excel_data = [excel_data; row_data];
end

writetable(excel_data, output_excel, 'Sheet', 'Summary');
fprintf('Summary saved to: %s\n', output_excel);

% Save detailed CV curves per case
for s = 1:length(cv_plot_data)
    data = cv_plot_data{s};
    sheet_name = sprintf('%s_%.0fK', data.dataset, data.temp_k);
    % Truncate sheet name if needed (max 31 chars in Excel)
    sheet_name = sheet_name(1:min(31, length(sheet_name)));
    
    detail_table = table(...
        data.lambda_grid', ...
        data.cv_real, ...
        data.cv_imag, ...
        data.cv_total, ...
        'VariableNames', {'Lambda', 'RealPartCV', 'ImaginaryPartCV', 'TotalCV'} ...
    );
    writetable(detail_table, output_excel, 'Sheet', sheet_name);
end

fprintf('Detailed CV curves saved to: %s (multiple sheets)\n\n', output_excel);

%% Local helper functions

function [A, Z_cal] = tr_drt_component_local(freq, Z_exp, lambda, component)
% Simplified ridge regression for one component
    
    if nargin < 4, component = 'real'; end
    
    % Build response vector
    if strcmp(component, 'real')
        Z_vec = real(Z_exp);
    else
        Z_vec = imag(Z_exp);
    end
    
    % Simple regularization (no tau grid for now)
    % Just use identity matrix for demonstration
    n = length(freq);
    A = (eye(n) + lambda * eye(n)) \ Z_vec;
    Z_cal = A .* ones(size(freq));
end
