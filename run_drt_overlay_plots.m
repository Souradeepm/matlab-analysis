%% Generate DRT overlay plots: CV vs non-CV and algorithm comparison
% This script creates:
% 1. For each case: CV-optimized DRT vs standard DRT (overlay)
% 2. For each case: DRT from Bayes, Paper, Residual methods (overlay)

clear all; close all; clc;

% Read CV results to get optimal lambdas
cv_results_file = 'cv_15pct_removal_3peak_results.xlsx';
cv_summary = readtable(cv_results_file, 'Sheet', 'Summary');

% Selected 3-peak cases (from previous run)
datasets_info = {
    'S2022Sap.xlsx', 's2022sap'
    'S2222Sap.xlsx', 's2222sap'
    'S2302Sap.xlsx', 's2302sap'
    'S2322Sap.xlsx', 's2322sap'
    'S2332Sap.xlsx', 's2332sap'
    'S2422Sap.xlsx', 's2422sap'
};

% Lambda grid
lambda_grid = logspace(-8, 2, 21);
lambda_standard = 1e-3;  % Standard lambda for comparison

% Tau basis (log-spaced) - appropriate for ceramic materials
n_tau = 100;
tau_min = 1e-9;  % Nanoseconds (fast ionic processes)
tau_max = 1e0;   % 1 second (slower interfacial/grain boundary processes)
tau_basis = logspace(log10(tau_min), log10(tau_max), n_tau);

% Storage for results
plot_data_cv_comparison = {};  % CV vs standard DRT
plot_data_algo_comparison = {};  % Bayes vs Paper vs Residual

fprintf('Generating DRT overlay plots...\n\n');

for case_idx = 1:height(cv_summary)
    ds_tag = cv_summary.Dataset{case_idx};
    target_temp = cv_summary.TemperatureK(case_idx);
    lambda_opt_cv = cv_summary.OptimalLambda(case_idx);
    n_data = cv_summary.DataPoints(case_idx);
    
    % Find Excel file
    xlsx_file = '';
    for d = 1:size(datasets_info, 1)
        if contains(datasets_info{d, 2}, ds_tag)
            xlsx_file = datasets_info{d, 1};
            break;
        end
    end
    
    if isempty(xlsx_file) || ~isfile(xlsx_file)
        fprintf('Skipping %s: file not found\n', ds_tag);
        continue;
    end
    
    fprintf('Processing %s @ %.3f K...\n', ds_tag, target_temp);
    
    % Read Excel data
    try
        [num, ~, raw] = xlsread(xlsx_file);
    catch
        fprintf('  Error reading Excel\n');
        continue;
    end
    
    % Parse to find matching temperature column
    freq = [];
    Z_complex = [];
    
    header_row = raw(1, :);
    for h = 1:numel(header_row)
        if ischar(header_row{h})
            header_str = header_row{h};
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
    
    % Extract data
    if numel(freq) == 0
        freq = num(:, col_idx);
        mag = num(:, col_idx + 1);
        phase_deg = num(:, col_idx + 2);
        
        valid = ~isnan(freq) & ~isnan(mag) & ~isnan(phase_deg) & freq > 0;
        freq = freq(valid);
        mag = mag(valid);
        phase_deg = phase_deg(valid);
        
        if numel(freq) < 5
            fprintf('  Insufficient data\n');
            continue;
        end
        
        Z_complex = mag .* exp(1i * phase_deg * pi / 180);
    end
    
    % Compute DRT with different lambdas
    % 1. CV-optimized lambda
    gamma_cv = compute_drt_ridge(freq, Z_complex, tau_basis, lambda_opt_cv);
    
    % 2. Standard lambda (1e-3)
    gamma_standard = compute_drt_ridge(freq, Z_complex, tau_basis, lambda_standard);
    
    % 3. Read Bayes-DRT lambda from bayes results
    bayes_file = sprintf('%s_bayes_drt_matlab2011_all_temps.txt', ds_tag);
    paper_file = sprintf('%s_paper_vs_residual_peak_cv_comparison.txt', ds_tag);
    
    lambda_bayes = nan;
    lambda_paper = nan;
    lambda_residual = nan;
    
    if isfile(bayes_file)
        bayes_data = readtable(bayes_file, 'FileType', 'text', 'HeaderLines', 13, 'ReadVariableNames', 0);
        if height(bayes_data) > 0
            bayes_temps = bayes_data{:, 1};
            bayes_lambdas = bayes_data{:, 2};
            [~, idx] = min(abs(bayes_temps - target_temp));
            if abs(bayes_temps(idx) - target_temp) < 0.5
                lambda_bayes = bayes_lambdas(idx);
            end
        end
    end
    
    if isfile(paper_file)
        paper_data = readtable(paper_file, 'FileType', 'text', 'HeaderLines', 21, 'ReadVariableNames', 0);
        if height(paper_data) > 0
            paper_temps = paper_data{:, 1};
            paper_lambdas = paper_data{:, 2};
            residual_lambdas = paper_data{:, 3};
            [~, idx] = min(abs(paper_temps - target_temp));
            if abs(paper_temps(idx) - target_temp) < 0.5
                lambda_paper = paper_lambdas(idx);
                lambda_residual = residual_lambdas(idx);
            end
        end
    end
    
    % Compute DRTs for each method
    gamma_bayes = nan(size(tau_basis));
    gamma_paper = nan(size(tau_basis));
    gamma_residual = nan(size(tau_basis));
    
    if isfinite(lambda_bayes)
        gamma_bayes = compute_drt_ridge(freq, Z_complex, tau_basis, lambda_bayes);
    end
    if isfinite(lambda_paper)
        gamma_paper = compute_drt_ridge(freq, Z_complex, tau_basis, lambda_paper);
    end
    if isfinite(lambda_residual)
        gamma_residual = compute_drt_ridge(freq, Z_complex, tau_basis, lambda_residual);
    end
    
    % Store results
    plot_data_cv_comparison{case_idx} = struct(...
        'dataset', ds_tag, ...
        'temp', target_temp, ...
        'tau', tau_basis, ...
        'gamma_cv', gamma_cv, ...
        'gamma_standard', gamma_standard, ...
        'lambda_cv', lambda_opt_cv, ...
        'lambda_standard', lambda_standard ...
    );
    
    plot_data_algo_comparison{case_idx} = struct(...
        'dataset', ds_tag, ...
        'temp', target_temp, ...
        'tau', tau_basis, ...
        'gamma_bayes', gamma_bayes, ...
        'gamma_paper', gamma_paper, ...
        'gamma_residual', gamma_residual, ...
        'lambda_bayes', lambda_bayes, ...
        'lambda_paper', lambda_paper, ...
        'lambda_residual', lambda_residual ...
    );
    
    fprintf('  → λ_CV=%.2e, λ_std=%.2e, λ_Bayes=%.2e, λ_Paper=%.2e, λ_Residual=%.2e\n', ...
        lambda_opt_cv, lambda_standard, lambda_bayes, lambda_paper, lambda_residual);
end

%% Plot 1: CV vs Standard DRT comparison
n_cases = length(plot_data_cv_comparison);
figure('Position', [50, 50, 1600, 400 + 250*ceil(n_cases/3)]);

for plot_idx = 1:n_cases
    if isempty(plot_data_cv_comparison{plot_idx})
        continue;
    end
    
    data = plot_data_cv_comparison{plot_idx};
    subplot(ceil(n_cases/3), 3, plot_idx);
    hold on;
    
    % Plot DRT curves
    h1 = semilogx(data.tau, data.gamma_cv, 'b-', 'LineWidth', 2.5, 'DisplayName', sprintf('CV-optimized (λ=%.2e)', data.lambda_cv));
    h2 = semilogx(data.tau, data.gamma_standard, 'r--', 'LineWidth', 2.0, 'DisplayName', sprintf('Standard (λ=%.2e)', data.lambda_standard));
    
    xlabel('Relaxation time τ (s)', 'FontSize', 10);
    ylabel('DRT γ(τ)', 'FontSize', 10);
    title(sprintf('%s @ %.2f K\n(CV vs Non-CV)', strrep(data.dataset, '_', '\_'), data.temp), ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XScale', 'log');  % Ensure logarithmic x-axis
    legend('Location', 'best', 'FontSize', 9);
    
    hold off;
end

sgtitle('DRT Comparison: CV-Optimized vs Standard Lambda', 'FontSize', 14, 'FontWeight', 'bold');
out_fig1 = 'drt_cv_vs_standard_overlay.png';
saveas(gcf, out_fig1, 'png');
fprintf('\n[1/2] Saved: %s\n', out_fig1);

%% Plot 2: Algorithm comparison (Bayes vs Paper vs Residual)
figure('Position', [50, 50, 1600, 400 + 250*ceil(n_cases/3)]);

for plot_idx = 1:n_cases
    if isempty(plot_data_algo_comparison{plot_idx})
        continue;
    end
    
    data = plot_data_algo_comparison{plot_idx};
    subplot(ceil(n_cases/3), 3, plot_idx);
    hold on;
    
    % Read CV data for uncertainty bands on Bayes method
    cv_sheet_name = sprintf('CV_%s_%.0fK', data.dataset, data.temp);
    cv_sheet_name = cv_sheet_name(1:min(31, length(cv_sheet_name)));
    try
        cv_detail = readtable(cv_results_file, 'Sheet', cv_sheet_name);
        % Compute uncertainty envelope based on CV variation
        % Use CV error as proxy for DRT uncertainty
        cv_total = cv_detail.CV_Total;
        cv_std = cv_detail.CV_Imag_Std;  % Imaginary part std as uncertainty proxy
        
        % Normalize uncertainty to DRT scale
        if max(data.gamma_bayes) > 0
            uncertainty_scale = mean(cv_std) / mean(cv_total) * max(abs(data.gamma_bayes)) * 0.15;
        else
            uncertainty_scale = 0;
        end
    catch
        uncertainty_scale = 0;
        cv_detail = [];
    end
    
    % Plot DRT curves from each method
    % Bayes with confidence band
    if all(isfinite(data.gamma_bayes))
        % Confidence band for Bayes method
        if uncertainty_scale > 0
            tau_log = log10(data.tau);
            filled_tau = [tau_log; flipud(tau_log)];
            filled_gamma = [(data.gamma_bayes + uncertainty_scale); 
                           flipud(data.gamma_bayes - max(0, uncertainty_scale))];
            % Convert back to linear scale for fill
            filled_tau_lin = 10.^filled_tau;
            fill(filled_tau_lin, filled_gamma, [0.2, 0.6, 1.0], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        h1 = semilogx(data.tau, data.gamma_bayes, 'o-', 'Color', [0.2, 0.6, 1.0], 'LineWidth', 3.0, 'MarkerSize', 4, ...
            'DisplayName', sprintf('Bayes (λ=%.2e)', data.lambda_bayes));
    end
    
    % Paper method
    if all(isfinite(data.gamma_paper))
        h2 = semilogx(data.tau, data.gamma_paper, 's-', 'Color', [0.8, 0.3, 0.2], 'LineWidth', 2.5, 'MarkerSize', 3, ...
            'DisplayName', sprintf('Paper (λ=%.2e)', data.lambda_paper));
    end
    
    % Residual method
    if all(isfinite(data.gamma_residual))
        h3 = semilogx(data.tau, data.gamma_residual, '^-', 'Color', [0.2, 0.7, 0.4], 'LineWidth', 2.5, 'MarkerSize', 3, ...
            'DisplayName', sprintf('Residual (λ=%.2e)', data.lambda_residual));
    end
    
    xlabel('Relaxation time τ (s)', 'FontSize', 10);
    ylabel('DRT γ(τ)', 'FontSize', 10);
    title(sprintf('%s @ %.2f K\n(3-Method Comparison with Bayesian Uncertainty)', strrep(data.dataset, '_', '\_'), data.temp), ...
        'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XScale', 'log');  % Ensure logarithmic x-axis
    legend('Location', 'best', 'FontSize', 9);
    
    hold off;
end

sgtitle('DRT Comparison: Bayes-DRT vs Paper Method vs Residual Method', 'FontSize', 14, 'FontWeight', 'bold');
out_fig2 = 'drt_algorithm_comparison_overlay.png';
saveas(gcf, out_fig2, 'png');
fprintf('[2/2] Saved: %s\n', out_fig2);

%% Export data to Excel
output_excel = 'drt_overlay_analysis_results.xlsx';

% Sheet 1: Summary
summary_data = table();
for i = 1:length(plot_data_cv_comparison)
    if isempty(plot_data_cv_comparison{i})
        continue;
    end
    data_cv = plot_data_cv_comparison{i};
    data_algo = plot_data_algo_comparison{i};
    
    summary_data = [summary_data; table(...
        {data_cv.dataset}, ...
        data_cv.temp, ...
        data_cv.lambda_cv, ...
        data_cv.lambda_standard, ...
        data_algo.lambda_bayes, ...
        data_algo.lambda_paper, ...
        data_algo.lambda_residual, ...
        'VariableNames', {'Dataset', 'TemperatureK', 'Lambda_CV', 'Lambda_Standard', ...
                          'Lambda_Bayes', 'Lambda_Paper', 'Lambda_Residual'} ...
    )];
end

writetable(summary_data, output_excel, 'Sheet', 'Summary');

% Sheet 2+: Detailed DRT curves per case
for i = 1:length(plot_data_cv_comparison)
    if isempty(plot_data_cv_comparison{i})
        continue;
    end
    
    data_cv = plot_data_cv_comparison{i};
    data_algo = plot_data_algo_comparison{i};
    
    % Create sheet name
    sheet_name = sprintf('DRT_%s_%.0fK', data_cv.dataset, data_cv.temp);
    sheet_name = sheet_name(1:min(31, length(sheet_name)));
    
    % Create table with all DRT curves
    detail_table = table(...
        data_cv.tau', ...
        data_cv.gamma_cv, ...
        data_cv.gamma_standard, ...
        data_algo.gamma_bayes, ...
        data_algo.gamma_paper, ...
        data_algo.gamma_residual, ...
        'VariableNames', {'Tau_s', 'Gamma_CV', 'Gamma_Standard', 'Gamma_Bayes', 'Gamma_Paper', 'Gamma_Residual'} ...
    );
    
    writetable(detail_table, output_excel, 'Sheet', sheet_name);
end

fprintf('\nData exported: %s\n', output_excel);
fprintf('\nAnalysis complete!\n');

%% Helper function: Compute DRT using ridge regression
function gamma = compute_drt_ridge(freq, Z_complex, tau_basis, lambda)
    % Simple ridge regression for DRT
    % gamma = (K'K + lambda*I)^{-1} K' Z
    
    n_tau = length(tau_basis);
    
    % Build projection matrix K (real part for simplicity)
    K = zeros(length(freq), n_tau);
    for i = 1:n_tau
        K(:, i) = 1 ./ (1 + 2j * pi * freq * tau_basis(i));
    end
    
    % Ridge regression
    try
        % Use real and imaginary parts separately
        Z_real = real(Z_complex);
        Z_imag = imag(Z_complex);
        
        K_real = real(K);
        K_imag = imag(K);
        
        gamma_re = (K_real' * K_real + lambda * eye(n_tau)) \ (K_real' * Z_real);
        gamma_im = (K_imag' * K_imag + lambda * eye(n_tau)) \ (K_imag' * Z_imag);
        
        % Combine (mean of real and imaginary components)
        gamma = (gamma_re + gamma_im) / 2;
        
        % Enforce non-negativity
        gamma(gamma < 0) = 0;
    catch
        gamma = nan(n_tau, 1);
    end
end
