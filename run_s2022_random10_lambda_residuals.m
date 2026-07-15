function run_s2022_random10_lambda_residuals()
% run_s2022_random10_lambda_residuals
% MATLAB 2011-compatible lambda study on 10 randomly selected S2022 temperatures.
% Compares residual-only selection against ChemElectroChem-style multi-criterion
% selection using RID, Re/Im cross-validation, and resampling variance.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_random10_lambda_residuals.xlsx');
out_summary_txt = fullfile(repo_root, 's2022_random10_lambda_residuals_summary.txt');
out_detail_csv = fullfile(repo_root, 's2022_random10_lambda_residuals_detail.csv');
out_aggregate_csv = fullfile(repo_root, 's2022_random10_lambda_residuals_aggregate.csv');
out_png = fullfile(repo_root, 's2022_random10_lambda_residuals.png');
out_drt_png = fullfile(repo_root, 's2022_random10_lambda_drt_compare.png');

rng_seed = 20260710;
n_select = 10;
lambda_values = logspace(-4, -1, 11).';  % Range [1e-4, 1e-1] constrained

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text] = read_excel_robust_local(input_path);
temperature = parse_temperature_labels_local(text, size(data, 2));
if numel(temperature) < n_select
    error('Requested %d temperatures, but only %d are available.', n_select, numel(temperature));
end

rand('twister', rng_seed);
randn('state', rng_seed + 1);
selected_idx = randperm(numel(temperature));
selected_idx = selected_idx(1:n_select);
[selected_temperature, order] = sort(temperature(selected_idx));
selected_idx = selected_idx(order);

n_lambda = numel(lambda_values);
mean_abs_residual = zeros(n_select, n_lambda);
rmse_residual = zeros(n_select, n_lambda);
rid_metric = zeros(n_select, n_lambda);
cv_real_metric = zeros(n_select, n_lambda);
cv_imag_metric = zeros(n_select, n_lambda);
cv_total_metric = zeros(n_select, n_lambda);
resample_var_metric = zeros(n_select, n_lambda);
score_metric = zeros(n_select, n_lambda);
R_inf_all = zeros(n_select, n_lambda);
n_points = zeros(n_select, 1);
residual_best_idx = zeros(n_select, 1);
paper_best_idx = zeros(n_select, 1);

fprintf('S2022 random 10-temperature lambda study (MATLAB 2011 compatible)\n');
fprintf('Random seed: %d\n', rng_seed);
fprintf('Selected temperatures (K):\n');
fprintf('  %.2f', selected_temperature(1));
for k = 2:n_select
    fprintf(', %.2f', selected_temperature(k));
end
fprintf('\n\n');

for t = 1:n_select
    src_cols = (3 * selected_idx(t) - 2):(3 * selected_idx(t));
    sel = data(:, src_cols);
    freq_vec = sel(:,1);
    mag = sel(:,2);
    phase_deg = sel(:,3);

    valid = isfinite(freq_vec) & isfinite(mag) & isfinite(phase_deg) & (freq_vec > 0);
    freq_vec = freq_vec(valid);
    mag = mag(valid);
    phase_deg = phase_deg(valid);

    Z_exp = mag .* exp(1i * phase_deg * pi / 180);
    [freq_vec, sort_order] = sort(freq_vec);
    Z_exp = Z_exp(sort_order);
    n_points(t) = numel(freq_vec);

    fprintf('Temperature %.2f K (%d points)\n', selected_temperature(t), n_points(t));
    for i = 1:n_lambda
        lam = lambda_values(i);

        [gamma_i, R_inf_i] = TR_DRT_local(freq_vec, Z_exp, lam);
        Z_fit_i = calculate_EIS_local(freq_vec, gamma_i, R_inf_i);
        abs_resid = abs(Z_exp - Z_fit_i);

        [gamma_re_only, R_re_only] = TR_DRT_component_local(freq_vec, Z_exp, lam, 'real');
        [gamma_im_only, ~] = TR_DRT_component_local(freq_vec, Z_exp, lam, 'imag');
        Z_from_re_only = calculate_EIS_local(freq_vec, gamma_re_only, R_re_only);
        R_from_im = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, gamma_im_only);
        Z_from_im_only = calculate_EIS_local(freq_vec, gamma_im_only, R_from_im);

        mean_abs_residual(t, i) = mean(abs_resid);
        rmse_residual(t, i) = sqrt(mean(abs_resid .^ 2));
        rid_metric(t, i) = mean((gamma_re_only - gamma_im_only).^2);
        cv_imag_metric(t, i) = mean((imag(Z_exp) - imag(Z_from_re_only)).^2);
        cv_real_metric(t, i) = mean((real(Z_exp) - real(Z_from_im_only)).^2);
        cv_total_metric(t, i) = cv_real_metric(t, i) + cv_imag_metric(t, i);
        resample_var_metric(t, i) = estimate_resample_variance_local(freq_vec, Z_exp, lam, 8);
        R_inf_all(t, i) = R_inf_i;

        fprintf('  lambda=%8.1e | resid=%10.6f | RID=%10.3e | CV=%10.3e | Var=%10.3e\n', ...
            lam, mean_abs_residual(t, i), rid_metric(t, i), cv_total_metric(t, i), resample_var_metric(t, i));
    end

    score_metric(t, :) = (normalize_metric_local(rid_metric(t, :).') + ...
        normalize_metric_local(cv_total_metric(t, :).') + ...
        normalize_metric_local(resample_var_metric(t, :).')).';
    [~, residual_best_idx(t)] = min(mean_abs_residual(t, :));
    paper_best_idx(t) = select_lambda_idx_local(rid_metric(t, :).', cv_total_metric(t, :).', resample_var_metric(t, :).');

    fprintf('  residual-only best : %.1e (mean abs residual %.6f)\n', ...
        lambda_values(residual_best_idx(t)), mean_abs_residual(t, residual_best_idx(t)));
    fprintf('  paper-style best   : %.1e (RID %.3e, CV %.3e, Var %.3e)\n\n', ...
        lambda_values(paper_best_idx(t)), rid_metric(t, paper_best_idx(t)), ...
        cv_total_metric(t, paper_best_idx(t)), resample_var_metric(t, paper_best_idx(t)));
end

aggregate_mean_residual = mean(mean_abs_residual, 1).';
aggregate_std_residual = std(mean_abs_residual, 0, 1).';
aggregate_rid = mean(rid_metric, 1).';
aggregate_cv = mean(cv_total_metric, 1).';
aggregate_var = mean(resample_var_metric, 1).';
aggregate_score = normalize_metric_local(aggregate_rid) + normalize_metric_local(aggregate_cv) + normalize_metric_local(aggregate_var);

[~, aggregate_residual_idx] = min(aggregate_mean_residual);
aggregate_paper_idx = select_lambda_idx_local(aggregate_rid, aggregate_cv, aggregate_var);

write_summary_local(out_summary_txt, rng_seed, selected_temperature, lambda_values, ...
    mean_abs_residual, rid_metric, cv_total_metric, resample_var_metric, ...
    residual_best_idx, paper_best_idx, aggregate_mean_residual, aggregate_std_residual, ...
    aggregate_rid, aggregate_cv, aggregate_var, aggregate_score, ...
    aggregate_residual_idx, aggregate_paper_idx);

write_detail_csv_local(out_detail_csv, selected_temperature, lambda_values, n_points, R_inf_all, ...
    mean_abs_residual, rmse_residual, rid_metric, cv_real_metric, cv_imag_metric, ...
    cv_total_metric, resample_var_metric, score_metric, residual_best_idx, paper_best_idx);

write_aggregate_csv_local(out_aggregate_csv, lambda_values, aggregate_mean_residual, ...
    aggregate_std_residual, aggregate_rid, aggregate_cv, aggregate_var, aggregate_score, ...
    aggregate_residual_idx, aggregate_paper_idx);

write_excel_best_effort_local(out_xlsx, selected_temperature, lambda_values, n_points, R_inf_all, ...
    mean_abs_residual, rmse_residual, rid_metric, cv_real_metric, cv_imag_metric, ...
    cv_total_metric, resample_var_metric, score_metric, residual_best_idx, paper_best_idx, ...
    aggregate_mean_residual, aggregate_std_residual, aggregate_rid, aggregate_cv, ...
    aggregate_var, aggregate_score, aggregate_residual_idx, aggregate_paper_idx, rng_seed);

plot_results_local(out_png, selected_temperature, lambda_values, mean_abs_residual, ...
    aggregate_mean_residual, aggregate_std_residual, aggregate_rid, aggregate_cv, ...
    aggregate_var, aggregate_score, aggregate_residual_idx, aggregate_paper_idx);

plot_drt_comparison_local(out_drt_png, data, selected_idx, selected_temperature, lambda_values, ...
    residual_best_idx, paper_best_idx);

fprintf('Wrote summary     : %s\n', out_summary_txt);
fprintf('Wrote detail CSV  : %s\n', out_detail_csv);
fprintf('Wrote aggregate CSV: %s\n', out_aggregate_csv);
fprintf('Wrote plot        : %s\n', out_png);
fprintf('Wrote DRT plot    : %s\n', out_drt_png);
fprintf('Aggregate residual-only best lambda: %.1e\n', lambda_values(aggregate_residual_idx));
fprintf('Aggregate paper-style best lambda : %.1e\n', lambda_values(aggregate_paper_idx));
end

function [data, text] = read_excel_robust_local(file_path)
try
    [data, text] = xlsread(file_path);
catch ME
    error('Failed to read Excel file %s with xlsread: %s', file_path, ME.message);
end
end

function temperature = parse_temperature_labels_local(text, n_cols)
if iscell(text)
    header = text(1, :);
else
    header = cellstr(char(text(1, :)));
end

temperature = [];
for k = 1:numel(header)
    tok = header{k};
    if ischar(tok)
        tok(tok == 'K') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end + 1, 1) = val; %#ok<AGROW>
        end
    end
end

n_sets = floor(n_cols / 3);
temperature = temperature(1:min(numel(temperature), n_sets));
end

function idx_best = select_lambda_idx_local(rid_vec, cv_vec, var_vec)
% Paper Method v2: Score-based selection
% Combines three metrics equally (normalized):
% - RID: Re/Im divergence (fit consistency)
% - CV: Cross-validation error (stability)
% - Var: Bootstrap resampling variance (robustness)
%
% Lambda is selected that minimizes: score = norm(RID) + norm(CV) + norm(Var)
% This is the published algorithm from ChemElectroChem 2019.
%
% Previous version had a RID-constrained path that could be overwritten by
% score-based selection, causing ambiguity. This version uses score-based
% consistently as the primary selection criterion.

norm_rid = normalize_metric_local(rid_vec);
norm_cv = normalize_metric_local(cv_vec);
norm_var = normalize_metric_local(var_vec);

score = norm_rid + norm_cv + norm_var;
[~, idx_best] = min(score);
end

function write_summary_local(out_path, rng_seed, selected_temperature, lambda_values, ...
    mean_abs_residual, rid_metric, cv_total_metric, resample_var_metric, ...
    residual_best_idx, paper_best_idx, aggregate_mean_residual, aggregate_std_residual, ...
    aggregate_rid, aggregate_cv, aggregate_var, aggregate_score, ...
    aggregate_residual_idx, aggregate_paper_idx)

fid = fopen(out_path, 'w');
if fid < 0
    error('Failed to open summary file for writing: %s', out_path);
end

fprintf(fid, 'S2022 random 10-temperature lambda study (MATLAB 2011 compatible)\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Random seed: %d\n', rng_seed);
fprintf(fid, 'Selected temperatures (K): ');
fprintf(fid, '%.6f ', selected_temperature);
fprintf(fid, '\n\n');
fprintf(fid, 'Aggregate residual-only best lambda: %.8e\n', lambda_values(aggregate_residual_idx));
fprintf(fid, 'Aggregate paper-style best lambda : %.8e\n\n', lambda_values(aggregate_paper_idx));

fprintf(fid, 'Per-temperature selections:\n');
for t = 1:numel(selected_temperature)
    fprintf(fid, '  %.6f K -> residual best %.8e (resid %.8f), paper best %.8e (RID %.8e, CV %.8e, Var %.8e)\n', ...
        selected_temperature(t), lambda_values(residual_best_idx(t)), mean_abs_residual(t, residual_best_idx(t)), ...
        lambda_values(paper_best_idx(t)), rid_metric(t, paper_best_idx(t)), ...
        cv_total_metric(t, paper_best_idx(t)), resample_var_metric(t, paper_best_idx(t)));
end

fprintf(fid, '\nAggregate trend:\n');
fprintf(fid, 'lambda,mean_abs_residual,std_abs_residual,mean_RID,mean_CV,mean_Var,score\n');
for i = 1:numel(lambda_values)
    fprintf(fid, '%.8e,%.8f,%.8f,%.8e,%.8e,%.8e,%.8f\n', ...
        lambda_values(i), aggregate_mean_residual(i), aggregate_std_residual(i), ...
        aggregate_rid(i), aggregate_cv(i), aggregate_var(i), aggregate_score(i));
end

fclose(fid);
end

function write_detail_csv_local(out_path, selected_temperature, lambda_values, n_points, R_inf_all, ...
    mean_abs_residual, rmse_residual, rid_metric, cv_real_metric, cv_imag_metric, ...
    cv_total_metric, resample_var_metric, score_metric, residual_best_idx, paper_best_idx)

fid = fopen(out_path, 'w');
if fid < 0
    error('Failed to open detail CSV for writing: %s', out_path);
end

fprintf(fid, 'TemperatureK,PointCount,Lambda,R_inf,MeanAbsResidual,RMSE,RID,CVReal,CVImag,CVTotal,ResampleVariance,Score,IsResidualBest,IsPaperBest\n');
for t = 1:numel(selected_temperature)
    for i = 1:numel(lambda_values)
        fprintf(fid, '%.6f,%d,%.8e,%.8f,%.8f,%.8f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8f,%d,%d\n', ...
            selected_temperature(t), n_points(t), lambda_values(i), R_inf_all(t, i), ...
            mean_abs_residual(t, i), rmse_residual(t, i), rid_metric(t, i), ...
            cv_real_metric(t, i), cv_imag_metric(t, i), cv_total_metric(t, i), ...
            resample_var_metric(t, i), score_metric(t, i), i == residual_best_idx(t), i == paper_best_idx(t));
    end
end

fclose(fid);
end

function write_aggregate_csv_local(out_path, lambda_values, aggregate_mean_residual, ...
    aggregate_std_residual, aggregate_rid, aggregate_cv, aggregate_var, aggregate_score, ...
    aggregate_residual_idx, aggregate_paper_idx)

fid = fopen(out_path, 'w');
if fid < 0
    error('Failed to open aggregate CSV for writing: %s', out_path);
end

fprintf(fid, 'Lambda,MeanAbsResidualAcrossTemps,StdAbsResidualAcrossTemps,MeanRID,MeanCV,MeanVar,Score,IsAggregateResidualBest,IsAggregatePaperBest\n');
for i = 1:numel(lambda_values)
    fprintf(fid, '%.8e,%.8f,%.8f,%.8e,%.8e,%.8e,%.8f,%d,%d\n', ...
        lambda_values(i), aggregate_mean_residual(i), aggregate_std_residual(i), ...
        aggregate_rid(i), aggregate_cv(i), aggregate_var(i), aggregate_score(i), ...
        i == aggregate_residual_idx, i == aggregate_paper_idx);
end

fclose(fid);
end

function write_excel_best_effort_local(out_xlsx, selected_temperature, lambda_values, n_points, R_inf_all, ...
    mean_abs_residual, rmse_residual, rid_metric, cv_real_metric, cv_imag_metric, ...
    cv_total_metric, resample_var_metric, score_metric, residual_best_idx, paper_best_idx, ...
    aggregate_mean_residual, aggregate_std_residual, aggregate_rid, aggregate_cv, ...
    aggregate_var, aggregate_score, aggregate_residual_idx, aggregate_paper_idx, rng_seed)

try
    meta_header = {'random_seed', 'n_temperatures', 'aggregate_residual_best_lambda', 'aggregate_paper_best_lambda'};
    meta_values = [rng_seed, numel(selected_temperature), lambda_values(aggregate_residual_idx), lambda_values(aggregate_paper_idx)];
    xlswrite(out_xlsx, meta_header, 'meta', 'A1');
    xlswrite(out_xlsx, meta_values, 'meta', 'A2');

    temp_header = {'TemperatureK', 'ResidualBestLambda', 'ResidualBestMeanAbsResidual', 'PaperBestLambda', 'PaperBestRID', 'PaperBestCV', 'PaperBestVar'};
    temp_values = zeros(numel(selected_temperature), 7);
    for t = 1:numel(selected_temperature)
        temp_values(t, :) = [selected_temperature(t), lambda_values(residual_best_idx(t)), ...
            mean_abs_residual(t, residual_best_idx(t)), lambda_values(paper_best_idx(t)), ...
            rid_metric(t, paper_best_idx(t)), cv_total_metric(t, paper_best_idx(t)), ...
            resample_var_metric(t, paper_best_idx(t))];
    end
    xlswrite(out_xlsx, temp_header, 'per_temperature_selection', 'A1');
    xlswrite(out_xlsx, temp_values, 'per_temperature_selection', 'A2');

    detail_header = {'TemperatureK', 'PointCount', 'Lambda', 'R_inf', 'MeanAbsResidual', 'RMSE', 'RID', 'CVReal', 'CVImag', 'CVTotal', 'ResampleVariance', 'Score', 'IsResidualBest', 'IsPaperBest'};
    detail_values = zeros(numel(selected_temperature) * numel(lambda_values), 14);
    row = 1;
    for t = 1:numel(selected_temperature)
        for i = 1:numel(lambda_values)
            detail_values(row, :) = [selected_temperature(t), n_points(t), lambda_values(i), R_inf_all(t, i), ...
                mean_abs_residual(t, i), rmse_residual(t, i), rid_metric(t, i), cv_real_metric(t, i), ...
                cv_imag_metric(t, i), cv_total_metric(t, i), resample_var_metric(t, i), score_metric(t, i), ...
                i == residual_best_idx(t), i == paper_best_idx(t)];
            row = row + 1;
        end
    end
    xlswrite(out_xlsx, detail_header, 'detail', 'A1');
    xlswrite(out_xlsx, detail_values, 'detail', 'A2');

    aggregate_header = {'Lambda', 'MeanAbsResidualAcrossTemps', 'StdAbsResidualAcrossTemps', 'MeanRID', 'MeanCV', 'MeanVar', 'Score', 'IsAggregateResidualBest', 'IsAggregatePaperBest'};
    aggregate_values = [lambda_values, aggregate_mean_residual, aggregate_std_residual, aggregate_rid, aggregate_cv, aggregate_var, aggregate_score, ...
        double((1:numel(lambda_values)).' == aggregate_residual_idx), double((1:numel(lambda_values)).' == aggregate_paper_idx)];
    xlswrite(out_xlsx, aggregate_header, 'aggregate', 'A1');
    xlswrite(out_xlsx, aggregate_values, 'aggregate', 'A2');
catch ME
    fprintf('Excel export skipped: %s\n', ME.message);
end
end

function plot_results_local(out_png, selected_temperature, lambda_values, mean_abs_residual, ...
    aggregate_mean_residual, aggregate_std_residual, aggregate_rid, aggregate_cv, ...
    aggregate_var, aggregate_score, aggregate_residual_idx, aggregate_paper_idx)

fig = figure('Color', 'w', 'Name', 'S2022 random 10 temperature lambda study');

subplot(3,1,1);
cm = lines(numel(selected_temperature));
for t = 1:numel(selected_temperature)
    semilogx(lambda_values, mean_abs_residual(t, :), 'o-', 'Color', cm(t, :), 'LineWidth', 1.1, 'MarkerSize', 4); hold on;
end
grid on;
xlabel('lambda');
ylabel('Mean abs residual');
title('Residual trend for selected temperatures');
legend_text = cell(numel(selected_temperature), 1);
for t = 1:numel(selected_temperature)
    legend_text{t} = sprintf('%.2f K', selected_temperature(t));
end
legend(legend_text, 'Location', 'EastOutside');

subplot(3,1,2);
semilogx(lambda_values, aggregate_mean_residual, 'ko-', 'LineWidth', 1.5, 'MarkerSize', 5); hold on;
semilogx(lambda_values, aggregate_mean_residual + aggregate_std_residual, 'k--', 'LineWidth', 1.0);
semilogx(lambda_values, aggregate_mean_residual - aggregate_std_residual, 'k--', 'LineWidth', 1.0);
yl = ylim;
plot([lambda_values(aggregate_residual_idx) lambda_values(aggregate_residual_idx)], yl, 'r--', 'LineWidth', 1.0);
plot([lambda_values(aggregate_paper_idx) lambda_values(aggregate_paper_idx)], yl, 'b--', 'LineWidth', 1.0);
ylim(yl);
grid on;
xlabel('lambda');
ylabel('Aggregate mean abs residual');
title('Aggregate residual trend');
legend('mean residual', 'mean + std', 'mean - std', 'residual best', 'paper best', 'Location', 'Best');

subplot(3,1,3);
semilogx(lambda_values, normalize_metric_local(aggregate_rid), 'r-o', 'LineWidth', 1.2, 'MarkerSize', 4); hold on;
semilogx(lambda_values, normalize_metric_local(aggregate_cv), 'g-s', 'LineWidth', 1.2, 'MarkerSize', 4);
semilogx(lambda_values, normalize_metric_local(aggregate_var), 'b-^', 'LineWidth', 1.2, 'MarkerSize', 4);
semilogx(lambda_values, aggregate_score, 'k-d', 'LineWidth', 1.3, 'MarkerSize', 4);
grid on;
xlabel('lambda');
ylabel('Normalized metric / score');
title('Paper-style selection metrics across selected temperatures');
legend('RID', 'CV', 'Variance', 'Score', 'Location', 'Best');

saveas(fig, out_png);
close(fig);
end

function plot_drt_comparison_local(out_png, data, selected_idx, selected_temperature, lambda_values, residual_best_idx, paper_best_idx)
fig = figure('Color', 'w', 'Name', 'S2022 random 10 DRT comparison');

for t = 1:numel(selected_temperature)
    src_cols = (3 * selected_idx(t) - 2):(3 * selected_idx(t));
    sel = data(:, src_cols);
    freq_vec = sel(:,1);
    mag = sel(:,2);
    phase_deg = sel(:,3);

    valid = isfinite(freq_vec) & isfinite(mag) & isfinite(phase_deg) & (freq_vec > 0);
    freq_vec = freq_vec(valid);
    mag = mag(valid);
    phase_deg = phase_deg(valid);

    Z_exp = mag .* exp(1i * phase_deg * pi / 180);
    [freq_vec, sort_order] = sort(freq_vec);
    Z_exp = Z_exp(sort_order);

    lambda_res = lambda_values(residual_best_idx(t));
    lambda_paper = lambda_values(paper_best_idx(t));

    [gamma_res, ~] = TR_DRT_local(freq_vec, Z_exp, lambda_res);
    [gamma_paper, ~] = TR_DRT_local(freq_vec, Z_exp, lambda_paper);
    tau = 1 ./ (2 * pi * freq_vec);

    subplot(5, 2, t);
    semilogx(tau, gamma_res, 'r-', 'LineWidth', 1.0); hold on;
    semilogx(tau, gamma_paper, 'b-', 'LineWidth', 1.0);
    set(gca, 'XDir', 'reverse');
    grid on;
    title(sprintf('T=%.2f K | res %.1e | paper %.1e', selected_temperature(t), lambda_res, lambda_paper), 'FontSize', 8);
    if mod(t, 2) == 1
        ylabel('gamma(tau)');
    end
    if t > 8
        xlabel('tau (s)');
    end
    if t == 1
        legend('Residual best', 'Paper best', 'Location', 'Best');
    end
end

saveas(fig, out_png);
close(fig);
end

function [gamma, R_inf] = TR_DRT_local(freq_vec, Z_exp, el)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];

ub = max(abs(Z_exp)) * ones(n + 1, 1);
if exist('lsqlin', 'file') == 2
    lb = zeros(n + 1, 1);
    options = optimset('Display', 'off');
    x = lsqlin(M, b, [], [], [], [], lb, ub, [], options);
else
    x = lsqnonneg(M, b);
    x = min(x, ub);
end

gamma = x(1:end-1);
R_inf = x(end);
end

function [gamma, R_inf] = TR_DRT_component_local(freq_vec, Z_exp, el, component)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

switch lower(component)
    case 'real'
        M = [A_re, ones(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
        b = [Z_re; zeros(n,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(el / 2) * eye(n), zeros(n,1)];
        b = [Z_im; zeros(n,1)];
    otherwise
        error('Unknown component mode: %s', component);
end

ub = max(abs(Z_exp)) * ones(n + 1, 1);
if exist('lsqlin', 'file') == 2
    lb = zeros(n + 1, 1);
    options = optimset('Display', 'off');
    x = lsqlin(M, b, [], [], [], [], lb, ub, [], options);
else
    x = lsqnonneg(M, b);
    x = min(x, ub);
end

gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calculate_EIS_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function R_inf = estimate_Rinf_for_gamma_local(freq_vec, Z_exp, gamma)
A_re = calc_A_re_local(freq_vec);
R_inf = mean(real(Z_exp(:)) - A_re * gamma(:));
if ~(R_inf > 0)
    R_inf = 0;
end
end

function mean_var = estimate_resample_variance_local(freq_vec, Z_exp, el, n_boot)
n = numel(freq_vec);
if nargin < 4 || isempty(n_boot)
    n_boot = 8;
end

sigma_re = max(0.005 * max(abs(real(Z_exp))), eps);
sigma_im = max(0.005 * max(abs(imag(Z_exp))), eps);

gamma_samples = zeros(n, n_boot);
for k = 1:n_boot
    noise = sigma_re * randn(n,1) + 1i * sigma_im * randn(n,1);
    Z_boot = Z_exp(:) + noise;
    [gk, ~] = TR_DRT_local(freq_vec, Z_boot, el);
    gamma_samples(:, k) = gk(:);
end

mean_var = mean(var(gamma_samples, 0, 2));
end

function nvec = normalize_metric_local(vec)
vec = vec(:);
vmin = min(vec);
vmax = max(vec);
if vmax > vmin
    nvec = (vec - vmin) / (vmax - vmin);
else
    nvec = zeros(size(vec));
end
end

function A_re = calc_A_re_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);

A_re = zeros(n, n);
for p = 1:n
    for q = 1:n
        if q == 1
            log_term = log(tau(2) / tau(1));
        elseif q == n
            log_term = log(tau(n) / tau(n-1));
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
            log_term = log(tau(2) / tau(1));
        elseif q == n
            log_term = log(tau(n) / tau(n-1));
        else
            log_term = log(tau(q+1) / tau(q-1));
        end
        A_im(p, q) = 0.5 * (omega(p) * tau(q)) / (1 + (omega(p) * tau(q))^2) * log_term;
    end
end
end