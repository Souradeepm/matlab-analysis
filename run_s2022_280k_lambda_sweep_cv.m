function run_s2022_280k_lambda_sweep_cv()
% run_s2022_280k_lambda_sweep_cv
% S2022 near-280 K lambda sweep with refined peak counting and 15% holdout CV.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_280k_lambda_sweep_cv.xlsx');
out_sweep_png = fullfile(repo_root, 's2022_280k_lambda_sweep_refined_peaks.png');
out_cv_png = fullfile(repo_root, 's2022_280k_cv_resampling.png');
out_cv_iter_png = fullfile(repo_root, 's2022_280k_cv_iterations_lambda_1e3.png');
out_cv_drt_overlay_png = fullfile(repo_root, 's2022_280k_cv_drt_overlay_lambda_1e3.png');
out_lambda1e3_refined_overlay_png = fullfile(repo_root, 's2022_280k_lambda_1e3_refined_peak_overlay.png');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text, raw] = xlsread(input_path); %#ok<ASGLU>
text = char(text(1,:));
[p1, ~] = size(text);
temperature = zeros(p1,1);
for p = 1:p1
    temp = text(p,:);
    temp(temp == 'K') = [];
    temp = str2double(temp);
    temperature(p) = temp;
end
temperature = temperature(~isnan(temperature));

n_sets = floor(size(data, 2) / 3);
temperature = temperature(1:min(numel(temperature), n_sets));
[~, idx_sel] = min(abs(temperature - 280));
selected_temperature = temperature(idx_sel);
selected_data = data(:, (3 * idx_sel - 2):(3 * idx_sel));

freq_vec = selected_data(:,1);
mag = selected_data(:,2);
phase_deg = selected_data(:,3);
Z_exp = mag .* exp(1i * phase_deg * pi / 180);
valid = freq_vec > 0;
freq_vec = freq_vec(valid);
Z_exp = Z_exp(valid);
[freq_vec, sort_idx] = sort(freq_vec);
Z_exp = Z_exp(sort_idx);
tau = 1 ./ (2 * pi * freq_vec);

lambda_values = [1e-4; 1e-3; 5e-3; 1e-2; 1e-1];
n_lambda = numel(lambda_values);

mean_resid = zeros(n_lambda, 1);
R_inf_all = zeros(n_lambda, 1);
refined_peak_count = zeros(n_lambda, 1);
base_peak_count = zeros(n_lambda, 1);
hidden_peak_count = zeros(n_lambda, 1);
fit_rmse = zeros(n_lambda, 1);
gamma_all = zeros(numel(freq_vec), n_lambda);
lambda_peak_overlay = 1e-3;
fit_peaks_1e3 = [];
peak_refine_1e3 = [];

fprintf('S2022 near 280 K lambda sweep with refined peak counting\n');
fprintf('Selected temperature: %.2f K\n', selected_temperature);
for i = 1:n_lambda
    lam = lambda_values(i);
    tau_basis = 1 ./ freq_vec;
    [gamma_i, R_inf_i] = tr_drt_general(freq_vec, Z_exp, lam, tau_basis);
    Z_fit_i = calc_eis_general(freq_vec, gamma_i, R_inf_i, tau_basis);
    [fit_peaks, peak_refine] = refine_rbf_peak_fitting_general(tau_basis, gamma_i);

    gamma_all(:, i) = gamma_i(:);
    mean_resid(i) = mean(abs(Z_exp - Z_fit_i));
    R_inf_all(i) = R_inf_i;
    refined_peak_count(i) = numel(fit_peaks);
    hidden_peak_count(i) = sum([fit_peaks.is_hidden]);
    base_peak_count(i) = refined_peak_count(i) - hidden_peak_count(i);
    fit_rmse(i) = peak_refine.fit_rmse;

    if abs(lam - lambda_peak_overlay) <= 1e-12
        fit_peaks_1e3 = fit_peaks;
        peak_refine_1e3 = peak_refine;
    end

    fprintf('  lambda=%8.1e | mean residual=%9.6f | refined peaks=%d (base=%d hidden=%d)\n', ...
        lam, mean_resid(i), refined_peak_count(i), base_peak_count(i), hidden_peak_count(i));
end

% If 1e-4 and 1e-3 are nearly identical in fit quality, prefer 1e-3.
resid_1e4 = mean_resid(1);
resid_1e3 = mean_resid(2);
resid_diff_pct = 100 * abs(resid_1e3 - resid_1e4) / max(resid_1e4, eps);
small_diff_threshold_pct = 2.0;
if resid_diff_pct <= small_diff_threshold_pct
    lambda_cv = 1e-3;
    lambda_cv_reason = 1;
else
    lambda_cv = 1e-4;
    lambda_cv_reason = 0;
end

fprintf('\nLambda choice for CV:\n');
fprintf('  Residual difference between 1e-4 and 1e-3: %.3f %%\n', resid_diff_pct);
fprintf('  Small-difference threshold                : %.3f %%\n', small_diff_threshold_pct);
fprintf('  Selected CV lambda                        : %.1e\n', lambda_cv);

% 15% holdout CV using the selected lambda.
rand('twister', 280);
randn('state', 281);
n = numel(freq_vec);
holdout_n = max(1, round(0.15 * n));
n_cv = 5;
train_resid = zeros(n_cv, 1);
holdout_resid = zeros(n_cv, 1);
train_peak_count = zeros(n_cv, 1);
lambda_cv_overlay = 1e-3;
gamma_cv_overlay = zeros(numel(tau), n_cv);

for r = 1:n_cv
    order = randperm(n);
    idx_hold = sort(order(1:holdout_n));
    idx_train = sort(order(holdout_n+1:end));

    freq_train = freq_vec(idx_train);
    Z_train = Z_exp(idx_train);
    freq_hold = freq_vec(idx_hold);
    Z_hold = Z_exp(idx_hold);

    tau_basis_train = 1 ./ freq_train;
    [gamma_train, R_inf_train] = tr_drt_general(freq_train, Z_train, lambda_cv, tau_basis_train);
    Z_train_fit = calc_eis_general(freq_train, gamma_train, R_inf_train, tau_basis_train);
    Z_hold_fit = calc_eis_general(freq_hold, gamma_train, R_inf_train, tau_basis_train);
    [fit_peaks_cv, peak_refine_cv] = refine_rbf_peak_fitting_general(tau_basis_train, gamma_train); %#ok<ASGLU>

    % Always store a fixed-lambda (1e-3) training DRT on the full tau grid for overlay plotting.
    tau_basis_overlay = tau;
    [gamma_train_overlay, ~] = tr_drt_general(freq_train, Z_train, lambda_cv_overlay, tau_basis_overlay);
    gamma_cv_overlay(:, r) = gamma_train_overlay(:);

    train_resid(r) = mean(abs(Z_train - Z_train_fit));
    holdout_resid(r) = mean(abs(Z_hold - Z_hold_fit));
    train_peak_count(r) = numel(fit_peaks_cv);
    fprintf('  CV run %d/%d | train resid=%.6f | holdout resid=%.6f | peaks=%d\n', ...
        r, n_cv, train_resid(r), holdout_resid(r), train_peak_count(r));
end

fprintf('\n15%% holdout CV summary at lambda %.1e:\n', lambda_cv);
fprintf('  Mean train residual   : %.6f\n', mean(train_resid));
fprintf('  Mean holdout residual : %.6f\n', mean(holdout_resid));
fprintf('  Std holdout residual  : %.6f\n', std(holdout_resid));
fprintf('  Peak count mode       : %d\n', mode(train_peak_count));

figure('Color', 'w', 'Name', 'S2022 280 K lambda sweep refined peaks');
cm = lines(n_lambda);
for i = 1:n_lambda
    semilogx(tau, gamma_all(:, i), 'Color', cm(i,:), 'LineWidth', 1.2); hold on;
end
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('S2022 %.2f K lambda sweep with refined peak counts', selected_temperature));
grid on;
legend_text = cell(n_lambda, 1);
for i = 1:n_lambda
    legend_text{i} = sprintf('\\lambda=%.1e | peaks=%d', lambda_values(i), refined_peak_count(i));
end
legend(legend_text, 'Location', 'best');
saveas(gcf, out_sweep_png);

if ~isempty(fit_peaks_1e3) && ~isempty(peak_refine_1e3)
    tau_dense_1e3 = peak_refine_1e3.tau_dense(:);
    x_dense_1e3 = log(tau_dense_1e3);
    gamma_rbf_1e3 = peak_refine_1e3.gamma_rbf(:);
    fit_curve_1e3 = peak_refine_1e3.fit_curve(:);
    n_comp = numel(fit_peaks_1e3);
    comp_curves = zeros(numel(tau_dense_1e3), n_comp);
    if n_comp > 0
        baseline_1e3 = fit_peaks_1e3(1).baseline;
    else
        baseline_1e3 = 0;
    end

    for k = 1:n_comp
        mu_k = log(fit_peaks_1e3(k).tau);
        sigma_k = fit_peaks_1e3(k).sigma;
        amp_k = fit_peaks_1e3(k).amplitude;
        comp_curves(:, k) = amp_k .* exp(-0.5 * ((x_dense_1e3 - mu_k) ./ sigma_k).^2);
    end

    figure('Color', 'w', 'Name', 'S2022 280 K lambda=1e-3 refined peaks overlay');
    semilogx(tau_dense_1e3, gamma_rbf_1e3, 'k-', 'LineWidth', 1.6); hold on;
    semilogx(tau_dense_1e3, fit_curve_1e3, 'r-', 'LineWidth', 1.6);
    cm_comp = lines(max(n_comp, 1));
    for k = 1:n_comp
        if fit_peaks_1e3(k).is_hidden
            ls_k = '--';
        else
            ls_k = '-';
        end
        semilogx(tau_dense_1e3, comp_curves(:, k), 'Color', cm_comp(k,:), 'LineStyle', ls_k, 'LineWidth', 1.1);
    end
    if baseline_1e3 > 0
        semilogx(tau_dense_1e3, baseline_1e3 * ones(size(tau_dense_1e3)), 'Color', [0.35 0.35 0.35], 'LineStyle', ':', 'LineWidth', 1.0);
    end
    set(gca, 'XDir', 'reverse');
    xlabel('Relaxation time tau (s)');
    ylabel('gamma(tau)');
    title('DRT and refined peak components at \lambda=1e-3');
    grid on;

    legend_items = cell(n_comp + 2 + (baseline_1e3 > 0), 1);
    legend_items{1} = 'RBF-smoothed DRT';
    legend_items{2} = 'Refined total fit';
    for k = 1:n_comp
        if fit_peaks_1e3(k).is_hidden
            legend_items{2 + k} = sprintf('Peak %d (hidden)', k);
        else
            legend_items{2 + k} = sprintf('Peak %d', k);
        end
    end
    if baseline_1e3 > 0
        legend_items{end} = 'Baseline';
    end
    legend(legend_items, 'Location', 'best');
    saveas(gcf, out_lambda1e3_refined_overlay_png);
end

figure('Color', 'w', 'Name', 'S2022 280 K holdout CV');
subplot(2,1,1);
plot(1:n_cv, train_resid, 'bo-'); hold on;
plot(1:n_cv, holdout_resid, 'rs-');
xlabel('Resample run');
ylabel('Mean residual');
title(sprintf('15%% holdout CV at \\lambda=%.1e', lambda_cv));
grid on;
legend('Train residual', 'Holdout residual', 'Location', 'best');

subplot(2,1,2);
stem(1:n_cv, train_peak_count, 'filled');
xlabel('Resample run');
ylabel('Refined peak count');
grid on;
legend('Peak count on training fit', 'Location', 'best');
saveas(gcf, out_cv_png);

% Dedicated per-iteration plot for lambda=1e-3-style CV view.
figure('Color', 'w', 'Name', 'S2022 280 K CV iterations detail');
cv_gap = holdout_resid - train_resid;
subplot(3,1,1);
bar(1:n_cv, [train_resid, holdout_resid], 0.9);
xlabel('CV iteration');
ylabel('Mean residual');
title(sprintf('CV iterations at selected \\lambda=%.1e', lambda_cv));
grid on;
legend('Train', 'Holdout', 'Location', 'best');

subplot(3,1,2);
plot(1:n_cv, cv_gap, 'k-o', 'LineWidth', 1.2, 'MarkerSize', 5);
xlabel('CV iteration');
ylabel('Holdout - Train residual');
grid on;

subplot(3,1,3);
stem(1:n_cv, train_peak_count, 'filled');
xlabel('CV iteration');
ylabel('Refined peak count');
grid on;
saveas(gcf, out_cv_iter_png);

figure('Color', 'w', 'Name', 'S2022 280 K CV DRT overlay (lambda 1e-3)');
cm_cv = lines(n_cv);
for r = 1:n_cv
    semilogx(tau, gamma_cv_overlay(:, r), 'Color', cm_cv(r,:), 'LineWidth', 1.2); hold on;
end
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title('Training DRT overlay across CV iterations at \lambda=1e-3');
grid on;
legend_text_cv = cell(n_cv, 1);
for r = 1:n_cv
    legend_text_cv{r} = sprintf('CV run %d', r);
end
legend(legend_text_cv, 'Location', 'best');
saveas(gcf, out_cv_drt_overlay_png);

xlswrite(out_xlsx, {'selected_temperature_K','resid_diff_pct_1e4_vs_1e3','small_diff_threshold_pct','selected_cv_lambda','used_1e3_flag'}, 'meta', 'A1');
xlswrite(out_xlsx, [selected_temperature, resid_diff_pct, small_diff_threshold_pct, lambda_cv, lambda_cv_reason], 'meta', 'A2');

xlswrite(out_xlsx, {'lambda','mean_abs_residual','R_inf','refined_peak_count','base_peak_count','hidden_peak_count','fit_rmse'}, 'lambda_sweep', 'A1');
xlswrite(out_xlsx, [lambda_values, mean_resid, R_inf_all, refined_peak_count, base_peak_count, hidden_peak_count, fit_rmse], 'lambda_sweep', 'A2');

xlswrite(out_xlsx, {'cv_run','train_mean_residual','holdout_mean_residual','train_refined_peak_count'}, 'cv_resampling', 'A1');
xlswrite(out_xlsx, [(1:n_cv).', train_resid, holdout_resid, train_peak_count], 'cv_resampling', 'A2');

curve_header = cell(1, n_lambda + 1);
curve_header{1} = 'tau_s';
for i = 1:n_lambda
    curve_header{i+1} = sprintf('gamma_lambda_%0.1e', lambda_values(i));
end
xlswrite(out_xlsx, curve_header, 'drt_curves', 'A1');
xlswrite(out_xlsx, [tau, gamma_all], 'drt_curves', 'A2');

cv_overlay_header = cell(1, n_cv + 1);
cv_overlay_header{1} = 'tau_s';
for r = 1:n_cv
    cv_overlay_header{r+1} = sprintf('gamma_cv_run_%d_lambda_1e3', r);
end
xlswrite(out_xlsx, cv_overlay_header, 'cv_drt_overlay_1e3', 'A1');
xlswrite(out_xlsx, [tau, gamma_cv_overlay], 'cv_drt_overlay_1e3', 'A2');

fprintf('\nSaved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_sweep_png);
fprintf('Saved: %s\n', out_cv_png);
fprintf('Saved: %s\n', out_cv_iter_png);
fprintf('Saved: %s\n', out_cv_drt_overlay_png);
if ~isempty(fit_peaks_1e3) && ~isempty(peak_refine_1e3)
    fprintf('Saved: %s\n', out_lambda1e3_refined_overlay_png);
end
end

function [gamma, R_inf] = tr_drt_general(freq_eval, Z_exp, lam, tau_basis)
A_re = calc_A_re_general(freq_eval, tau_basis);
A_im = calc_A_im_general(freq_eval, tau_basis);

Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n_eval = numel(freq_eval);
n_basis = numel(tau_basis);

M = [A_re, ones(n_eval,1); A_im, zeros(n_eval,1); sqrt(lam/2) * eye(n_basis), zeros(n_basis,1)];
b = [Z_re; Z_im; zeros(n_basis,1)];

x = lsqnonneg(M, b);
gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calc_eis_general(freq_eval, gamma, R_inf, tau_basis)
A_re = calc_A_re_general(freq_eval, tau_basis);
A_im = calc_A_im_general(freq_eval, tau_basis);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function A_re = calc_A_re_general(freq_eval, tau_basis)
omega = 2 * pi * freq_eval(:);
tau_basis = tau_basis(:);
n_eval = numel(freq_eval);
n_basis = numel(tau_basis);
A_re = zeros(n_eval, n_basis);
for p = 1:n_eval
    for q = 1:n_basis
        if q == 1
            if n_basis == 1
                log_term = 1;
            else
                log_term = log(tau_basis(q+1) / tau_basis(q));
            end
        elseif q == n_basis
            log_term = log(tau_basis(q) / tau_basis(q-1));
        else
            log_term = log(tau_basis(q+1) / tau_basis(q-1));
        end
        A_re(p, q) = -0.5 / (1 + (omega(p) * tau_basis(q))^2) * log_term;
    end
end
end

function A_im = calc_A_im_general(freq_eval, tau_basis)
omega = 2 * pi * freq_eval(:);
tau_basis = tau_basis(:);
n_eval = numel(freq_eval);
n_basis = numel(tau_basis);
A_im = zeros(n_eval, n_basis);
for p = 1:n_eval
    for q = 1:n_basis
        if q == 1
            if n_basis == 1
                log_term = 1;
            else
                log_term = log(tau_basis(q+1) / tau_basis(q));
            end
        elseif q == n_basis
            log_term = log(tau_basis(q) / tau_basis(q-1));
        else
            log_term = log(tau_basis(q+1) / tau_basis(q-1));
        end
        A_im(p, q) = 0.5 * (omega(p) * tau_basis(q)) / (1 + (omega(p) * tau_basis(q))^2) * log_term;
    end
end
end

function [fitted_peaks, peak_refine] = refine_rbf_peak_fitting_general(tau_vals, gamma_vals)
tau_vals = tau_vals(:);
gamma_vals = gamma_vals(:);
tau_log = log(tau_vals);

rbf_model = create_rbf_model_2011(tau_log, gamma_vals);
tau_dense = logspace(log10(min(tau_vals)), log10(max(tau_vals)), 1000).';
tau_dense_log = log(tau_dense);
gamma_rbf = rbf_eval_2011(rbf_model, tau_dense_log);
gamma_rbf(gamma_rbf < 0) = 0;

base_thr = max(0.05 * max(gamma_rbf), eps);
[base_peak_idx, ~] = find_peaks_simple(gamma_rbf, base_thr, 20);
hidden_peak_idx = detect_hidden_peaks_second_derivative_2011(tau_dense_log, gamma_rbf, base_peak_idx);
all_peak_idx = unique([base_peak_idx(:); hidden_peak_idx(:)]);
all_peak_idx = sort(all_peak_idx);

hidden_flag = zeros(size(all_peak_idx));
for k = 1:numel(all_peak_idx)
    if any(abs(hidden_peak_idx - all_peak_idx(k)) <= 2)
        hidden_flag(k) = 1;
    end
end

[fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2011(tau_dense, gamma_rbf, all_peak_idx, hidden_flag);
peak_refine = struct('tau_dense', tau_dense, 'gamma_rbf', gamma_rbf, 'base_peak_idx', base_peak_idx, ...
    'hidden_peak_idx', hidden_peak_idx, 'all_peak_idx', all_peak_idx, 'fit_curve', fit_curve, 'fit_rmse', fit_rmse, ...
    'rbf_epsilon', rbf_model.epsilon);
end

function hidden_peak_idx = detect_hidden_peaks_second_derivative_2011(tau_log, gamma_rbf, base_peak_idx)
tau_log = tau_log(:);
gamma_rbf = gamma_rbf(:);
n = numel(gamma_rbf);
if n < 5
    hidden_peak_idx = [];
    return;
end

dx = mean(diff(tau_log));
if ~(dx > 0)
    dx = 1;
end

d2 = zeros(n,1);
for i = 2:n-1
    d2(i) = (gamma_rbf(i+1) - 2 * gamma_rbf(i) + gamma_rbf(i-1)) / (dx * dx);
end
curvature_signal = -d2;
curvature_signal(curvature_signal < 0) = 0;
curvature_signal = smooth_sgolay_2011(curvature_signal, 7, 2);
curv_thr = max(0.30 * max(curvature_signal), 3.2 * std(curvature_signal));
[curv_idx, ~] = find_peaks_simple(curvature_signal, curv_thr, 28);

hidden_peak_idx = [];
for i = 1:numel(curv_idx)
    idx = curv_idx(i);
    near_base = any(abs(base_peak_idx - idx) <= 12);
    amp_ok = gamma_rbf(idx) >= 0.05 * max(gamma_rbf);
    j1 = max(1, idx - 18);
    j2 = min(n, idx + 18);
    local_floor = min(gamma_rbf(j1:j2));
    shoulder_depth_ok = (gamma_rbf(idx) - local_floor) >= 0.02 * max(gamma_rbf);

    l1 = max(1, idx - 25);
    l2 = idx;
    r1 = idx;
    r2 = min(n, idx + 25);
    left_min = min(gamma_rbf(l1:l2));
    right_min = min(gamma_rbf(r1:r2));
    prominence_ok = (gamma_rbf(idx) - max(left_min, right_min)) >= 0.015 * max(gamma_rbf);

    if ~near_base && amp_ok && shoulder_depth_ok && prominence_ok && idx > 5 && idx < n - 5
        hidden_peak_idx(end+1,1) = idx; %#ok<AGROW>
    end
end

% Mid-gap rescue: recover one weak central shoulder between neighboring base peaks.
if numel(base_peak_idx) >= 2
    for k = 1:(numel(base_peak_idx) - 1)
        i1 = base_peak_idx(k);
        i2 = base_peak_idx(k + 1);
        if i2 - i1 < 24
            continue;
        end
        seg = (i1 + 3):(i2 - 3);
        if isempty(seg)
            continue;
        end

        [seg_max, rel_idx] = max(gamma_rbf(seg));
        cand = seg(rel_idx);
        near_any = any(abs([base_peak_idx(:); hidden_peak_idx(:)] - cand) <= 10);
        if near_any
            continue;
        end

        valley_left = min(gamma_rbf(i1:cand));
        valley_right = min(gamma_rbf(cand:i2));
        local_prom = seg_max - max(valley_left, valley_right);
        if seg_max >= 0.04 * max(gamma_rbf) && local_prom >= 0.012 * max(gamma_rbf)
            hidden_peak_idx(end+1,1) = cand; %#ok<AGROW>
        end
    end
end

if ~isempty(hidden_peak_idx)
    hidden_peak_idx = unique(sort(hidden_peak_idx));
end

max_hidden = max(0, floor(numel(base_peak_idx) / 2));
if numel(hidden_peak_idx) > max_hidden
    [~, ord_h] = sort(gamma_rbf(hidden_peak_idx), 'descend');
    hidden_peak_idx = hidden_peak_idx(ord_h(1:max_hidden));
    hidden_peak_idx = sort(hidden_peak_idx);
end
end

function [fitted_peaks, fit_curve, fit_rmse] = fit_rbf_gaussian_peaks_2011(tau_dense, gamma_rbf, peak_idx, hidden_flag)
if isempty(peak_idx)
    fitted_peaks = [];
    fit_curve = zeros(size(gamma_rbf));
    fit_rmse = 0;
    return;
end

tau_dense = tau_dense(:);
gamma_rbf = gamma_rbf(:);
xall = log(tau_dense);
centers = xall(peak_idx(:));
[fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);

% Residual-based rescue: if the Gaussian-sum underfits between two peaks, add one central hidden basis.
if numel(centers) >= 2
    residual_pos = gamma_rbf - fit_curve;
    added_centers = [];
    for k = 1:(numel(centers) - 1)
        seg = find(xall > centers(k) & xall < centers(k+1));
        if numel(seg) < 12
            continue;
        end
        [res_max, rel_idx] = max(residual_pos(seg));
        cand = seg(rel_idx);
        if gamma_rbf(cand) < 0.035 * max(gamma_rbf)
            continue;
        end
        if res_max < max(0.02 * max(gamma_rbf), 1.2 * sqrt(mean((fit_curve - gamma_rbf).^2)))
            continue;
        end
        if any(abs([centers(:); added_centers(:)] - xall(cand)) < 0.06)
            continue;
        end
        added_centers(end+1,1) = xall(cand); %#ok<AGROW>
    end

    if ~isempty(added_centers)
        centers = sort([centers(:); added_centers(:)]);
        hidden_flag = [hidden_flag(:); ones(numel(added_centers), 1)];
        [~, ord_h] = sort([centers(:)]);
        hidden_flag = hidden_flag(ord_h);
        [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);
    end
end

peak_area = amplitudes .* sigma_vec;
amp_cutoff = 0.05 * max(amplitudes);
area_cutoff = 0.05 * max(peak_area);
keep = (amplitudes >= amp_cutoff) & (peak_area >= area_cutoff);
if any(~keep) && any(keep)
    centers = centers(keep);
    hidden_flag = hidden_flag(keep);
    [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers);
end

fit_rmse = sqrt(mean((fit_curve - gamma_rbf).^2));
m = numel(centers);
fitted_peaks = repmat(struct('peak_id', 0, 'tau', 0, 'amplitude', 0, 'sigma', 0, 'R', 0, 'C', 0, 'n', 0, 'is_hidden', 0, 'baseline', 0), m, 1);
for k = 1:m
    tau_peak = exp(centers(k));
    amp = amplitudes(k);
    sig = sigma_vec(k);
    R_est = amp * sig * sqrt(2 * pi);
    C_est = tau_peak / max(R_est, eps);
    n_est = 1.144 / (1 + sig);
    fitted_peaks(k).peak_id = k;
    fitted_peaks(k).tau = tau_peak;
    fitted_peaks(k).amplitude = amp;
    fitted_peaks(k).sigma = sig;
    fitted_peaks(k).R = R_est;
    fitted_peaks(k).C = C_est;
    fitted_peaks(k).n = n_est;
    fitted_peaks(k).is_hidden = hidden_flag(k);
    fitted_peaks(k).baseline = baseline;
end
all_tau = [fitted_peaks.tau].';
[~, order] = sort(all_tau, 'descend');
fitted_peaks = fitted_peaks(order);
for k = 1:numel(fitted_peaks)
    fitted_peaks(k).peak_id = k;
end
end

function [fit_curve, amplitudes, sigma_vec, baseline] = solve_gaussian_rbf_sum_fit_2011(xall, gamma_rbf, centers)
m = numel(centers);
sigma0 = zeros(m,1);
for k = 1:m
    if m == 1
        neighbor_span = max(xall) - min(xall);
    elseif k == 1
        neighbor_span = abs(centers(k+1) - centers(k));
    elseif k == m
        neighbor_span = abs(centers(k) - centers(k-1));
    else
        neighbor_span = min(abs(centers(k) - centers(k-1)), abs(centers(k+1) - centers(k)));
    end
    sigma0(k) = min(max(0.12, 0.35 * neighbor_span), 0.75);
end
p0 = [log(sigma0); log(max(min(gamma_rbf), eps))];
opts = optimset('Display', 'off', 'MaxFunEvals', 6000, 'MaxIter', 2500);
obj = @(p) gaussian_rbf_sum_objective_2011(p, xall, gamma_rbf, centers);
p = fminsearch(obj, p0, opts);
[fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2011(p, xall, gamma_rbf, centers);
end

function sse = gaussian_rbf_sum_objective_2011(p, xall, yall, centers)
[fit_curve, amplitudes, sigma_vec] = evaluate_gaussian_rbf_sum_2011(p, xall, yall, centers);
resid = fit_curve - yall;
sse = sum(resid .^ 2);

sigma_penalty = 0;
for k = 1:numel(sigma_vec)
    if sigma_vec(k) < 0.04
        sigma_penalty = sigma_penalty + (0.04 - sigma_vec(k)) ^ 2;
    elseif sigma_vec(k) > 0.90
        sigma_penalty = sigma_penalty + (sigma_vec(k) - 0.90) ^ 2;
    end
end
amp_penalty = 0;
for k = 1:numel(amplitudes)
    if amplitudes(k) < 1e-8 * max(yall)
        amp_penalty = amp_penalty + (1e-8 * max(yall) - amplitudes(k)) ^ 2;
    end
end
sse = sse + 1e3 * sigma_penalty + 1e-2 * amp_penalty;
end

function [fit_curve, amplitudes, sigma_vec, baseline] = evaluate_gaussian_rbf_sum_2011(p, xall, yall, centers)
sigma_vec = exp(p(1:numel(centers)));
Phi = zeros(numel(xall), numel(centers));
for k = 1:numel(centers)
    Phi(:, k) = exp(-0.5 * ((xall - centers(k)) / sigma_vec(k)).^2);
end
Phi_aug = [Phi, ones(numel(xall), 1)];
coef = lsqnonneg(Phi_aug, yall);
amplitudes = coef(1:numel(centers));
baseline = coef(end);
fit_curve = Phi_aug * coef;
end

function model = create_rbf_model_2011(x, y)
x = x(:);
y = y(:);
n = numel(x);
if n > 1
    ds = abs(diff(sort(x)));
    med_ds = median(ds);
    if med_ds <= eps
        epsilon = 1;
    else
        epsilon = 1 / med_ds;
    end
else
    epsilon = 1;
end
lambda = 1e-8;
DX = x(:, ones(1,n)) - x(:, ones(1,n)).';
Phi = exp(-(epsilon * DX).^2);
w = (Phi + lambda * eye(n)) \ y;
model.x = x;
model.w = w;
model.epsilon = epsilon;
end

function yq = rbf_eval_2011(model, xq)
xq = xq(:);
mx = model.x(:).';
DX = xq(:, ones(1,numel(mx))) - mx(ones(numel(xq),1), :);
Phi_q = exp(-(model.epsilon * DX).^2);
yq = Phi_q * model.w;
end

function [peak_indices, peak_values] = find_peaks_simple(y, min_height, min_distance)
y = y(:);
n = numel(y);
cand_idx = [];
for i = 2:n-1
    if y(i) >= y(i-1) && y(i) > y(i+1) && y(i) >= min_height
        cand_idx(end+1,1) = i; %#ok<AGROW>
    end
end
if isempty(cand_idx)
    peak_indices = [];
    peak_values = [];
    return;
end
cand_val = y(cand_idx);
[~, order] = sort(cand_val, 'descend');
selected = [];
for k = 1:numel(order)
    idx = cand_idx(order(k));
    if isempty(selected)
        selected(end+1,1) = idx; %#ok<AGROW>
    else
        if all(abs(selected - idx) >= min_distance)
            selected(end+1,1) = idx; %#ok<AGROW>
        end
    end
end
selected = sort(selected);
peak_indices = selected;
peak_values = y(selected);
end

function ysm = smooth_mavg(y, w)
y = y(:);
n = numel(y);
if w <= 1
    ysm = y;
    return;
end
half = floor(w / 2);
ysm = zeros(n,1);
for i = 1:n
    i1 = max(1, i - half);
    i2 = min(n, i + half);
    ysm(i) = mean(y(i1:i2));
end
end

function ysm = smooth_sgolay_2011(y, frame_len, poly_order)
% smooth_sgolay_2011 - Lightweight Savitzky-Golay smoothing without toolbox dependency

y = y(:);
n = numel(y);

if n < 5
    ysm = y;
    return;
end

if mod(frame_len, 2) == 0
    frame_len = frame_len + 1;
end
frame_len = min(frame_len, n);
if mod(frame_len, 2) == 0
    frame_len = frame_len - 1;
end

poly_order = min(poly_order, frame_len - 1);
if poly_order < 1 || frame_len < 3
    ysm = smooth_mavg(y, 3);
    return;
end

half = floor(frame_len / 2);
ysm = zeros(n, 1);
for i = 1:n
    i1 = max(1, i - half);
    i2 = min(n, i + half);
    xi = (i1:i2).' - i;
    yi = y(i1:i2);

    V = zeros(numel(xi), poly_order + 1);
    for p = 0:poly_order
        V(:, p + 1) = xi .^ p;
    end

    c = V \ yi;
    ysm(i) = c(1);
end
end
