function run_s2022_lambda_sweep_overlay_3peak_near100K()
% run_s2022_lambda_sweep_overlay_3peak_near100K
% Perform the same lambda sweep on S2022 at a temperature with 3 peaks at lambda=1e-3,
% chosen as the closest such temperature to 100 K.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_3peak_near100K_lambda_sweep_overlay.xlsx');
out_png = fullfile(repo_root, 's2022_3peak_near100K_lambda_sweep_overlay.png');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text, ~] = xlsread(input_path);
if isempty(data) || size(data,2) < 3
    error('S2022 data appears empty or malformed in: %s', input_path);
end

n_sets = floor(size(data,2) / 3);
if iscell(text) && ~isempty(text)
    header_row = text(1,:);
else
    header_row = {};
end

temps = nan(1, numel(header_row));
for i = 1:numel(header_row)
    token = header_row{i};
    if ischar(token)
        token_u = upper(strtrim(token));
        token_u(token_u == 'K') = [];
        temps(i) = str2double(token_u);
    end
end

temps = temps(~isnan(temps));
temps = temps(1:min(numel(temps), n_sets));
if isempty(temps)
    error('Could not parse temperature labels from S2022 header row.');
end

% First find temperatures with 3 peaks at lambda = 1e-3.
lambda_ref = 1e-3;
peak_height_thr = 0.05;
peak_count_ref = zeros(numel(temps),1);

for t = 1:numel(temps)
    cols = (3*t-2):(3*t);
    slice = data(:, cols);

    freq_vec_t = slice(:,1);
    mag_t = slice(:,2);
    phase_t = slice(:,3);

    Z_t = mag_t .* exp(1i * phase_t * pi / 180);
    valid_t = freq_vec_t > 0;
    freq_vec_t = freq_vec_t(valid_t);
    Z_t = Z_t(valid_t);
    [freq_vec_t, idx_t] = sort(freq_vec_t);
    Z_t = Z_t(idx_t);

    [gamma_t, ~] = tr_drt_local(freq_vec_t, Z_t, lambda_ref);
    thr_t = peak_height_thr * max(gamma_t);
    [pk_idx_t, ~] = find_peaks_simple_local(gamma_t, thr_t, 10);
    peak_count_ref(t) = numel(pk_idx_t);
end

idx3 = find(peak_count_ref == 3);
if isempty(idx3)
    error('No S2022 temperature has exactly 3 peaks at lambda=1e-3 with this threshold.');
end

target_K = 100;
[~, k_local] = min(abs(temps(idx3) - target_K));
t_idx = idx3(k_local);
selected_temperature = temps(t_idx);

cols = (3*t_idx-2):(3*t_idx);
slice = data(:, cols);

freq_vec = slice(:,1);
mag = slice(:,2);
phase_deg = slice(:,3);

Z_exp = mag .* exp(1i * phase_deg * pi / 180);
valid = freq_vec > 0;
freq_vec = freq_vec(valid);
Z_exp = Z_exp(valid);
[freq_vec, idx] = sort(freq_vec);
Z_exp = Z_exp(idx);

tau = 1 ./ (2 * pi * freq_vec);

lambda_values = [1e-4, 1e-3, 5e-3, 1e-2, 1e-2, 1e-1];
nl = numel(lambda_values);
nf = numel(freq_vec);

gamma_all = zeros(nf, nl);
R_inf_all = zeros(nl, 1);
mean_resid_all = zeros(nl, 1);
peak_count_all = zeros(nl, 1);

fprintf('S2022 same sweep at 3-peak temperature closest to 100 K\n');
fprintf('Selected temperature: %.2f K\n', selected_temperature);

for i = 1:nl
    lam = lambda_values(i);
    [gamma_i, rinf_i] = tr_drt_local(freq_vec, Z_exp, lam);
    z_fit_i = calc_eis_local(freq_vec, gamma_i, rinf_i);

    gamma_all(:, i) = gamma_i(:);
    R_inf_all(i) = rinf_i;
    mean_resid_all(i) = mean(abs(Z_exp - z_fit_i));

    thr_i = peak_height_thr * max(gamma_i);
    [pk_idx, ~] = find_peaks_simple_local(gamma_i, thr_i, 10);
    peak_count_all(i) = numel(pk_idx);

    fprintf('  lambda=%8.1e | mean residual=%10.6f | peaks=%d\n', lam, mean_resid_all(i), peak_count_all(i));
end

figure('Color', 'w', 'Name', sprintf('S2022 DRT overlay (%.2f K)', selected_temperature));
cm = lines(nl);
for i = 1:nl
    semilogx(tau, gamma_all(:,i), 'LineWidth', 1.3, 'Color', cm(i,:)); hold on;
end
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('S2022 DRT overlay at %.2f K (3 peaks at \\lambda=1e-3)', selected_temperature));
grid on;

legend_text = cell(nl,1);
for i = 1:nl
    legend_text{i} = sprintf('\\lambda=%.1e (peaks=%d)', lambda_values(i), peak_count_all(i));
end
legend(legend_text, 'Location', 'best');
saveas(gcf, out_png);

xlswrite(out_xlsx, {'target_temperature_K','selected_temperature_K','selected_index','lambda_ref_for_3peak','peak_threshold_ratio'}, 'meta', 'A1');
xlswrite(out_xlsx, [target_K, selected_temperature, t_idx, lambda_ref, peak_height_thr], 'meta', 'A2');

xlswrite(out_xlsx, {'lambda','R_inf','mean_abs_residual','peak_count_5pct'}, 'summary', 'A1');
xlswrite(out_xlsx, [lambda_values(:), R_inf_all, mean_resid_all, peak_count_all], 'summary', 'A2');

detail_header = cell(1, nl + 1);
detail_header{1} = 'tau_s';
for i = 1:nl
    detail_header{i+1} = sprintf('gamma_lambda_%0.1e', lambda_values(i));
end
xlswrite(out_xlsx, detail_header, 'drt_curves', 'A1');
xlswrite(out_xlsx, [tau, gamma_all], 'drt_curves', 'A2');

fprintf('Saved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_png);
end

function [gamma, R_inf] = tr_drt_local(freq_vec, Z_exp, lam)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);

Z_re = real(Z_exp(:));
Z_im = imag(Z_exp(:));
n = numel(freq_vec);

M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(lam/2) * eye(n), zeros(n,1)];
b = [Z_re; Z_im; zeros(n,1)];

x = lsqnonneg(M, b);
gamma = x(1:end-1);
R_inf = x(end);
end

function Z_cal = calc_eis_local(freq_vec, gamma, R_inf)
A_re = calc_A_re_local(freq_vec);
A_im = calc_A_im_local(freq_vec);
Z_cal = R_inf + A_re * gamma + 1i * (A_im * gamma);
end

function A_re = calc_A_re_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);
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

function A_im = calc_A_im_local(freq)
omega = 2 * pi * freq(:);
tau = 1 ./ freq(:);
n = numel(freq);
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

function [peak_indices, peak_values] = find_peaks_simple_local(y, min_height, min_distance)
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
