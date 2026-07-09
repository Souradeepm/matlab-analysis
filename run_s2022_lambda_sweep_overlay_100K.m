function run_s2022_lambda_sweep_overlay_100K()
% run_s2022_lambda_sweep_overlay_100K
% Run the same lambda sweep overlay on S2022 at the temperature nearest 100 K.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_100K_lambda_sweep_overlay.xlsx');
out_png = fullfile(repo_root, 's2022_100K_lambda_sweep_overlay.png');

target_temperature_K = 100;

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text, ~] = xlsread(input_path);
if isempty(data) || size(data,2) < 3
    error('S2022 data appears empty or malformed in: %s', input_path);
end

n_sets = floor(size(data,2) / 3);

% Parse available temperatures from first header row (same style as existing workflow).
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
if isempty(temps)
    error('Could not parse temperature labels from S2022 header row.');
end

% Keep only as many temperatures as available triplet sets.
temps = temps(1:min(numel(temps), n_sets));
if isempty(temps)
    error('No temperature entries matched numeric triplet sets.');
end

[~, t_idx] = min(abs(temps - target_temperature_K));
selected_temperature = temps(t_idx);

col_start = 3 * t_idx - 2;
col_end = 3 * t_idx;
if col_end > size(data,2)
    error('Selected temperature index exceeds data columns.');
end

slice = data(:, col_start:col_end);
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

% Same lambda list as requested previously (includes duplicate 1e-2 entry).
lambda_values = [1e-4, 1e-3, 5e-3, 1e-2, 1e-2, 1e-1];
nl = numel(lambda_values);
nf = numel(freq_vec);

gamma_all = zeros(nf, nl);
R_inf_all = zeros(nl, 1);
mean_resid_all = zeros(nl, 1);
peak_count_all = zeros(nl, 1);
peak_height_thr = 0.05;

fprintf('S2022 lambda sweep near %.1f K (selected %.1f K)\n', target_temperature_K, selected_temperature);
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

figure('Color', 'w', 'Name', sprintf('S2022 DRT overlay (%.1f K)', selected_temperature));
cm = lines(nl);
for i = 1:nl
    semilogx(tau, gamma_all(:,i), 'LineWidth', 1.3, 'Color', cm(i,:)); hold on;
end
set(gca, 'XDir', 'reverse');
xlabel('Relaxation time tau (s)');
ylabel('gamma(tau)');
title(sprintf('S2022 DRT overlay near 100 K (selected %.1f K)', selected_temperature));
grid on;

legend_text = cell(nl,1);
for i = 1:nl
    legend_text{i} = sprintf('\\lambda=%.1e (peaks=%d)', lambda_values(i), peak_count_all(i));
end
legend(legend_text, 'Location', 'best');
saveas(gcf, out_png);

xlswrite(out_xlsx, {'target_temperature_K','selected_temperature_K','temperature_index','n_temperature_sets'}, 'meta', 'A1');
xlswrite(out_xlsx, [target_temperature_K, selected_temperature, t_idx, n_sets], 'meta', 'A2');

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
