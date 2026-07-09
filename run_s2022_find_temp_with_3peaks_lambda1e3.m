function run_s2022_find_temp_with_3peaks_lambda1e3()
% run_s2022_find_temp_with_3peaks_lambda1e3
% Scan S2022 temperatures and find where lambda=1e-3 yields 3 distinct DRT peaks.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_lambda1e3_peakcount_scan.xlsx');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text, ~] = xlsread(input_path);
if isempty(data) || size(data,2) < 3
    error('S2022 data appears empty or malformed in: %s', input_path);
end

n_sets = floor(size(data,2) / 3);

% Parse temperature labels from first header row.
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
    error('Could not parse temperature labels for S2022 sets.');
end

lambda = 1e-3;
peak_height_thr = 0.05;

n_eval = numel(temps);
peak_count = zeros(n_eval,1);
mean_resid = zeros(n_eval,1);
R_inf_all = zeros(n_eval,1);

fprintf('Scanning S2022 sets for lambda=%.1e with 5%% peak threshold...\n', lambda);

for t = 1:n_eval
    cols = (3*t-2):(3*t);
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

    [gamma_i, rinf_i] = tr_drt_local(freq_vec, Z_exp, lambda);
    z_fit_i = calc_eis_local(freq_vec, gamma_i, rinf_i);

    thr_i = peak_height_thr * max(gamma_i);
    [pk_idx, ~] = find_peaks_simple_local(gamma_i, thr_i, 10);

    peak_count(t) = numel(pk_idx);
    mean_resid(t) = mean(abs(Z_exp - z_fit_i));
    R_inf_all(t) = rinf_i;

    fprintf('  T=%7.2f K | peaks=%d | mean residual=%10.6f\n', temps(t), peak_count(t), mean_resid(t));
end

header = {'temperature_K','peak_count_lambda_1e_3','mean_abs_residual','R_inf'};
xlswrite(out_xlsx, header, 'scan', 'A1');
xlswrite(out_xlsx, [temps(:), peak_count, mean_resid, R_inf_all], 'scan', 'A2');

idx3 = find(peak_count == 3);
if isempty(idx3)
    fprintf('\nNo temperature produced exactly 3 peaks at lambda=1e-3.\n');
else
    fprintf('\nTemperatures with 3 peaks at lambda=1e-3:\n');
    for k = 1:numel(idx3)
        fprintf('  %.2f K\n', temps(idx3(k)));
    end
    xlswrite(out_xlsx, {'temp_with_3_peaks_K'}, 'three_peaks', 'A1');
    xlswrite(out_xlsx, temps(idx3(:)), 'three_peaks', 'A2');
end

fprintf('Saved: %s\n', out_xlsx);
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
