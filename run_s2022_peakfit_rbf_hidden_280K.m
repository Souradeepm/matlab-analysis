function run_s2022_peakfit_rbf_hidden_280K()
% run_s2022_peakfit_rbf_hidden_280K
% Peak fitting on S2022 near 280 K using RBF smoothing and hidden-peak
% detection from second derivative, then count peaks per lambda.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_xlsx = fullfile(repo_root, 's2022_280K_rbf_hidden_peakfit.xlsx');
out_png = fullfile(repo_root, 's2022_280K_rbf_hidden_peakfit.png');

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

temps = parse_temperatures(header_row);
temps = temps(1:min(numel(temps), n_sets));
if isempty(temps)
    error('Could not parse temperature labels from S2022 header row.');
end

% Use same lambda sweep list as before.
lambda_values = [1e-4, 1e-3, 5e-3, 1e-2, 1e-2, 1e-1];

% Select the nearest temperature to 280 K.
[~, t_idx] = min(abs(temps - 280));
sel_T = temps(t_idx);

cols = (3*t_idx-2):(3*t_idx);
slice = data(:, cols);
freq_vec = slice(:,1);
mag = slice(:,2);
phase_deg = slice(:,3);

Z_exp = mag .* exp(1i * phase_deg * pi / 180);
valid = freq_vec > 0;
freq_vec = freq_vec(valid);
Z_exp = Z_exp(valid);
[freq_vec, ord] = sort(freq_vec);
Z_exp = Z_exp(ord);

tau = 1 ./ (2 * pi * freq_vec);
xlog = log(tau(:));

nl = numel(lambda_values);
peak_count = zeros(nl,1);
base_count = zeros(nl,1);
hidden_count = zeros(nl,1);
mean_resid = zeros(nl,1);
R_inf_all = zeros(nl,1);

all_peak_rows = [];

figure('Color', 'w', 'Name', sprintf('S2022 %.2f K RBF + 2nd derivative hidden peaks', sel_T));

fprintf('RBF + hidden-peak (2nd derivative) analysis at %.2f K\n', sel_T);
for i = 1:nl
    lam = lambda_values(i);

    [gamma_i, rinf_i] = tr_drt_local(freq_vec, Z_exp, lam);
    z_fit = calc_eis_local(freq_vec, gamma_i, rinf_i);

    mean_resid(i) = mean(abs(Z_exp - z_fit));
    R_inf_all(i) = rinf_i;

    % RBF smoothing on log(tau) domain.
    rbf_model = create_rbf_model_2011(xlog, gamma_i);
    tau_dense = logspace(log10(min(tau)), log10(max(tau)), 1200).';
    x_dense = log(tau_dense);
    gamma_rbf = rbf_eval_2011(rbf_model, x_dense);
    gamma_rbf(gamma_rbf < 0) = 0;

    [idx_base, idx_hidden, idx_all] = detect_peaks_with_second_derivative(x_dense, gamma_rbf);

    base_count(i) = numel(idx_base);
    hidden_count(i) = numel(idx_hidden);
    peak_count(i) = numel(idx_all);

    fprintf('  lambda=%8.1e | peaks total=%d (base=%d, hidden=%d) | mean residual=%9.6f\n', ...
        lam, peak_count(i), base_count(i), hidden_count(i), mean_resid(i));

    % Local Gaussian parameterization around each detected peak.
    fit_struct = fit_peaks_gaussian_local(x_dense, gamma_rbf, idx_all, idx_hidden);

    for p = 1:numel(fit_struct)
        all_peak_rows(end+1, :) = [i, lam, p, fit_struct(p).tau, fit_struct(p).amplitude, ...
            fit_struct(p).sigma, fit_struct(p).R, fit_struct(p).C, fit_struct(p).n, fit_struct(p).is_hidden]; %#ok<AGROW>
    end

    subplot(3,2,i);
    semilogx(tau_dense, gamma_rbf, 'b-', 'LineWidth', 1.1); hold on;
    if ~isempty(idx_base)
        plot(tau_dense(idx_base), gamma_rbf(idx_base), 'go', 'MarkerSize', 5, 'LineWidth', 1.1);
    end
    if ~isempty(idx_hidden)
        plot(tau_dense(idx_hidden), gamma_rbf(idx_hidden), 'mo', 'MarkerSize', 6, 'LineWidth', 1.1);
    end
    set(gca, 'XDir', 'reverse');
    title(sprintf('\\lambda=%.1e | total=%d', lam, peak_count(i)));
    xlabel('tau (s)');
    ylabel('gamma');
    grid on;
end

saveas(gcf, out_png);

% Export summary and peak tables.
xlswrite(out_xlsx, {'selected_temperature_K','temperature_index','target_temperature_K'}, 'meta', 'A1');
xlswrite(out_xlsx, [sel_T, t_idx, 280], 'meta', 'A2');

xlswrite(out_xlsx, {'lambda','R_inf','mean_abs_residual','base_peak_count','hidden_peak_count','total_peak_count'}, 'summary', 'A1');
xlswrite(out_xlsx, [lambda_values(:), R_inf_all, mean_resid, base_count, hidden_count, peak_count], 'summary', 'A2');

xlswrite(out_xlsx, {'lambda_index','lambda','peak_id','tau_peak_s','amplitude','sigma_logtau','R_est','C_est','n_est','is_hidden'}, 'peak_table', 'A1');
if ~isempty(all_peak_rows)
    xlswrite(out_xlsx, all_peak_rows, 'peak_table', 'A2');
end

fprintf('Saved: %s\n', out_xlsx);
fprintf('Saved: %s\n', out_png);
end

function temps = parse_temperatures(header_row)
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
end

function [idx_base, idx_hidden, idx_all] = detect_peaks_with_second_derivative(x, y)
% Base peaks from RBF curve + hidden peaks from 2nd derivative channel.

x = x(:);
y = y(:);

% Base peaks from amplitude threshold.
base_thr = max(0.05 * max(y), eps);
[idx_base, ~] = find_peaks_simple_local(y, base_thr, 20);

% Numerical second derivative on log(tau) axis.
dx = mean(diff(x));
if ~(dx > 0)
    dx = 1;
end

d2 = zeros(size(y));
for i = 2:numel(y)-1
    d2(i) = (y(i+1) - 2*y(i) + y(i-1)) / (dx*dx);
end
kappa = -d2; % peaks/shoulders become positive in this channel
kappa(kappa < 0) = 0;

% Detect hidden candidates from curvature peaks.
curv_thr = max(0.25 * max(kappa), 3.0 * std(kappa));
[idx_curv, ~] = find_peaks_simple_local(kappa, curv_thr, 20);

% Keep only curvature peaks not already represented by base peaks.
idx_hidden = [];
for i = 1:numel(idx_curv)
    ii = idx_curv(i);
    near_base = any(abs(idx_base - ii) <= 10);
    amp_ok = y(ii) >= 0.03 * max(y);

    j1 = max(1, ii-18);
    j2 = min(numel(y), ii+18);
    local_floor = min(y(j1:j2));
    shoulder_depth_ok = (y(ii) - local_floor) >= 0.01 * max(y);

    if ~near_base && amp_ok && shoulder_depth_ok
        idx_hidden(end+1,1) = ii; %#ok<AGROW>
    end
end

idx_all = unique([idx_base(:); idx_hidden(:)]);
idx_all = sort(idx_all);

% Remove boundary artifacts.
valid = idx_all > 5 & idx_all < (numel(y)-5);
idx_all = idx_all(valid);

if isempty(idx_hidden)
    return;
end
keep_hidden = idx_hidden > 5 & idx_hidden < (numel(y)-5);
idx_hidden = idx_hidden(keep_hidden);
end

function fit_struct = fit_peaks_gaussian_local(x, y, idx_all, idx_hidden)
if isempty(idx_all)
    fit_struct = [];
    return;
end

n = numel(y);
fit_struct = repmat(struct('tau',0,'amplitude',0,'sigma',0,'R',0,'C',0,'n',0,'is_hidden',0), numel(idx_all), 1);

for i = 1:numel(idx_all)
    c = idx_all(i);
    w = max(24, round(n/80));
    i1 = max(1, c-w);
    i2 = min(n, c+w);

    xx = x(i1:i2);
    yy = y(i1:i2);

    a0 = max(y(c), eps);
    m0 = x(c);
    s0 = max(0.10, 0.18 * (max(xx)-min(xx) + eps));
    b0 = max(min(yy), 0);

    p0 = [log(a0); m0; log(s0); b0];
    obj = @(p) sum((exp(p(1)) * exp(-0.5 * ((xx - p(2)) / exp(p(3))).^2) + p(4) - yy).^2);
    opts = optimset('Display','off','MaxFunEvals',3000,'MaxIter',1500);
    p = fminsearch(obj, p0, opts);

    amp = max(exp(p(1)), 0);
    mu = p(2);
    sig = max(exp(p(3)), 1e-6);

    tau_peak = exp(mu);
    R_est = amp * sig * sqrt(2*pi);
    C_est = tau_peak / max(R_est, eps);
    n_est = 1.144 / (1 + sig);

    fit_struct(i).tau = tau_peak;
    fit_struct(i).amplitude = amp;
    fit_struct(i).sigma = sig;
    fit_struct(i).R = R_est;
    fit_struct(i).C = C_est;
    fit_struct(i).n = n_est;
    fit_struct(i).is_hidden = double(any(abs(idx_hidden - c) <= 2));
end

all_tau = [fit_struct.tau].';
[~, ord] = sort(all_tau, 'descend');
fit_struct = fit_struct(ord);
end

function [idx, vals] = find_peaks_simple_local(y, min_height, min_distance)
y = y(:);
n = numel(y);

cand = [];
for i = 2:n-1
    if y(i) >= y(i-1) && y(i) > y(i+1) && y(i) >= min_height
        cand(end+1,1) = i; %#ok<AGROW>
    end
end

if isempty(cand)
    idx = [];
    vals = [];
    return;
end

v = y(cand);
[~, ord] = sort(v, 'descend');
sel = [];
for k = 1:numel(ord)
    ii = cand(ord(k));
    if isempty(sel) || all(abs(sel - ii) >= min_distance)
        sel(end+1,1) = ii; %#ok<AGROW>
    end
end

sel = sort(sel);
idx = sel;
vals = y(sel);
end

function model = create_rbf_model_2011(x, y)
x = x(:);
y = y(:);

if numel(x) ~= numel(y)
    error('RBF input x and y must be same length');
end

n = numel(x);
if n < 3
    error('Need at least 3 points for RBF fitting');
end

x_sorted = sort(x);
dx = diff(x_sorted);
dx = dx(dx > 0);
if isempty(dx)
    scale = 1;
else
    scale = median(dx);
end
epsilon = 1 / max(scale, eps);

X = x(:, ones(1,n));
DX = X - X.';
Phi = exp(-(epsilon * DX).^2);
reg = 1e-8 * eye(n);
w = (Phi + reg) \ y;

model = struct('x', x, 'w', w, 'epsilon', epsilon);
end

function yq = rbf_eval_2011(model, xq)
xq = xq(:);
mx = model.x(:).';
DX = xq(:, ones(1,numel(mx))) - mx(ones(numel(xq),1), :);
Phi_q = exp(-(model.epsilon * DX).^2);
yq = Phi_q * model.w;
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
