function run_chemelectrochem_test()
% run_chemelectrochem_test
% Self-contained ChemElectroChem DRT optimizer test.
% Avoids xlsread/COM/figures - writes all output to a log file.

repo_root = fileparts(mfilename('fullpath'));
log_file  = fullfile(repo_root, 'chemelectrochem_test_output.txt');
input_path = fullfile(repo_root, 'S2022Sap.xlsx');

fid = fopen(log_file, 'w');
if fid < 0
    error('Cannot open log file: %s', log_file);
end

lg = @(varargin) fprintf(fid, varargin{:});
lgs = @(varargin) fprintf(varargin{:}); % also print to stdout

function dual_print(fmt, varargin)
    fprintf(fid, fmt, varargin{:});
    fprintf(fmt, varargin{:});
end

dual_print('=== ChemElectroChem DRT Optimizer Test ===\n');
dual_print('Date: %s\n', datestr(now));
dual_print('Input: %s\n\n', input_path);

% ---- Load data using readmatrix (no COM/Excel dependency) ----
try
    data = readmatrix(input_path);
    hdr  = readcell(input_path, 'Range', '1:1');
    dual_print('Data loaded via readmatrix: %d rows x %d cols\n', size(data,1), size(data,2));
catch ME
    dual_print('readmatrix failed: %s\nFalling back to xlsread...\n', ME.message);
    try
        [data, txt] = xlsread(input_path);
        hdr = txt(1,:);
        dual_print('xlsread succeeded.\n');
    catch ME2
        dual_print('xlsread also failed: %s\n', ME2.message);
        fclose(fid);
        error('Cannot load input data.');
    end
end

% ---- Parse temperatures from header ----
temperature = [];
for k = 1:numel(hdr)
    tok = hdr{k};
    if ischar(tok) || isstring(tok)
        tok = char(tok);
        tok(tok == 'K') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end+1) = val; %#ok<AGROW>
        end
    end
end
n_sets = floor(size(data,2) / 3);
temperature = temperature(1:min(numel(temperature), n_sets));

[~, idx_sel] = min(abs(temperature - 280));
selected_temperature = temperature(idx_sel);
dual_print('Selected temperature: %.4f K\n\n', selected_temperature);

sel = data(:, (3*idx_sel-2):(3*idx_sel));
freq_vec  = sel(:,1);
mag       = sel(:,2);
phase_deg = sel(:,3);
Z_exp     = mag .* exp(1i * phase_deg * pi/180);

valid = isfinite(freq_vec) & freq_vec > 0 & isfinite(mag) & isfinite(phase_deg);
freq_vec = freq_vec(valid);
Z_exp    = Z_exp(valid);
[freq_vec, si] = sort(freq_vec);
Z_exp = Z_exp(si);
dual_print('Frequency points: %d (%.2f Hz to %.2f Hz)\n\n', ...
    numel(freq_vec), min(freq_vec), max(freq_vec));

% ---- Lambda sweep with ChemElectroChem metrics ----
lambda_values = [1e-4; 1e-3; 5e-3; 1e-2; 1e-1];
n_lam = numel(lambda_values);
tau_basis = 1 ./ freq_vec;

mean_resid        = zeros(n_lam, 1);
rid_metric        = zeros(n_lam, 1);
cv_real_metric    = zeros(n_lam, 1);
cv_imag_metric    = zeros(n_lam, 1);
resample_var_metric = zeros(n_lam, 1);

dual_print('--- Lambda sweep ---\n');
for i = 1:n_lam
    lam = lambda_values(i);

    % Full DRT
    [gamma_i, R_inf_i] = tr_drt_local(freq_vec, Z_exp, lam, tau_basis);
    Z_fit_i = calc_eis_local(freq_vec, gamma_i, R_inf_i, tau_basis);
    mean_resid(i) = mean(abs(Z_exp - Z_fit_i));

    % Real-only inversion
    [g_re, R_re] = tr_drt_component_local(freq_vec, Z_exp, lam, tau_basis, 'real');
    % Imag-only inversion
    [g_im, ~]   = tr_drt_component_local(freq_vec, Z_exp, lam, tau_basis, 'imag');

    % RID (Eq 4 in paper): ||DRT_re - DRT_im||^2
    rid_metric(i) = mean((g_re - g_im).^2);

    % Cross-validation (Eq 5): predict Im from Re-DRT and vice versa
    Z_from_re = calc_eis_local(freq_vec, g_re, R_re, tau_basis);
    R_im_est  = max(0, mean(real(Z_exp) - calc_A_re_local(freq_vec, tau_basis) * g_im));
    Z_from_im = calc_eis_local(freq_vec, g_im, R_im_est, tau_basis);
    cv_imag_metric(i) = mean((imag(Z_exp) - imag(Z_from_re)).^2);
    cv_real_metric(i) = mean((real(Z_exp) - real(Z_from_im)).^2);

    % Resampling variance (paper's repetitive-measurement idea)
    resample_var_metric(i) = compute_resample_var(freq_vec, Z_exp, lam, tau_basis, 8);

    dual_print('  lam=%8.1e | resid=%9.4f | RID=%9.3e | CV=%9.3e | Var=%9.3e\n', ...
        lam, mean_resid(i), rid_metric(i), ...
        cv_real_metric(i)+cv_imag_metric(i), resample_var_metric(i));
end

% ---- ChemElectroChem selection: boundary from CV+Var, min RID inside ----
cv_total = cv_real_metric + cv_imag_metric;
[~, idx_cv]  = min(cv_total);
[~, idx_var] = min(resample_var_metric);
idx_low  = min(idx_cv, idx_var);
idx_high = max(idx_cv, idx_var);
if idx_low == idx_high
    cand = (1:n_lam).';
else
    cand = (idx_low:idx_high).';
end
[~, rr] = min(rid_metric(cand));
idx_best = cand(rr);

% Normalized combined score
sc = norm01(rid_metric) + norm01(cv_total) + norm01(resample_var_metric);
[~, idx_score] = min(sc);
if idx_score ~= idx_best
    idx_best = idx_best; % prefer RID-in-boundary choice
end
lambda_selected = lambda_values(idx_best);

dual_print('\n--- ChemElectroChem Selection ---\n');
dual_print('  CV boundary index range: [%d, %d]\n', idx_low, idx_high);
dual_print('  Best index (min RID in boundary): %d\n', idx_best);
dual_print('  Selected lambda: %.4e\n', lambda_selected);
dual_print('  RID  at selected: %.4e\n', rid_metric(idx_best));
dual_print('  CV   at selected: %.4e\n', cv_total(idx_best));
dual_print('  Var  at selected: %.4e\n', resample_var_metric(idx_best));
dual_print('  Mean residual   : %.6f\n', mean_resid(idx_best));

dual_print('\n--- Full score table ---\n');
dual_print('  %-10s %-12s %-12s %-12s %-12s %-10s %-10s\n', ...
    'lambda','mean_resid','RID','CV','Var','score','chosen');
for i = 1:n_lam
    ch = ' ';
    if i == idx_best; ch = '*'; end
    dual_print('  %-10.2e %-12.4e %-12.4e %-12.4e %-12.4e %-10.4f %s\n', ...
        lambda_values(i), mean_resid(i), rid_metric(i), cv_total(i), ...
        resample_var_metric(i), sc(i), ch);
end

fclose(fid);
fprintf('\nLog written to: %s\n', log_file);
end

% ============================================================
%  LOCAL HELPER FUNCTIONS (no external dependencies)
% ============================================================

function [gamma, R_inf] = tr_drt_local(freq, Z, lam, tau)
A_re = calc_A_re_local(freq, tau);
A_im = calc_A_im_local(freq, tau);
n = numel(freq);  nb = numel(tau);
M = [A_re, ones(n,1); A_im, zeros(n,1); sqrt(lam/2)*eye(nb), zeros(nb,1)];
b = [real(Z(:)); imag(Z(:)); zeros(nb,1)];
x = lsqnonneg(M, b);
gamma = x(1:end-1);  R_inf = x(end);
end

function [gamma, R_inf] = tr_drt_component_local(freq, Z, lam, tau, comp)
A_re = calc_A_re_local(freq, tau);
A_im = calc_A_im_local(freq, tau);
n = numel(freq);  nb = numel(tau);
switch lower(comp)
    case 'real'
        M = [A_re, ones(n,1); sqrt(lam/2)*eye(nb), zeros(nb,1)];
        b = [real(Z(:)); zeros(nb,1)];
    case 'imag'
        M = [A_im, zeros(n,1); sqrt(lam/2)*eye(nb), zeros(nb,1)];
        b = [imag(Z(:)); zeros(nb,1)];
end
x = lsqnonneg(M, b);
gamma = x(1:end-1);  R_inf = x(end);
end

function Z_cal = calc_eis_local(freq, gamma, R_inf, tau)
A_re = calc_A_re_local(freq, tau);
A_im = calc_A_im_local(freq, tau);
Z_cal = R_inf + A_re*gamma + 1i*(A_im*gamma);
end

function A_re = calc_A_re_local(freq, tau)
omega = 2*pi*freq(:);
tau   = tau(:);
n = numel(freq);  nb = numel(tau);
A_re = zeros(n, nb);
for p = 1:n
    for q = 1:nb
        if q == 1
            lt = (nb>1)*log(tau(2)/tau(1)); if lt==0; lt=1; end
        elseif q == nb
            lt = log(tau(nb)/tau(nb-1));
        else
            lt = log(tau(q+1)/tau(q-1));
        end
        A_re(p,q) = -0.5/(1+(omega(p)*tau(q))^2)*lt;
    end
end
end

function A_im = calc_A_im_local(freq, tau)
omega = 2*pi*freq(:);
tau   = tau(:);
n = numel(freq);  nb = numel(tau);
A_im = zeros(n, nb);
for p = 1:n
    for q = 1:nb
        if q == 1
            lt = (nb>1)*log(tau(2)/tau(1)); if lt==0; lt=1; end
        elseif q == nb
            lt = log(tau(nb)/tau(nb-1));
        else
            lt = log(tau(q+1)/tau(q-1));
        end
        A_im(p,q) = 0.5*(omega(p)*tau(q))/(1+(omega(p)*tau(q))^2)*lt;
    end
end
end

function mv = compute_resample_var(freq, Z, lam, tau, n_boot)
n = numel(freq);  nb = numel(tau);
sr = max(0.005*max(abs(real(Z))), eps);
si = max(0.005*max(abs(imag(Z))), eps);
gs = zeros(nb, n_boot);
for k = 1:n_boot
    noise = sr*randn(n,1) + 1i*si*randn(n,1);
    [gk,~] = tr_drt_local(freq, Z+noise, lam, tau);
    gs(:,k) = gk(:);
end
mv = mean(var(gs, 0, 2));
end

function v = norm01(x)
x = x(:);
mn = min(x);  mx = max(x);
if mx > mn;  v = (x-mn)/(mx-mn);
else;        v = zeros(size(x));
end
end
