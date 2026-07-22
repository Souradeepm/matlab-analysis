function results = compare_n_formula_residuals(repo_root)
% compare_n_formula_residuals - Compare Nyquist residuals for two n formulas.
%
% Formula A (current workflow): n = 1.144 / (1 + sigma)
% Formula B (FWHM-derived): solve sqrt(2*log(2))*sigma*n = acosh(2 + cos(pi*n))
%
% Uses previously exported S2022 workbooks:
%   - s2022sap_peak_detail_result.xlsx
%   - s2022sap_temperature_plot_data.xlsx
%   - s2022sap_kk_result.xlsx
%
% Returns struct array and writes a text summary.

if nargin < 1 || isempty(repo_root)
    repo_root = fileparts(mfilename('fullpath'));
end

peak_file = fullfile(repo_root, 's2022sap_peak_detail_result.xlsx');
plot_file = fullfile(repo_root, 's2022sap_temperature_plot_data.xlsx');
kk_file = fullfile(repo_root, 's2022sap_kk_result.xlsx');
out_file = fullfile(repo_root, 'n_formula_residual_comparison_s2022.txt');

if exist(peak_file, 'file') ~= 2
    error('Missing file: %s', peak_file);
end
if exist(plot_file, 'file') ~= 2
    error('Missing file: %s', plot_file);
end
if exist(kk_file, 'file') ~= 2
    error('Missing file: %s', kk_file);
end

[~, peak_sheets] = xlsfinfo(peak_file);
if isempty(peak_sheets)
    error('No sheets found in %s', peak_file);
end

[kk_num, ~] = xlsread(kk_file, 'summary');
if isempty(kk_num)
    error('KK summary sheet is empty in %s', kk_file);
end

results = struct('temperature', {}, 'n_peaks', {}, ...
    'mean_abs_resid_heur', {}, 'max_abs_resid_heur', {}, ...
    'mean_abs_resid_fwhm', {}, 'max_abs_resid_fwhm', {}, ...
    'best_formula', {}, 'improvement_pct', {});

fid = fopen(out_file, 'w');
if fid < 0
    error('Unable to open output file: %s', out_file);
end

fprintf(fid, 'n-formula Nyquist residual comparison for S2022\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now, 31));
fprintf(fid, 'Formula A (heuristic): n = 1.144 / (1 + sigma)\n');
fprintf(fid, 'Formula B (FWHM-derived): solve sqrt(2*log(2))*sigma*n = acosh(2 + cos(pi*n))\n\n');
fprintf(fid, 'Temp(K)\tNpeaks\tMeanResA\tMeanResB\tBest\tImprove(%%)\n');

ridx = 0;
for i = 1:numel(peak_sheets)
    sh = peak_sheets{i};
    if ~strncmpi(sh, 'PEAK_', 5)
        continue;
    end

    T = parse_temp_from_sheet_2011(sh);
    if isnan(T)
        continue;
    end

    peak_num = xlsread(peak_file, sh);
    if isempty(peak_num)
        continue;
    end

    % PEAK columns: [id, tau, amp, sigma, fwhm, R, C, n_est, is_hidden]
    sigma = peak_num(:, 4);
    Rvec = peak_num(:, 6);
    Cvec = peak_num(:, 7);

    nyq_sheet = make_sheet_name_2011_local('NYQRCN', T);
    nyq_num = xlsread(plot_file, nyq_sheet);
    if isempty(nyq_num) || size(nyq_num,2) < 3
        continue;
    end

    freq = nyq_num(:, 1);
    Z_exp = nyq_num(:, 2) + 1i * nyq_num(:, 3);

    kk_row = find(abs(kk_num(:,1) - T) < 1e-6, 1);
    if isempty(kk_row)
        continue;
    end
    R_inf = kk_num(kk_row, 3);

    n_heur = 1.144 ./ (1 + sigma);
    n_fwhm = zeros(size(sigma));
    for k = 1:numel(sigma)
        n_fwhm(k) = solve_n_from_sigma_fwhm_2011(sigma(k));
    end

    Z_heur = rebuild_series_zarc_2011(freq, R_inf, Rvec, Cvec, n_heur);
    Z_fwhm = rebuild_series_zarc_2011(freq, R_inf, Rvec, Cvec, n_fwhm);

    resid_heur = abs(Z_exp - Z_heur);
    resid_fwhm = abs(Z_exp - Z_fwhm);

    mean_heur = mean(resid_heur);
    max_heur = max(resid_heur);
    mean_fwhm = mean(resid_fwhm);
    max_fwhm = max(resid_fwhm);

    if mean_heur <= mean_fwhm
        best = 'heuristic';
        improve_pct = 100 * (mean_fwhm - mean_heur) / max(mean_fwhm, eps);
    else
        best = 'fwhm';
        improve_pct = 100 * (mean_heur - mean_fwhm) / max(mean_heur, eps);
    end

    ridx = ridx + 1;
    results(ridx).temperature = T; %#ok<AGROW>
    results(ridx).n_peaks = numel(Rvec); %#ok<AGROW>
    results(ridx).mean_abs_resid_heur = mean_heur; %#ok<AGROW>
    results(ridx).max_abs_resid_heur = max_heur; %#ok<AGROW>
    results(ridx).mean_abs_resid_fwhm = mean_fwhm; %#ok<AGROW>
    results(ridx).max_abs_resid_fwhm = max_fwhm; %#ok<AGROW>
    results(ridx).best_formula = best; %#ok<AGROW>
    results(ridx).improvement_pct = improve_pct; %#ok<AGROW>

    fprintf(fid, '%.1f\t%d\t%.6g\t%.6g\t%s\t%.2f\n', T, numel(Rvec), mean_heur, mean_fwhm, best, improve_pct);
end

if isempty(results)
    fprintf(fid, '\nNo comparable PEAK/NYQRCN temperature sheets were found.\n');
    fclose(fid);
    error('No comparable temperature sheets found.');
end

mean_all_heur = mean([results.mean_abs_resid_heur]);
mean_all_fwhm = mean([results.mean_abs_resid_fwhm]);

if mean_all_heur <= mean_all_fwhm
    global_best = 'heuristic';
    global_improve = 100 * (mean_all_fwhm - mean_all_heur) / max(mean_all_fwhm, eps);
else
    global_best = 'fwhm';
    global_improve = 100 * (mean_all_heur - mean_all_fwhm) / max(mean_all_heur, eps);
end

fprintf(fid, '\nGlobal mean residual (heuristic): %.6g\n', mean_all_heur);
fprintf(fid, 'Global mean residual (fwhm): %.6g\n', mean_all_fwhm);
fprintf(fid, 'Global best formula: %s (improvement %.2f%%)\n', global_best, global_improve);

fclose(fid);

fprintf('\nSaved comparison report: %s\n', out_file);
fprintf('Global mean residual heuristic: %.6g\n', mean_all_heur);
fprintf('Global mean residual fwhm     : %.6g\n', mean_all_fwhm);
fprintf('Global best formula           : %s\n', global_best);

end

function Z = rebuild_series_zarc_2011(freq, R_inf, Rvec, Cvec, nvec)
omega = 2 * pi * freq(:);
Z = R_inf + zeros(size(omega));
for k = 1:numel(Rvec)
    Rk = Rvec(k);
    Ck = Cvec(k);
    nk = nvec(k);
    if ~(Rk > 0) || ~(Ck > 0) || ~(nk > 0)
        continue;
    end
    Z = Z + Rk ./ (1 + (1i * omega * Rk * Ck) .^ nk);
end
end

function n = solve_n_from_sigma_fwhm_2011(sig)
if ~(sig > 0)
    n = 1;
    return;
end

obj = @(x) (sqrt(2*log(2))*sig*x - acosh(2 + cos(pi*x)))^2;
n = fminbnd(obj, 0.05, 0.999);
end

function T = parse_temp_from_sheet_2011(sh)
T = NaN;
if numel(sh) < 6
    return;
end
u = upper(sh);
if ~strcmp(u(1:5), 'PEAK_')
    return;
end
tok = sh(6:end);
tok = strrep(tok, 'K', '');
tok = strrep(tok, 'k', '');
tok = strrep(tok, 'p', '.');
tok = strrep(tok, 'm', '-');
val = str2double(tok);
if ~isnan(val)
    T = val;
end
end

function sheet_name = make_sheet_name_2011_local(prefix, temperature_value)
sheet_name = sprintf('%s_%0.1fK', prefix, temperature_value);
sheet_name = strrep(sheet_name, '.', 'p');
sheet_name = strrep(sheet_name, '-', 'm');
if numel(sheet_name) > 31
    sheet_name = sheet_name(1:31);
end
end
