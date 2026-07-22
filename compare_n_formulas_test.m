%% compare_n_formulas_test.m
% Direct comparison of two n-estimation methods on S2022 peak data.
% Tests: (1) heuristic n = 1.144 / (1 + sigma), (2) FWHM-based implicit.

clear all; close all;

% Load generated peak results from prior S2022 run.
peak_file = 'S2022sap_peak_detail_result.xlsx';
plot_data_file = 'S2022sap_temperature_plot_data.xlsx';

if ~isfile(peak_file) || ~isfile(plot_data_file)
    error('Missing input files: %s or %s', peak_file, plot_data_file);
end

try
    [peak_numeric, peak_text] = xlsread(peak_file);
    [plot_data_numeric, plot_data_text] = xlsread(plot_data_file);
catch ME
    fprintf('Error reading Excel files: %s\n', ME.message);
    return;
end

% Extract peak data: assume columns after header are
% [peak_id, tau_s, amplitude, sigma_logtau, fwhm_logtau, R_ohm, C_F, n_est, is_hidden]
if isempty(peak_numeric) || size(peak_numeric, 2) < 8
    error('Peak data has unexpected structure.');
end

sigma_vals = peak_numeric(:, 4);  % sigma_logtau from Gaussian fit
tau_vals = peak_numeric(:, 2);    % tau_s from peak
amplitude_vals = peak_numeric(:, 3);
fwhm_logtau_vals = peak_numeric(:, 5);
R_original = peak_numeric(:, 6);
C_original = peak_numeric(:, 7);
n_heuristic = peak_numeric(:, 8);  % Current formula

fprintf('\n=== N-Formula Comparison Test ===\n');
fprintf('Peak data loaded: %d peaks, %d sigma values\n', size(peak_numeric, 1), numel(sigma_vals));

% Compute alternative n values from FWHM-based implicit formula.
n_fwhm_based = zeros(size(sigma_vals));
for i = 1:numel(sigma_vals)
    sig = sigma_vals(i);
    if sig > 0
        % FWHM_logtau = 2*sqrt(2*log(2))*sigma
        % FWHM_logtau = (2/n)*acosh(2 + cos(pi*n))
        % Rearrange: sqrt(2*log(2))*sigma*n = acosh(2 + cos(pi*n))
        % Solve numerically.
        fwhm_tgt = 2 * sqrt(2 * log(2)) * sig;
        n_fwhm_based(i) = solve_n_from_fwhm_implicit(fwhm_tgt);
    else
        n_fwhm_based(i) = 1.0;  % Default if sigma is zero
    end
end

% Display side-by-side comparison.
fprintf('\nPeak-by-peak n estimates:\n');
fprintf('Peak |  Sigma  | N_heur  | N_fwhm  | Tau (s) | Amp (Ohm) |\n');
fprintf('----------------------------------------------------------\n');
for i = 1:min(numel(sigma_vals), 10)
    fprintf('%4d | %7.4f | %7.4f | %7.4f | %.2e | %.2e |\n', ...
        i, sigma_vals(i), n_heuristic(i), n_fwhm_based(i), tau_vals(i), amplitude_vals(i));
end

% Load frequency and impedance for Nyquist calculation.
% Assume plot_data has sheets like NYQ_280pK containing [freq_Hz, Zre_exp, Zim_exp, ...]
% For simplicity, extract first temperature sheet.
[nyq_data_numeric, nyq_text] = xlsread(plot_data_file, 'NYQ_280p0K');
if isempty(nyq_data_numeric)
    fprintf('Could not load NYQ_280p0K sheet; using available sheet list.\n');
    % Try to find any NYQ sheet.
    [sheet_names] = xlsfinfo(plot_data_file);
    nyq_sheets = sheet_names(startsWith(sheet_names, 'NYQ'));
    if ~isempty(nyq_sheets)
        nyq_sheet = nyq_sheets{1};
        [nyq_data_numeric] = xlsread(plot_data_file, nyq_sheet);
        fprintf('Using sheet: %s\n', nyq_sheet);
    else
        fprintf('No NYQ sheets found in workbook.\n');
        return;
    end
end

freq_vec = nyq_data_numeric(:, 1);
Z_re_exp = nyq_data_numeric(:, 2);
Z_im_exp = nyq_data_numeric(:, 3);
Z_exp = Z_re_exp + 1i * Z_im_exp;

fprintf('\nNyquist data loaded: %d frequencies\n', numel(freq_vec));

% Reconstruct Nyquist with heuristic n.
Z_rcn_heur = reconstruct_rcn_impedance(freq_vec, R_original, C_original, n_heuristic);
resid_heur = abs(Z_exp - Z_rcn_heur);
residual_mean_heur = mean(resid_heur);
residual_max_heur = max(resid_heur);
residual_rms_heur = sqrt(mean(resid_heur.^2));
residual_rel_pct_heur = 100 * mean(resid_heur ./ max(abs(Z_exp), eps));

% Reconstruct Nyquist with FWHM-based n.
Z_rcn_fwhm = reconstruct_rcn_impedance(freq_vec, R_original, C_original, n_fwhm_based);
resid_fwhm = abs(Z_exp - Z_rcn_fwhm);
residual_mean_fwhm = mean(resid_fwhm);
residual_max_fwhm = max(resid_fwhm);
residual_rms_fwhm = sqrt(mean(resid_fwhm.^2));
residual_rel_pct_fwhm = 100 * mean(resid_fwhm ./ max(abs(Z_exp), eps));

fprintf('\n=== Nyquist Residual Comparison ===\n');
fprintf('Metric                   | Heuristic     | FWHM-based    | Winner\n');
fprintf('-------------------------------------------------------------\n');

if residual_mean_heur < residual_mean_fwhm
    winner_mean = 'Heuristic';
else
    winner_mean = 'FWHM';
end
fprintf('Mean residual (Ohm)      | %.6e | %.6e | %s\n', residual_mean_heur, residual_mean_fwhm, winner_mean);

if residual_max_heur < residual_max_fwhm
    winner_max = 'Heuristic';
else
    winner_max = 'FWHM';
end
fprintf('Max residual (Ohm)       | %.6e | %.6e | %s\n', residual_max_heur, residual_max_fwhm, winner_max);

if residual_rms_heur < residual_rms_fwhm
    winner_rms = 'Heuristic';
else
    winner_rms = 'FWHM';
end
fprintf('RMS residual (Ohm)       | %.6e | %.6e | %s\n', residual_rms_heur, residual_rms_fwhm, winner_rms);

if residual_rel_pct_heur < residual_rel_pct_fwhm
    winner_rel = 'Heuristic';
else
    winner_rel = 'FWHM';
end
fprintf('Mean relative (percent)  | %.2f %%        | %.2f %%        | %s\n', residual_rel_pct_heur, residual_rel_pct_fwhm, winner_rel);

fprintf('\nDelta (Heuristic - FWHM):\n');
fprintf('  Mean:  %.6e Ohm\n', residual_mean_heur - residual_mean_fwhm);
fprintf('  Max:   %.6e Ohm\n', residual_max_heur - residual_max_fwhm);
fprintf('  RMS:   %.6e Ohm\n', residual_rms_heur - residual_rms_fwhm);
fprintf('  Rel:   %.2f %%\n', residual_rel_pct_heur - residual_rel_pct_fwhm);

% Plot comparison.
figure('Color', 'w', 'Name', 'N-formula residual comparison');
subplot(2,2,1);
plot(abs(Z_exp), resid_heur, 'r.', 'MarkerSize', 4); hold on;
plot(abs(Z_exp), resid_fwhm, 'b.', 'MarkerSize', 4);
xlabel('|Z_exp| (Ohm)');
ylabel('Residual (Ohm)');
legend('Heuristic', 'FWHM-based');
grid on;
title('Residual vs impedance magnitude');

subplot(2,2,2);
semilogx(freq_vec, resid_heur, 'r.-'); hold on;
semilogx(freq_vec, resid_fwhm, 'b.-');
xlabel('Frequency (Hz)');
ylabel('Residual (Ohm)');
legend('Heuristic', 'FWHM-based');
grid on;
title('Residual vs frequency');

subplot(2,2,3);
plot(real(Z_exp), -imag(Z_exp), 'k.', 'MarkerSize', 3); hold on;
plot(real(Z_rcn_heur), -imag(Z_rcn_heur), 'r--', 'LineWidth', 1.2);
plot(real(Z_rcn_fwhm), -imag(Z_rcn_fwhm), 'b--', 'LineWidth', 1.2);
axis equal;
xlabel('Re(Z) (Ohm)');
ylabel('-Im(Z) (Ohm)');
legend('Measured', 'Heuristic', 'FWHM-based');
grid on;
title('Nyquist overlay');

subplot(2,2,4);
bar([residual_mean_heur, residual_mean_fwhm; residual_rms_heur, residual_rms_fwhm]);
set(gca, 'XTickLabel', {'Mean', 'RMS'});
legend('Heuristic', 'FWHM-based');
ylabel('Residual (Ohm)');
title('Residual metrics comparison');
grid on;

saveas(gcf, 'n_formula_comparison_plot.png');
fprintf('\nPlot saved as n_formula_comparison_plot.png\n');

fprintf('\n=== Conclusion ===\n');
fprintf('Based on %d frequency points:\n', numel(freq_vec));
fprintf('  Heuristic formula (n = 1.144/(1+sigma)) mean residual: %.6e Ohm\n', residual_mean_heur);
fprintf('  FWHM-based formula mean residual:                     %.6e Ohm\n', residual_mean_fwhm);
if residual_mean_heur < residual_mean_fwhm
    pct_improvement = 100 * (residual_mean_fwhm - residual_mean_heur) / residual_mean_fwhm;
    fprintf('  Winner: HEURISTIC (%.1f%% better)\n', pct_improvement);
else
    pct_improvement = 100 * (residual_mean_heur - residual_mean_fwhm) / residual_mean_heur;
    fprintf('  Winner: FWHM-BASED (%.1f%% better)\n', pct_improvement);
end

end

function n_val = solve_n_from_fwhm_implicit(fwhm_target)
% Solve: sqrt(2*log(2))*sigma*n = acosh(2 + cos(pi*n)) for n in (0, 1].
% Equivalently: (2/n)*acosh(2 + cos(pi*n)) = fwhm_target.

    objective = @(n) abs((2/n) * acosh(2 + cos(pi*n)) - fwhm_target);
    
    % Search in valid range (0, 1).
    n_candidates = linspace(0.01, 0.99, 100);
    errors = arrayfun(objective, n_candidates);
    [~, best_idx] = min(errors);
    n_init = n_candidates(best_idx);
    
    % Refine with fminsearch.
    opts = optimset('Display', 'off', 'TolX', 1e-8, 'MaxIter', 500);
    n_val = fminsearch(objective, n_init, opts);
    
    % Clamp to valid range.
    n_val = max(0.01, min(n_val, 0.99));
end

function Z_rcn = reconstruct_rcn_impedance(freq_vec, R_vec, C_vec, n_vec)
% Reconstruct series sum of R||CPE branches: Z = sum_k R_k / (1 + (j*omega*R_k*C_k)^n_k).

    omega = 2 * pi * freq_vec(:);
    Z_rcn = zeros(size(omega));
    
    for k = 1:numel(R_vec)
        Rk = R_vec(k);
        Ck = C_vec(k);
        nk = n_vec(k);
        
        if Rk > 0 && Ck > 0 && nk > 0
            Zk = Rk ./ (1 + (1i * omega * Rk * Ck) .^ nk);
            Z_rcn = Z_rcn + Zk;
        end
    end
end
