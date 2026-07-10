% Plot tau grid from the random-10 report output (MATLAB 2011 compatible style).

clear;
clc;

csv_path = fullfile(pwd, 's2022_random10_tau_grid.csv');
if ~exist(csv_path, 'file')
    error('Tau CSV not found: %s. Run run_s2022_tau_grid_report first.', csv_path);
end

M = dlmread(csv_path, ',', 1, 0); % skip header row
if isempty(M) || size(M, 2) < 2
    error('Tau CSV has unexpected format: %s', csv_path);
end

tau = M(:, 2);
idx = (1:numel(tau)).';

h = figure('Color', 'w', 'Position', [100 100 900 520]);
semilogy(idx, tau, '-o', 'LineWidth', 1.2, 'MarkerSize', 3);
grid on;
xlabel('Tau Index');
ylabel('Tau (s)');
title('S2022 Random-10 Tau Grid');

png_path = fullfile(pwd, 's2022_random10_tau_grid_plot.png');
saveas(h, png_path);

% Fallback for systems where saveas can fail silently for PNG.
if ~exist(png_path, 'file')
    print(h, '-dpng', '-r200', png_path);
end

fprintf('Wrote tau grid plot: %s\n', png_path);
