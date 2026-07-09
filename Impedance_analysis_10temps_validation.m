clc;
clear all;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
[data, text, raw] = xlsread(input_path); %#ok<ASGLU>
text = char(text(1,:));
[p1, q] = size(text); %#ok<NASGU>
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

sample_targets = linspace(2, 280, 10);
selected_idx = zeros(numel(sample_targets), 1);
used = zeros(numel(temperature), 1);
for k = 1:numel(sample_targets)
    diff_vec = abs(temperature(:) - sample_targets(k));
    diff_vec(logical(used)) = inf;
    [~, idx] = min(diff_vec);
    selected_idx(k) = idx;
    used(idx) = 1;
end
selected_idx = sort(selected_idx);
selected_temperature = temperature(selected_idx);

selected_data = zeros(size(data, 1), 3 * numel(selected_idx));
for k = 1:numel(selected_idx)
    src_cols = (3 * selected_idx(k) - 2):(3 * selected_idx(k));
    dst_cols = (3 * k - 2):(3 * k);
    selected_data(:, dst_cols) = data(:, src_cols);
end

fprintf('Running 10-temperature validation on selected temperatures (K):\n');
fprintf('  %.2f', selected_temperature(1));
for k = 2:numel(selected_temperature)
    fprintf(', %.2f', selected_temperature(k));
end
fprintf('\n');

drt_input_analysis_matlab2011(repo_root, input_path, selected_temperature, selected_data);
