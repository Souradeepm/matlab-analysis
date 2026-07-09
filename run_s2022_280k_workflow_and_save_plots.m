function run_s2022_280k_workflow_and_save_plots()
% Run the full DRT workflow for the S2022 slice nearest 280 K and save plots.

clc;
clear all;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
plot_dir = fullfile(repo_root, 'plots_2011', 's2022_280k');

if exist(plot_dir, 'dir') ~= 7
    mkdir(plot_dir);
end

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

[~, idx] = min(abs(temperature - 280));
selected_temperature = temperature(idx);
selected_data = data(:, (3 * idx - 2):(3 * idx));

fprintf('Running S2022 workflow at selected temperature %.2f K\n', selected_temperature);
drt_input_analysis_matlab2011(repo_root, input_path, selected_temperature, selected_data);

figs = sort(findall(0, 'Type', 'figure'));
for k = 1:numel(figs)
    fig = figs(k);
    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        fig_name = sprintf('figure_%02d', k);
    end
    fig_name = strrep(fig_name, ' ', '_');
    fig_name = strrep(fig_name, '.', 'p');
    fig_name = strrep(fig_name, '/', '_');
    fig_name = strrep(fig_name, '\\', '_');
    out_png = fullfile(plot_dir, sprintf('%02d_%s.png', k, fig_name));
    saveas(fig, out_png);
    fprintf('Saved plot: %s\n', out_png);
end
end
