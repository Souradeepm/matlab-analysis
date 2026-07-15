%% Simplified Peak Analysis and RGB Visualization
% Extract peaks from all datasets and create color-coded visualizations

fprintf('Peak Frequency Analysis and RGB Visualization\n');
fprintf('=============================================\n\n');

% Define datasets
datasets = {
    'S2022Sap', 's2022sap';
    'S2022Al', 's2022al';
    'S2222Sap', 's2222sap';
    'S2222Al', 's2222al';
    'S2302Sap', 's2302sap';
    'S2302Al', 's2302al';
    'S2322Sap', 's2322sap';
    'S2322Al', 's2322al';
    'S2332Sap', 's2332sap';
    'S2422Sap', 's2422sap';
    'S2422Al', 's2422al'
};

methods = {'Bayes-DRT', 'Paper Method', 'Residual Method'};

% Create output directories
if ~isfolder('peak_visualizations')
    mkdir('peak_visualizations');
end

% Color palette for peak groups (RGB)
color_palette = [
    1.0, 0.0, 0.0;   % Red - 1 peak
    0.0, 1.0, 0.0;   % Green - 2 peaks
    0.0, 0.0, 1.0;   % Blue - 3 peaks
    1.0, 1.0, 0.0;   % Yellow - 4 peaks
    1.0, 0.0, 1.0;   % Magenta - 5 peaks
    0.0, 1.0, 1.0;   % Cyan - 6 peaks
    1.0, 0.5, 0.0;   % Orange - 7+ peaks
];

% Collect peak statistics
peak_stats = {};
idx = 1;

fprintf('Collecting peak data from all batch files...\n\n');

% Scan all datasets
for d = 1:size(datasets, 1)
    ds_name = datasets{d, 1};
    ds_tag = datasets{d, 2};
    
    % Process Bayes-DRT
    bayes_files = dir([ds_tag '_bayes_drt_matlab2011_b*.txt']);
    for f = 1:numel(bayes_files)
        [temp_arr, peak_arr] = extract_temps_and_peaks(bayes_files(f).name);
        if ~isempty(temp_arr)
            peak_stats{idx, 1} = ds_name;
            peak_stats{idx, 2} = 'Bayes-DRT';
            peak_stats{idx, 3} = temp_arr;
            peak_stats{idx, 4} = peak_arr;
            idx = idx + 1;
        end
    end
    
    % Process Paper Method
    paper_files = dir([ds_tag '_paper_method_b*.txt']);
    for f = 1:numel(paper_files)
        [temp_arr, peak_arr] = extract_temps_and_peaks(paper_files(f).name);
        if ~isempty(temp_arr)
            peak_stats{idx, 1} = ds_name;
            peak_stats{idx, 2} = 'Paper Method';
            peak_stats{idx, 3} = temp_arr;
            peak_stats{idx, 4} = peak_arr;
            idx = idx + 1;
        end
    end
    
    % Process Residual Method
    resid_files = dir([ds_tag '_residual_method_b*.txt']);
    for f = 1:numel(resid_files)
        [temp_arr, peak_arr] = extract_temps_and_peaks(resid_files(f).name);
        if ~isempty(temp_arr)
            peak_stats{idx, 1} = ds_name;
            peak_stats{idx, 2} = 'Residual Method';
            peak_stats{idx, 3} = temp_arr;
            peak_stats{idx, 4} = peak_arr;
            idx = idx + 1;
        end
    end
end

% Identify peak groups
all_peak_counts = [];
for i = 1:size(peak_stats, 1)
    all_peak_counts = [all_peak_counts; peak_stats{i, 4}(:)]; %#ok<AGROW>
end
unique_peak_counts = sort(unique(all_peak_counts));

fprintf('Peak Groups Found:\n');
fprintf('------------------\n');
for i = 1:numel(unique_peak_counts)
    pc = unique_peak_counts(i);
    count = sum(all_peak_counts == pc);
    if i <= size(color_palette, 1)
        rgb = color_palette(i, :);
        fprintf('Group %d: %d peaks (RGB: %.1f, %.1f, %.1f) - %d instances\n', ...
            i, pc, rgb(1), rgb(2), rgb(3), count);
    else
        fprintf('Group %d: %d peaks - %d instances\n', i, pc, count);
    end
end

% Generate visualizations for each dataset-method
fprintf('\n\nGenerating visualizations...\n');

for i = 1:size(peak_stats, 1)
    ds_name = peak_stats{i, 1};
    method = peak_stats{i, 2};
    temps = peak_stats{i, 3};
    peaks = peak_stats{i, 4};
    
    fprintf('  %s - %s: %d temperatures\n', ds_name, method, numel(temps));
    
    % Create RGB image
    img_height = numel(peaks);
    img_width = numel(unique_peak_counts) + 2;
    rgb_img = ones(img_height, img_width, 3) * 0.95;  % Light background
    
    % Color each row by peak group
    for t = 1:img_height
        pc = peaks(t);
        group_idx = find(unique_peak_counts == pc);
        
        if ~isempty(group_idx) && group_idx(1) <= size(color_palette, 1)
            rgb = color_palette(group_idx(1), :);
            rgb_img(t, :, 1) = rgb(1);
            rgb_img(t, :, 2) = rgb(2);
            rgb_img(t, :, 3) = rgb(3);
        end
    end
    
    % Save image
    fig = figure('Visible', 'off');
    imagesc(rgb_img);
    axis image;
    xlabel('Peak Group Channel');
    ylabel(sprintf('Temperature (%d)', numel(temps)));
    title(sprintf('%s - %s\nPeak Distribution', ds_name, method));
    
    filename = sprintf('peak_visualizations/%s_%s_peak_dist.png', ...
        lower(strrep(ds_name, ' ', '_')), lower(strrep(method, ' ', '_')));
    saveas(fig, filename);
    close(fig);
end

% Create summary statistics file
fid = fopen('peak_analysis_summary.txt', 'w');
fprintf(fid, 'PEAK FREQUENCY ANALYSIS AND RGB VISUALIZATION\n');
fprintf(fid, '=============================================\n\n');
fprintf(fid, 'Peak Groups and RGB Mapping:\n');
fprintf(fid, '---------------------------\n');

for i = 1:numel(unique_peak_counts)
    pc = unique_peak_counts(i);
    count = sum(all_peak_counts == pc);
    if i <= size(color_palette, 1)
        rgb = color_palette(i, :);
        fprintf(fid, 'Group %d: %d-Peak Process\n', i, pc);
        fprintf(fid, '  RGB Color: (%.1f, %.1f, %.1f)\n', rgb(1)*100, rgb(2)*100, rgb(3)*100);
        fprintf(fid, '  Occurrences: %d\n', count);
    else
        fprintf(fid, 'Group %d: %d-Peak Process (occurrences: %d)\n', i, pc, count);
    end
    fprintf(fid, '\n');
end

fprintf(fid, '\n\nDataset-Method Peak Profiles:\n');
fprintf(fid, '----------------------------\n');

unique_datasets = unique(peak_stats(:, 1));
for d = 1:numel(unique_datasets)
    ds_name = unique_datasets{d};
    fprintf(fid, '\n%s:\n', ds_name);
    
    for m = 1:numel(methods)
        method = methods{m};
        mask = strcmp(peak_stats(:, 1), ds_name) & strcmp(peak_stats(:, 2), method);
        rows = find(mask);
        
        if ~isempty(rows)
            fprintf(fid, '  %s:\n', method);
            for r = rows(:)'
                peak_dist = peak_stats{r, 4};
                unique_peaks = unique(peak_dist);
                fprintf(fid, '    Peak distribution: %s\n', strjoin(string(unique_peaks), ', '));
            end
        end
    end
end

fprintf(fid, '\n\nVisualization Guide:\n');
fprintf(fid, '-------------------\n');
fprintf(fid, '- Each row represents one temperature measurement\n');
fprintf(fid, '- Row color indicates the number of detected relaxation processes\n');
fprintf(fid, '- Same color = same peak count\n');
fprintf(fid, '- RGB channels encode peak group information\n');
fprintf(fid, '- Visualizations saved to: peak_visualizations/\n');

fclose(fid);

fprintf('\n\nAnalysis complete!\n');
fprintf('Summary saved to: peak_analysis_summary.txt\n');
fprintf('Visualizations saved to: peak_visualizations/\n');

exit;

%% Helper function to extract temperatures and peak counts
function [temp_arr, peak_arr] = extract_temps_and_peaks(fname)
    temp_arr = [];
    peak_arr = [];
    
    try
        fid = fopen(fname, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        in_detail = false;
        for l = 1:numel(lines{1})
            line = lines{1}{l};
            
            if contains(line, 'TemperatureK,') || contains(line, 'Per-temperature')
                in_detail = true;
                continue;
            end
            
            if in_detail && ~isempty(line) && ~contains(line, 'Temperature')
                parts = strsplit(line, ',');
                if numel(parts) >= 7
                    try
                        t = str2double(parts{1});
                        p = str2double(parts{7});
                        if ~isnan(t) && ~isnan(p) && t > 0 && p > 0
                            temp_arr = [temp_arr; t]; %#ok<AGROW>
                            peak_arr = [peak_arr; p]; %#ok<AGROW>
                        end
                    catch
                    end
                end
            end
        end
    catch
    end
end
