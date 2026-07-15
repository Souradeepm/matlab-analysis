%% Peak Frequency Extraction and RGB Visualization
% Extracts peak frequencies from all datasets/methods, groups them, and creates RGB images

fprintf('Peak Frequency Analysis and RGB Visualization\n');
fprintf('=============================================\n\n');

% Define datasets
datasets = {
    'S2022Sap', 's2022sap', 69;
    'S2022Al', 's2022al', 77;
    'S2222Sap', 's2222sap', 61;
    'S2222Al', 's2222al', 112;
    'S2302Sap', 's2302sap', 86;
    'S2302Al', 's2302al', 104;
    'S2322Sap', 's2322sap', 69;
    'S2322Al', 's2322al', 63;
    'S2332Sap', 's2332sap', 220;
    'S2422Sap', 's2422sap', 114;
    'S2422Al', 's2422al', 106
};

methods = {'Bayes-DRT', 'Paper Method', 'Residual Method'};

% Initialize peak database
all_peaks = {};
peak_idx = 1;

fprintf('Extracting peak data from all batch files...\n\n');

% Process each dataset and method
for d = 1:size(datasets, 1)
    ds_name = datasets{d, 1};
    ds_tag = datasets{d, 2};
    n_temps = datasets{d, 3};
    
    fprintf('%s: ', ds_name);
    
    % Extract Bayes-DRT peaks
    bayes_files = dir([ds_tag '_bayes_drt_matlab2011_b*.txt']);
    for f = 1:numel(bayes_files)
        fname = bayes_files(f).name;
        fid = fopen(fname, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        % Count peaks from each temperature
        n_peaks_found = 0;
        for l = 1:numel(lines{1})
            line = lines{1}{l};
            if contains(line, 'Peak count:') && ~contains(line, 'Mean')
                parts = strsplit(line, ':');
                peak_count = str2double(strtrim(parts{2}));
                if ~isnan(peak_count)
                    n_peaks_found = n_peaks_found + peak_count;
                end
            end
        end
        
        if n_peaks_found > 0
            all_peaks{peak_idx, 1} = ds_name;
            all_peaks{peak_idx, 2} = 'Bayes-DRT';
            all_peaks{peak_idx, 3} = n_peaks_found;
            all_peaks{peak_idx, 4} = fname;
            peak_idx = peak_idx + 1;
        end
    end
    
    % Extract Paper Method peaks
    paper_files = dir([ds_tag '_paper_method_b*.txt']);
    for f = 1:numel(paper_files)
        fname = paper_files(f).name;
        fid = fopen(fname, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        n_peaks_found = 0;
        for l = 1:numel(lines{1})
            line = lines{1}{l};
            if contains(line, 'Peak count:') && ~contains(line, 'Mean')
                parts = strsplit(line, ':');
                peak_count = str2double(strtrim(parts{2}));
                if ~isnan(peak_count)
                    n_peaks_found = n_peaks_found + peak_count;
                end
            end
        end
        
        if n_peaks_found > 0
            all_peaks{peak_idx, 1} = ds_name;
            all_peaks{peak_idx, 2} = 'Paper Method';
            all_peaks{peak_idx, 3} = n_peaks_found;
            all_peaks{peak_idx, 4} = fname;
            peak_idx = peak_idx + 1;
        end
    end
    
    % Extract Residual Method peaks
    residual_files = dir([ds_tag '_residual_method_b*.txt']);
    for f = 1:numel(residual_files)
        fname = residual_files(f).name;
        fid = fopen(fname, 'r');
        lines = textscan(fid, '%s', 'Delimiter', '\n');
        fclose(fid);
        
        n_peaks_found = 0;
        for l = 1:numel(lines{1})
            line = lines{1}{l};
            if contains(line, 'Peak count:') && ~contains(line, 'Mean')
                parts = strsplit(line, ':');
                peak_count = str2double(strtrim(parts{2}));
                if ~isnan(peak_count)
                    n_peaks_found = n_peaks_found + peak_count;
                end
            end
        end
        
        if n_peaks_found > 0
            all_peaks{peak_idx, 1} = ds_name;
            all_peaks{peak_idx, 2} = 'Residual Method';
            all_peaks{peak_idx, 3} = n_peaks_found;
            all_peaks{peak_idx, 4} = fname;
            peak_idx = peak_idx + 1;
        end
    end
    
    fprintf('OK\n');
end

% Analyze peak distribution
fprintf('\n\nPeak Distribution Analysis\n');
fprintf('===========================\n\n');

% Find unique peak groups
unique_peak_counts = unique([all_peaks{:, 3}]);
fprintf('Peak count distribution:\n');
for pc = unique_peak_counts
    count = sum([all_peaks{:, 3}] == pc);
    fprintf('  %d peaks: %d instances\n', pc, count);
end

% Create peak grouping
fprintf('\n\nPeak Grouping (by frequency count):\n');
fprintf('===================================\n\n');

peak_groups = {};
group_idx = 1;

for pc = sort(unique_peak_counts)
    indices = find([all_peaks{:, 3}] == pc);
    fprintf('Group %d: %d peaks (%d instances)\n', group_idx, pc, numel(indices));
    
    datasets_in_group = unique(all_peaks(indices, 1));
    methods_in_group = unique(all_peaks(indices, 2));
    
    fprintf('  Datasets: %s\n', strjoin(datasets_in_group, ', '));
    fprintf('  Methods: %s\n\n', strjoin(methods_in_group, ', '));
    
    peak_groups{group_idx, 1} = pc;
    peak_groups{group_idx, 2} = datasets_in_group;
    peak_groups{group_idx, 3} = methods_in_group;
    group_idx = group_idx + 1;
end

% Create RGB mapping
fprintf('\n\nRGB Channel Assignment\n');
fprintf('======================\n\n');

n_groups = size(peak_groups, 1);
colors = distinguishable_colors(n_groups);

fprintf('Peak Group -> Color Mapping:\n');
for g = 1:n_groups
    n_peaks = peak_groups{g, 1};
    rgb = colors(g, :);
    fprintf('  Group %d (%d peaks): R=%.2f, G=%.2f, B=%.2f\n', ...
        g, n_peaks, rgb(1), rgb(2), rgb(3));
end

% Create summary file
fid = fopen('peak_frequency_analysis.txt', 'w');
fprintf(fid, 'PEAK FREQUENCY ANALYSIS AND GROUPING\n');
fprintf(fid, '====================================\n\n');
fprintf(fid, 'Total Datasets: %d\n', size(datasets, 1));
fprintf(fid, 'Total Methods: %d\n', numel(methods));
fprintf(fid, 'Total Peak Records: %d\n\n', size(all_peaks, 1));

fprintf(fid, 'PEAK GROUPS\n');
fprintf(fid, '===========\n\n');
for g = 1:n_groups
    n_peaks = peak_groups{g, 1};
    fprintf(fid, 'Group %d: %d-Peak Relaxation Processes\n', g, n_peaks);
    fprintf(fid, '  Count: %d instances\n', sum([all_peaks{:, 3}] == n_peaks));
    fprintf(fid, '  Color (RGB): (%.3f, %.3f, %.3f)\n', colors(g, 1), colors(g, 2), colors(g, 3));
    fprintf(fid, '\n');
end

fprintf(fid, '\n\nDATASET-METHOD PEAK PROFILES\n');
fprintf(fid, '============================\n\n');
for d = 1:size(datasets, 1)
    ds_name = datasets{d, 1};
    fprintf(fid, '%s:\n', ds_name);
    
    for m = 1:numel(methods)
        method = methods{m};
        mask = strcmp(all_peaks(:, 1), ds_name) & strcmp(all_peaks(:, 2), method);
        entries = all_peaks(mask, :);
        
        if ~isempty(entries)
            peak_counts = unique([entries{:, 3}]);
            fprintf(fid, '  %s: %s peaks\n', method, strjoin(string(peak_counts), ', '));
        else
            fprintf(fid, '  %s: No data\n', method);
        end
    end
    fprintf(fid, '\n');
end

fclose(fid);

fprintf('\n\nAnalysis saved to: peak_frequency_analysis.txt\n');

% Save RGB mapping for visualization
rgb_map_file = 'peak_rgb_mapping.mat';
save(rgb_map_file, 'peak_groups', 'colors', 'all_peaks', 'datasets', 'methods');
fprintf('RGB mapping saved to: %s\n', rgb_map_file);

fprintf('\nPeak analysis complete!\n');
exit;

%% Helper function to create distinguishable colors
function colors = distinguishable_colors(N)
    % Create N distinguishable colors using predefined palette
    if N <= 1
        colors = [0, 0, 0];
        return;
    end
    
    % Predefined colors for peak groups
    base_colors = [
        1.0, 0.0, 0.0;   % Red
        0.0, 1.0, 0.0;   % Green
        0.0, 0.0, 1.0;   % Blue
        1.0, 1.0, 0.0;   % Yellow
        1.0, 0.0, 1.0;   % Magenta
        0.0, 1.0, 1.0;   % Cyan
        1.0, 0.5, 0.0;   % Orange
        0.5, 0.0, 1.0;   % Purple
        1.0, 0.0, 0.5;   % Pink
        0.0, 0.5, 1.0;   % Sky Blue
    ];
    
    if N <= size(base_colors, 1)
        colors = base_colors(1:N, :);
    else
        % Generate additional colors using HSV
        colors = zeros(N, 3);
        for i = 1:min(N, size(base_colors, 1))
            colors(i, :) = base_colors(i, :);
        end
        
        for i = (size(base_colors, 1) + 1):N
            h = (i - size(base_colors, 1)) / (N - size(base_colors, 1));
            colors(i, :) = [sin(h*pi), cos(h*pi/2), 0.7];
            colors(i, :) = abs(colors(i, :));
        end
    end
end
