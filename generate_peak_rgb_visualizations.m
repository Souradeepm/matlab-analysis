%% Generate RGB Visualizations for Peak Frequencies
% Creates color-coded images showing peak distributions across temperatures

fprintf('Generating RGB Peak Visualizations\n');
fprintf('===================================\n\n');

% Load peak mapping data
load('peak_rgb_mapping.mat', 'peak_groups', 'colors', 'all_peaks', 'datasets', 'methods');

% Create output directory
if ~isfolder('peak_visualizations')
    mkdir('peak_visualizations');
end

n_groups = size(peak_groups, 1);

% Process each dataset
for d = 1:size(datasets, 1)
    ds_name = datasets{d, 1};
    ds_tag = datasets{d, 2};
    n_temps = datasets{d, 3};
    
    fprintf('Processing %s (%d temperatures)...\n', ds_name, n_temps);
    
    % Create figure for this dataset showing all three methods
    fig = figure('Position', [100, 100, 1200, 400]);
    
    % Process each method
    for m = 1:numel(methods)
        method = methods{m};
        
        % Find batch files for this dataset-method combination
        if strcmp(method, 'Bayes-DRT')
            files = dir([ds_tag '_bayes_drt_matlab2011_b*.txt']);
        elseif strcmp(method, 'Paper Method')
            files = dir([ds_tag '_paper_method_b*.txt']);
        else  % Residual Method
            files = dir([ds_tag '_residual_method_b*.txt']);
        end
        
        if isempty(files)
            continue;
        end
        
        % Read temperature and peak data
        temps = [];
        peak_counts = [];
        
        for f = 1:numel(files)
            fname = files(f).name;
            fid = fopen(fname, 'r');
            lines = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
            
            % Parse per-temperature data
            in_detail_section = false;
            for l = 1:numel(lines{1})
                line = lines{1}{l};
                
                if contains(line, 'Per-temperature detail') || ...
                   contains(line, 'TemperatureK,Selected')
                    in_detail_section = true;
                    continue;
                end
                
                if in_detail_section && ~isempty(line) && ...
                   ~contains(line, 'Temperature') && ~contains(line, 'Mean')
                    parts = strsplit(line, ',');
                    if numel(parts) >= 7
                        try
                            temp = str2double(parts{1});
                            peak_count = str2double(parts{7});
                            if ~isnan(temp) && ~isnan(peak_count)
                                temps = [temps; temp]; %#ok<AGROW>
                                peak_counts = [peak_counts; peak_count]; %#ok<AGROW>
                            end
                        catch
                        end
                    end
                end
            end
        end
        
        if isempty(peak_counts)
            continue;
        end
        
        % Create RGB image: each row is a temperature, colored by peak group
        n_temps_found = numel(peak_counts);
        img_width = n_groups + 1;  % Width = number of peak groups + padding
        img_height = n_temps_found;
        rgb_img = ones(img_height, img_width, 3);  % White background
        
        % Map peak counts to RGB
        for t = 1:n_temps_found
            pc = peak_counts(t);
            
            % Find which group this peak count belongs to
            group_idx = find([peak_groups{:, 1}] == pc);
            
            if ~isempty(group_idx)
                group_idx = group_idx(1);
                % Assign RGB color to this row
                rgb_img(t, :, 1) = colors(group_idx, 1);
                rgb_img(t, :, 2) = colors(group_idx, 2);
                rgb_img(t, :, 3) = colors(group_idx, 3);
            else
                % Gray for unknown peak count
                rgb_img(t, :, 1) = 0.5;
                rgb_img(t, :, 2) = 0.5;
                rgb_img(t, :, 3) = 0.5;
            end
        end
        
        % Plot in subplot
        subplot(1, 3, m);
        imagesc(rgb_img);
        axis image;
        set(gca, 'YTick', 1:min(n_temps_found, 10));
        xlabel('Peak Group');
        ylabel('Temperature Index');
        title(sprintf('%s\n%s', ds_name, method));
        
        % Save individual visualization
        fig_method = figure('Position', [100, 100, 400, 600]);
        imagesc(rgb_img);
        axis image;
        ylabel('Temperature Index (sorted)');
        title(sprintf('%s - %s\nPeak Distribution', ds_name, method));
        colorbar;
        
        % Add peak count legend
        legend_text = {};
        for g = 1:n_groups
            legend_text{g} = sprintf('Group %d: %d peaks', g, peak_groups{g, 1});
        end
        
        saveas(fig_method, sprintf('peak_visualizations/%s_%s_peak_distribution.png', ...
            lower(strrep(ds_name, ' ', '_')), lower(strrep(method, ' ', '_'))));
        close(fig_method);
        
        fprintf('  %s: %d temperatures mapped\n', method, n_temps_found);
    end
    
    % Save combined figure
    saveas(fig, sprintf('peak_visualizations/%s_all_methods_peaks.png', ...
        lower(strrep(ds_name, ' ', '_'))));
    close(fig);
    
    fprintf('  Saved: %s_all_methods_peaks.png\n\n', lower(strrep(ds_name, ' ', '_')));
end

% Create master legend figure
fig_legend = figure('Position', [100, 100, 400, 300]);
hold on;

for g = 1:n_groups
    n_peaks = peak_groups{g, 1};
    rgb = colors(g, :);
    
    % Create color patch
    patch([0.1, 0.3, 0.3, 0.1], [g-0.4, g-0.4, g+0.4, g+0.4], ...
        reshape(rgb, 1, 1, 3), 'EdgeColor', 'black');
    
    % Add label
    text(0.4, g, sprintf('Group %d: %d-Peak Processes', g, n_peaks), ...
        'VerticalAlignment', 'middle', 'FontSize', 10);
end

axis([0, 1, 0.5, n_groups+0.5]);
axis off;
title('Peak Group Color Legend');

saveas(fig_legend, 'peak_visualizations/peak_group_legend.png');
close(fig_legend);

fprintf('\nRGB visualizations complete!\n');
fprintf('All visualizations saved to: peak_visualizations/\n');

% Create summary report
fid = fopen('peak_rgb_visualization_summary.txt', 'w');
fprintf(fid, 'PEAK FREQUENCY RGB VISUALIZATION SUMMARY\n');
fprintf(fid, '========================================\n\n');
fprintf(fid, 'Color Mapping (RGB Channels):\n');
fprintf(fid, '----------------------------\n');
for g = 1:n_groups
    n_peaks = peak_groups{g, 1};
    rgb = colors(g, :);
    fprintf(fid, 'Group %d (%d peaks): RGB(%.1f%%, %.1f%%, %.1f%%)\n', ...
        g, n_peaks, rgb(1)*100, rgb(2)*100, rgb(3)*100);
end

fprintf(fid, '\n\nVisualization Format:\n');
fprintf(fid, '-------------------\n');
fprintf(fid, 'Each image row represents one temperature\n');
fprintf(fid, 'Each row is colored according to the peak group\n');
fprintf(fid, 'Color intensity corresponds to RGB channel assignments\n');
fprintf(fid, 'Rows with same color have same number of detected peaks\n\n');

fprintf(fid, 'Datasets Processed: %d\n', size(datasets, 1));
fprintf(fid, 'Methods per Dataset: %d\n', numel(methods));
fprintf(fid, 'Total Peak Groups: %d\n', n_groups);

fclose(fid);

fprintf('Summary report saved to: peak_rgb_visualization_summary.txt\n');

exit;
