function run_s2022_tau_grid_report()
% run_s2022_tau_grid_report
% Report tau grids for the random-selected S2022 temperatures used in the workflow.

clc;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2022Sap.xlsx');
out_txt = fullfile(repo_root, 's2022_random10_tau_grid_report.txt');
out_csv = fullfile(repo_root, 's2022_random10_tau_grid.csv');

rng_seed = 20260710;
n_select = 10;

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text] = xlsread(input_path);
temperature = parse_temperature_labels_local(text, size(data, 2));

rand('twister', rng_seed);
randn('state', rng_seed + 1);
selected_idx = randperm(numel(temperature));
selected_idx = selected_idx(1:n_select);
[selected_temperature, order] = sort(temperature(selected_idx));
selected_idx = selected_idx(order);

tau_ref = [];
is_same_grid = true;

fid_txt = fopen(out_txt, 'w');
if fid_txt < 0
    error('Cannot open report file: %s', out_txt);
end

fid_csv = fopen(out_csv, 'w');
if fid_csv < 0
    fclose(fid_txt);
    error('Cannot open CSV file: %s', out_csv);
end

fprintf(fid_txt, 'S2022 random-10 tau grid report\n');
fprintf(fid_txt, 'Generated: %s\n', datestr(now));
fprintf(fid_txt, 'Random seed: %d\n\n', rng_seed);
fprintf(fid_csv, 'TemperatureK,PointIndex,Tau_s\n');

for t = 1:n_select
    src_cols = (3 * selected_idx(t) - 2):(3 * selected_idx(t));
    sel = data(:, src_cols);

    freq_vec = sel(:,1);
    mag = sel(:,2);
    phase_deg = sel(:,3);
    valid = isfinite(freq_vec) & isfinite(mag) & isfinite(phase_deg) & (freq_vec > 0);
    freq_vec = freq_vec(valid);

    [freq_vec, ~] = sort(freq_vec);
    tau = 1 ./ (2 * pi * freq_vec);

    if isempty(tau_ref)
        tau_ref = tau;
    else
        if numel(tau_ref) ~= numel(tau) || any(abs(tau_ref - tau) > 1e-12)
            is_same_grid = false;
        end
    end

    fprintf(fid_txt, 'Temperature %.6f K\n', selected_temperature(t));
    fprintf(fid_txt, '  N points: %d\n', numel(tau));
    fprintf(fid_txt, '  tau min : %.8e s\n', min(tau));
    fprintf(fid_txt, '  tau max : %.8e s\n', max(tau));
    fprintf(fid_txt, '  tau[1:5]: %.8e, %.8e, %.8e, %.8e, %.8e\n\n', tau(1), tau(2), tau(3), tau(4), tau(5));

    for i = 1:numel(tau)
        fprintf(fid_csv, '%.6f,%d,%.12e\n', selected_temperature(t), i, tau(i));
    end
end

fprintf(fid_txt, 'All selected temperatures share identical tau grid: %d\n', is_same_grid);
if is_same_grid
    fprintf(fid_txt, '\nCommon tau grid preview (first 10 points):\n');
    for i = 1:min(10, numel(tau_ref))
        fprintf(fid_txt, '  %3d: %.12e\n', i, tau_ref(i));
    end
end

fclose(fid_txt);
fclose(fid_csv);

fprintf('Wrote tau report: %s\n', out_txt);
fprintf('Wrote tau CSV   : %s\n', out_csv);
end

function temperature = parse_temperature_labels_local(text, n_cols)
if iscell(text)
    header = text(1, :);
else
    header = cellstr(char(text(1, :)));
end

temperature = [];
for k = 1:numel(header)
    tok = header{k};
    if ischar(tok)
        tok(tok == 'K') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end + 1, 1) = val; %#ok<AGROW>
        end
    end
end

n_sets = floor(n_cols / 3);
temperature = temperature(1:min(numel(temperature), n_sets));
end
