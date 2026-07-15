% run_all_datasets_comprehensive.m
% Runs all three lambda selection methods on all datasets
% in batches of 10 temperatures per MATLAB session

clear all;

% ---- Configuration ----
batch_size = 10;
repo_root = pwd;
setenv('MATLAB_ANALYSIS_REPO_ROOT', repo_root);

datasets = {
  'S2022Sap.xlsx'
  'S2022Al.xlsx'
  'S2222Sap.xlsx'
  'S2222Al.xlsx'
  'S2302Sap.xlsx'
  'S2302Al.xlsx'
  'S2322Sap.xlsx'
  'S2322Al.xlsx'
  'S2332Sap.xlsx'
  'S2422Sap.xlsx'
  'S2422Al.xlsx'
};

output_file = fullfile(repo_root, 'comprehensive_batch_status.txt');
fid_status = fopen(output_file, 'w');
fprintf(fid_status, 'Comprehensive Batch Processing\n');
fprintf(fid_status, 'Started: %s\n', datetime('now'));
fprintf(fid_status, '=====================================\n\n');

total_batches = 0;
success_batches = 0;
failed_batches = 0;

%% ---- Process each dataset ----
for d = 1:numel(datasets)
    dataset_file = datasets{d};
    fprintf('Processing %s...\n', dataset_file);
    
    if ~isfile(dataset_file)
        fprintf('  SKIP: File not found\n');
        fprintf(fid_status, 'SKIP: %s (file not found)\n', dataset_file);
        continue;
    end
    
    % Count temperatures using xlsread
    try
        [data_tmp, ~, raw_tmp] = xlsread(dataset_file); %#ok<XLSRD>
        n_cols = size(data_tmp, 2);
        n_sets = floor(n_cols / 3);
        temperature = parse_temperature_from_raw(raw_tmp, n_sets);
        
        if isempty(temperature)
            fprintf('  SKIP: Unable to parse temperatures\n');
            fprintf(fid_status, 'SKIP: %s (parse error)\n', dataset_file);
            continue;
        end
        
        n_temps = numel(temperature);
        fprintf('  Found %d temperatures\n', n_temps);
        fprintf(fid_status, 'Dataset: %s (%d temps)\n', dataset_file, n_temps);
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        fprintf(fid_status, 'ERROR: %s - %s\n', dataset_file, ME.message);
        continue;
    end
    
    % Process in batches
    n_batch = ceil(n_temps / batch_size);
    fprintf('  Batches: %d\n', n_batch);
    
    for b = 1:n_batch
        start_idx = (b-1)*batch_size + 1;
        end_idx = min(b*batch_size, n_temps);
        
        total_batches = total_batches + 1;
        fprintf('    Batch %d [%d-%d]... ', b, start_idx, end_idx);
        
        try
            tic;
            run_all_three_methods_batch(dataset_file, start_idx, end_idx);
            elapsed = toc;
            fprintf('OK (%.0fs)\n', elapsed);
            success_batches = success_batches + 1;
            fprintf(fid_status, '  Batch %d: OK\n', b);
        catch ME
            fprintf('FAILED (%s)\n', ME.message);
            failed_batches = failed_batches + 1;
            fprintf(fid_status, '  Batch %d: FAILED - %s\n', b, ME.message);
        end
    end
    
    fprintf('\n');
end

%% ---- Summary ----
fprintf('\n============================================\n');
fprintf('Processing Complete\n');
fprintf('============================================\n');
fprintf('Total batches: %d\n', total_batches);
fprintf('Successful: %d\n', success_batches);
fprintf('Failed: %d\n', failed_batches);

fprintf(fid_status, '\n=====================================\n');
fprintf(fid_status, 'Summary\n');
fprintf(fid_status, 'Total batches: %d\n', total_batches);
fprintf(fid_status, 'Successful: %d\n', success_batches);
fprintf(fid_status, 'Failed: %d\n', failed_batches);
fprintf(fid_status, 'Ended: %s\n', datetime('now'));
fclose(fid_status);

fprintf('\nStatus file: %s\n', output_file);
exit;

%% ---- Helper: Parse temperatures from xlsread raw cell array ----
function temperature = parse_temperature_from_raw(raw, n_sets)
temperature = [];
if isempty(raw) || size(raw, 1) < 1
    return;
end

header = raw(1, :);
for k = 1:numel(header)
    tok = header{k};
    if ischar(tok)
        tok(tok == 'K' | tok == 'k') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end + 1, 1) = val; %#ok<AGROW>
        end
    end
end

if numel(temperature) >= n_sets
    temperature = temperature(1:n_sets);
end
end
