% Quick test: Run all three methods on first 10 temperatures of S2022Sap
setenv('MATLAB_ANALYSIS_REPO_ROOT', pwd);

fid = fopen('test_three_methods_output.txt', 'w');
fprintf(fid, 'Testing all three methods on S2022Sap (first 10 temps)...\n');
fprintf(fid, 'Working directory: %s\n', pwd);

dataset = 'S2022Sap.xlsx';
fprintf(fid, 'Dataset path: %s\n', dataset);
fprintf(fid, 'Dataset exists: %d\n', isfile(dataset));

if ~isfile(dataset)
    fprintf(fid, '\n✗ Dataset not found!\n');
    fclose(fid);
    exit(1);
end

% Direct test of one temperature to check data extraction
fprintf(fid, '\n--- Direct data test (t=1) ---\n');
try
    [data, ~, raw] = xlsread(dataset); %#ok<XLSRD>
    fprintf(fid, 'xlsread OK: data %dx%d, raw %dx%d\n', size(data,1), size(data,2), size(raw,1), size(raw,2));
    block = data(:, 1:3);
    valid = all(isfinite(block), 2) & block(:,1) > 0;
    block = block(valid, :);
    fprintf(fid, 'Valid rows: %d\n', size(block,1));
    if size(block,1) >= 5
        freq_vec = block(:,1);
        Z_exp = block(:,2) .* exp(1i * block(:,3) * pi / 180);
        fprintf(fid, 'freq range: %.2f to %.2f Hz\n', min(freq_vec), max(freq_vec));
        fprintf(fid, '|Z| range: %.4e to %.4e Ohm\n', min(abs(Z_exp)), max(abs(Z_exp)));
        fprintf(fid, 'Data extraction: OK\n');
    else
        fprintf(fid, 'Not enough valid rows!\n');
    end
catch ME
    fprintf(fid, 'Direct test error: %s\n', ME.message);
end

fprintf(fid, '\nCalling run_all_three_methods_batch (t=1 only)...\n');
try
    run_all_three_methods_batch(dataset, 1, 1);
    fprintf(fid, '✓ Batch function returned OK\n');
catch ME
    fprintf(fid, '✗ Batch function error: %s\n', ME.message);
    fprintf(fid, 'Stack: %s\n', ME.getReport());
end

% Check output file
bayes_file = 's2022sap_bayes_drt_matlab2011_b1.txt';
if isfile(bayes_file)
    content = fileread(bayes_file);
    fprintf(fid, '\nBayes output file (%d bytes):\n%s\n', numel(content), content);
else
    fprintf(fid, '\nBayes output file NOT created!\n');
end

fclose(fid);
exit(0);
