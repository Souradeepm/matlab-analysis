% Debug batch to trace exactly why temps = 0
setenv('MATLAB_ANALYSIS_REPO_ROOT', pwd);

dataset = 'S2022Sap.xlsx';
fid = fopen('debug_batch_trace.txt', 'w');

fprintf(fid, 'Reading %s\n', dataset);

[data, ~, raw] = xlsread(dataset);
fprintf(fid, 'data size: %d x %d\n', size(data,1), size(data,2));

n_cols = size(data, 2);
n_sets = floor(n_cols / 3);
fprintf(fid, 'n_sets: %d\n', n_sets);

% Parse temperature
temperature = [];
header = raw(1, :);
for k = 1:numel(header)
    tok = header{k};
    if ischar(tok)
        tok(tok == 'K' | tok == 'k') = [];
        val = str2double(strtrim(tok));
        if ~isnan(val)
            temperature(end+1, 1) = val;
        end
    end
end
fprintf(fid, 'temperatures found: %d\n', numel(temperature));
if numel(temperature) > 0
    fprintf(fid, 'first 3 temps: %.4f %.4f %.4f\n', temperature(1), temperature(2), temperature(3));
end

% Test column access for t=1
t = 1;
c1 = (t-1)*3 + 1;
c3 = (t-1)*3 + 3;
fprintf(fid, 'Temp 1 columns: %d to %d\n', c1, c3);

block = data(:, c1:c3);
fprintf(fid, 'block size: %d x %d\n', size(block,1), size(block,2));
fprintf(fid, 'block isfinite rows: %d / %d\n', sum(all(isfinite(block),2)), size(block,1));
fprintf(fid, 'first 3 rows of block:\n');
for r = 1:min(3, size(block,1))
    fprintf(fid, '  row %d: %.6e  %.6e  %.6e\n', r, block(r,1), block(r,2), block(r,3));
end

valid = all(isfinite(block), 2);
block = block(valid, :);
fprintf(fid, 'valid rows after filter: %d\n', size(block,1));

fclose(fid);
exit;
