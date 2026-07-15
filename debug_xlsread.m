% Debug xlsread
dataset = 'S2022Sap.xlsx';

fprintf('Testing xlsread on %s...\n', dataset);

try
    [data, text, raw] = xlsread(dataset);
    
    fprintf('Data size: %d × %d\n', size(data, 1), size(data, 2));
    fprintf('Text size: %d × %d\n', size(text, 1), size(text, 2));
    fprintf('Raw size: %d × %d\n', size(raw, 1), size(raw, 2));
    
    fprintf('\nFirst row of raw (headers):\n');
    for i = 1:min(9, size(raw, 2))
        v = raw{1, i};
        if ischar(v)
            fprintf('  [%d]: "%s"\n', i, v);
        else
            fprintf('  [%d]: %s\n', i, class(v));
        end
    end
    
    fprintf('\nFirst row of data:\n');
    for i = 1:min(9, size(data, 2))
        fprintf('  [%d]: %.6e\n', i, data(1, i));
    end
    
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

exit;
