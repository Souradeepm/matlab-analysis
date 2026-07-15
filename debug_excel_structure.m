% Debug: Check Excel structure
dataset = 'S2022Sap.xlsx';

fprintf('Reading %s...\n', dataset);
try
    data = readtable(dataset, 'Sheet', 1);
    fprintf('Dimensions: %d rows × %d cols\n', height(data), width(data));
    fprintf('Variable names:\n');
    var_names = data.Properties.VariableNames;
    for i = 1:min(9, numel(var_names))
        fprintf('  %d: %s\n', i, var_names{i});
    end
    
    % Try to parse temperatures
    fprintf('\nParsing temperatures...\n');
    temperature = [];
    for k = 1:numel(var_names)
        name_str = var_names{k};
        fprintf('  Parsing "%s"... ', name_str);
        name_str_clean = name_str;
        name_str_clean(name_str_clean == 'K') = [];
        val = str2double(strtrim(name_str_clean));
        fprintf('→ %.6f\n', val);
        if ~isnan(val)
            temperature(end + 1, 1) = val;
        end
    end
    
    fprintf('\nFound %d temperatures\n', numel(temperature));
    if numel(temperature) > 0
        fprintf('First 3: %.6f, %.6f, %.6f\n', temperature(1), temperature(2), temperature(3));
    end
    
catch ME
    fprintf('ERROR: %s\n', ME.message);
end

exit;
