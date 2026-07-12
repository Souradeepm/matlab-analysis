function n = get_dataset_temp_count(dataset_file)
% Return number of temperature sets available in dataset workbook.

if nargin < 1 || isempty(dataset_file)
    dataset_file = getenv('COMPARE_DATASET');
end
if isempty(dataset_file)
    error('Dataset file is required.');
end

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, dataset_file);
if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

[data, text, raw] = xlsread(input_path); %#ok<ASGLU>
n_sets = detect_n_sets_local(text, raw, size(data,2));
if n_sets < 1
    n_sets = floor(size(data,2) / 3);
end

temperature = parse_temperature_labels_local(text, raw, n_sets);
if isempty(temperature)
    temperature = (1:n_sets).';
end
n = numel(temperature(1:min(numel(temperature), n_sets)));
end

function n_sets = detect_n_sets_local(text, raw, n_data_cols)
n_sets = 0;
if ~isempty(text) && iscell(text)
    n_sets = count_header_triplets_local(text(1,:));
end
if n_sets < 1 && ~isempty(raw)
    n_sets = count_header_triplets_local(raw(1,:));
end
if n_sets < 1
    n_sets = floor(n_data_cols / 3);
end
end

function n_triplets = count_header_triplets_local(row_cells)
n_triplets = 0;
if isempty(row_cells), return; end
n_cols = numel(row_cells); max_start = n_cols - 2;
for c = 1:3:max_start
    a = row_cells{c}; b = row_cells{c+1}; d = row_cells{c+2};
    if ischar(a) && ischar(b) && ischar(d)
        if strcmpi(strtrim(a), 'freq') && strcmpi(strtrim(b), 'z') && strcmpi(strtrim(d), 'theta')
            n_triplets = n_triplets + 1;
        end
    end
end
end

function temperature = parse_temperature_labels_local(text, raw, n_sets)
temperature = [];
if ~isempty(text)
    if iscell(text), header = text(1, :); else, header = cellstr(char(text(1, :))); end
    for k = 1:numel(header)
        tok = header{k};
        if ischar(tok)
            tok(tok == 'K' | tok == 'k') = [];
            val = str2double(strtrim(tok));
            if ~isnan(val)
                temperature(end+1,1) = val; %#ok<AGROW>
            end
        end
    end
end
if numel(temperature) >= n_sets
    temperature = temperature(1:n_sets);
    return;
end
temperature = [];
if isempty(raw), return; end
scan_rows = min(size(raw,1), 12);
best_row = 0; best_count = 0; best_offset = 0;
for r = 1:scan_rows
    for off = 0:2
        count = 0;
        for c = (1+off):3:(3*n_sets)
            v = raw{r,c};
            if ischar(v) && (~isempty(strfind(v, 'K')) || ~isempty(strfind(v, 'k')))
                count = count + 1;
            end
        end
        if count > best_count
            best_count = count; best_row = r; best_offset = off;
        end
    end
end
if best_row == 0, return; end
temperature = nan(n_sets,1);
for i = 1:n_sets
    cols3 = (3*i-2):(3*i); found = NaN;
    c0 = 3*(i-1) + (1 + best_offset);
    if c0 >= 1 && c0 <= size(raw,2)
        v0 = raw{best_row,c0};
        if ischar(v0)
            tok0 = strtrim(v0);
            tok0(tok0 == 'K' | tok0 == 'k') = [];
            val0 = str2double(tok0);
            if ~isnan(val0), found = val0; end
        end
    end
    for c = cols3
        if ~isnan(found), break; end
        v = raw{best_row,c};
        if ischar(v)
            tok = strtrim(v);
            tok(tok == 'K' | tok == 'k') = [];
            val = str2double(tok);
            if ~isnan(val)
                found = val;
                break;
            end
        end
    end
    temperature(i) = found;
end
temperature = temperature(~isnan(temperature));
temperature = temperature(1:min(numel(temperature), n_sets));
end
