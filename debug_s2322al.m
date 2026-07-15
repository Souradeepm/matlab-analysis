% Debug S2322Al header structure
dataset = 'S2322Al.xlsx';
fid = fopen('debug_s2322al_output.txt', 'w');

fprintf(fid, 'Debugging S2322Al.xlsx\n');
fprintf(fid, '======================\n\n');

try
    [num, txt, raw] = xlsread(dataset);
    fprintf(fid, 'xlsread OK\n');
    fprintf(fid, 'num: %dx%d\n', size(num,1), size(num,2));
    fprintf(fid, 'txt: %dx%d\n', size(txt,1), size(txt,2));
    fprintf(fid, 'raw: %dx%d\n\n', size(raw,1), size(raw,2));

    fprintf(fid, '--- First 3 rows of raw (first 12 cols) ---\n');
    for r = 1:min(3, size(raw,1))
        fprintf(fid, 'Row %d:\n', r);
        for c = 1:min(12, size(raw,2))
            v = raw{r,c};
            if ischar(v)
                fprintf(fid, '  [%d]: (str) "%s"\n', c, v);
            elseif isnumeric(v)
                fprintf(fid, '  [%d]: (num) %.6g\n', c, v);
            else
                fprintf(fid, '  [%d]: (%s)\n', c, class(v));
            end
        end
    end

    fprintf(fid, '\n--- txt first row (first 12 cols) ---\n');
    if ~isempty(txt)
        for c = 1:min(12, size(txt,2))
            fprintf(fid, '  [%d]: "%s"\n', c, txt{1,c});
        end
    else
        fprintf(fid, '  txt is empty\n');
    end

    fprintf(fid, '\n--- Detecting n_sets ---\n');
    n_cols = size(num,2);
    n_sets = floor(n_cols / 3);
    fprintf(fid, 'n_cols: %d, n_sets: %d\n', n_cols, n_sets);

    fprintf(fid, '\n--- Temperature parse attempt ---\n');
    if ~isempty(raw)
        header_row = raw(1, :);
        for h = 1:min(9, numel(header_row))
            v = header_row{h};
            if ischar(v)
                val_str = regexprep(v, '[^\d.]', '');
                fprintf(fid, '  [%d]: str="%s" stripped="%s"\n', h, v, val_str);
            else
                fprintf(fid, '  [%d]: type=%s val=%s\n', h, class(v), mat2str(v));
            end
        end
    end

catch ME
    fprintf(fid, 'ERROR: %s\n', ME.message);
end

fclose(fid);
exit;
