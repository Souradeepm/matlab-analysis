% Quick K-K validation rerun for S2322Al and S2422Al only
setenv('MATLAB_ANALYSIS_REPO_ROOT', pwd);

fid = fopen('kk_rerun_s2322al_s2422al.txt', 'w');
fprintf(fid, 'K-K Validation Rerun (S2322Al and S2422Al)\n');
fprintf(fid, 'Started: %s\n', datetime('now'));
fprintf(fid, '===============================================\n\n');

datasets = {'S2322Al.xlsx', 's2322al'; 'S2422Al.xlsx', 's2422al'};

for d = 1:size(datasets, 1)
    xlsx_file = datasets{d,1};
    ds_tag = datasets{d,2};
    
    fprintf(fid, '%s:\n', ds_tag);
    
    try
        [num, ~, raw] = xlsread(xlsx_file); %#ok<XLSRD>
        
        n_sets = floor(size(num,2) / 3);
        
        % FIXED parsing - handle both numeric and string headers
        temps = [];
        if ~isempty(raw) && size(raw,1) >= 1
            header_row = raw(1, :);
            for h = 1:min(numel(header_row), 3*n_sets)
                v = header_row{h};
                if isnumeric(v) && isscalar(v) && isfinite(v) && v > 0 && v < 1000
                    temps = [temps; v]; %#ok<AGROW>
                elseif ischar(v)
                    val_str = regexprep(v, '[^\d.]', '');
                    if ~isempty(val_str)
                        vv = str2double(val_str);
                        if isfinite(vv) && vv > 0 && vv < 1000
                            temps = [temps; vv]; %#ok<AGROW>
                        end
                    end
                end
            end
        end
        
        fprintf(fid, '  Temperatures found: %d\n', numel(temps));
        
        if numel(temps) > 0
            fprintf(fid, '  Sample temps: %.1f K, %.1f K, %.1f K\n', temps(1), temps(min(2,end)), temps(end));
        end
        
    catch ME
        fprintf(fid, '  ERROR: %s\n', ME.message);
    end
    fprintf(fid, '\n');
end

fprintf(fid, '===============================================\n');
fprintf(fid, 'Completed: %s\n', datetime('now'));
fclose(fid);

fprintf('Done. Check kk_rerun_s2322al_s2422al.txt\n');
exit;
