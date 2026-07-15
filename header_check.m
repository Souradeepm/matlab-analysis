% Check actual header structure for S2322Al
fid = fopen('header_check.txt', 'w');

datasets = {'S2322Al.xlsx', 's2322al'; 'S2422Al.xlsx', 's2422al'};

for d = 1:size(datasets, 1)
    dataset = datasets{d,1};
    tag = datasets{d,2};
    
    fprintf(fid, '\n%s\n', tag);
    fprintf(fid, '================\n');
    
    try
        [num, txt, raw] = xlsread(dataset);
        fprintf(fid, 'num size: %dx%d\n', size(num,1), size(num,2));
        fprintf(fid, 'txt size: %dx%d\n', size(txt,1), size(txt,2));
        fprintf(fid, 'raw size: %dx%d\n', size(raw,1), size(raw,2));
        
        fprintf(fid, '\nFirst 9 header cells (raw):\n');
        for k = 1:min(9, size(raw,2))
            v = raw{1,k};
            fprintf(fid, '  [%d] class=%s ', k, class(v));
            if ischar(v)
                fprintf(fid, 'str="%s"\n', v);
            elseif isnumeric(v)
                fprintf(fid, 'num=%.6g scalar=%d\n', v, isscalar(v));
            elseif iscell(v)
                fprintf(fid, 'cell size=%s\n', mat2str(size(v)));
            else
                fprintf(fid, 'val=%s\n', mat2str(v));
            end
        end
        
        fprintf(fid, '\nFirst 9 header cells (txt):\n');
        for k = 1:min(9, size(txt,2))
            v = txt{1,k};
            fprintf(fid, '  [%d] "%s"\n', k, v);
        end
        
    catch ME
        fprintf(fid, 'ERROR: %s\n', ME.message);
    end
end

fclose(fid);
exit;
