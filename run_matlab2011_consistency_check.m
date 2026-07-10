function run_matlab2011_consistency_check()
% run_matlab2011_consistency_check
% Static consistency check report for MATLAB 2011-compatible workflow files.

repo_root = fileparts(mfilename('fullpath'));
files = {
    fullfile(repo_root, 'Impedance_analysis.m');
    fullfile(repo_root, 'drt_input_analysis_matlab2011.m')
};

out_path = fullfile(repo_root, 'matlab2011_consistency_report.txt');
fid = fopen(out_path, 'w');
if fid < 0
    error('Failed to open report file: %s', out_path);
end

fprintf(fid, 'MATLAB checkcode consistency report\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

for i = 1:numel(files)
    f = files{i};
    if exist(f, 'file') ~= 2
        fprintf(fid, 'File: %s\nIssues: N/A (missing file)\n\n', f);
        continue;
    end

    msgs = checkcode(f, '-id');
    fprintf(fid, 'File: %s\n', f);
    fprintf(fid, 'Issues: %d\n', numel(msgs));
    for k = 1:numel(msgs)
        fprintf(fid, '  Line %d, Col %d, ID %s: %s\n', ...
            msgs(k).line, msgs(k).column, msgs(k).id, msgs(k).message);
    end
    fprintf(fid, '\n');
end

fclose(fid);
fprintf('Wrote MATLAB consistency report: %s\n', out_path);
end
