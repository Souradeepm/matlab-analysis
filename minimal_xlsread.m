% Minimal xlsread diagnostic
fprintf('Checking S2322Al.xlsx...\n');
try
    [num, txt, raw] = xlsread('S2322Al.xlsx');
    fprintf('  num: %dx%d\n', size(num,1), size(num,2));
    fprintf('  txt: %dx%d\n', size(txt,1), size(txt,2));
    fprintf('  raw: %dx%d\n', size(raw,1), size(raw,2));
    if size(raw,1) >= 3 && size(raw,2) >= 3
        fprintf('  raw{1,1:3}: %s, %s, %s\n', char(raw{1,1}), char(raw{1,2}), char(raw{1,3}));
        fprintf('  raw{2,1:3}: %s, %s, %s\n', char(raw{2,1}), char(raw{2,2}), char(raw{2,3}));
        fprintf('  raw{3,1:3}: %s, %s, %s\n', char(raw{3,1}), char(raw{3,2}), char(raw{3,3}));
    end
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

fprintf('\nChecking S2422Al.xlsx...\n');
try
    [num, txt, raw] = xlsread('S2422Al.xlsx');
    fprintf('  num: %dx%d\n', size(num,1), size(num,2));
    fprintf('  txt: %dx%d\n', size(txt,1), size(txt,2));
    fprintf('  raw: %dx%d\n', size(raw,1), size(raw,2));
    if size(raw,1) >= 3 && size(raw,2) >= 3
        fprintf('  raw{1,1:3}: %s, %s, %s\n', char(raw{1,1}), char(raw{1,2}), char(raw{1,3}));
        fprintf('  raw{2,1:3}: %s, %s, %s\n', char(raw{2,1}), char(raw{2,2}), char(raw{2,3}));
        fprintf('  raw{3,1:3}: %s, %s, %s\n', char(raw{3,1}), char(raw{3,2}), char(raw{3,3}));
    end
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

exit;
