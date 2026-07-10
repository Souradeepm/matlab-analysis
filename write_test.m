fid = fopen('batch_write_test.txt','w');
if fid < 0
    error('open failed');
end
fprintf(fid, 'hello\n');
fclose(fid);
