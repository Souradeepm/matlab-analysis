% Run comprehensive batch with all three methods (FULL VERSION)
% This calls run_all_three_methods_batch which now has robust parsing

fprintf('\n====== FULL Comprehensive Batch Run (Robust Parsing) ======\n');
fprintf('Started: %s\n\n', datetime('now'));

setenv('MATLAB_ANALYSIS_REPO_ROOT', pwd);

datasets = {
    'S2022Sap.xlsx', 'S2022Sap';
    'S2022Al.xlsx', 'S2022Al';
    'S2222Sap.xlsx', 'S2222Sap';
    'S2222Al.xlsx', 'S2222Al';
    'S2302Sap.xlsx', 'S2302Sap';
    'S2302Al.xlsx', 'S2302Al';
    'S2322Sap.xlsx', 'S2322Sap';
    'S2322Al.xlsx', 'S2322Al';
    'S2332Sap.xlsx', 'S2332Sap';
    'S2422Sap.xlsx', 'S2422Sap';
    'S2422Al.xlsx', 'S2422Al'
};

output_file = 'comprehensive_batch_full_robust_status.txt';
fid_status = fopen(output_file, 'w');
fprintf(fid_status, 'FULL Comprehensive Batch (All 3 Methods, Robust Parsing)\n');
fprintf(fid_status, 'Started: %s\n', datetime('now'));
fprintf(fid_status, '================================================\n\n');

total_batches = 0;
total_success = 0;
total_failed = 0;

for d = 1:size(datasets, 1)
    dataset_file = datasets{d, 1};
    ds_tag = datasets{d, 2};
    
    fprintf('Processing %s...\n', ds_tag);
    fprintf(fid_status, 'Dataset: %s (%d)\n', dataset_file, d);
    
    try
        % Run batch processing for this dataset
        run_all_three_methods_batch(dataset_file);
        fprintf('  OK\n');
        fprintf(fid_status, '  Status: OK\n');
        total_success = total_success + 1;
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        fprintf(fid_status, '  ERROR: %s\n', ME.message);
        total_failed = total_failed + 1;
    end
    
    fprintf(fid_status, '\n');
end

fprintf(fid_status, '================================================\n');
fprintf(fid_status, 'Summary\n');
fprintf(fid_status, 'Total datasets: %d\n', size(datasets,1));
fprintf(fid_status, 'Successful: %d\n', total_success);
fprintf(fid_status, 'Failed: %d\n', total_failed);
fprintf(fid_status, 'Ended: %s\n', datetime('now'));
fclose(fid_status);

fprintf('\n%s created\n', output_file);
fprintf('Summary: %d/%d datasets processed\n', total_success, size(datasets,1));

exit;
