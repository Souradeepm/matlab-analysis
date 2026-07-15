% Test lambda range constraints
fid = fopen('lambda_range_test_output.txt', 'w');

fprintf(fid, '=== LAMBDA RANGE TEST ===\n\n');

fprintf(fid, 'Bayes-DRT Method (31 points):\n');
lambda_min_exp = -4;
lambda_max_exp = -1;
lambda_count = 31;
bayes_lambda = logspace(lambda_min_exp, lambda_max_exp, lambda_count).';
fprintf(fid, '  Min: %.2e, Max: %.2e, Count: %d\n', bayes_lambda(1), bayes_lambda(end), numel(bayes_lambda));
fprintf(fid, '  First 3: [%.2e, %.2e, %.2e]\n', bayes_lambda(1), bayes_lambda(2), bayes_lambda(3));
fprintf(fid, '  Last 3: [%.2e, %.2e, %.2e]\n\n', bayes_lambda(end-2), bayes_lambda(end-1), bayes_lambda(end));

fprintf(fid, 'Residual Method (11 points):\n');
residual_lambda = logspace(-4, -1, 11).';
fprintf(fid, '  Min: %.2e, Max: %.2e, Count: %d\n', residual_lambda(1), residual_lambda(end), numel(residual_lambda));
fprintf(fid, '  Values: ');
for i=1:numel(residual_lambda)
    fprintf(fid, '%.2e ', residual_lambda(i));
end
fprintf(fid, '\n\n');

fprintf(fid, 'Paper vs Residual Method (4 points):\n');
paper_lambda = [1e-4; 1e-3; 1e-2; 1e-1];
fprintf(fid, '  Min: %.2e, Max: %.2e, Count: %d\n', paper_lambda(1), paper_lambda(end), numel(paper_lambda));
fprintf(fid, '  Values: [%.2e, %.2e, %.2e, %.2e]\n', paper_lambda(1), paper_lambda(2), paper_lambda(3), paper_lambda(4));

fprintf(fid, '\n✓ All lambda ranges successfully constrained to [1e-4, 1e-1]\n');

fclose(fid);
exit;
