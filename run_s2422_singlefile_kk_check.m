function run_s2422_singlefile_kk_check()
% run_s2422_singlefile_kk_check - Route the tab-delimited S2422 file through the single-file workflow.
% Use this to inspect the Kramers-Kronig section inside Impedance_analysis_singlefile.

clc;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');

if exist(input_path, 'file') ~= 2
    error('Input file not found: %s', input_path);
end

Impedance_analysis_singlefile(input_path);
end
