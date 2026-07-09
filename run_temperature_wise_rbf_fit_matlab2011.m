% run_temperature_wise_rbf_fit_matlab2011
% Simple runner for the MATLAB 2011-compatible Gaussian RBF fitting code.
%
% Edit the two Excel file paths below and run this script.

clc;
clear;
close all;

repo_root = fileparts(mfilename('fullpath'));
curve_excel_file = fullfile(repo_root, 'curve_input.xlsx');
peak_excel_file = fullfile(repo_root, 'peak_input.xlsx');
selected_temperature = [];
output_dir = fullfile(repo_root, 'rbf_fit_outputs_2011');

results = temperature_wise_rbf_fit_matlab2011( ...
    curve_excel_file, peak_excel_file, selected_temperature, output_dir);

disp('RBF fitting complete.');
disp(results);
