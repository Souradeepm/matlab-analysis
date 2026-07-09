clc;
clear;
close all;

repo_root = fileparts(mfilename('fullpath'));
input_path = fullfile(repo_root, 'S2422-SapIMPEDANCE 29 October 2024 2.26 .xls');
data = dlmread(input_path, '\t');

% Single dataset in this file; temperature inferred from filename tag "S2422".
temperature = 242.2;

drt_input_analysis_matlab2011(repo_root, input_path, temperature, data);
