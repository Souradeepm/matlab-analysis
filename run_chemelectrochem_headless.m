function run_chemelectrochem_headless()
% run_chemelectrochem_headless
% Runs the ChemElectroChem DRT lambda optimization in headless mode:
% no figure windows, no Excel writes - outputs only printed metrics
% and a plain-text summary file.

setenv('CHEMELECTROCHEM_HEADLESS', '1');
run_s2022_280k_lambda_sweep_cv();
setenv('CHEMELECTROCHEM_HEADLESS', '');
end
