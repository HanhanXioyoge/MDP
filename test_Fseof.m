% test_Fseof.m  Single-algorithm debug script for FSEOF (no workspace)
addpath(genpath('src'));

model = loadModel('input_models/iHM1533/iHM1533_heparosan.mat');
result = algFseof(model, 'EX_heparosan_e', ...
    'BIOMASS_EcN_iHM1533_core_59p80M', ...
    'Iterations', 21, 'Coefficient', 0.99);

fprintf('\n=== FSEOF Test Results ===\n');
fprintf('WT biomass:  %g\n', result.diagnostics.wtBiomass);
fprintf('Target max:  %g\n', result.diagnostics.targetMaxFlux);
fprintf('Candidates:  %d\n', height(result.candidates));
fprintf('Fallback:    %d\n', result.fallbackUsed);
if ~result.fallbackUsed
    disp(result.candidates(1:min(10, height(result.candidates)), :));
end
