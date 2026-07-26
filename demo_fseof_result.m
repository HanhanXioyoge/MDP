% demo_fseof_result.m  Minimal smoke test: load model + run fseof + capture into `result`.

addpath(genpath('src'));

modelPath  = 'input_models/iHM1533/iHM1533_heparosan.mat';
targetRxn  = 'EX_heparosan_e';
biomassRxn = 'BIOMASS_EcN_iHM1533_core_59p80M';

%% --- 1) Load model ---
fprintf('--- Loading model ---\n');
model = loadModel(modelPath);
fprintf('Model loaded: %d reactions (modelType=%s)\n', ...
        numel(model.rxns), model.modelType);

%% --- 2) Run fseof ---
fprintf('\n--- Running algFseof() — Iterations=21, Coefficient=0.99 ---\n');
fseof_result = algFseof(model, biomassRxn, targetRxn, 21, 0.99);

%% --- 3) Capture into `result` ---
result.model      = 'iHM1533_heparosan.mat';
result.biomassRxn = biomassRxn;
result.targetRxn  = targetRxn;
result.fseof      = fseof_result.rows;

%% --- 4) Show what we got ---
fprintf('\n--- result ---\n');
disp(result);
fprintf('\nnumel(result.fseof) = %d  (one struct per FSEOF-hit reaction)\n', ...
        numel(result.fseof));
fprintf('Tab-delimited text (same data) written to:\n  %s\n', ...
        fseof_result.outputFile);
