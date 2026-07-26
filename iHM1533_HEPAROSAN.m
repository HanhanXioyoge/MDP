% iHM1533_HEPAROSAN.m  Main test script: batch strain design with workspace

modelPath = 'input_models/iHM1533/iHM1533_heparosan.mat';
targetRxn = 'EX_heparosan_e';
biomassRxn = 'BIOMASS_EcN_iHM1533_core_59p80M';

%% --- 1) Load model ---
fprintf('--- Loading model ---\n');
model = loadModel(modelPath);
model=setParam(model,'obj',biomassRxn,1);
fprintf('Model loaded: %d reactions (modelType=%s)\n', ...
        numel(model.rxns), model.modelType);

%% --- 2) Run algorithms ---
fprintf('\n--- Running algFseof() — Iterations=21, Coefficient=0.99 ---\n');
fseof_result = algFseof(model, biomassRxn, targetRxn, 21, 0.99);


fprintf('\n--- Running algOptknock() — NumDel=1, MaxSolutions=200, minGrowth=0.5*WT ---\n');

% Filter to reactions with GPR (gene-protein-reaction) rule —
% only these can be genetically knocked out. Excludes '()' / '[]' / empty rules.
selectedRxnList = rxnsWithGpr(model);
fprintf('KO candidates: %d reactions (with GPR) out of %d total\n', ...
        numel(selectedRxnList), numel(model.rxns));

optknock_result = algOptknock(model, biomassRxn, targetRxn, selectedRxnList, ...
                              2, 10, 0.5, 'G', true);
%{
fprintf('\n--- Running algOptforce() — K=3, NSets=2 ---\n');
optforce_result = algOptforce(model, biomassRxn, targetRxn, ...
                           'K', 3, 'NSets', 2);
%}
%% --- 3) Wrap into a single `combined` struct (saved to workspace) ---
fprintf('\n--- Wrapping all results into one struct ---\n');
%{
combined = combineResults(model, biomassRxn, targetRxn, ...
                          'FSEOF',    fseof_result, ...
                          'OptKnock', optknock_result, ...
                          'optForce', optforce_result);
%}
combined = combineResults(model, biomassRxn, targetRxn, ...
                          'FSEOF',    fseof_result, ...
                          'OptKnock', optknock_result);
%% --- 4) Show ---
fprintf('\n--- combined struct ---\n');
disp(combined);
fprintf('\nFSEOF hit count        : %d\n', numel(combined.FSEOF));
%{
fprintf('OptKnock candidates   : %d\n', ...
        combined.OptKnock.fallbackUsed + height(combined.OptKnock.candidates));  % see note
fprintf('optForce candidates   : %d\n', ...
        combined.optForce.fallbackUsed + height(combined.optForce.candidates));
fprintf('Saved to              : %s\n', combined.combinedFile);
%}


% === Local helper ===

function rxns = rxnsWithGpr(model)
% rxnsWithGpr  Return reaction IDs that have a non-empty GPR rule.
%   Empty markers '()' and '[]' are treated as no-GPR.
%   Falls back to all reactions if model.grRules is missing.

    if ~isfield(model, 'grRules')
        warning('rxnsWithGpr:NoGrRules', ...
                'model.grRules missing; using all reactions.');
        rxns = model.rxns;
        return;
    end

    nRxns = numel(model.rxns);
    mask  = false(nRxns, 1);
    for i = 1:nRxns
        rule = model.grRules{i};
        if isempty(rule), continue; end
        rule = strtrim(rule);
        if isempty(rule), continue; end
        if strcmp(rule, '()') || strcmp(rule, '[]'), continue; end
        mask(i) = true;
    end

    rxns = model.rxns(mask);
end
