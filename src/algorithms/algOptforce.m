function result = algOptforce(model, biomassRxn, targetRxn, CsourceRxn, outputFile)
% algOptforce  Find forced flux changes using optForce (Mendoza et al.).
%   result = algOptforce(model, biomassRxn, targetRxn, CsourceRxn, outputFile)
%
%   Mirrors the optForce procedure:
%     1. FBA on WT (max biomass) and target (max target)
%     2. FVA on WT and mutant (constrained) strains
%     3. findMustU / findMustL (first-order must sets)
%     4. findMustUU / findMustLL / findMustUL (second-order must sets)
%     5. optForce (forced intervention sets, K=1 then K=2)
%
%   All mustU / mustL / mustUU / mustLL / mustUL sets are stored in result.targets.
%   Re-run overwrites previous optforce.* files in the workspace subfolder.
%   Per-set interventions are grouped together in the CSV; sets are separated
%   by a blank line.
%
%   Inputs:
%     model       struct, COBRA-compatible model
%     biomassRxn  char, biomass reaction ID
%     targetRxn   char, target reaction ID
%     CsourceRxn  char, carbon source uptake reaction ID
%     outputFile  char, output CSV path (empty -> stdout; default -> workspace)
%
%   Output:
%     result - struct (algFseof-style 7 fields):
%         config     - actual params used
%         biomassRxn - biomass reaction ID
%         targetRxn  - target reaction ID
%         outputFile - CSV output path (or '' if stdout)
%         matFile    - path to optforce_result.mat
%         rows       - per-intervention rows (grouped by force set)
%         targets    - all must sets + force sets + FVA + diagnostics

    % --- Type squash ---
    biomassRxn = char(biomassRxn);
    targetRxn  = char(targetRxn);
    CsourceRxn = char(CsourceRxn);

    % --- Solver setup ---
    changeCobraSolver('gurobi', 'ALL');

    % --- outputFile default (mirror algFseof/algOptknock) ---
    defaultCsv = defaultOutputFile(model, targetRxn);
    algoDir    = fileparts(defaultCsv);
    matFile    = fullfile(algoDir, 'optforce_result.mat');

    if nargin < 5
        outputFile = defaultCsv;
    end
    output = ~isempty(outputFile);

    % --- Re-run cleanup (mirror algFseof/algOptknock) ---
    for oldName = {'optforce.txt', 'optforce.csv', 'optforce_result.mat'}
        p = fullfile(algoDir, oldName{1});
        if exist(p, 'file')
            try
                delete(p);
            catch
            end
        end
    end

    % --- Step 1: max growth + max target ---
    model = setParam(model, 'obj', biomassRxn, 1);
    growthRate_sol = solveLP(model);
    fprintf('[algOptforce] The maximum growth rate is %1.2f\n', growthRate_sol.f);

    model = setParam(model, 'obj', targetRxn, 1);
    maxTarget_sol = solveLP(model);
    fprintf('[algOptforce] The maximum production rate of Target is %1.2f\n', maxTarget_sol.f);

    wtBiomass     = growthRate_sol.f;
    targetMaxFlux = maxTarget_sol.f;

    % --- Step 2: WT / mutant constraints ---
    constrWT = struct('rxnList', {{biomassRxn}}, 'rxnValues', 0.95*wtBiomass, 'rxnBoundType', 'b');

    constrMT = struct('rxnList', {{biomassRxn, targetRxn}}, 'rxnValues', [0, 0.95*targetMaxFlux], ...
                      'rxnBoundType', 'bb');

    % --- Step 3: FVA on both strains ---
    [minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model, constrWT, constrMT);

    % --- Step 4a: First-order must sets ---
    runID = 'TestOptForceM';
    constrOpt = struct('rxnList', {{CsourceRxn, biomassRxn, targetRxn}}, 'values', [-10, 0, 0.95*targetMaxFlux]');

    [mustLSet, pos_mustL] = findMustL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                      'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
                                      'outputFileName', 'MustL' , 'printExcel', 1, 'printText', 1, ...
                                      'printReport', 1, 'keepInputs', 1, 'verbose', 0);

    [mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                                      'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
                                      'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, ...
                                      'printReport', 1, 'keepInputs', 1, 'verbose', 0);

    % --- Step 4b: Second-order must sets ---
    exchangeRxns = model.rxns(cellfun(@isempty, strfind(model.rxns, 'EX_')) == 0);
    excludedRxns = unique([mustUSet; mustLSet; exchangeRxns]);

    [mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = ...
        findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                   'excludedRxns', excludedRxns,'runID', runID, ...
                   'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', ...
                   'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
                   'verbose', 1);

    [mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = ...
        findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
                   'excludedRxns', excludedRxns,'runID', runID, ...
                   'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', ...
                   'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
                   'verbose', 1);

    [mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = ...
        findMustUL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, ...
               'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', ...
               'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, ...
               'verbose', 1);

    % --- Aggregate must sets ---
    mustU = unique(union(mustUSet, mustUU));
    mustL = unique(union(mustLSet, mustLL));

    % --- Step 5: optForce (K=1, then K=2) ---
    k = 1;
    nSets = 1;
    constrOpt_force = struct('rxnList', {{CsourceRxn, biomassRxn}}, 'values', [-100, 0]);

    [optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
        optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
                 minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
                 'k', k, 'nSets', nSets, 'constrOpt', constrOpt_force, ...
                 'runID', runID, 'outputFolder', 'OutputsOptForce', ...
                 'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
                 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

    % --- Second optForce call (K=2, nSets=20) ---
    k = 2;
    nSets = 20;
    runID2 = 'TestOptForceM2';
    excludedRxns2 = struct('rxnList', {{optForceSets}}, 'typeReg','U');
    [optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = ...
        optForce(model, targetRxn, biomassRxn, mustU, mustL, ...
                 minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
                 'k', k, 'nSets', nSets, 'constrOpt', constrOpt_force, ...
                 'excludedRxns', excludedRxns2, ...
                 'runID', runID2, 'outputFolder', 'OutputsOptForce', ...
                 'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
                 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

    % --- Build result struct (algFseof-style 7 fields) ---
    result = struct( ...
        'config',     struct('K', 2, 'NSets', nSets, ...
                             'WTGrowthFrac', 0.95, 'MTGrowthFrac', 0.0, ...
                             'CSourceBound', -10, ...
                             'CsourceRxn', CsourceRxn, ...
                             'RunID', runID2), ...
        'biomassRxn', biomassRxn, ...
        'targetRxn',  targetRxn, ...
        'outputFile', outputFile, ...
        'matFile',    '', ...
        'rows',       struct('setID', {}, 'interventionType', {}, 'rxnID', {}, ...
                             'rxnName', {}, 'subsystems', {}, 'grRule', {}, ...
                             'postFlux', {}), ...
        'targets',    struct( ...
            'mustUSet',      {{mustUSet}}, ...
            'mustLSet',      {{mustLSet}}, ...
            'mustUU',        {{mustUU}}, ...
            'mustLL',        {{mustLL}}, ...
            'mustUL',        {{mustUL}}, ...
            'mustU',         {{mustU}}, ...
            'mustL',         {{mustL}}, ...
            'forceSets',     {{optForceSets}}, ...
            'typeReg',       {{typeRegOptForceSets}}, ...
            'fluxes',        {flux_optForceSets}, ...
            'minFluxesW',    minFluxesW, ...
            'maxFluxesW',    maxFluxesW, ...
            'minFluxesM',    minFluxesM, ...
            'maxFluxesM',    maxFluxesM, ...
            'wtBiomass',     wtBiomass, ...
            'targetMaxFlux', targetMaxFlux));

    % --- Open CSV (file or stdout) and write header ---
    if output
        outputFile = char(outputFile);
        fid = fopen(outputFile, 'w');
        fprintf(fid, 'SetID,InterventionType,RxnID,RxnName,Subsystems,GrRule,PostFlux\n');
    else
        fprintf('SetID,InterventionType,RxnID,RxnName,Subsystems,GrRule,PostFlux\n');
    end

    % --- Build result.rows and write CSV (per-set, grouped, separated) ---
    nSetsOut = numel(optForceSets);
    for k = 1:nSetsOut
        rxnsInSet   = optForceSets{k};
        typesInSet  = typeRegOptForceSets{k};
        fluxesInSet = flux_optForceSets{k};
        for j = 1:numel(rxnsInSet)
            % --- Build row data ---
            rowFlux = 0;
            if ~isempty(fluxesInSet) && numel(fluxesInSet) >= j
                rowFlux = fluxesInSet(j);
            end
            rxnName    = '';
            subsystems = '';
            grRule     = '';
            pos = find(strcmp(model.rxns, rxnsInSet{j}), 1);
            if ~isempty(pos)
                if isfield(model, 'rxnNames')
                    rxnName = char(model.rxnNames{pos});
                end
                if isfield(model, 'subSystems')
                    ss = model.subSystems{pos};
                    if ~iscell(ss), ss = {ss}; end
                    subsystems = strjoin(ss, ';');
                end
                if isfield(model, 'grRules')
                    grRule = char(model.grRules{pos});
                end
            end

            % --- Append to result.rows ---
            result.rows(end+1) = struct( ...
                'setID',            k, ...
                'interventionType', typesInSet{j}, ...
                'rxnID',            rxnsInSet{j}, ...
                'rxnName',          rxnName, ...
                'subsystems',       subsystems, ...
                'grRule',           grRule, ...
                'postFlux',         rowFlux);

            % --- Write to CSV ---
            if output
                fprintf(fid, '%d,%s,%s,%s,%s,%s,%g\n', ...
                    k, csvEscape(typesInSet{j}), csvEscape(rxnsInSet{j}), ...
                    csvEscape(rxnName), csvEscape(subsystems), csvEscape(grRule), rowFlux);
            else
                fprintf('%d,%s,%s,%s,%s,%s,%g\n', ...
                    k, csvEscape(typesInSet{j}), csvEscape(rxnsInSet{j}), ...
                    csvEscape(rxnName), csvEscape(subsystems), csvEscape(grRule), rowFlux);
            end
        end
        % Separator between sets
        if output
            fprintf(fid, '\n');
        else
            fprintf('\n');
        end
    end

    if output
        fclose(fid);
    end

    % --- Persist result.mat ---
    if ~isempty(algoDir) && ~exist(algoDir, 'dir')
        mkdir(algoDir);
    end
    result.matFile = matFile;
    save(result.matFile, 'result');
    fprintf('[algOptforce] Saved `result` to %s\n', result.matFile);
end


function s = csvEscape(s)
% csvEscape  Quote a string for CSV output if it contains a comma, double
%   quote, or newline. Embedded double quotes are doubled per RFC 4180.
    if isempty(s), return; end
    if any(s == ',') || any(s == '"') || any(s == sprintf('\n')) || any(s == sprintf('\r'))
        s = strrep(s, '"', '""');
        s = ['"' s '"'];
    end
end


function outFile = defaultOutputFile(model, targetRxn)
% defaultOutputFile  Build a default output path inside the optForce subfolder
%   of the MDP workspaces tree.
%   workspaces/{safeModelBase}_{safeTargetRxn}/optForce/optforce.csv
    candidate = '';
    for fld = {'id', 'name', 'modelID', 'modelName'}
        if isfield(model, fld{1}) && ~isempty(model.(fld{1}))
            v = model.(fld{1});
            if ischar(v) || (isstring(v) && isscalar(v))
                candidate = char(v);
                break;
            end
        end
    end
    if isempty(candidate)
        candidate = 'model';
    end
    safeModel  = regexprep(candidate, '[^\w\-]', '_');
    safeTarget = regexprep(targetRxn,  '[^\w\-]', '_');
    wsDir   = fullfile('workspaces', [safeModel '_' safeTarget]);
    algoDir = fullfile(wsDir, 'optForce');
    if ~exist(algoDir, 'dir')
        mkdir(algoDir);
    end
    outFile = fullfile(algoDir, 'optforce.csv');
end
