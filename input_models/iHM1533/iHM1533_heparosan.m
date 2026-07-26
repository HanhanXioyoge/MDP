% =========================================================================
% Script to add Heparin precursor (Heparosan) synthesis and transport
% reactions to the iHM1533 genome-scale metabolic model.
% Requirements: COBRA Toolbox initialized in MATLAB.
% =========================================================================

fprintf('\n');
fprintf('############################################################\n');
fprintf('#  iHM1533 Heparosan Production Modeling Pipeline           \n');
fprintf('#  Started: %s\n', datestr(now));
fprintf('############################################################\n\n');

% 1. Initialize COBRA Toolbox (Uncomment if not already initialized)
% initCobraToolbox(false);

% 2. Load the original iHM1533 model
% Note: Change 'iHM1533.xml' to your actual file path and format if different
model = readCbModel('iHM1533.xml');

% 3. Define Metabolite IDs
% IMPORTANT: You must check iHM1533 to ensure these IDs match the model's namespace.
% Commonly used BiGG IDs are used here as placeholders.
udp_glcnac = 'uacgam[c]';    % UDP-N-acetyl-D-glucosamine
udp_glca = 'udpglcur[c]';    % UDP-D-glucuronate
udp = 'udp[c]';              % UDP
h_c = 'h[c]';                % Proton
heparosan_c = 'heparosan[c]'; % Heparosan (cytosol)
heparosan_e = 'heparosan[e]'; % Heparosan (extracellular)

% 3.1 Add heparosan metabolites if they don't already exist in the model
% Heparosan disaccharide repeat unit: GlcNAc + GlcA, formula C14H21NO11
if ~any(strcmp(model.mets, heparosan_c))
    model = addMetabolite(model, heparosan_c, ...
        'metName', 'Heparosan (cytosolic disaccharide repeat unit)', ...
        'metFormula', 'C14H21NO11');
    fprintf('[OK] Added metabolite: %s\n', heparosan_c);
else
    fprintf('[INFO] Metabolite already present: %s\n', heparosan_c);
end
if ~any(strcmp(model.mets, heparosan_e))
    model = addMetabolite(model, heparosan_e, ...
        'metName', 'Heparosan (extracellular disaccharide repeat unit)', ...
        'metFormula', 'C14H21NO11');
    fprintf('[OK] Added metabolite: %s\n', heparosan_e);
else
    fprintf('[INFO] Metabolite already present: %s\n', heparosan_e);
end

% 3.2 Validate that required metabolites exist in the model
requiredMets = {udp_glcnac, udp_glca, udp, h_c, heparosan_c, heparosan_e};
missingMets = setdiff(requiredMets, model.mets);
if ~isempty(missingMets)
    fprintf('[ERROR] Missing metabolites: %s\n', strjoin(missingMets, ', '));
    error('iHM1533_heparosan:MissingMetabolites', ...
        ['The following metabolites were not found in the model. ', ...
         'Please update the IDs in section 3 to match your model namespace:\n%s'], ...
        strjoin(missingMets, ', '));
end
fprintf('[OK] All %d required metabolites are present in the model.\n', length(requiredMets));

% 4. Add Heparosan Synthesis Reaction (Cytosol)
% Formula: UDP-GlcNAc + UDP-GlcA -> Heparosan[c] + 2 UDP + 2 H+
% Note: Stoichiometry might vary based on your specific polymerization assumptions.
% Gene cluster kfiABC encodes the heparosan synthase complex.
model = addReaction(model, 'HEPAROSAN_SYN', ...
    'reactionFormula', sprintf('%s + %s -> %s + 2 %s + 2 %s', udp_glcnac, udp_glca, heparosan_c, udp, h_c), ...
    'geneRule', 'kfiA and kfiB and kfiC', ... % Heparosan synthase gene cluster
    'lowerBound', 0, ...                       % Irreversible reaction
    'upperBound', 1000, ...
    'subSystem', 'Heparosan biosynthesis', ...
    'reactionName', 'Heparosan synthase');

% 5. Add Heparosan Transport Reaction (Cytosol to Extracellular)
% Formula: Heparosan[c] -> Heparosan[e]
% kpsD, kpsM, kpsT encode the Kps transport system for surface polysaccharide export.
model = addReaction(model, 'HEPAROSAN_t', ...
    'reactionFormula', sprintf('%s -> %s', heparosan_c, heparosan_e), ...
    'geneRule', 'kpsD and kpsM and kpsT', ... % Kps polysaccharide transport system
    'lowerBound', 0, ...
    'upperBound', 1000, ...
    'subSystem', 'Transport, extracellular', ...
    'reactionName', 'Heparosan transport');

% 6. Add Heparosan Exchange Reaction (System Boundary)
% Formula: Heparosan[e] <=>
% Allows the product to leave the simulated system space.
model = addReaction(model, 'EX_heparosan_e', ...
    'reactionFormula', sprintf('%s -> ', heparosan_e), ...
    'lowerBound', 0, ...             % 0 means it can only be secreted, not consumed
    'upperBound', 1000, ...
    'subSystem', 'Exchange', ...
    'reactionName', 'Heparosan exchange');

% =========================================================================
% Optional: Test the modification by running FBA
% model = changeObjective(model, 'EX_heparosan_e');
% fbaResult = optimizeCbModel(model, 'max');
% fprintf('Maximum theoretical production rate of Heparosan: %f\n', fbaResult.f);
% =========================================================================

% =========================================================================
% 8. Flux Balance Analysis: Maximum Heparosan Production
% =========================================================================
% Constrain glucose uptake to 10 mmol/gDW/h
glcUptakeBound = -10;  % Negative indicates uptake in COBRA convention
glcExchangeIDs = {'EX_glc__D_e', 'EX_glc_e', 'EX_glc_D_e'};
glcConstrained = false;
for i = 1:length(glcExchangeIDs)
    if any(strcmp(model.rxns, glcExchangeIDs{i}))
        model = changeRxnBounds(model, glcExchangeIDs{i}, glcUptakeBound, 'l');
        fprintf('[OK] Constrained glucose uptake: %s lowerBound = %.1f (uptake = %.1f mmol/gDW/h)\n', ...
            glcExchangeIDs{i}, glcUptakeBound, abs(glcUptakeBound));
        glcConstrained = true;
        break;
    end
end
if ~glcConstrained
    fprintf('[WARN] Glucose exchange reaction not found. Tried: %s\n', ...
        strjoin(glcExchangeIDs, ', '));
    warning('iHM1533_heparosan:GlucoseNotFound', ...
        'Glucose exchange reaction not found. Tried: %s', strjoin(glcExchangeIDs, ', '));
end

% Set objective to maximize heparosan secretion
model = changeObjective(model, 'EX_heparosan_e');

% Run FBA
fbaSolution = optimizeCbModel(model, 'max');

% Initialize results struct
results = struct();
results.timestamp = datestr(now);
results.fba = struct();
results.fva = struct();
results.yield = struct();

fprintf('\n============================================================\n');
fprintf('  FBA: Maximum Heparosan Production\n');
fprintf('============================================================\n');

if fbaSolution.stat == 1
    fprintf('[OK] Solver status: Optimal\n');
    fprintf('     Maximum heparosan production rate: %.4f mmol/gDW/h\n', fbaSolution.f);
    results.fba.maxProduction = fbaSolution.f;
    results.fba.status = 'optimal';
else
    fprintf('[ERROR] Solver status: %d (not optimal). Check model feasibility.\n', fbaSolution.stat);
    fprintf('        Heparan production flux: %.4f mmol/gDW/h\n', fbaSolution.f);
    results.fba.maxProduction = fbaSolution.f;
    results.fba.status = 'not_optimal';
end

% =========================================================================
% 9. Theoretical Yield Calculations
% =========================================================================
% Approximate molecular weights (g/mol)
mwHeparosan_unit = 379.3;  % Heparosan disaccharide repeat unit (GlcNAc + GlcA)
mwGlucose = 180.16;        % D-Glucose

% Find glucose uptake exchange reaction (try common BiGG IDs)
glcExchangeIDs = {'EX_glc__D_e', 'EX_glc_e', 'EX_glc_D_e'};
glcUptakeRate = 0;
glcExchangeRxn = '';

for i = 1:length(glcExchangeIDs)
    idx = find(strcmp(model.rxns, glcExchangeIDs{i}), 1);
    if ~isempty(idx)
        glcExchangeRxn = glcExchangeIDs{i};
        % Glucose uptake is typically a negative flux in the exchange reaction
        glcUptakeRate = abs(min(0, fbaSolution.x(idx)));
        break;
    end
end

fprintf('\n============================================================\n');
fprintf('  Theoretical Yield Analysis\n');
fprintf('============================================================\n');

if glcUptakeRate > 0
    yieldMolar = fbaSolution.f / glcUptakeRate;
    yieldMass = (fbaSolution.f * mwHeparosan_unit) / (glcUptakeRate * mwGlucose);

    fprintf('[OK] Glucose exchange reaction: %s\n', glcExchangeRxn);
    fprintf('     Glucose uptake rate:       %.4f mmol/gDW/h\n', glcUptakeRate);
    fprintf('     Molar yield (mol/mol):     %.4f\n', yieldMolar);
    fprintf('     Mass yield (g/g):          %.4f (%.2f%%)\n', yieldMass, yieldMass * 100);
    fprintf('     Theoretical max (stoich):  ~0.5 mol heparosan / mol glucose\n');

    results.yield.glucoseUptakeRate = glcUptakeRate;
    results.yield.glucoseExchangeRxn = glcExchangeRxn;
    results.yield.molarYield_molPerMol = yieldMolar;
    results.yield.massYield_gPerg = yieldMass;
else
    fprintf('[WARN] Glucose exchange reaction not found. Tried: %s\n', strjoin(glcExchangeIDs, ', '));
    fprintf('        Please manually specify the glucose uptake reaction for yield calculations.\n');
    results.yield.glucoseUptakeRate = NaN;
end

% =========================================================================
% 10. Flux Variability Analysis (FVA) at Optimal Production
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  FVA: Flux Variability at Optimal Heparosan Production\n');
fprintf('============================================================\n');
if fbaSolution.stat == 1 && fbaSolution.f > 0
    % Run FVA within 99% of the optimal solution (Use percentage: 99)
    optPercentage = 99;

    % Ensure rxnNameList is a valid cellstr
    rxnsForFVA = cellstr(model.rxns);

    % Call standard COBRA fluxVariability which returns two arrays
    [minFlux, maxFlux] = fluxVariability(model, optPercentage, 'max', rxnsForFVA);

    % Store into a struct to maintain compatibility with the rest of your script
    fvaSolution_struct = struct();
    fvaSolution_struct.min = minFlux;
    fvaSolution_struct.max = maxFlux;

    % Identify fixed-flux reactions (essential or blocked)
    rangeFlux = fvaSolution_struct.max - fvaSolution_struct.min;
    blockedReactions = model.rxns(abs(fvaSolution_struct.max) < 1e-8 & abs(fvaSolution_struct.min) < 1e-8);
    fixedReactions = model.rxns(rangeFlux < 1e-8 & ~(abs(fvaSolution_struct.max) < 1e-8 & abs(fvaSolution_struct.min) < 1e-8));

    fprintf('[OK] FVA completed (%d%% of optimum)\n', optPercentage);
    fprintf('     Blocked reactions:                %d\n', length(blockedReactions));
    fprintf('     Fixed-flux (essential) reactions: %d\n', length(fixedReactions));

    % Show FVA range for key heparosan pathway reactions
    keyRxns = {'HEPAROSAN_SYN', 'HEPAROSAN_t', 'EX_heparosan_e'};
    fprintf('\n     Key reaction flux ranges:\n');
    fprintf('     %-22s %14s %14s %14s\n', 'Reaction', 'Min flux', 'Max flux', 'Range');
    fprintf('     %s\n', repmat('-', 1, 68));
    for i = 1:length(keyRxns)
        idx = find(strcmp(model.rxns, keyRxns{i}), 1);
        if ~isempty(idx)
            fprintf('     %-22s %14.4f %14.4f %14.4f\n', ...
                keyRxns{i}, fvaSolution_struct.min(idx), fvaSolution_struct.max(idx), rangeFlux(idx));
        end
    end

    results.fva.solution = fvaSolution_struct;
    results.fva.optFraction = optPercentage / 100;
    results.fva.numBlocked = length(blockedReactions);
    results.fva.numFixed = length(fixedReactions);
else
    fprintf('[WARN] Skipping FVA: FBA solution was not optimal or has zero production.\n');
    results.fva.status = 'skipped';
end
% =========================================================================
% 11. Save Analysis Results
% =========================================================================
fprintf('\n============================================================\n');
fprintf('  Saving Analysis Results\n');
fprintf('============================================================\n');

% Save the model and analysis results together
save('iHM1533_heparosan_analysis.mat', ...
    'model', 'fbaSolution', 'results', '-v7.3');
fprintf('[OK] Analysis results saved to: iHM1533_heparosan_analysis.mat\n');

% Also export FBA flux distribution as a CSV table for easy inspection
fluxTable = table(model.rxns, model.rxnNames, fbaSolution.x, ...
    'VariableNames', {'Reaction', 'ReactionName', 'Flux'});
writetable(fluxTable, 'iHM1533_heparosan_fluxes.csv');
fprintf('[OK] FBA flux distribution saved to: iHM1533_heparosan_fluxes.csv\n');

% Save a human-readable summary
fid = fopen('iHM1533_heparosan_summary.txt', 'w');
fprintf(fid, 'iHM1533 Heparosan Production Analysis Summary\n');
fprintf(fid, 'Generated: %s\n\n', results.timestamp);
fprintf(fid, 'FBA Maximum Production: %.4f mmol/gDW/h\n', results.fba.maxProduction);
if isfield(results.yield, 'molarYield_molPerMol') && ~isnan(results.yield.molarYield_molPerMol)
    fprintf(fid, 'Molar Yield (mol heparosan / mol glucose): %.4f\n', results.yield.molarYield_molPerMol);
    fprintf(fid, 'Mass Yield (g heparosan / g glucose): %.4f\n', results.yield.massYield_gPerg);
end
if isfield(results.fva, 'numBlocked')
    fprintf(fid, 'FVA blocked reactions: %d\n', results.fva.numBlocked);
    fprintf(fid, 'FVA fixed-flux reactions: %d\n', results.fva.numFixed);
end
fclose(fid);
fprintf('[OK] Text summary saved to: iHM1533_heparosan_summary.txt\n');

% 7. Save and export the modified model into multiple formats

% 7.0 Ensure gene metadata fields are consistent with model.genes
% (Required because addReaction with geneRule adds new genes (kfiA/B/C, kpsD/M/T)
%  without initializing their metadata, causing verifyModel to fail.)
nGenes = length(model.genes);
geneMetadataFields = {'geneNames', 'proteins', 'geneEntrezID', 'geneUniprotID', ...
                      'geneEcoGeneID', 'geneASAPID'};

for i = 1:length(geneMetadataFields)
    field = geneMetadataFields{i};
    if ~isfield(model, field)
        % Field missing entirely: create with appropriate default type
        model.(field) = repmat({''}, nGenes, 1);
    elseif length(model.(field)) ~= nGenes
        % Size mismatch: resize while preserving type
        nKeep = min(length(model.(field)), nGenes);
        if isnumeric(model.(field))
            newField = zeros(nGenes, 1);
            newField(1:nKeep) = model.(field)(1:nKeep);
        else
            newField = repmat({''}, nGenes, 1);
            for j = 1:nKeep
                if iscell(model.(field))
                    newField{j} = model.(field){j};
                elseif ischar(model.(field))
                    newField{j} = model.(field)(j, :);
                end
            end
        end
        model.(field) = newField;
    end
end
fprintf('[OK] Gene metadata fields aligned (nGenes = %d)\n', nGenes);

% 7.1 Export as SBML (XML format) - Standard for metabolic models
try
    writeCbModel(model, 'format', 'sbml', 'fileName', 'iHM1533_heparosan.xml');
    fprintf('[OK] Exported: iHM1533_heparosan.xml\n');
catch ME
    fprintf('[ERROR] Failed to export SBML: %s\n', ME.message);
end

% 7.2 Export as MATLAB native format (.mat)
try
    writeCbModel(model, 'format', 'mat', 'fileName', 'iHM1533_heparosan.mat');
    fprintf('[OK] Exported: iHM1533_heparosan.mat\n');
catch ME
    fprintf('[ERROR] Failed to export MAT: %s\n', ME.message);
end

% 7.3 Export as Excel format (.xls / .xlsx)
try
    writeCbModel(model, 'format', 'xls', 'fileName', 'iHM1533_heparosan.xlsx');
    fprintf('[OK] Exported: iHM1533_heparosan.xlsx\n');
catch ME
    fprintf('[ERROR] Failed to export XLS: %s\n', ME.message);
end

fprintf('\n############################################################\n');
fprintf('#  Pipeline completed: %s\n', datestr(now));
fprintf('#  Output files:\n');
fprintf('#    - iHM1533_heparosan.xml / .mat / .xlsx (model)\n');
fprintf('#    - iHM1533_heparosan_analysis.mat         (results)\n');
fprintf('#    - iHM1533_heparosan_fluxes.csv           (fluxes)\n');
fprintf('#    - iHM1533_heparosan_summary.txt          (summary)\n');
fprintf('############################################################\n');