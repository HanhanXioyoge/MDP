# MetabolicDesigner Pro MVP — Design Spec

**Date:** 2026-07-24
**Status:** Approved
**Scope:** MVP rebirth — 5 MATLAB functions, 0 classes, 3 algorithms

---

## 1. Overview

MetabolicDesigner Pro (MDP) is a MATLAB-based platform for constraint-based metabolic network in silico strain design. The MVP provides three classic algorithms (FSEOF, OptKnock, optForce) as standalone functions with a unified result structure and workspace-based caching.

**Key decisions:**
- Pure COBRA/RAVEN wrapper approach (no custom solver fallback in MVP)
- 5 function files, 0 classes, 0 +package folders
- English code, Chinese README
- Result structure is the consumer contract — future filter/report/GUI read from it

---

## 2. Directory Structure

```
<project_root>/
├── test_Fseof.m              # Single-algo debug script
├── iHM1533_HEPAROSAN.m       # Main test script (batch + workspace)
├── README.md                  # Chinese README
├── .gitignore                 # Exclude workspaces/, run-*
├── input_models/iHM1533/
│   └── iHM1533_heparosan.mat  # Test model (user provides)
├── workspaces/                # Runtime output (gitignored)
└── src/
    ├── loadModel.m            # Load model + ec-model sniff
    ├── strainDesign.m         # Batch entry + workspace management
    └── algorithms/
        ├── fseof.m            # FSEOF algorithm
        ├── optknock.m         # OptKnock algorithm
        └── optforce.m        # optForce algorithm
```

- `algorithms/` is a path group, NOT a namespace
- `addpath(genpath('src'))` includes all subdirectories
- Call as `fseof()` or `import algorithms.fseof` then `fseof()`

---

## 3. Result Structure Contract

### 3.1 Inner: algoResult (per-algorithm output)

```matlab
algoResult = struct( ...
    'algorithm',         'FSEOF', ...           % char
    'candidates',        candidatesTable, ...    % table, 11 columns
    'config',            struct(...), ...        % actual params used
    'diagnostics',       struct(...), ...        % wtBiomass, targetMaxFlux
    'fallbackUsed',      false, ...             % logical
    'failureReason',     '', ...                % char, empty on success
    'algorithmSpecific', struct(...) ...         % algo-private fields
);
```

### 3.2 Candidates Table (11 columns, shared framework)

| Column | Type | Description |
|--------|------|-------------|
| reaction | cellstr | Reaction ID (required) |
| name | cellstr | Full name |
| equation | cellstr | Formatted equation |
| direction | double | +1 / -1 |
| intervention | cellstr | 'OE' / 'KD' / 'KO' |
| score | double | Algorithm-specific score |
| scoreLabel | cellstr | Score meaning description |
| genes | cellstr | Gene IDs (comma-separated) |
| geneNames | cellstr | Standard gene names (comma-separated) |
| ncbiId | cellstr | NCBI IDs (comma-separated) |

Missing columns are filled with empty values. The 11-column schema is a framework — each algorithm fills what it can. `algorithmSpecific` preserves raw algorithm-specific data.

**Note:** This schema may be adjusted after seeing actual algorithm outputs.

### 3.3 algorithmSpecific (per-algorithm)

**FSEOF:**
```matlab
struct('slopes', [...], 'fluxLevels', [...], 'fluxMatrix', [...], 'wtSlope', [...])
```

**OptKnock:**
```matlab
struct('optKnockSol', solStructOrEmpty, 'selectedRxnList', {{...}})
```

**optForce:**
```matlab
struct('optForceSets', {{...}}, 'typeRegOptForceSets', {{...}}, 'fluxOptForceSets', {{...}})
```

### 3.4 Outer: combined (workspace-persisted)

```matlab
combined = struct( ...
    'modelName',  'iHM1533_heparosan', ...
    'biomassRxn', 'BIOMASS_EcN_iHM1533_core_59p80M', ...
    'targetRxn',  'EX_heparosan_e', ...
    'timestamp',  '2026-07-24T15:30:00', ...
    'algorithms', struct('FSEOF', ..., 'OptKnock', ..., 'optForce', ...) ...
);
```

4 fixed fields + 1 algorithms sub-struct.

---

## 4. Function Signatures

```matlab
function model = loadModel(path)
function result = fseof(model, targetRxn, biomassRxn, varargin)
function result = optknock(model, targetRxn, biomassRxn, varargin)
function result = optforce(model, targetRxn, biomassRxn, varargin)
function combined = strainDesign(modelPath, targetRxn, biomassRxn, algoSpec)
```

Positional args (model/targetRxn/biomassRxn) = "what to do"; optional varargin = "how to do it".

### 4.1 Algorithm Parameters

**fseof:**
- `'Iterations'` (default: 10) — scan steps
- `'Coefficient'` (default: 0.9) — target max flux ratio (must be < 1)

**optknock:**
- `'MaxCandidates'` (default: 200)
- `'NumDel'` (default: 5)
- `'MinGrowthFraction'` (default: 0.1)
- `'VMax'` (default: 1000)

**optforce:**
- `'K'` (default: 2)
- `'NSets'` (default: 1)
- `'MaxCandidates'` (default: 500)

---

## 5. Workspace Structure

```
workspaces/
└── {safeModelBase}_{safeTargetRxn}/
    ├── combined_result.mat        ← All algorithms combined
    ├── FSEOF/
    │   └── result.mat             ← FSEOF per-algorithm result
    ├── OptKnock/
    │   └── result.mat             ← OptKnock per-algorithm result
    └── optForce/
        └── result.mat             ← optForce per-algorithm result
```

- `safeModelBase` = `regexprep(modelBase, '[^\w\-]', '_')`
- `safeTargetRxn` = `regexprep(targetRxn, '[^\w\-]', '_')`
- Same modelBase + targetRxn share the workspace folder
- Re-running an algorithm: overwrite its subfolder `result.mat` + update `combined_result.mat`
- Other algorithms' results preserved
- First run: initialize combined (algorithms = struct()), create subfolder

---

## 6. Algorithm Implementation Strategy

### 6.1 fseof.m

1. Parse options (Iterations, Coefficient)
2. Validate: target/biomass in model, coefficient < 1
3. Inline FBA via `optimizeCbModel` → wtBiomass, targetMaxFlux
4. Scan: enforce target flux at `linspace(0, targetMax*coefficient, iterations)` levels, maximize biomass → fluxMatrix
5. If RAVEN `FSEOF()` available → use for candidate detection; else → self-decode from fluxMatrix slopes
6. Build candidates table + result struct

### 6.2 optknock.m

1. Parse options (MaxCandidates, NumDel, MinGrowthFraction, VMax)
2. Validate: target/biomass in model
3. Call COBRA `OptKnock()` with appropriate parameters
4. Map COBRA output → candidates table (KO reactions, direction=-1, intervention='KO')
5. Build result struct

### 6.3 optforce.m

1. Parse options (K, NSets, MaxCandidates)
2. Validate: target/biomass in model
3. Call COBRA `optForce()` or equivalent
4. Map output → candidates table (with OE/KD/KO intervention types)
5. Build result struct

### 6.4 loadModel.m

1. Detect extension (.mat / .xml / .yml / .json)
2. Load with appropriate loader
3. Ensure `rev` field exists
4. Validate `rxns` and `S` fields
5. Sniff ec-model type → set `model.modelType`

### 6.5 strainDesign.m

1. Load model via `loadModel`
2. Compute workspace dir from modelBase + targetRxn
3. Load existing combined or initialize new
4. Loop through algoSpec, run each algorithm
5. Save per-algorithm result to subfolder `result.mat`
6. Update combined, save to `combined_result.mat`

---

## 7. Error Handling

- Algorithm failure → return `fallbackUsed=true` + `failureReason`, candidates as empty table
- Does NOT block the strainDesign pipeline — other algorithms continue
- Model validation failure → `failureResult()` immediate return
- Workspace model mismatch → warning, not blocking
- All error identifiers: three-segment format `<file>:<Reason>`

---

## 8. Test Scripts

### test_Fseof.m (single-algo debug, no workspace)

```matlab
addpath(genpath('src'));
model = loadModel('input_models/iHM1533/iHM1533_heparosan.mat');
result = fseof(model, 'EX_heparosan_e', ...
    'BIOMASS_EcN_iHM1533_core_59p80M', ...
    'Iterations', 21, 'Coefficient', 0.99);
fprintf('WT biomass: %g\n', result.diagnostics.wtBiomass);
fprintf('Target max: %g\n', result.diagnostics.targetMaxFlux);
fprintf('Candidates: %d\n', height(result.candidates));
disp(result.candidates(1:min(10,end), :));
```

### iHM1533_HEPAROSAN.m (batch + workspace)

```matlab
addpath(genpath('src'));
spec = { ...
    struct('algo','fseof','params',{{'Iterations',21,'Coefficient',0.99}}), ...
    struct('algo','optknock','params',{{'MaxCandidates',150,'NumDel',3}}), ...
    struct('algo','optforce','params',{{'K',3,'NSets',2}}) ...
};
combined = strainDesign( ...
    'input_models/iHM1533/iHM1533_heparosan.mat', ...
    'EX_heparosan_e', ...
    'BIOMASS_EcN_iHM1533_core_59p80M', ...
    spec);
```

---

## 9. Coding Conventions

- Function names: lowercase (fseof / optknock / optforce / loadModel / strainDesign)
- File name = function name
- `arguments` block compatible with R2019b
- Comments/logs/errors in English; README in Chinese
- Error identifiers: three-segment `<file>:<Reason>`
- Local functions at file end (no separate utility files)
- UTF-8 encoding, LF line endings

---

## 10. Out of Scope (MVP)

- Any classdef
- Any +xxx package
- Separate utility files
- Filtering, combination, Pareto, reporting
- PipelineManager / ConfigBuilder / CLI_Runner / AlgorithmRegistry / workspace_manager
- Sweep mode, checkpoint system
- EC-model algorithm variants

---

## 11. Future Interface Slots (not implemented)

| File | Interface | When |
|------|-----------|------|
| filterCandidates.m | f(algoResult) → algoResult | Add filtering |
| combineCandidates.m | f(filtered results) → Pareto table | Add combination |
| writeReport.m | f(combined) → HTML | Add reporting |
| algorithms/fseof_ec.m | Enzyme-constrained FSEOF | Add EC path |
| app/DesignerApp.m | GUI class | Add GUI |
| clearWorkspace.m | f(workspace) → clear | Manual cleanup |
