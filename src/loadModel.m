function model = loadModel(path)
% loadModel  Load a constraint-based metabolic model from file.
%   model = loadModel(path) loads a model from the given file path.
%   Supported formats: .mat, .xml (SBML), .yml/.yaml/.json.
%   The returned model struct always has fields: rxns, S, rev, modelType.
%
%   Input:
%     path  - char vector, file path to the model file
%
%   Output:
%     model - struct, COBRA-compatible model with guaranteed fields

    arguments
        path  (1,:) char
    end

    [~, ~, ext] = fileparts(path);
    switch lower(ext)
        case '.mat'
            S = load(path);
            % Find the model variable: prefer 'model', else take first struct
            model = [];
            if isfield(S, 'model')
                model = S.model;
            else
                fnames = fieldnames(S);
                for i = 1:numel(fnames)
                    if isstruct(S.(fnames{i}))
                        model = S.(fnames{i});
                        break;
                    end
                end
            end
            if isempty(model)
                error('loadModel:NoModel', ...
                      'No model struct found in %s', path);
            end
        case '.xml'
            if exist('importSbml', 'file') == 2
                model = importSbml(path);
            else
                error('loadModel:NoSBML', ...
                      'importSbml not available. Install COBRA Toolbox.');
            end
        case {'.yml', '.yaml', '.json'}
            if exist('importYmlJson', 'file') == 2
                model = importYmlJson(path);
            elseif exist('readYaml', 'file') == 2
                model = readYaml(path);
            else
                error('loadModel:NoYaml', ...
                      'No YAML/JSON loader available for %s', path);
            end
        otherwise
            error('loadModel:BadExt', ...
                  'Unsupported file extension: %s', ext);
    end

    % Ensure rev field exists
    if ~isfield(model, 'rev')
        model.rev = double(model.lb < 0 & model.ub > 0);
    end

    % Validate essential fields
    if ~isfield(model, 'rxns') || ~isfield(model, 'S')
        error('loadModel:InvalidModel', ...
              'Model missing required fields: rxns or S');
    end

    % Sniff enzyme-constrained model type
    model.modelType = sniffEcModel(model);

    fprintf('[loadModel] Loaded %s (%d rxns, type=%s)\n', ...
            path, numel(model.rxns), model.modelType);
end

% === Local functions ===

function t = sniffEcModel(model)
% sniffEcModel  Detect whether a model is enzyme-constrained (ecGEM).
%   Score-based heuristic: enzymes, MWs, pathways, prot_pool_exchange.

    score = 0;
    if isfield(model, 'enzymes'),            score = score + 1; end
    if isfield(model, 'MWs'),                score = score + 1; end
    if isfield(model, 'pathways'),           score = score + 1; end
    if isfield(model, 'prot_pool_exchange'), score = score + 2; end
    if score >= 4
        t = 'enzyme_constrained';
    else
        t = 'stoichiometric';
    end
end
