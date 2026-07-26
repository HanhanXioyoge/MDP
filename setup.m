function setup(option)
%SETUP Add MetabolicDesigner Pro (MDP) to MATLAB path.
%   setup               Add src/ (and subdirs) to the current session.
%   setup('session')    Same as calling setup() with no argument.
%   setup('permanent')  Also persist to pathdef.m so the addition survives
%                       MATLAB restarts.
%
%   Locates the project root from this file's own location, so it works
%   no matter what the current working directory is when you call it.
%   Idempotent — running setup() twice does not duplicate entries.
%
%   After setup, you can call directly:
%       model    = loadModel('input_models/.../foo.mat')
%       result   = algFseof(model, target, biomass, 'Iterations', 21)
%       combined = strainDesign(modelPath, target, biomass, spec)
%
%   Use uninstall to remove. See README.md.

    arguments
        option {mustBeMember(option, {'session','permanent'})} = 'session'
    end

    % --- 1. Locate project root ---
    thisFile   = mfilename('fullpath');
    projectDir = fileparts(thisFile);
    srcDir     = fullfile(projectDir, 'src');

    if exist(srcDir, 'dir') ~= 7
        error('setup:NoSrc', ...
              'Cannot find src/ directory at "%s". Is setup.m in the project root?', ...
              srcDir);
    end

    % --- 2. Add src/ and its subdirs to current session ---
    newPath = genpath(srcDir);   % src/ + src/algorithms/
    addpath(newPath, '-begin');  % '-begin' => shadow conflicting names (defensive)
    fprintf('[setup] Added to current session:\n  %s\n', ...
            strrep(newPath, pathsep, sprintf('\n  ')));

    % --- 3. Optionally persist to pathdef.m ---
    if strcmpi(option, 'permanent')
        try
            savepath;
            fprintf('[setup] Saved to pathdef.m — survives MATLAB restart.\n');
        catch ME
            warning('setup:SavePathFailed', ...
                    'Could not save pathdef.m (%s). ' + ...
                    'Path is added for this session only. ' + ...
                    'Run MATLAB as administrator, or copy this setup.m to a writable location.', ...
                    ME.message);
        end
    end

    % --- 4. Friendliness ---
    fprintf(newline + '[setup] MDP ready. Callable now: ' + ...
            'loadModel, fseof, optknock, optforce, strainDesign.' + newline);
end
