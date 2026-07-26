function uninstall
%UNINSTALL Remove MetabolicDesigner Pro (MDP) from MATLAB path.
%   Removes src/ (and subdirs) from the current session and saves pathdef.m
%   so the removal survives MATLAB restarts.
%
%   Use this when:
%     - You no longer need MDP on your MATLAB path
%     - The project was moved or deleted (run uninstall from a clone)
%     - You want to roll back a prior setup('permanent')
%
%   Locates the project root from this file's own location, so it works
%   from any current working directory.

    % --- 1. Locate project root ---
    thisFile   = mfilename('fullpath');
    projectDir = fileparts(thisFile);
    srcDir     = fullfile(projectDir, 'src');

    if exist(srcDir, 'dir') ~= 7
        error('uninstall:NoSrc', ...
              'Cannot find src/ directory at "%s". ' + ...
              'If MDP has been moved, run uninstall from any clone ' + ...
              'at the same path, or remove it manually via pathtool.', ...
              srcDir);
    end

    % --- 2. Remove from current session ---
    oldPath = genpath(srcDir);
    rmpath(oldPath);
    fprintf('[uninstall] Removed from current session:\n  %s\n', ...
            strrep(oldPath, pathsep, sprintf('\n  ')));

    % --- 3. Persist removal to pathdef.m ---
    try
        savepath;
        fprintf('[uninstall] Saved pathdef.m — removal survives MATLAB restart.\n');
    catch ME
        warning('uninstall:SavePathFailed', ...
                'Could not save pathdef.m (%s). ' + ...
                'Run uninstall from a MATLAB session started with write ' + ...
                'permission to pathdef.m, or remove entries via pathtool.', ...
                ME.message);
    end

    fprintf(newline + '[uninstall] Done.' + newline);
end
