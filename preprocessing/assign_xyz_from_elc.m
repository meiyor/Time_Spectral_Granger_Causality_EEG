function EEG = assign_xyz_from_elc(EEG, elcfile)
%% Assign X/Y/Z coordinates to EEG.chanlocs by matching labels to an ASA .elc file
%% Description:
%%   Assigns X/Y/Z coordinates to EEG.chanlocs using an ASA .elc file.
%%
%% Syntax:
%%   EEG = assign_xyz_from_elc(EEG, elcfile)
%%
%% Inputs:
%%   EEG     - EEGLAB struct with channel labels.
%%   elcfile - Path to .elc coordinate file.
%%
%% Outputs:
%%   EEG     - Updated struct with coordinates.
%%
%% Notes:
%%   Matching is case-insensitive. Missing labels are reported.
%% -------------------------------------------------------------------------
    % Read whole file
    txt = fileread(elcfile);
    lines = regexp(txt, '\r\n|\n|\r', 'split');

    % Find key sections
    idx_num = find(startsWith(strtrim(lines), 'NumberPositions='), 1);
    idx_pos = find(strcmp(strtrim(lines), 'Positions'), 1);
    idx_lab = find(strcmp(strtrim(lines), 'Labels'), 1);

    if isempty(idx_num) || isempty(idx_pos) || isempty(idx_lab)
        error('Could not find NumberPositions / Positions / Labels sections in %s', elcfile);
    end

    % Number of positions
    num_str = strtrim(strrep(lines{idx_num}, 'NumberPositions=', ''));
    npos = str2double(num_str);

    if isnan(npos) || npos <= 0
        error('Could not parse NumberPositions in %s', elcfile);
    end

    % Parse XYZ rows
    XYZ = nan(npos, 3);
    for k = 1:npos
        row = strtrim(lines{idx_pos + k});
        vals = sscanf(row, '%f %f %f');
        if numel(vals) ~= 3
            error('Could not parse XYZ on row %d after Positions', k);
        end
        XYZ(k, :) = vals(:).';
    end

    % Parse labels
    tmplabels = cell(npos,1);
    for k = 1:npos
        tmplabels{k} = strtrim(lines{idx_lab + k});
    end

    % Match and assign
    eeglabels = {EEG.chanlocs.labels};
    for k = 1:numel(eeglabels)
        ix = find(strcmpi(strtrim(eeglabels{k}), tmplabels), 1);
        if ~isempty(ix)
            EEG.chanlocs(k).X = XYZ(ix,1);
            EEG.chanlocs(k).Y = XYZ(ix,2);
            EEG.chanlocs(k).Z = XYZ(ix,3);
        else
            fprintf('No template match for channel: %s\n', eeglabels{k});
        end
    end
end