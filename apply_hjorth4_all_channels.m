function EEG = apply_hjorth4_all_channels(EEG)
%% APPLY_HJORTH4_ALL_CHANNELS
%% Hjorth / local Laplacian rereference for each EEG channel:
%%
%%   X_hjorth(ch,:) = X(ch,:) - mean(X(nearest 4 neighbors,:), 1)
%%
%% Nearest neighbors are computed from EEG.chanlocs XYZ coordinates.

    EEG = eeg_checkset(EEG);

    labels = string({EEG.chanlocs.labels});
    nChan  = EEG.nbchan;

    %% Get XYZ coordinates
    xyz = nan(nChan,3);
    for ch = 1:nChan
        xyz(ch,:) = [EEG.chanlocs(ch).X, EEG.chanlocs(ch).Y, EEG.chanlocs(ch).Z];
    end

    data_old = double(EEG.data);
    data_new = nan(size(data_old));

    %% get the neighbor table here
    neighbor_table = struct([]);

    for ch = 1:nChan

        this_label = labels(ch);

        valid = true(nChan,1);
        valid(ch) = false;
        valid(any(isnan(xyz),2)) = false;

        idx_valid = find(valid);

        if numel(idx_valid) < 4
            warning('Channel %s has fewer than 4 valid neighbors. Keeping original signal.', this_label);
            data_new(ch,:) = data_old(ch,:);
            continue
        end

        %% Find 4 nearest channels in 3D electrode space
        d = sqrt(sum((xyz(idx_valid,:) - xyz(ch,:)).^2, 2));
        [~, ord] = sort(d, 'ascend');

        neigh = idx_valid(ord(1:4));
        fprintf('Implementing Hjorth subtraction for channel %s using channels %s\n', this_label, strjoin(cellstr(labels(neigh)), ', '));

        %% do the Hjorth subtraction here
        data_new(ch,:) = data_old(ch,:) - mean(data_old(neigh,:), 1);

        neighbor_table(ch).target = char(this_label);
        neighbor_table(ch).neighbors = cellstr(labels(neigh));
    end

    EEG.data = data_new;
    EEG.setname = [EEG.setname '_Hjorth4All'];

    EEG.etc.hjorth4_all.reference_type = ...
        'Hjorth local reference: each channel minus average of 4 nearest neighbors';
    EEG.etc.hjorth4_all.neighbor_table = neighbor_table;

    %% After this transform, ICA matrices no longer correspond to data, this will be done at the end of the preprocessing so this won't be necessary
    EEG.icaweights = [];
    EEG.icasphere = [];
    EEG.icawinv = [];
    EEG.icachansind = [];

    EEG = eeg_checkset(EEG);
end