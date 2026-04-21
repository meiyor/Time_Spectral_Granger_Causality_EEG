function X_eeg_recombined = recombine_cleaned_stage_segments(X_eeg_filtered, X_eeg_clean_s, do_taper, taper_ms)
% Reinsert stage-wise cleaned EEG excerpts back into the original
% continuous EEG structure using the original event latencies.
%
% Inputs
%   X_eeg_filtered : original continuous EEGLAB struct
%   X_eeg_clean_s  : cell array of cleaned stage EEGLAB structs
%                    each cell corresponds to one event/stage segment
%   do_taper       : true/false, whether to smooth seams
%   taper_ms       : taper window in milliseconds, e.g. 100
%
% Output
%   X_eeg_recombined : continuous EEGLAB struct with cleaned segments reinserted
%
% Notes
%   - Assumes each X_eeg_clean_s{k} was extracted from the original
%     X_eeg_filtered using the corresponding event latency/duration.
%   - Assumes same channel order and sampling rate.
%   - Uses the first N events, where N = numel(X_eeg_clean_s).

    if nargin < 3 || isempty(do_taper)
        do_taper = true;
    end
    if nargin < 4 || isempty(taper_ms)
        taper_ms = 100; % 100 ms default
    end

    % Start from original continuous signal
    X_eeg_recombined = X_eeg_filtered;

    fs = X_eeg_filtered.srate;
    taper_n = round((taper_ms / 1000) * fs);

    nSeg = numel(X_eeg_clean_s);
    nEvt = numel(X_eeg_filtered.event);

    if nSeg > nEvt
        error('More cleaned segments than events in X_eeg_filtered.event.');
    end

    for k = 1:nSeg
        EEG_k = X_eeg_clean_s{k};

        if isempty(EEG_k) || ~isstruct(EEG_k)
            warning('Skipping segment %d: empty or invalid struct.', k);
            continue;
        end

        evt = X_eeg_filtered.event(k);

        % Sample indices in continuous signal
        s1 = round(evt.latency);
        s2 = round(evt.latency + evt.duration - 1);

        % Basic checks
        if s1 < 1 || s2 > X_eeg_recombined.pnts || s2 < s1
            warning('Skipping segment %d: invalid sample bounds [%d %d].', k, s1, s2);
            continue;
        end

        seg_len_cont = s2 - s1 + 1;
        seg_len_clean = size(EEG_k.data, 2);

        if size(EEG_k.data, 1) ~= X_eeg_recombined.nbchan
            warning(['Skipping segment %d: channel count mismatch. ' ...
                     'Cleaned=%d, continuous=%d.'], ...
                     k, size(EEG_k.data, 1), X_eeg_recombined.nbchan);
            continue;
        end

        if seg_len_clean ~= seg_len_cont
            warning(['Segment %d length mismatch. ' ...
                     'Continuous=%d, cleaned=%d. Cropping to min length.'], ...
                     k, seg_len_cont, seg_len_clean);

            L = min(seg_len_cont, seg_len_clean);
            s2_use = s1 + L - 1;
            clean_data = EEG_k.data(:, 1:L);
        else
            s2_use = s2;
            clean_data = EEG_k.data;
            L = seg_len_cont;
        end

        % Original segment
        orig_data = X_eeg_recombined.data(:, s1:s2_use);

        % Hard replacement or tapered replacement
        if do_taper && L > 2 * taper_n && taper_n > 0
            % Build weights
            w = ones(1, L, 'single');

            ramp_up   = linspace(0, 1, taper_n);
            ramp_down = linspace(1, 0, taper_n);

            w(1:taper_n) = ramp_up;
            w(end-taper_n+1:end) = min(w(end-taper_n+1:end), ramp_down);

            % Expand to channel dimension
            W = repmat(w, size(clean_data, 1), 1);

            blended = (1 - W) .* orig_data + W .* clean_data;
            X_eeg_recombined.data(:, s1:s2_use) = blended;
        else
            X_eeg_recombined.data(:, s1:s2_use) = clean_data;
        end

        fprintf('Reinserted segment %d (%s): samples [%d %d]\n', ...
            k, char(string(evt.type)), s1, s2_use);
    end

    X_eeg_recombined.setname = [X_eeg_filtered.setname '_stageECGrecombined'];
    X_eeg_recombined = eeg_checkset(X_eeg_recombined);

end