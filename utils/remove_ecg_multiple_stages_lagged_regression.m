function [EEG_stage_corr_all, EEG_stage_raw_all, beta_all_stages, yhat_all_stages, stage_info] = remove_ecg_multiple_stages_lagged_regression(EEG, ECG_signal, stage_labels, cfg)
%% Remove ECG contamination from multiple annotated stages using lagged regression
%
% INPUTS
%   EEG         - preprocessed EEGLAB struct
%   ECG_signal  - preprocessed ECG vector, same length as EEG.data time axis
%   stage_labels- cell array / string array of stage labels, e.g. {'S4','S6'}
%   cfg         - optional struct
%
% OPTIONAL cfg fields
%   cfg.maxlag_ms        : lag window in ms (default 100)
%   cfg.lambda           : ridge penalty (default 10)
%   cfg.method           : 'ridge' or 'ols' (default 'ridge')
%   cfg.use_derivative   : include dECG lags (default true)
%   cfg.demean_channels  : demean EEG/ECG in segment (default true)
%   cfg.exclude_labels   : labels excluded from corrected output
%                          default {'EKG1','ECG','EKG'}
%   cfg.match_mode       : 'contains' or 'exact' (default 'contains')
%   cfg.use_all_matches  : if true, process all matching events for each label
%                          if false, process only first match (default false)
%
% OUTPUTS
%   EEG_stage_corr_all - cell array of corrected EEGLAB stage structures
%   beta_all_stages    - cell array of regression coefficient matrices
%   yhat_all_stages    - cell array of predicted ECG contamination
%   stage_info         - table summarizing processed stages

    if nargin < 4, cfg = struct; end
    if ~isfield(cfg,'maxlag_ms'),        cfg.maxlag_ms = 100; end
    if ~isfield(cfg,'lambda'),           cfg.lambda = 10; end
    if ~isfield(cfg,'method'),           cfg.method = 'ridge'; end
    if ~isfield(cfg,'use_derivative'),   cfg.use_derivative = true; end
    if ~isfield(cfg,'demean_channels'),  cfg.demean_channels = true; end
    if ~isfield(cfg,'exclude_labels'),   cfg.exclude_labels = {'EKG1','ECG','EKG'}; end
    if ~isfield(cfg,'match_mode'),       cfg.match_mode = 'contains'; end
    if ~isfield(cfg,'use_all_matches'),  cfg.use_all_matches = false; end

    if ischar(stage_labels) || isstring(stage_labels)
        stage_labels = cellstr(string(stage_labels));
    end

    if ~iscell(stage_labels) || isempty(stage_labels)
        error('stage_labels must be a non-empty cell array or string array.');
    end

    if ~isfield(EEG,'data') || ~isfield(EEG,'srate') || ~isfield(EEG,'event')
        error('EEG must contain data, srate, and event fields.');
    end

    ECG_signal = double(ECG_signal(:)');
    if numel(ECG_signal) ~= size(EEG.data,2)
        error('ECG_signal length (%d) must match EEG.data time length (%d).', ...
            numel(ECG_signal), size(EEG.data,2));
    end

    % First collect all requested event indices
    jobs = struct('stage_label', {}, 'event_idx', {});
    for s = 1:numel(stage_labels)
        this_label = stage_labels{s};
        idx = find_stage_event_local(EEG, this_label, cfg.match_mode);

        if isempty(idx)
            warning('Stage "%s" not found. Skipping.', char(string(this_label)));
            continue;
        end

        if ~cfg.use_all_matches
            idx = idx(1);
        end

        for k = 1:numel(idx)
            jobs(end+1).stage_label = char(string(this_label)); %#ok<AGROW>
            jobs(end).event_idx = idx(k);
        end
    end

    nJobs = numel(jobs);
    EEG_stage_corr_all = cell(nJobs,1);
    EEG_stage_raw_all = cell(nJobs,1);
    beta_all_stages    = cell(nJobs,1);
    yhat_all_stages    = cell(nJobs,1);

    info_stage_label = cell(nJobs,1);
    info_event_idx   = nan(nJobs,1);
    info_start_idx   = nan(nJobs,1);
    info_end_idx     = nan(nJobs,1);
    info_duration_s  = nan(nJobs,1);
    info_nchan       = nan(nJobs,1);

    for j = 1:nJobs
        [EEG_stage_corr_all{j}, EEG_stage_raw_all{j}, beta_all_stages{j}, yhat_all_stages{j}, meta] = ...
            remove_single_stage_by_eventidx_local(EEG, ECG_signal, jobs(j).stage_label, jobs(j).event_idx, cfg);

        info_stage_label{j} = meta.stage_label;
        info_event_idx(j)   = meta.event_idx;
        info_start_idx(j)   = meta.start_idx;
        info_end_idx(j)     = meta.end_idx;
        info_duration_s(j)  = meta.duration_s;
        info_nchan(j)       = meta.nchan;
    end

    stage_info = table(info_stage_label, info_event_idx, info_start_idx, info_end_idx, ...
        info_duration_s, info_nchan, ...
        'VariableNames', {'StageLabel','EventIdx','StartSample','EndSample','DurationSec','NChannels'});
end


function [EEG_stage_corr, EEG_stage_raw, beta_all, yhat_all, meta] = ...
    remove_single_stage_by_eventidx_local(EEG, ECG_signal, stage_label, event_idx, cfg)

    fs = EEG.srate;

    start_idx = round(EEG.event(event_idx).latency);

    dur = 0;
    if isfield(EEG.event(event_idx),'duration') && ~isempty(EEG.event(event_idx).duration)
        dur = round(EEG.event(event_idx).duration);
    end
    if dur <= 0
        error('Stage "%s" at event %d has duration <= 0.', char(string(stage_label)), event_idx);
    end

    end_idx = start_idx + dur - 1;
    start_idx = max(1, start_idx);
    end_idx   = min(size(EEG.data,2), end_idx);

    if end_idx <= start_idx
        error('Invalid bounds for stage "%s" at event %d.', char(string(stage_label)), event_idx);
    end

    % Channel selection
    chanlabels_all = get_chanlabels_local(EEG);
    keep_idx = true(1, size(EEG.data,1));

    for k = 1:numel(chanlabels_all)
        lbl = string(chanlabels_all{k});
        if any(strcmpi(lbl, string(cfg.exclude_labels)))
            keep_idx(k) = false;
        end
    end

    EEG_stage_raw = double(EEG.data(keep_idx, start_idx:end_idx));
    ECG_stage     = double(ECG_signal(start_idx:end_idx));
    chanlabels    = chanlabels_all(keep_idx);

    % Demean
    if cfg.demean_channels
        EEG_stage_raw = EEG_stage_raw - mean(EEG_stage_raw, 2);
        ECG_stage     = ECG_stage     - mean(ECG_stage, 2);
    end

    % Robust scale ECG
    ECG_stage = ECG_stage ./ max(mad(ECG_stage,1), eps);

    % Lagged design matrix
    maxlag = round(cfg.maxlag_ms * fs / 1000);
    lags = -maxlag:maxlag;

    X_ecg = build_lagged_matrix_local(ECG_stage, lags);
    X_parts = {X_ecg};

    if cfg.use_derivative
        dECG = [0 diff(ECG_stage)];
        dECG = dECG ./ max(mad(dECG,1), eps);
        X_d = build_lagged_matrix_local(dECG, lags);
        X_parts{end+1} = X_d;
    end

    X = cat(2, X_parts{:});
    X = [ones(size(X,1),1), X];

    % Regression
    nChan = size(EEG_stage_raw,1);
    T     = size(EEG_stage_raw,2);
    nReg  = size(X,2);

    EEG_stage_corr_data = zeros(size(EEG_stage_raw));
    yhat_all            = zeros(size(EEG_stage_raw));
    beta_all            = zeros(nReg, nChan);

    switch lower(cfg.method)
        case 'ridge'
            I = eye(nReg);
            I(1,1) = 0;
            XtX = X' * X;
            Xt  = X';

            for ch = 1:nChan
                y = EEG_stage_raw(ch,:).';
                beta = (XtX + cfg.lambda * I) \ (Xt * y);
                yhat = X * beta;

                beta_all(:,ch)            = beta;
                yhat_all(ch,:)            = yhat.';
                EEG_stage_corr_data(ch,:) = (y - yhat).';
            end

        case 'ols'
            for ch = 1:nChan
                y = EEG_stage_raw(ch,:).';
                beta = X \ y;
                yhat = X * beta;

                beta_all(:,ch)            = beta;
                yhat_all(ch,:)            = yhat.';
                EEG_stage_corr_data(ch,:) = (y - yhat).';
            end

        otherwise
            error('Unknown cfg.method "%s". Use "ridge" or "ols".', cfg.method);
    end

    % Build new stage EEG structure
    EEG_stage_corr = EEG;
    EEG_stage_corr.data   = EEG_stage_corr_data;
    EEG_stage_corr.nbchan = size(EEG_stage_corr_data,1);
    EEG_stage_corr.pnts   = size(EEG_stage_corr_data,2);
    EEG_stage_corr.trials = 1;
    EEG_stage_corr.srate  = fs;
    EEG_stage_corr.xmin   = 0;
    EEG_stage_corr.xmax   = (EEG_stage_corr.pnts - 1) / fs;
    EEG_stage_corr.times  = (0:EEG_stage_corr.pnts-1) / fs * 1000;

    if isfield(EEG_stage_corr,'chanlocs') && ~isempty(EEG_stage_corr.chanlocs)
        EEG_stage_corr.chanlocs = EEG_stage_corr.chanlocs(keep_idx);
        for c = 1:min(numel(EEG_stage_corr.chanlocs), numel(chanlabels))
            EEG_stage_corr.chanlocs(c).labels = chanlabels{c};
        end
    else
        EEG_stage_corr.chanlocs = struct('labels','');
        for c = 1:numel(chanlabels)
            EEG_stage_corr.chanlocs(c).labels = chanlabels{c};
        end
    end

    if isfield(EEG_stage_corr,'icachansind'), EEG_stage_corr.icachansind = []; end
    if isfield(EEG_stage_corr,'icaact'),      EEG_stage_corr.icaact = []; end
    if isfield(EEG_stage_corr,'icawinv'),     EEG_stage_corr.icawinv = []; end
    if isfield(EEG_stage_corr,'icasphere'),   EEG_stage_corr.icasphere = []; end
    if isfield(EEG_stage_corr,'icaweights'),  EEG_stage_corr.icaweights = []; end
    if isfield(EEG_stage_corr,'epoch'),       EEG_stage_corr.epoch = []; end
    if isfield(EEG_stage_corr,'specdata'),    EEG_stage_corr.specdata = []; end
    if isfield(EEG_stage_corr,'specicaact'),  EEG_stage_corr.specicaact = []; end
    if isfield(EEG_stage_corr,'stats'),       EEG_stage_corr.stats = []; end
    if isfield(EEG_stage_corr,'reject'),      EEG_stage_corr.reject = []; end

    EEG_stage_corr.event = struct('type', {}, 'latency', {}, 'duration', {}, 'urevent', {});
    EEG_stage_corr.urevent = struct('type', {}, 'latency', {}, 'duration', {});

    EEG_stage_corr.event(1).type     = char(string(stage_label));
    EEG_stage_corr.event(1).latency  = 1;
    EEG_stage_corr.event(1).duration = EEG_stage_corr.pnts;
    EEG_stage_corr.event(1).urevent  = 1;

    EEG_stage_corr.urevent(1).type     = char(string(stage_label));
    EEG_stage_corr.urevent(1).latency  = 1;
    EEG_stage_corr.urevent(1).duration = EEG_stage_corr.pnts;

    if isfield(EEG_stage_corr,'setname') && ~isempty(EEG_stage_corr.setname)
        EEG_stage_corr.setname = sprintf('%s_%s_evt%d_ecgLagClean', ...
            EEG_stage_corr.setname, char(string(stage_label)), event_idx);
    else
        EEG_stage_corr.setname = sprintf('%s_evt%d_ecgLagClean', ...
            char(string(stage_label)), event_idx);
    end

    if exist('eeg_checkset', 'file') == 2
        EEG_stage_corr = eeg_checkset(EEG_stage_corr, 'eventconsistency');
        EEG_stage_corr = eeg_checkset(EEG_stage_corr);
    end

    meta = struct;
    meta.stage_label = char(string(stage_label));
    meta.event_idx   = event_idx;
    meta.start_idx   = start_idx;
    meta.end_idx     = end_idx;
    meta.duration_s  = T / fs;
    meta.nchan       = nChan;
end


function idx = find_stage_event_local(EEG, stage_label, match_mode)
    idx = [];
    target = char(string(stage_label));

    for k = 1:numel(EEG.event)
        if ~isfield(EEG.event(k),'type') || isempty(EEG.event(k).type)
            continue;
        end

        evt = char(string(EEG.event(k).type));

        switch lower(match_mode)
            case 'exact'
                if strcmpi(strtrim(evt), strtrim(target))
                    idx(end+1) = k; %#ok<AGROW>
                end
            otherwise
                if contains(lower(strtrim(evt)), lower(strtrim(target)))
                    idx(end+1) = k; %#ok<AGROW>
                end
        end
    end
end


function labels = get_chanlabels_local(EEG)
    labels = cell(1, size(EEG.data,1));

    if isfield(EEG,'chanlocs') && ~isempty(EEG.chanlocs)
        for k = 1:min(numel(EEG.chanlocs), numel(labels))
            if isfield(EEG.chanlocs(k),'labels') && ~isempty(EEG.chanlocs(k).labels)
                labels{k} = EEG.chanlocs(k).labels;
            else
                labels{k} = sprintf('Ch%d', k);
            end
        end
    else
        for k = 1:numel(labels)
            labels{k} = sprintf('Ch%d', k);
        end
    end
end


function X = build_lagged_matrix_local(x, lags)
    x = double(x(:)');
    T = numel(x);
    nL = numel(lags);
    X = zeros(T, nL);

    for k = 1:nL
        X(:,k) = shift_with_zeros_local(x, lags(k)).';
    end
end


function y = shift_with_zeros_local(x, lag)
    x = x(:)';
    T = numel(x);
    y = zeros(1,T);

    if lag > 0
        y((lag+1):end) = x(1:end-lag);
    elseif lag < 0
        lag2 = abs(lag);
        y(1:end-lag2) = x((lag2+1):end);
    else
        y = x;
    end
end