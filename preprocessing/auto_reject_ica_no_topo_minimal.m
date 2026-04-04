function [bad_ic, feat] = auto_reject_ica_no_topo_minimal(EEG, max_removal, cfg)
    %% minimal ICA artifact removal allowing N worst ICs
    %% same flags, improved cardiac + EMG specificity
    %% Description:
    %%   Conservative ICA rejection returning limited number of worst components.
    %%
    %% Syntax:
    %%   [bad_ic, feat] = auto_reject_ica_no_topo_minimal(EEG, max_removal, cfg)
    %%
    %% Inputs:
    %%   EEG         - EEGLAB struct.
    %%   max_removal - Max ICs to reject.
    %%   cfg         - Optional config struct.
    %%
    %% Outputs:
    %%   bad_ic - Selected IC indices.
    %%   feat   - Feature table.
    %%
    %% Notes:
    %%   Uses robust scaling and spectral heuristics. IC severity is measured by flags of
    %%   flags of IC artifactual severity, min_removal and min_flags is defined by the user
    %% -------------------------------------------------------------------------

    if nargin < 3, cfg = struct; end
    if ~isfield(cfg,'kurt_thr'),       cfg.kurt_thr = 20; end
    if ~isfield(cfg,'hf_ratio_thr'),   cfg.hf_ratio_thr = 0.5; end
    if ~isfield(cfg,'line_ratio_thr'), cfg.line_ratio_thr = 0.3; end
    if ~isfield(cfg,'card_ac_thr'),    cfg.card_ac_thr = 0.5; end
    if ~isfield(cfg,'max_remove'),     cfg.max_remove = max_removal; end
    if ~isfield(cfg,'min_flags'),      cfg.min_flags = 1; end

    EEG = eeg_checkset(EEG, 'ica');

    W = EEG.icaweights * EEG.icasphere;
    icaact = W * double(EEG.data);
    fs = EEG.srate;
    nIC = size(icaact,1);

    kurt_v     = zeros(nIC,1);
    hf_ratio   = zeros(nIC,1);
    line_ratio = zeros(nIC,1);
    card_ac    = zeros(nIC,1);

    for i = 1:nIC
        x = double(icaact(i,:));
        x = x - mean(x);

        %% robust scaling to reduce effect of huge bursts dominating PSD
        xz = x / max(mad(x,1), eps);

        %% kurtosis
        kurt_v(i) = kurtosis(xz);

        %% PSD-based features
        [Pxx,F] = pwelch(xz, round(2*fs), round(fs), [], fs);
        Pxx = Pxx(:);
        F   = F(:);

        total_band = bandpower_from_psd(Pxx,F,1,80);
        hf_band    = bandpower_from_psd(Pxx,F,30,80);
        line_band  = bandpower_from_psd(Pxx,F,59,61);

        hf_ratio(i)   = hf_band   / max(total_band, eps);
        line_ratio(i) = line_band / max(total_band, eps);

        %% cardiac periodicity feature
        try
            x_card = bandpass(xz,[5 25],fs);
        catch
            x_card = xz;
        end

        maxLag = round(1.5*fs);
        ac = xcorr(x_card, maxLag, 'coeff');
        ac = ac(maxLag+1:end);

        %% plausible heart beat interval range: ~50-130 bpm
        lag0 = round(0.46*fs);
        lag1 = round(1.20*fs);
        lag1 = min(lag1, numel(ac));

        if lag0 < lag1
            seg = ac(lag0:lag1);
            [pks,~] = findpeaks(seg);
            if isempty(pks)
                card_ac(i) = 0;
            else
                card_ac(i) = max(pks);
            end
        else
            card_ac(i) = 0;
        end
    end

    %% same four flags
    flags = zeros(nIC,1);
    flags = flags + (kurt_v     > cfg.kurt_thr);
    flags = flags + (hf_ratio   > cfg.hf_ratio_thr);
    flags = flags + (line_ratio > cfg.line_ratio_thr);
    flags = flags + (card_ac    > cfg.card_ac_thr);

    %% weighted severity, same metrics but better prioritization
    sev = 0.5*max(0, (kurt_v     - cfg.kurt_thr)      / max(cfg.kurt_thr,eps)) + ...
          1.5*max(0, (hf_ratio   - cfg.hf_ratio_thr)  / max(cfg.hf_ratio_thr,eps)) + ...
          0.5*max(0, (line_ratio - cfg.line_ratio_thr)/ max(cfg.line_ratio_thr,eps)) + ...
          1.5*max(0, (card_ac    - cfg.card_ac_thr)   / max(cfg.card_ac_thr,eps));

    cand = find(flags >= cfg.min_flags);

    [~,ord] = sort(sev(cand), 'descend');
    cand = cand(ord);

    bad_ic = cand(1:min(cfg.max_remove, numel(cand)));

    disp("Removing ")
    disp(bad_ic)
    disp("ICs based on severity across the entire estimation and all the criteria")

    feat = table((1:nIC)', kurt_v, hf_ratio, line_ratio, card_ac, flags, sev, ...
        'VariableNames', {'IC','Kurtosis','HF_Ratio','Line60_Ratio','Cardiac_AC','NFlags','Severity'});
end

function bp = bandpower_from_psd(Pxx,F,f1,f2)
    idx = (F >= f1) & (F <= f2);
    if ~any(idx)
        bp = 0;
    else
        bp = trapz(F(idx), Pxx(idx));
    end
end