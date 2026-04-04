function [bad_ic, feat] = auto_reject_ica_no_topo(EEG, cfg)
%% Automatic IC rejection without topography/readlocs.
%% Description:
%%   Automatically flags ICA components using signal-based heuristics
%%   (kurtosis, HF power, line noise, cardiac autocorrelation).
%%
%% Syntax:
%%   [bad_ic, feat] = auto_reject_ica_no_topo(EEG, cfg)
%%
%% Inputs:
%%   EEG - EEGLAB struct with ICA computed.
%%   cfg - Optional threshold struct.
%%
%% Outputs:
%%   bad_ic - Indices of rejected ICs.
%%   feat   - Table of extracted features.
%%
%% Notes:
%%   No topography required. Uses ICA weights. This is differs on ADJUST on that matter
%% -------------------------------------------------------------------------

    if nargin < 2, cfg = struct; end
    if ~isfield(cfg,'kurt_thr'),      cfg.kurt_thr = 50; end
    if ~isfield(cfg,'hf_ratio_thr'),  cfg.hf_ratio_thr = 0.4; end
    if ~isfield(cfg,'line_ratio_thr'),cfg.line_ratio_thr = 0.25; end
    if ~isfield(cfg,'card_ac_thr'),   cfg.card_ac_thr = 0.5; end

    EEG = eeg_checkset(EEG, 'ica');

    W = EEG.icaweights * EEG.icasphere;
    icaact = W * double(EEG.data);   % nIC x nTime
    fs = EEG.srate;
    nIC = size(icaact,1);

    kurt_v     = zeros(nIC,1);
    hf_ratio   = zeros(nIC,1);
    line_ratio = zeros(nIC,1);
    card_ac    = zeros(nIC,1);

    for i = 1:nIC
        x = double(icaact(i,:));
        x = x - mean(x);

        % Kurtosis
        kurt_v(i) = kurtosis(x);

        % Welch PSD
        [Pxx,F] = pwelch(x, [], [], [], fs);
        Pxx = Pxx(:); F = F(:);

        total_band = bandpower_from_psd(Pxx,F,1,100);
        hf_band    = bandpower_from_psd(Pxx,F,20,100);
        line_band  = bandpower_from_psd(Pxx,F,58,62);

        hf_ratio(i)   = hf_band   / max(total_band, eps);
        line_ratio(i) = line_band / max(total_band, eps);

        % Cardiac-like periodicity: autocorr peak around ~0.4 to 1.5 s
        maxLag = round(1.5*fs);
        ac = xcorr(x, maxLag, 'coeff');
        ac = ac(maxLag+1:end); % nonnegative lags
        lag0 = round(0.4*fs);
        lag1 = round(1.5*fs);
        if lag1 > numel(ac), lag1 = numel(ac); end
        if lag0 < lag1
            card_ac(i) = max(ac(lag0:lag1));
        else
            card_ac(i) = 0;
        end
    end

    is_bad = (kurt_v > cfg.kurt_thr) | ...
             (hf_ratio > cfg.hf_ratio_thr) | ...
             (line_ratio > cfg.line_ratio_thr) | ...
             (card_ac > cfg.card_ac_thr);

    bad_ic = find(is_bad);

    feat = table((1:nIC)', kurt_v, hf_ratio, line_ratio, card_ac, is_bad, ...
        'VariableNames', {'IC','Kurtosis','HF_Ratio','Line60_Ratio','Cardiac_AC','Reject'});
end

function bp = bandpower_from_psd(Pxx,F,f1,f2)
    idx = (F >= f1) & (F <= f2);
    if ~any(idx)
        bp = 0;
    else
        bp = trapz(F(idx), Pxx(idx));
    end
end