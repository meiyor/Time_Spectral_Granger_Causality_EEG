function [EEG_clean, beta_all, info] = remove_ecg_regression_lagged(EEG_data, ECG, fs, max_lag_sec, lambda)
%% remove_ecg_regression_lagged_stronger
%% Remove ECG contamination from EEG using lagged regression with ECG + dECG
%%
%% INPUTS:
%%   EEG_data    : [n_channels x n_samples]
%%   ECG         : [1 x n_samples] or [n_samples x 1]
%%   fs          : sampling rate (Hz)
%%   max_lag_sec : max lag in seconds (default 0.12)
%%   lambda      : ridge regularization (default 0; e.g. 1e-2 or 1e-1)
%%
%% OUTPUTS:
%%   EEG_clean   : cleaned EEG [n_channels x n_samples]
%%   beta_all    : regression coefficients per channel
%%   info        : struct with lags and design info
%% Description:
%%   Removes ECG artifacts using lagged regression with optional ridge penalty.
%%
%% Syntax:
%%   [EEG_clean, beta_all, info] = remove_ecg_regression_lagged(EEG_data, ECG, fs, max_lag_sec, lambda)
%%
%% Inputs:
%%   EEG_data    - EEG matrix.
%%   ECG         - ECG vector.
%%   fs          - Sampling rate.
%%   max_lag_sec - Max lag (seconds).
%%   lambda      - Ridge parameter.
%%
%% Outputs:
%%   EEG_clean - Cleaned EEG.
%%   beta_all  - Regression coefficients.
%%   info      - Metadata struct.
%%
%% Notes:
%%   Uses lagged ECG and derivative regressors. Minimal for now**
%% -------------------------------------------------------------------------

    if nargin < 4 || isempty(max_lag_sec)
        max_lag_sec = 0.12; % ±120 ms
    end
    if nargin < 5 || isempty(lambda)
        lambda = 0; % 0 = ordinary least squares
    end

    [nCh, nT] = size(EEG_data);

    ecg = double(ECG(:));
    if length(ecg) ~= nT
        error('ECG length must match EEG time dimension');
    end

    % -------- 1) preprocess ECG reference --------
    ecg = ecg - median(ecg);

    % Mild bandpass for ECG morphology / leakage structure
    try
        ecg_f = bandpass(ecg, [1 25], fs);
    catch
        ecg_f = ecg;
    end

    % z-score for numerical stability
    ecg_f = ecg_f / max(std(ecg_f), eps);

    % derivative reference
    decg = [0; diff(ecg_f)];
    decg = decg / max(std(decg), eps);

    % -------- 2) define lags --------
    max_lag = round(max_lag_sec * fs);
    lags = -max_lag:max_lag;
    nL = numel(lags);

    ECG_lag  = zeros(nT, nL);
    DECG_lag = zeros(nT, nL);

    % zero-padded lagging (no circular wrap)
    for i = 1:nL
        l = lags(i);

        if l < 0
            ECG_lag(1:end+l, i)  = ecg_f(1-l:end);
            DECG_lag(1:end+l, i) = decg(1-l:end);
        elseif l > 0
            ECG_lag(1+l:end, i)  = ecg_f(1:end-l);
            DECG_lag(1+l:end, i) = decg(1:end-l);
        else
            ECG_lag(:, i)  = ecg_f;
            DECG_lag(:, i) = decg;
        end
    end

    % -------- 3) build design matrix --------
    X = [ECG_lag, DECG_lag, ones(nT,1)];
    nReg = size(X,2);

    EEG_clean = zeros(nCh, nT);
    beta_all  = zeros(nReg, nCh);

    % Ridge matrix (do not penalize intercept)
    Ireg = eye(nReg);
    Ireg(end,end) = 0;

    % -------- 4) solve per channel --------
    for ch = 1:nCh
        x = double(EEG_data(ch,:)');
        x = x - mean(x);

        if lambda > 0
            beta = (X' * X + lambda * Ireg) \ (X' * x);
        else
            beta = X \ x;
        end

        beta_all(:,ch) = beta;

        ecg_est = X * beta;
        EEG_clean(ch,:) = (x - ecg_est)';
    end

    info.lags = lags;
    info.max_lag_samples = max_lag;
    info.lambda = lambda;
    info.n_regressors = nReg;
end