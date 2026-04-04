function MVGC_application(X_eeg_clean, title_path, type_eval, p_order, freq)
%% ============================================================
%% MVGC ANALYSIS ON CLEANED EEG: X_eeg_clean
%% Add this block AFTER:
%% X_eeg_clean = pop_subcomp(X_eeg_clean, bad_ic, 0);
%% type_eval = string with the type of evaluation (channel-wise)
%% p_order = int defining the maximump value or maxorder for GC inference
%% freq = cell with the frequencies given in the paper for calculation
%% ============================================================
%% Description:
%%   Runs multivariate Granger causality (MVGC) analysis on cleaned EEG data.
%%   Extracts channels, preprocesses, segments into pseudo-trials, fits VAR,
%%   computes time-domain and spectral GC, and saves results.
%%
%% Syntax:
%%   MVGC_application(X_eeg_clean, title_path, type_eval, p_order, freq, process_single)
%%
%% Inputs:
%%   X_eeg_clean   - EEGLAB struct with cleaned EEG data.
%%   title_path    - String prefix for saved outputs.
%%   type_eval     - 'all' or 'custom' channel selection.
%%   p_order       - VAR model order.
%%   freq          - Cell array of frequency bands.
%%
%% Outputs:
%%   None (results saved to disk).
%%
%% Notes:
%%   Requires EEGLAB + MVGC toolbox. Assumes preprocessed EEG.

disp('--- Starting MVGC analysis on X_eeg_clean ---');
% ----------------------------
% 1) Extract cleaned data
% ----------------------------
all_labels = {X_eeg_clean.chanlocs.labels};

if strcmp(type_eval, 'custom') 
    tpo_labels = {'T3','T4','T5','T6','P3','P4','O1','O2','Pz','C3','C4','Cz','F3','F4','Fz'};
    idx_tpo = find(ismember(all_labels, tpo_labels));
    data_cont = double(X_eeg_clean.data(idx_tpo, :));
    chnames   = all_labels(idx_tpo);
end
if strcmp(type_eval, 'all') 
   data_cont = double(X_eeg_clean.data(:, :));
   chnames   = all_labels;
end
fs = X_eeg_clean.srate;

[nvars, nsamples_total] = size(data_cont);

fprintf('MVGC input continuous data size: %d channels x %d samples\n', nvars, nsamples_total);

% ----------------------------
% 2) Basic preprocessing for GC
% ----------------------------
% IMPORTANT:
% GC assumes covariance-stationary data.
% For continuous EEG, a common practical step is to:
%   - detrend each channel
%   - z-score per channel
%   - segment into short windows

for ch = 1:nvars
    x = data_cont(ch,:);
    x = detrend(x);                         % remove linear trend
    x = (x - mean(x)) / (std(x) + eps);     % z-score
    data_cont(ch,:) = x;
end

% ----------------------------
% 3) Segment continuous EEG into pseudo-trials
% ----------------------------
% You can tune these.
win_sec = 5.0;                       %% 5-second windows. This can be another important hyperparameter**
step_sec = 2.0;                      %% non-overlapping windows
win_len  = round(win_sec * fs);
step_len = round(step_sec * fs);

if win_len < 20
    error('Window too short. Increase win_sec.');
end

start_idx = 1:step_len:(nsamples_total - win_len + 1);
ntrials   = numel(start_idx);

if ntrials < 5
    error('Too few windows/trials for MVGC. Need more usable data or shorter windows.');
end

X = zeros(nvars, win_len, ntrials);

for k = 1:ntrials
    idx = start_idx(k):(start_idx(k) + win_len - 1);
    X(:,:,k) = data_cont(:, idx);
end

fprintf('Segmented into %d pseudo-trials of %d samples each (%.2f sec)\n', ...
    ntrials, win_len, win_sec);

% ----------------------------
% 4) Check for bad values
% ----------------------------
if any(~isfinite(X(:)))
    error('MVGC input contains NaN or Inf values.');
end

% ----------------------------
% 5) Choose model order
% ----------------------------
% Use information criteria to pick a VAR order.
% For EEG with 2 s windows, a modest order is usually safer.
% You may adjust momax depending on fs and data length.

%regmode = 'LWR';
%icregmode = 'LWR';
regmode = 'OLS';      % ordinary least squares
icregmode = 'LWR';    % recommended in MVGC demos for order selection
momax = p_order;           % try 10-30 depending on sampling rate/window length

fprintf('Estimating information criteria up to model order %d...\n', momax);

[AIC, BIC, moAIC, moBIC] = tsdata_to_infocrit(X, momax, icregmode, false);

% Pick BIC by default (more conservative)
morder  = moBIC; %min(moBIC, 3); % moBIC

% Safety fallback
if isempty(morder) || ~isfinite(morder) || morder < 1
    warning('Invalid model order from BIC. Falling back to order = 5');
    morder = 5;
end

fprintf('Selected model order (BIC): %d\n', morder);
fprintf('AIC suggested order: %d\n', moAIC);
fprintf('BIC suggested order: %d\n', moBIC);

% ----------------------------
% 6) Fit VAR model
% ----------------------------
fprintf('Fitting VAR model...\n');

[A, SIG] = tsdata_to_var(X, morder, regmode);

assert(~isbad(A),   'VAR estimation failed: bad A');
assert(~isbad(SIG), 'VAR estimation failed: bad SIG');

% ----------------------------
% 7) Autocovariance sequence
% ----------------------------
fprintf('Computing autocovariance sequence...\n');

[G, info] = var_to_autocov(A, SIG);
disp(info)

% var_to_autocov(A, SIG);

assert(~isbad(G), 'Autocovariance calculation failed');

if info.error
    warning('var_to_autocov reported error code: %d', info.error);
end

if info.aclags < info.acminlags
    warning('Autocovariance may not have converged well: aclags=%d, acminlags=%d', ...
        info.aclags, info.acminlags);
end

% ----------------------------
% 8) Spectral Radius analysis
% ----------------------------

rho = spectral_radius_var(A);

fprintf('Spectral radius: %.6f\n', rho);

if rho < 1
    disp('VAR is stable');
else
    warning('VAR is UNSTABLE');
end

% ----------------------------
% 9) Pairwise-conditional time-domain GC
% ----------------------------
fprintf('Computing pairwise-conditional GC...\n');

tStart = tic;

F = autocov_to_pwcgc(G);

elapsed_time = toc(tStart);
fprintf('Elapsed time in pairwise-conditional GC calculation: %.3f seconds\n', elapsed_time);

F = real(F);
F(F < 0) = 0;   % clip small negatives
Fmin = min(F(:));
Fmax = max(F(:));

F_prev = F; %% this for saving purposes
F = (F - Fmin) / (Fmax - Fmin + 1e-8);
F_norm = F;

assert(~isbad(F, false), 'GC calculation failed');

% F(i,j) = GC j -> i
% row = target, column = source

% Zero diagonal for display
F(1:nvars+1:end) = 0;

% ----------------------------
% 9) Statistical significance
% ----------------------------
% Significance test for pairwise-conditional GC
alpha  = 0.05;
mhtc   = 'Holm';   % false discovery rate (dependent) Holm more strict

fprintf('Computing significance thresholds...\n');

pval = mvgc_pval(F, morder, win_len, ntrials, 1, 1, nvars-2, 'F');
%% just to inspect the p-values here, most of the time they survive the correction in the pairwise comparison
pval
pval_corr = significance(pval, alpha, mhtc);
signif = 1-pval_corr < 0.05;

% ----------------------------
% 10) Plot time-domain GC matrix
% ----------------------------
fig = figure('Color','w', ...
    'ToolBar','none', ...
    'MenuBar','none', ...
    'Units','pixels', ...
    'Position',[100 100 600 500]); 

ax = axes('Parent', fig);
hold(ax, 'on');

imagesc(ax, F);
axis(ax, 'square');
colorbar(ax, 'FontSize', 14);

title(ax, [sprintf('Pairwise GC (p=%d)', morder), ' ', title_path], ...
    'Interpreter','none', 'FontSize',16);

xlabel(ax, 'Source', 'FontSize',16);
ylabel(ax, 'Target', 'FontSize',16);

% ticks
set(ax, ...
    'XTick', 1:nvars, ...
    'YTick', 1:nvars, ...
    'XTickLabel', chnames, ...
    'YTickLabel', chnames, ...
    'XTickLabelRotation', 90, ...
    'FontSize',14);  

savefig(fig, [pwd, '/result_plots/', title_path , '_p', num2str(momax), '_GC_results.fig']);
saveas(gcf, [pwd, '/result_plots/', title_path , '_p', num2str(momax), '_GC_results.jpg']); 

% ----------------------------
% 11) Pairwise-conditional spectral GC
% ----------------------------
fprintf('Computing pairwise-conditional spectral GC...\n');
% User-defined frequency band [Hz]
fres  = 2048;       % spectral resolution

tStart = tic;

%% use this measure for all the frequency bands given
Fspec = autocov_to_spwcgc(G, fres);

elapsed_time = toc(tStart);
fprintf('Elapsed time in pairwise-conditional spectral GC calculation: %.3f seconds %s \n', elapsed_time);

assert(~isbad(Fspec, false), 'Spectral GC calculation failed');

% Frequency axis
freqs = linspace(0, fs/2, fres+1);

%% plot across all the possible frequencies given by entry
for index_freq = 1:length(freq)
    %% define the frequecies here
    if strcmp(freq{index_freq}, 'delta')
        fband = [0,4];
        band_tex = '\delta';
    elseif strcmp(freq{index_freq}, 'theta')
        fband = [4,8];
        band_tex = '\theta';
    elseif strcmp(freq{index_freq}, 'alpha')
        fband = [8,13];
        band_tex = '\alpha';
    elseif strcmp(freq{index_freq}, 'beta')
        fband = [13,25];
        band_tex = '\beta';
    elseif strcmp(freq{index_freq}, 'gamma1')
        fband = [25,55];
        band_tex = '\gamma_{1}';
    elseif strcmp(freq{index_freq}, 'gamma2')
        fband = [80,100];
        band_tex = '\gamma_{2}';
    end;
    % Select desired frequency range
    fidx = (freqs >= fband(1)) & (freqs <= fband(2));
    
    if ~any(fidx)
        error('No frequency bins found in the requested band [%g %g] Hz.', fband(1), fband(2));
    end
    
    % Average spectral GC inside the selected band
    F_band = mean(Fspec(:,:,fidx), 3);
    
    % Optional: integrate instead of mean
    % F_band = trapz(freqs(fidx), Fspec(:,:,fidx), 3);
    
    F_band = real(F_band);
    F_band(F_band < 0) = 0;   % clip small negatives
    
    Fmin = min(F_band(:));
    Fmax = max(F_band(:));
    
    F_band_norm = (F_band - Fmin) / (Fmax - Fmin + 1e-8);
    
    % Zero diagonal for display
    F_band(1:nvars+1:end) = 0;
    F_band_norm(1:nvars+1:end) = 0;
    
    fprintf('Computed band-averaged spectral GC in [%g %g] Hz using %d bins.\n', ...
        fband(1), fband(2), nnz(fidx));

    % ----------------------------
    % 12) Plot time-domain spectral GC matrix inside the inner fold for
    % avoiding recalculation
    % ----------------------------

    fig = figure('Color','w', ...
        'ToolBar','none', ...
        'MenuBar','none', ...
        'Units','pixels', ...
        'Position',[100 100 600 500]); 
    
    ax = axes('Parent', fig);
    hold(ax, 'on');
    
    imagesc(ax, F_band_norm);
    axis(ax, 'square');
    colorbar(ax, 'FontSize', 14);
    
    title_path_tex = strrep(title_path, '_', '\_');
    
    title(ax, sprintf('Pairwise spectral GC (p=%d) - $%s$ band - %s', ...
        morder, band_tex, title_path_tex), ...
        'Interpreter','latex', 'FontSize',16);
    
    xlabel(ax, 'Source', 'FontSize',16);
    ylabel(ax, 'Target', 'FontSize',16);
    
    % ticks
    set(ax, ...
        'XTick', 1:nvars, ...
        'YTick', 1:nvars, ...
        'XTickLabel', chnames, ...
        'YTickLabel', chnames, ...
        'XTickLabelRotation', 90, ...
        'FontSize',14);  
    
    %% save the results and the plots
    save([pwd, '/GC_predictions/', title_path, '_', band_tex, '_p', num2str(momax), '_gc_band_results.mat'], 'F_band', 'F_band_norm');
    savefig(fig, [pwd, '/result_plots/', title_path , '_', freq{index_freq}, '_p', num2str(momax), '_spectral_GC_results.fig']);
    saveas(gcf, [pwd, '/result_plots/', title_path , '_', freq{index_freq}, '_p', num2str(momax), '_spectral_GC_results.jpg']);
end
save([pwd, '/GC_predictions/', title_path, '_p', num2str(momax), '_time_gc_results.mat'], 'F_prev', 'F_norm');