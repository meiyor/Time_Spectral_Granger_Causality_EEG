function MVGC_application_regularized(X_eeg_clean, title_path, type_eval)
%% ============================================================
%% MVGC ANALYSIS ON CLEANED EEG: X_eeg_clean
%% Add this block AFTER:
%% X_eeg_clean = pop_subcomp(X_eeg_clean, bad_ic, 0);
%% type_eval = string with the type of evaluation (channel-wise)
%% ============================================================

disp('--- Starting MVGC analysis on X_eeg_clean ---');

% ----------------------------
% 1) Extract cleaned data
% ----------------------------
all_labels = {X_eeg_clean.chanlocs.labels};
if strcmp(type_eval, 'custom') 
    tpo_labels = {'T3','T4','T5','T6','P3','P4','O1','O2','Pz','C3','C4'};
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
win_sec = 5.0;                       % 5-second windows
step_sec = 2.0;                      % overlapping windows
%win_sec = 5.0;                       % 10-second windows
%step_sec = 2.0;                      % non-overlapping windows
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
% 5) Choose model order and do cusrtom VAR regression
% ----------------------------
% Use information criteria to pick a VAR order.
% For EEG with 2 s windows, a modest order is usually safer.
% You may adjust momax depending on fs and data length.

lambda = 0.05;   % tune this
p=2; %% choose p in the value you consider
[A, SIG] = ridge_var_regularized(X, p, lambda);

assert(~isbad(A),   'VAR estimation failed: bad A');
assert(~isbad(SIG), 'VAR estimation failed: bad SIG');

% ----------------------------
% 6) Autocovariance sequence
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
% 7) Spectral Radius analysis
% ----------------------------

rho = spectral_radius_var(A);

fprintf('Spectral radius: %.6f\n', rho);

if rho < 1
    disp('VAR is stable');
else
    warning('VAR is UNSTABLE');
end

% ----------------------------
% 8) Pairwise-conditional time-domain GC
% ----------------------------
fprintf('Computing pairwise-conditional GC...\n');

F = autocov_to_pwcgc(G);

F = real(F);
F(F < 0) = 0;   % clip small negatives
Fmin = min(F(:));
Fmax = max(F(:));

F = (F - Fmin) / (Fmax - Fmin + 1e-8);

assert(~isbad(F, false), 'GC calculation failed');

% F(i,j) = GC j -> i
% row = target, column = source

% Zero diagonal for display
F(1:nvars+1:end) = 0;

% ----------------------------
% 10) Statistical significance
% ----------------------------
% Significance test for pairwise-conditional GC
alpha  = 0.05;
mhtc   = 'Holm';   % false discovery rate (dependent) Holm more strict

fprintf('Computing significance thresholds...\n');

pval = mvgc_pval(F, p, win_len, ntrials, 1, 1, nvars-2, 'F');
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

title(ax, [sprintf('Pairwise GC (m=%d)', p), ' ', title_path], ...
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

savefig(fig, [pwd, '/result_plots/', title_path , '_GC_results.fig']);
saveas(gcf, [pwd, '/result_plots/', title_path , '_GC_results.jpg']); 

%% DON'T DO THIS FOR NOW**
% ==== overlay significance ====
%[row_sig, col_sig] = find(signif);
%plot(ax, col_sig, row_sig, 'k*', 'MarkerSize', 8);
%hold(ax, 'off');
X=1;