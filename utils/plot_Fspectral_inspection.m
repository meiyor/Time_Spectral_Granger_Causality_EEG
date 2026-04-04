function plot_Fspectral_inspection(path_patient, list_stages, frequencies, clim_vals)
%% plot the imagesc for inspection without normalization and check activation
%% Description:
%%   Loads spectral GC results and plots matrices and directed connectivity.
%%
%% Syntax:
%%   plot_Fspectral_inspection(path_patient, list_stages, frequencies, clim_vals)
%%
%% Inputs:
%%   path_patient - File prefix.
%%   list_stages  - Cell array of stages.
%%   frequencies  - Cell array of frequency bands.
%%   clim_vals    - Color limits [min max].
%%
%% Outputs:
%%   None (plots generated).
%%
%% Notes:
%%   Requires saved GC results and channel coordinates.

close all
restoredefaultpath
rehash toolboxcache
clear gcd
which gcd -all

%% loading the paths
addpath(genpath([pwd, '/eeglab_current/']));
addpath(genpath([pwd '/MVGC1-master/']));

%% remove freemat conflicts
allp = strsplit(path, pathsep);
badmask = contains(allp, 'freemat3.5_DISABLED');
badpaths = allp(badmask);

rmpath(genpath('/home/jmm/Borjigin_Lab/MVGC1-master/utils/legacy/rng/rng.m'));
rmpath(genpath('/home/jmm/Borjigin_Lab/MVGC1-master/deprecated/utils/newline.m'));

if ~isempty(badpaths)
    rmpath(strjoin(badpaths, pathsep));
end

rehash toolboxcache
clear gcd

%% add data paths
addpath(genpath([pwd '/Data/']));
addpath(genpath([pwd '/preprocessed_save/']));

%% base file for chanlocs
X_base = load('JF_20250225_clean_and_splitted_S1.mat');
all_labels = {X_base.EEG_k.chanlocs.labels};

%% channel subset
tpo_labels = {'T3','T4','T5','T6','P3','P4','O1','O2','Pz','C3','C4','Cz','F3','F4','Fz'};
idx_tpo = find(ismember(all_labels, tpo_labels));
chnames = all_labels(idx_tpo);
chanlocs_tpo = X_base.EEG_k.chanlocs(idx_tpo);
nvars = numel(idx_tpo);

%% threshold
thr = clim_vals(1) + 0.5*diff(clim_vals);

%% =========================================================
%% FIXED 2D ELECTRODE COORDINATES FOR SCALP PLOT
%% use theta/radius if available; otherwise derive from X/Y
%% =========================================================
xy = zeros(nvars,2);

for k = 1:nvars
    hasThetaRadius = isfield(chanlocs_tpo(k),'theta') && ~isempty(chanlocs_tpo(k).theta) && ...
                     isfield(chanlocs_tpo(k),'radius') && ~isempty(chanlocs_tpo(k).radius);

    hasXY = isfield(chanlocs_tpo(k),'X') && ~isempty(chanlocs_tpo(k).X) && ...
            isfield(chanlocs_tpo(k),'Y') && ~isempty(chanlocs_tpo(k).Y);

    if hasThetaRadius
        % EEGLAB topoplot-like 2D coordinates
        th = deg2rad(chanlocs_tpo(k).theta);
        rr = chanlocs_tpo(k).radius;
        xy(k,1) = rr * sin(th);   % left-right
        xy(k,2) = rr * cos(th);   % posterior-anterior
    elseif hasXY
        % fallback: use X/Y and normalize
        xy(k,1) = chanlocs_tpo(k).X;
        xy(k,2) = chanlocs_tpo(k).Y;
    else
        error('Missing coordinates for channel %s.', chnames{k});
    end
end

%% center and scale to fit inside scalp
xy = xy - mean(xy,1);

rxy = sqrt(sum(xy.^2,2));
rmax = max(rxy);

if rmax < eps
    error('Electrode coordinates collapsed to zero.');
end

% fit comfortably inside the head circle
xy = 0.46 * xy / rmax;

%% =========================================================
%% ROTATE 90 DEGREES CLOCKWISE
%% so Fs go to top, Os to bottom, T4 right, T3 left
%% =========================================================
x_old = xy(:,1);
y_old = xy(:,2);

xy(:,1) = y_old;
xy(:,2) = -x_old;

% push all electrodes slightly posterior / downward - just slightly
xy(:,2) = xy(:,2) - 0.05;

for index_status = 1:length(list_stages)
    for index_freq = 1:length(frequencies)

        %% load GC result
        X_eval = load([pwd, '/GC_predictions/', path_patient, '_', ...
            list_stages{index_status}, '_', frequencies{index_freq}, ...
            '_gc_band_results.mat']);

        F_band = X_eval.F_band;
        F_band(1:nvars+1:end) = 0;

        %% ======================
        %% 1) MATRIX PLOT (UNCHANGED)
        %% ======================
        fig = figure('Color','w', ...
            'ToolBar','none', ...
            'MenuBar','none', ...
            'Units','pixels', ...
            'Position',[100 100 600 500]);

        ax = axes('Parent', fig);
        hold(ax, 'on');

        imagesc(ax, F_band);
        axis(ax, 'square');
        colorbar(ax, 'FontSize', 14);
        clim(ax, clim_vals);
        colormap(ax, parula);

        title_path_tex = strrep([path_patient, '_', list_stages{index_status}], '_', '\_');

        title(ax, sprintf('Pairwise spectral GC (p=%d) - $%s$ band - %s', ...
            20, frequencies{index_freq}, title_path_tex), ...
            'Interpreter','latex', 'FontSize',16);

        xlabel(ax, 'Source', 'FontSize',16);
        ylabel(ax, 'Target', 'FontSize',16);

        set(ax, ...
            'XTick', 1:nvars, ...
            'YTick', 1:nvars, ...
            'XTickLabel', chnames, ...
            'YTickLabel', chnames, ...
            'XTickLabelRotation', 90, ...
            'FontSize',14);

        %% ============================================
        %% 2) DARK TOPOPLOT-STYLE DIRECTED GC ARROW PLOT
        %% ============================================
        fig2 = figure('Color',[0.02 0.05 0.18], ...
            'ToolBar','none', ...
            'MenuBar','none', ...
            'Units','pixels', ...
            'Position',[800 100 700 650]);

        ax2 = axes('Parent', fig2);
        hold(ax2, 'on');
        axis(ax2, 'equal');
        axis(ax2, 'off');
        set(ax2, 'Color', [0.02 0.05 0.18]);

        draw_dark_head(ax2);

        cmapGY = parula(256);

        idx_CENTRAL = find(ismember(chnames, {'C3','C4','Cz'}));
        idx_TPO = find(ismember(chnames, {'T3','T4', 'T6', 'T5', 'P3', 'P4', 'Pz', 'O2', 'O1'}));
        idx_all = [idx_CENTRAL, idx_TPO];
        %% draw only strong arrows
        for itarget = idx_all
            for jsource = idx_all

                if jsource == itarget
                    continue;
                end

                val = F_band(itarget, jsource);

                if val >= thr
                    % source -> target
                    x1 = xy(jsource,1);
                    y1 = xy(jsource,2);
                    x2 = xy(itarget,1);
                    y2 = xy(itarget,2);

                    % shorten edges so arrowheads do not cover markers
                    alpha1 = 0.02;
                    alpha2 = 0.05;

                    xs = x1 + alpha1*(x2-x1);
                    ys = y1 + alpha1*(y2-y1);
                    xe = x2 - alpha2*(x2-x1);
                    ye = y2 - alpha2*(y2-y1);

                    cidx = round(1 + (size(cmapGY,1)-1) * ...
                        (val - clim_vals(1)) / max(eps, diff(clim_vals)));
                    cidx = min(max(cidx,1), size(cmapGY,1));
                    thisColor = cmapGY(cidx,:);

                    %% adjust line width depending on the intensity of the GC
                    lw = 2.0 + 3.0*(val - thr)/max(eps,(clim_vals(2)-thr));
                    lw = min(max(lw,1.5),4.5);

                    plot(ax2, [xs xe], [ys ye], '-', ...
                        'Color', thisColor, 'LineWidth', lw);

                    draw_arrowhead_filled(ax2, [xs ys], [xe ye], thisColor, 0.035, 0.015);
                end
            end
        end

        %% electrode markers
        scatter(ax2, xy(:,1), xy(:,2), 70, ...
            'MarkerFaceColor', [0.78 0.88 1.00], ...
            'MarkerEdgeColor', [1.00 1.00 1.00], ...
            'LineWidth', 2.0);

        %% labels
        for k = 1:nvars
            text(ax2, xy(k,1), xy(k,2)+0.038, chnames{k}, ...
                'Color', [0.95 0.98 1.00], ...
                'FontSize', 14, ...
                'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center');
        end

        title(ax2, sprintf('$%s$ - %s', ...
            frequencies{index_freq}, list_stages{index_status}), ...
            'Interpreter','latex', 'FontSize',18, 'Color', [0.95 0.98 1.00]);

        %% colorbar
        colormap(ax2, parula);
        clim(ax2, clim_vals);
        cb2 = colorbar(ax2, 'eastoutside');
        cb2.Color = [0.95 0.98 1.00];
        cb2.FontSize = 12;
        ylabel(cb2, 'F_{spectral}', 'Color', [0.95 0.98 1.00], 'FontSize', 13);

        %% fixed limits
        xlim(ax2, [-0.62 0.62]);
        ylim(ax2, [-0.62 0.68]);
    end
end
end

%% ==========================================================
%% helper: dark scalp outline
%% ==========================================================
function draw_dark_head(ax)

th = linspace(0, 2*pi, 400);
r  = 0.54;
headColor = [0.92 0.96 1.00];

% head circle
plot(ax, r*cos(th), r*sin(th), '-', ...
    'Color', headColor, 'LineWidth', 3.0);

% nose
plot(ax, [0.00 -0.06 0.00 0.06 0.00], ...
         [r    0.585 0.635 0.585 r], '-', ...
         'Color', headColor, 'LineWidth', 3.0);

% left ear
plot(ax, [-r -0.60 -0.585 -r], ...
         [0.12 0.07 -0.07 -0.12], '-', ...
         'Color', headColor, 'LineWidth', 3.0);

% right ear
plot(ax, [r 0.60 0.585 r], ...
         [0.12 0.07 -0.07 -0.12], '-', ...
         'Color', headColor, 'LineWidth', 3.0);
end

%% ==========================================================
%% helper: filled arrowhead
%% ==========================================================
function draw_arrowhead_filled(ax, p1, p2, col, ahlen, ahwid)

v = p2 - p1;
nv = norm(v);
if nv < eps
    return;
end
v = v / nv;

vp = [-v(2), v(1)];

tip   = p2;
base1 = p2 - ahlen*v + ahwid*vp;
base2 = p2 - ahlen*v - ahwid*vp;

patch(ax, ...
    [tip(1) base1(1) base2(1)], ...
    [tip(2) base1(2) base2(2)], ...
    col, ...
    'EdgeColor', col, ...
    'LineWidth', 1.5);
end