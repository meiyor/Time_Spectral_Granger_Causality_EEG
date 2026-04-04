function stacked_eeg_plot_struct_list(EEG_list, chan_idx, seg_list_names, ...
    offset, uv_bar, title_path)
% stacked_eeg_plot_struct_list
% Plot multiple EEGLAB structs sequentially on one stacked EEG figure.
%
% INPUTS:
%   EEG_list   : cell array of EEGLAB structs, e.g. {EEG1, EEG2, EEG3}
%   chan_idx   : channel indices to plot (default: all channels from first EEG)
%   seg_list_names: cells with strings of all the states for adding in plot
%   offset     : vertical spacing between channels
%   uv_bar     : amplitude scale bar in microvolts
%   title_path : filename prefix for saving
%
% NOTES:
% - Assumes all EEG structs have the same channels/order and same sampling rate.
% - Draws dashed vertical lines between consecutive structures.
% - Adds text labels EEG1, EEG2, ... above each segment.

    if ~iscell(EEG_list) || isempty(EEG_list)
        error('EEG_list must be a non-empty cell array of EEGLAB structs.');
    end

    nSeg = numel(EEG_list);
    EEG0 = EEG_list{1};

    if nargin < 2 || isempty(chan_idx)
        chan_idx = 1:EEG0.nbchan;
    end

    if nargin < 5 || isempty(title_path)
        title_path = 'stacked_structs';
    end

    fs = EEG0.srate;
    nChan = numel(chan_idx);

    % Check consistency
    for k = 1:nSeg
        EEG = EEG_list{k};
        if EEG.srate ~= fs
            error('All EEG structs must have the same sampling rate.');
        end
        if max(chan_idx) > EEG.nbchan
            error('chan_idx exceeds nbchan in EEG_list{%d}.', k);
        end
    end

    % Build concatenated matrix
    total_pnts = 0;
    for k = 1:nSeg
        total_pnts = total_pnts + EEG_list{k}.pnts;
    end

    X = zeros(nChan, total_pnts, 'double');
    boundaries = zeros(1, nSeg+1);
    labels_mid = zeros(1, nSeg);

    cursor = 1;
    boundaries(1) = 1;

    for k = 1:nSeg
        EEG = EEG_list{k};
        pnts_k = EEG.pnts;

        X(:, cursor:cursor+pnts_k-1) = double(EEG.data(chan_idx, :));

        labels_mid(k) = cursor + floor((pnts_k-1)/2);

        cursor = cursor + pnts_k;
        boundaries(k+1) = cursor;
    end

    t = (0:size(X,2)-1) / fs;

    if nargin < 3 || isempty(offset)
        offset = max(abs(X(:))) * 1.5;
        if offset == 0
            offset = 1;
        end
    end

    if nargin < 4 || isempty(uv_bar)
        uv_bar = round(max(abs(X(:))) / 2);
        if uv_bar == 0
            uv_bar = 1;
        end
    end

    fig = figure('Color','w','ToolBar','none','MenuBar','none');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    % Plot channels
    for ch = 1:nChan
        plot(ax, t, X(ch,:) + (ch-1)*offset, 'LineWidth', 1.2);
    end

    yticks(ax, (0:nChan-1)*offset);
    yticklabels(ax, {EEG0.chanlocs(chan_idx).labels});

    xlabel(ax, 'Time [s]', 'FontSize', 16);
    ylabel(ax, 'Channels', 'FontSize', 16);
    set(ax, 'FontSize', 16);
    box(ax, 'off');
    grid(ax, 'on');

    ylim_now = [(0)*offset - 0.5*offset, (nChan-1)*offset + 0.5*offset];
    ylim(ax, ylim_now);

    % Draw boundaries between structures
    y_text = ylim_now(2) - 0.04*(ylim_now(2)-ylim_now(1));

    for k = 2:nSeg
        x_boundary = (boundaries(k)-1) / fs;
        xline(ax, x_boundary, '--k', 'LineWidth', 1.5);
    end

    % Labels per segment
    for k = 1:nSeg
        x_mid = (labels_mid(k)-1) / fs;

        if isfield(EEG_list{k}, 'setname') && ~isempty(EEG_list{k}.setname)
            seg_name = seg_list_names{k};
        else
            seg_name = sprintf('EEG%d', k);
        end

        text(ax, x_mid, y_text, seg_name, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top', ...
            'FontSize', 12, ...
            'FontWeight', 'bold', ...
            'BackgroundColor', 'w', ...
            'Margin', 1);
    end

    % Amplitude scale bar
    xlim_now = xlim(ax);
    ylim_now = ylim(ax);

    x_range = xlim_now(2) - xlim_now(1);
    y_range = ylim_now(2) - ylim_now(1);

    x_bar = xlim_now(2) - 0.03 * x_range;
    y_bar_bottom = ylim_now(2) - 0.25 * y_range;
    y_bar_top = y_bar_bottom + uv_bar;

    line(ax, [x_bar x_bar], [y_bar_bottom y_bar_top], 'Color', 'k', 'LineWidth', 2);

    cap = 0.01 * x_range;
    line(ax, [x_bar-cap x_bar+cap], [y_bar_bottom y_bar_bottom], 'Color', 'k', 'LineWidth', 4);
    line(ax, [x_bar-cap x_bar+cap], [y_bar_top y_bar_top], 'Color', 'k', 'LineWidth', 4);

    text(ax, x_bar + 0.01*x_range, ...
        (y_bar_bottom + y_bar_top)/2, ...
        sprintf('%g \\muV', uv_bar), ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'Interpreter', 'tex');

    % Save
    outdir = fullfile(pwd, 'result_plots');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    savefig(fig, fullfile(outdir, [title_path '_stacked_time_series_events.fig']));
    saveas(fig, fullfile(outdir, [title_path '_stacked_time_series_events.jpg']));
end