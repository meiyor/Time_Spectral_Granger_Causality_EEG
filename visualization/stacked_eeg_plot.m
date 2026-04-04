function stacked_eeg_plot(EEG, chan_idx, sample_range, offset, uv_bar, title_path)
%% Description:
%%   Plots stacked EEG channels over time.
%%
%% Syntax:
%%   stacked_eeg_plot(EEG, chan_idx, sample_range, offset, uv_bar, title_path)
%%
%% Inputs:
%%   EEG          - EEGLAB struct.
%%   chan_idx     - Channels to plot.
%%   sample_range - Time window.
%%   offset       - Vertical spacing.
%%   uv_bar       - Scale bar.
%%   title_path   - Plot title.
%%
%% Outputs:
%%   None.
%%
%% Notes:
%%   Useful for visual inspection adjusting the sample_range as desired.
%% -------------------------------------------------------------------------
%% stacked_eeg_plot - plot EEG channels stacked with spacing
%%
%% EEG           : EEGLAB struct
%% chan_idx      : channels to plot (e.g., 1:19)
%% sample_range  : samples (e.g., 1:20000)
%% offset        : spacing
%% uv_bar        : amplitude bar in microvolts

    if nargin < 2 || isempty(chan_idx)
        chan_idx = 1:EEG.nbchan;
    end

    if nargin < 3 || isempty(sample_range)
        sample_range = 1:EEG.pnts;
    end

    % force row vector and keep valid indices only
    sample_range = sample_range(:)';
    sample_range = sample_range(sample_range >= 1 & sample_range <= EEG.pnts);

    if isempty(sample_range)
        error('sample_range is empty after bounds checking.');
    end

    % Optional protection against giant plots
    max_plot_points = 50000;
    if numel(sample_range) > max_plot_points
        warning('Too many samples requested (%d). Truncating to first %d samples for plotting.', ...
            numel(sample_range), max_plot_points);
        sample_range = sample_range(1:max_plot_points);
    end

    data = double(EEG.data(chan_idx, sample_range));
    fs = EEG.srate;

    t = (sample_range - sample_range(1)) / fs;

    if nargin < 4 || isempty(offset)
        offset = max(abs(data(:))) * 1.5;
    end

    if nargin < 5 || isempty(uv_bar)
        uv_bar = round(max(abs(data(:))) / 2);
    end

    fig = figure('Color','w','ToolBar','none','MenuBar','none');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    for ch = 1:size(data,1)
        plot(ax, t, data(ch,:) + (ch-1)*offset, 'LineWidth', 1.5);
    end

    yticks(ax, (0:length(chan_idx)-1)*offset);
    yticklabels(ax, {EEG.chanlocs(chan_idx).labels});

    xlabel(ax, 'Time [s]', 'Fontsize', 16);
    ylabel(ax, 'Channels', 'Fontsize', 16);
    box(ax, 'off');

    % ===== Add amplitude scale bar =====
    xlim_now = xlim(ax);
    ylim_now = ylim(ax);

    x_range = xlim_now(2) - xlim_now(1);
    y_range = ylim_now(2) - ylim_now(1);

    x_bar = xlim_now(2) - 0.03 * x_range;
    y_bar_bottom = ylim_now(2) - 0.25 * y_range;
    y_bar_top = y_bar_bottom + uv_bar;

    line(ax, [x_bar x_bar], [y_bar_bottom y_bar_top], 'Color','k','LineWidth',2);

    cap = 0.01 * x_range;
    line(ax, [x_bar-cap x_bar+cap], [y_bar_bottom y_bar_bottom], 'Color','k','LineWidth',4);
    line(ax, [x_bar-cap x_bar+cap], [y_bar_top y_bar_top], 'Color','k','LineWidth',4);
    set(ax, 'FontSize', 16)

    text(ax, x_bar + 0.01*x_range, ...
        (y_bar_bottom + y_bar_top)/2, ...
        sprintf('%g \\muV', uv_bar), ...
        'VerticalAlignment','middle', ...
        'FontSize',14, ...
        'Interpreter','tex');
    grid(ax, 'on')

    %% saving the plots here
    savefig(fig, [pwd, '/result_plots/', title_path , '_stacked_time_series.fig']);
    saveas(gcf, [pwd, '/result_plots/', title_path , '_stacked_time_series.jpg']); 

end