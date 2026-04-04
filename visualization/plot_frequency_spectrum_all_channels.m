function plot_frequency_spectrum_all_channels(X_eeg, title_path)
    %% Description:
    %%   Computes and plots PSD of all EEG channels using Welch method.
    %%
    %% Syntax:
    %%   plot_frequency_spectrum_all_channels(X_eeg, title_path)
    %%
    %% Inputs:
    %%   X_eeg      - EEGLAB struct.
    %%   title_path - Output prefix.
    %%
    %% Outputs:
    %%   None (figures saved).
    %%
    %% Notes:
    %%   Uses first ~20000 samples per channel. Jus for the purpose of visualization

    fig = figure; %figure('Color','w', 'ToolBar','none', 'MenuBar','none');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    nchan = size(X_eeg.data,1);

    for ch = 1:nchan
        x = double(X_eeg.data(ch,1:min(end,20000)));
        [pxx, f] = pwelch(x, [], [], [], X_eeg.srate);
        %ydb = 10*log10(pxx);
        semilogy(ax, f, 10*log10(pxx), 'LineWidth', 4);
    end

    restoredefaultpath
    rehash toolboxcache
    legend(X_eeg.chanlocs.labels)
    xlabel(ax, 'Frequency (Hz)', 'FontSize', 16);
    ylabel(ax, 'Power [dB]', 'FontSize', 16);
    title(ax, ['Power Spectrum ', title_path],'FontSize', 16);
    xlim(ax, [1 100]);
    grid(ax, 'on');

    % optional: only show legend if not too many channels
    set(fig, 'Visible', 'on');
    set(ax, 'FontSize',14)
    drawnow;
    %if nchan <= 19
    %    legend(ax, legtxt, 'Interpreter', 'none', 'Location', 'best')
    %end

%% saving the plots here
savefig(fig, [pwd, '/result_plots/', title_path , '_spectrum_results.fig']);
saveas(gcf, [pwd, '/result_plots/', title_path , '_spectrum_results.jpg']); 

%% return back to the addpaths
addpath(genpath([pwd, '/eeglab_current/']));
addpath(genpath([pwd '/MVGC1-master/']));
%% remove the freemat here
allp = strsplit(path, pathsep);
badmask = contains(allp, 'freemat3.5_DISABLED');
badpaths = allp(badmask);
rmpath(genpath(['/home/jmm/Borjigin_Lab/MVGC1-master/utils/legacy/rng/rng.m']))
rmpath((genpath(['/home/jmm/Borjigin_Lab/MVGC1-master/deprecated/utils/newline.m'])));


if ~isempty(badpaths)
    rmpath(strjoin(badpaths, pathsep));
end

end

