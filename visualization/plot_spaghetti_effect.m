function plot_spaghetti_effect(VALUES_GC, VALUES_GC_SD, cond_labels, suffix_save, ylab)
%% ---- collect values from cells ----
%% Plot paired spaghetti/effect plots of GC connection strength across two
%% conditions, separately for each frequency band.
%%
%% This function takes mean GC values and corresponding standard deviations
%% stored in cell arrays, where each cell represents one condition. For each
%% frequency band, it plots the four directional/network connection classes
%% as paired lines across the two conditions, with SD error bars.
%%
%% INPUTS
%% ------
%% VALUES_GC : 1 x 2 cell array
%%     Each cell contains a numeric matrix of size:
%%
%%         [Nfreq x 4]
%%
%%     where rows are frequency bands and columns are connection classes:
%%
%%         1 = FC -> TPO
%%         2 = TPO -> FC
%%         3 = FC -> FC
%%         4 = TPO -> TPO
%%
%% VALUES_GC_SD : 1 x 2 cell array
%%     Same size and structure as VALUES_GC, containing standard deviations
%%     associated with each GC value.
%%
%% cond_labels : 1 x 2 cell array of char/string
%%     Labels for the two conditions shown on the x-axis.
%%
%% suffix_save : char/string
%%     Prefix used when saving output figures.
%%
%% ylab : char/string
%%     Label for the y-axis, e.g. 'Mean GC magnitude'.
%%
%% OUTPUTS
%% -------
%% None.
rmpath(genpath('/home/jmm/Borjigin_Lab/MVGC1-master/deprecated/utils/'));
Y  = zeros(2,size(VALUES_GC{1},1),size(VALUES_GC{1},2));
SD = zeros(2,size(VALUES_GC{1},1),size(VALUES_GC{1},2));
frequencies={'\theta', '\alpha', '\beta', '\gamma_{1}', '\gamma_{2}'};

for i = 1:2
    Y(i,:,:)  = VALUES_GC{i};
    SD(i,:,:) = VALUES_GC_SD{i};
end

% ---- labels ----
conn_labels = {'FC$\rightarrow$TPO', ...
               'TPO$\rightarrow$FC', ...
               'FC$\rightarrow$FC', ...
               'TPO$\rightarrow$TPO'};

% ---- plot ----
for index_s = 1:size(VALUES_GC{1},1)
    %figure; hold on; box on;
    fig2 = figure;
    ax2 = axes('Parent', fig2);
    hold(ax2, 'on');

    x = 1:2;
    
    for c = 1:size(VALUES_GC{1},2)
        errorbar(x, Y(:,index_s,c), SD(:,index_s,c), ...
            '-o', ...
            'LineWidth', 4, ...
            'MarkerSize', 8, ...
            'CapSize', 8);
    end
    
    xticks(x);
    xticklabels(cond_labels);
    xtickangle(45);
    
    ylabel(ylab);
    
    legend(conn_labels, ...
        'Interpreter', 'latex', ...
        'Location', 'best');
    
    set(ax2, ...
        'FontSize', 16, ...
        'TickLabelInterpreter', 'latex');
    title(ax2, sprintf('$%s$', ...
            frequencies{index_s}), ...
            'Interpreter','latex', 'FontSize',18);
    
    grid on;
    savefig(fig2, [pwd, '/comparison_plots/', suffix_save, '_', frequencies{index_s}, '_results.fig']);
    saveas(gcf, [pwd,  '/comparison_plots/', suffix_save, '_', erase(frequencies{index_s},'\'), '_results.jpg']);
end