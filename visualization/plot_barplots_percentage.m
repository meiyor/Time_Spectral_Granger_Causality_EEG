function plot_barplots_percentage(counts_number_matrix, freqs, path_patient)
%%
%%   INPUTS:
%%       counts_number_matrix : Numeric 3D array of size [N x F x 4], where:
%%           - N: Number of samples/observations/subjects.
%%           - F: Number of frequency bands (must match length(FREQS)).
%%           - 3rd dimension indices connection types:
%%               1: TPO → FC
%%               2: FC → TPO
%%               3: FC → FC
%%               4: TPO → TPO
%%       freqs              : 1xF cell array of character vectors or strings 
%%                            containing frequency band names (e.g., 
%%                            {'delta', 'theta', 'alpha'}). 
%%       path_patient       : Character vector or string used as a filename 
%%                            prefix (typically a patient or subject ID).
%%
%%   OUTPUTS:
%%       None. Figures are automatically saved to:
%%           <current_directory>/plot_topo_no_norm/
%%
%%   NOTES:
%%       • Requires MATLAB R2018b or later (due to sum(...,'all') syntax).
%%       • The y-axis is hard-coded to [0, 70]. Adjust if percentages exceed 70%.
%%       • The target directory './plot_topo_no_norm/' must exist prior to execution.
%%       • Known label mismatch: Column 1 corresponds to TPO→FC but is 
%%         labeled 'FC→TPO' on the x-axis, and vice versa for column 2.
%%
%%   SEE ALSO: bar, savefig, saveas
num_freqs = length(freqs);
close all;
perc = zeros(num_freqs,4); % [TPO->FC, FC->TPO, FC->FC, TPO->TPO]
M = counts_number_matrix;
for freqs_index = 1:num_freqs
    % Sum connections
    tpo_fc = sum(M(:, freqs_index, 1), 'all');
    fc_tpo = sum(M(:, freqs_index, 2), 'all');
    fc_fc  = sum(M(:, freqs_index, 3), 'all');
    tpo_tpo= sum(M(:, freqs_index, 4), 'all');
    
    total = tpo_fc + fc_tpo + fc_fc + tpo_tpo;
    
    perc(freqs_index,:) = [tpo_fc, fc_tpo, fc_fc, tpo_tpo] / total * 100;

    fig = figure;
    b = bar(perc(freqs_index,:));
    b.FaceColor = 'flat';

    % One color per bar
    b.CData = [
        0.2 0.6 1.0  % TPO->FC
        1.0 0.4 0.4  % FC->TPO
        0.4 0.8 0.4  % FC->FC
        0.7 0.3 0.9  % TPO->TPO
    ];

    xticks(1:4);
    xticklabels({'FC→TPO', 'TPO→FC', 'FC→FC','TPO→TPO'});
    ylabel('Occurence [%]', 'FontSize', 18);
    set(gca, 'FontSize', 16);
    ylim([0 70]);
    grid on;

    
    title(sprintf('%s $%s$', ...
            'Connection ROI contribution in ', freqs{freqs_index}), 'Interpreter','latex', 'FontSize', 18);

    savefig(fig, [pwd, '/plot_topo_no_norm/', path_patient, '_', erase(freqs{freqs_index}, '\'), '_', '_p20', '_contribution_results.fig']);
    saveas(gcf, [pwd, '/plot_topo_no_norm/', path_patient, '_', erase(freqs{freqs_index}, '\'), '_', '_p20', '_contribution_results.jpg']);
    grid on;
end

