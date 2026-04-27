function plot_comparison_spaghetti_plots(path_patient, list_stages, frequencies, clim_vals, suffixes)
%% this plots the spaghetti plot comparison across S2-S11 stages
%% comparing amplitudes and also amount of connections adding spaghetti effective plots
%% across different regions in the same plot, add error/SD and put the point in the mean
%% this function visualizes true changes for (1) GC values and (2) GC occurrences 
%% across FC->TPO, TPO->FC, FC->FC, TPO->TPO directions. 

%% Inputs:
%%   path_patient - File prefix.
%%   list_stages  - Cell array of stages.
%%   frequencies  - Cell array of frequency bands.
%%   clim_vals    - Cell of Color limits [min max] for each pair
%%   suffixes - suffixes cell of suffixe to compare
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
rmpath(genpath('/home/jmm/Borjigin_Lab/eeglab_current/eeglab2026.0.0/plugins/Fieldtrip-lite250523'));

if ~isempty(badpaths)
    rmpath(strjoin(badpaths, pathsep));
end

rehash toolboxcache
clear gcd

%% add data paths
addpath(genpath([pwd '/Data/']));
addpath(genpath([pwd '/preprocessed_save/']));
%% add the auxiliary function paths here
addpath(genpath([pwd '/utils/']));
addpath(genpath([pwd '/visualization/']));
addpath(genpath([pwd '/core/']));
addpath(genpath([pwd '/preprocessing/']));

%% base file for chanlocs % run this always if you want to read the pipeline from the beginning
%X_base = load(['JF_20250225_clean_and_splitted_S1', suffix ,'.mat']);
%% do this just in order to read the base file without suffix in the beginning
X_base = load(['JF_20250225_clean_and_splited_average_reref_ASR_ICA_ECG_rem.mat']);
all_labels = {X_base.X_eeg_clean.chanlocs.labels};

%% channel subset
%%tpo_labels = {'T3','T4','T5','T6','P3','P4','O1','O2','Pz','C3','C4','Cz','F3','F4','Fz'};
tpo_labels = {'T3','T4','T5','T6','P3','P4','O1','O2','C3','C4','F3','F4','F7','F8'};
idx_tpo = find(ismember(all_labels, tpo_labels));
chnames = all_labels(idx_tpo);
chanlocs_tpo = X_base.X_eeg_clean.chanlocs(idx_tpo);
nvars = numel(idx_tpo);


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

number_connections = zeros(length(list_stages), length(frequencies), 4);
value_connections = zeros(length(list_stages), length(frequencies), 4);
matrix_connections = cell(1,length(suffixes));
value_connections_m = cell(1,length(suffixes));
VALUES_GC = cell(1,length(suffixes));
VALUES_GC_SD = cell(1,length(suffixes));
OCCURS_GC = cell(1,length(suffixes));
OCCURS_GC_SD = cell(1,length(suffixes));

for index_suffix = 1:length(suffixes)
  suffix=suffixes{index_suffix};
  for index_freq = 1:length(frequencies)
    close all;
    for index_status = 1:length(list_stages)
        %% threshold for each band
        thr = clim_vals{index_freq}{index_suffix}(1) + 0.5*diff(clim_vals{index_freq}{index_suffix});
        %% initialize the count per frequency
        fc_to_tpo = 0;
        tpo_to_fc = 0;
        tpo_to_tpo = 0;
        fc_to_fc = 0;

        %% get the values of those connections
        val_fc_to_tpo = 0;
        val_tpo_to_fc = 0;
        val_tpo_to_tpo = 0;
        val_fc_to_fc = 0;
        %% load GC result
        X_eval = load([pwd, '/GC_predictions/', path_patient, '_', ...
            list_stages{index_status}, suffix, suffix, '_', frequencies{index_freq}, ...
            '_p20_gc_band_results.mat']);

        F_band = X_eval.F_band;
        F_band(1:nvars+1:end) = 0;

        idx_CENTRAL = find(ismember(chnames, {'F3', 'F4', 'F7', 'F8', 'C3','C4','Cz'}));
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
                    %% evaluates the amount of times a connection is made
                    %% between FC and TPO or visceversa
                    %% get here the direction of the GC score across ROIs
                    if any(idx_CENTRAL == jsource) && any(idx_TPO == itarget)
                        %% measure here the number of connections going from FC to TPO
                        fc_to_tpo = fc_to_tpo + 1;
                        val_fc_to_tpo = val_fc_to_tpo + val;
                    elseif any(idx_TPO == jsource) && any(idx_CENTRAL == itarget)
                        %% measure the number of connections TPO to FC
                        tpo_to_fc = tpo_to_fc + 1;
                        val_tpo_to_fc = val_tpo_to_fc + val;
                    elseif any(idx_TPO == jsource) && any(idx_TPO == itarget)
                        %% measure the number of connections TPO to TPO
                        tpo_to_tpo = tpo_to_tpo + 1;
                        val_tpo_to_tpo = val_tpo_to_tpo + val;
                    elseif any(idx_CENTRAL == jsource) && any(idx_CENTRAL == itarget) 
                        %% measure the number of connections FC to FC
                        fc_to_fc = fc_to_fc + 1;
                        val_fc_to_fc = val_fc_to_fc + val;
                    end
                end
            end
        end

        %% fill up the matrix here
        number_connections(index_status, index_freq,:) = [fc_to_tpo, tpo_to_fc, fc_to_fc, tpo_to_tpo];
        value_connections(index_status, index_freq,:) = [val_fc_to_tpo/fc_to_tpo, val_tpo_to_fc/tpo_to_fc, val_fc_to_fc/fc_to_fc, val_tpo_to_tpo/tpo_to_tpo];
    end
  end
  matrix_connections{index_suffix} = number_connections;
  value_connections_m{index_suffix} = value_connections;
  %% get the magnitude of the TPO->TPO, FC->TPO, FC->FC, and TPO->FC connections for each frequency range 
  FC_TPO = mean(value_connections_m{index_suffix}(2:end, :, 1), 1);
  TPO_FC = mean(value_connections_m{index_suffix}(2:end, :, 2), 1);
  FC_FC  = mean(value_connections_m{index_suffix}(2:end, :, 3), 1);
  TPO_TPO = mean(value_connections_m{index_suffix}(2:end, :, 4), 1);

  SD_FC_TPO = std(value_connections_m{index_suffix}(2:end, :, 1), 1);
  SD_TPO_FC = std(value_connections_m{index_suffix}(2:end, :, 2), 1);
  SD_FC_FC  = std(value_connections_m{index_suffix}(2:end, :, 3), 1);
  SD_TPO_TPO = std(value_connections_m{index_suffix}(2:end, :, 4), 1);

  VALUES_GC{index_suffix} = [FC_TPO', TPO_FC', FC_FC', TPO_TPO'];
  VALUES_GC_SD{index_suffix} = [SD_FC_TPO', SD_TPO_FC', SD_FC_FC', SD_TPO_TPO'];

  %% get the occurence percentage of the TPO->TPO, FC->TPO, FC->FC, and TPO->FC connections for each frequency range 
  FC_TPO_o = sum(matrix_connections{index_suffix}(2:end, :, 1), 1);
  TPO_FC_o = sum(matrix_connections{index_suffix}(2:end, :, 2), 1);
  FC_FC_o  = sum(matrix_connections{index_suffix}(2:end, :, 3), 1);
  TPO_TPO_o = sum(matrix_connections{index_suffix}(2:end, :, 4), 1);

  TOTAL = FC_TPO_o + TPO_FC_o + FC_FC_o + TPO_TPO_o;

  FC_TPO_o = FC_TPO_o ./ TOTAL .* 100;
  TPO_FC_o = TPO_FC_o ./ TOTAL .* 100;
  FC_FC_o  = FC_FC_o ./ TOTAL .* 100;
  TPO_TPO_o = TPO_TPO_o ./ TOTAL .* 100;

  SD_TPO_FC_o = std(matrix_connections{index_suffix}(2:end, :, 1), 1);
  SD_FC_TPO_o = std(matrix_connections{index_suffix}(2:end, :, 2), 1);
  SD_FC_FC_o  = std(matrix_connections{index_suffix}(2:end, :, 3), 1);
  SD_TPO_TPO_o = std(matrix_connections{index_suffix}(2:end, :, 4), 1);

  TOTAL_SD = SD_FC_TPO_o + SD_TPO_FC_o + SD_FC_FC_o + SD_TPO_TPO_o;

  SD_FC_TPO_o = SD_FC_TPO_o ./ TOTAL_SD;
  SD_TPO_FC_o = SD_TPO_FC_o ./ TOTAL_SD;
  SD_FC_FC_o  = SD_FC_FC_o ./ TOTAL_SD;
  SD_TPO_TPO_o = SD_TPO_TPO_o ./ TOTAL_SD;


  VALUES_GC{index_suffix} = [FC_TPO', TPO_FC', FC_FC', TPO_TPO'];
  VALUES_GC_SD{index_suffix} = [SD_FC_TPO', SD_TPO_FC', SD_FC_FC', SD_TPO_TPO'];
  OCCURS_GC{index_suffix} = [FC_TPO_o', TPO_FC_o', FC_FC_o', TPO_TPO_o'];
  OCCURS_GC_SD{index_suffix} = [SD_FC_TPO_o', SD_TPO_FC_o', SD_FC_FC_o', SD_TPO_TPO_o'];
end


idx_1 = find(suffixes{1} == '_');  
third_pos_1 = idx_1(3);
idx_2 = find(suffixes{2} == '_');  
third_pos_2 = idx_2(3);

plot_spaghetti_effect(VALUES_GC, VALUES_GC_SD, {suffixes{1}(third_pos_1+1:end), suffixes{2}(third_pos_2+1:end)}, [suffixes{1}(2:end), '_vs_', suffixes{2}(2:end)], 'GC [a.u]');
plot_spaghetti_effect(OCCURS_GC, OCCURS_GC_SD, {suffixes{1}(third_pos_1+1:end), suffixes{2}(third_pos_2+1:end)}, [suffixes{1}(2:end), '_vs_', suffixes{2}(2:end), '_occur_percentage'], '%');