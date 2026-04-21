function reading_eeg_saved_MVGC(patient_path, suffix_path, suffix)
%% READING THE PREPROCESSED SAVED SUBSTRUCTURES FOR EACH STAGE S1-S11
%% Description:
%%   Loads segmented EEG stages, plots them, and runs MVGC analysis.
%%
%% Syntax:
%%   reading_eeg_saved_MVGC(patient_path)
%%
%% Inputs:
%%   patient_path - File prefix.
%%   suffix_path - string suffix for reading the filename 
%%   suffix - string suffix for reading the filename 
%%   similar to the previous preprocessing step
%%
%% Outputs:
%%   None (plots + MVGC results generated).
%%
%% Notes:
%%   Requires pre-saved stage files.
%% load the variables again here
%% set the random seed here
close all;
rng(42)  % set seed for reproducibility

%% restoring the path first
restoredefaultpath
rehash toolboxcache
clear gcd
which gcd -all

%% loading the paths here
addpath(genpath([pwd, '/eeglab_current/']));
addpath(genpath([pwd '/MVGC1-master/']));
%% remove the freemat here
allp = strsplit(path, pathsep);
badmask = contains(allp, 'freemat3.5_DISABLED');
badpaths = allp(badmask);
rmpath(genpath(['/home/jmm/Borjigin_Lab/MVGC1-master/utils/legacy/rng/rng.m']));
rmpath((genpath(['/home/jmm/Borjigin_Lab/MVGC1-master/deprecated/utils/newline.m'])));
rmpath(genpath('/home/jmm/Borjigin_Lab/eeglab_current/eeglab2026.0.0/plugins/Fieldtrip-lite250523'));

if ~isempty(badpaths)
    rmpath(strjoin(badpaths, pathsep));
end

%% restore the gcd here NOTE: avoid plotting the signals using eegplot for conflict with datevec and datenum functions
%% in Matlab original functions
rehash toolboxcache
clear gcd

%% add the datapath here
addpath(genpath([pwd '/Data/']));
addpath(genpath([pwd '/preprocessed_save/']));
%% add the auxiliary function paths here
addpath(genpath([pwd '/utils/']));
addpath(genpath([pwd '/visualization/']));
addpath(genpath([pwd '/core/']));
addpath(genpath([pwd '/preprocessing/']));

%% plot the stacked 
EEG_list={};
stage_names={};
for stages = 1:11
    % X_sub_eeg = load([patient_path, '_clean_and_splitted_S', num2str(stages), suffix, '.mat']);
    X_sub_eeg = load([patient_path, '_clean_and_splitted_S', num2str(stages), suffix_path, '.mat']);
    EEG_list{end+1} = X_sub_eeg.EEG_k;
    stage_names{end+1} = ['S', num2str(stages)];

end
stacked_eeg_plot_struct_list(EEG_list, 1:size(X_sub_eeg.EEG_k.data,1), stage_names, 50, 50, [patient_path, suffix, '_stacked']);

%% continue here reading the structiures if necessary
freq_vector = {'theta', 'alpha', 'beta', 'gamma1', 'gamma2'};
for stages = 1:11 %% generate the GC estimation and outputs for all the stages S1-S11 
    % X_sub_eeg = load([patient_path, '_clean_and_splitted_S', num2str(stages), suffix,'.mat']);
    X_sub_eeg = load([patient_path, '_clean_and_splitted_S', num2str(stages), suffix_path, '.mat']);
    %% don't filter the signal if this is really necessary**
    %X_sub_eeg.EEG_k = pop_eegfiltnew(X_sub_eeg.EEG_k, 25, 55);
    stacked_eeg_plot(X_sub_eeg.EEG_k, 1:size(X_sub_eeg.EEG_k.data,1), 1:min(5120, size(X_sub_eeg.EEG_k.data,2)), 20, 20, [patient_path, '_filter_', num2str(stages), suffix]);
    %% do it first for a value of p=20 %% generate the GC interim predictions depends on the p value you select
    MVGC_application(X_sub_eeg.EEG_k, [patient_path, '_S', num2str(stages), suffix, suffix_path], 'custom', 20, freq_vector);
end
