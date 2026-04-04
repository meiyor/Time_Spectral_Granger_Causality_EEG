function time_frequency_gc_single_edf(path_patient, suffix)
%% get here the loaded environments to run Granger Causality following the MVGC Barnett's toolbox for multi-channel and multi-trial time-series.
%% This snippet of code will be tested in a tested .edf files from the Borjigin Data folder downloaded from Zenodo and representing a complement
%% for the biomarkers calculated in [Xu et al 2023]
%% Description:
%%   Full preprocessing + GC pipeline for a single EDF file.
%%
%% Syntax:
%%   time_frequency_gc_single_edf(path_patient, suffix)
%%
%% Inputs:
%%   path_patient - EDF file prefix.
%%   suffix       - Processing variant (e.g., reference type).
%%
%% Outputs:
%%   None (results saved).
%%
%% Notes:
%%   Uses EEGLAB + MVGC. Handles referencing, filtering, ICA, as preprocessing stages, and time/frequency GC
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

process_status_separated = 0;
%% add the datapath here for reading the edfs as it comes - with the current file name
if strcmp(path_patient, 'Pt2_S2-end')
    X_eeg_raw = pop_biosig([path_patient '.edf'], 'blockrange', [0 1000]);
elseif strcmp(path_patient, 'JF_20250225') %% processing here the large files using edfread here**
    %% the whole process using the whole file can take a lot of time**
    X_eeg_raw = reading_large_file([path_patient '.edf']);
    process_status_separated = 1;
else
    X_eeg_raw = pop_biosig([path_patient '.edf']);
end

%% get the ECG signal here..
ECG_signal = X_eeg_raw.data(20,:);
old_srate = X_eeg_raw.srate;

%% remove ECG before ICA application
idx_keep = find(~strcmpi({X_eeg_raw.chanlocs.labels}, 'EKG1'));
X_eeg_raw = pop_select(X_eeg_raw, 'channel', idx_keep);

%% re-reference to the average of the whole signals
if contains(suffix,'Fz')
   X_eeg_ref = pop_reref(X_eeg_raw, 'Fz');
else
   X_eeg_ref = pop_reref(X_eeg_raw, []);
end

%% downsample to 256
X_eeg_ref = pop_resample(X_eeg_ref, 256);
%% apply the notch filter on 60Hz first..
X_eeg_filter_notch = pop_eegfiltnew(X_eeg_ref, 59, 61, [], 1);
%% filter the signal between 0.1 and 100Hz
X_eeg_filtered = pop_eegfiltnew(X_eeg_filter_notch, 0.1, 100);
%% use cleanraw data pluging here to eliminate EMGs, use ASR in a very light way...
if process_status_separated == 0
   if strcmp(path_patient, 'Pt1_S2-end')
       X_eeg_clean=clean_rawdata(X_eeg_filtered,5,[0.25 0.75],0.85,-1,200,-1);
   else
       X_eeg_clean=clean_rawdata(X_eeg_filtered,5,[0.25 0.75],0.85,-1,50,-1);
   end
else
   %% do not remove any channel in the last preprocessing
   X_eeg_clean=clean_rawdata(X_eeg_filtered,5,[0.25 0.75],-1,-1,5,-1);     
end

%% apply ICA decomp here
rng(42)  % set seed for reproducibility
% DON'T SET LEARNING RATE OR DATA RANK FOR GUARANTEEING MINIMUM
% REPLICABILITY AND COLINEARITY..THIS VERY IMPORTANT FOR ASSURING MATRIX
% NON-SINGULARITY IN THE SUBSEQUENT MVGC GC ESTIMATION PROCESSES
[X_eeg_clean.icaweights,X_eeg_clean.icasphere]=runica(X_eeg_clean.data(:,:),'sphering','on', 'extended', 1, 'maxsteps', 1000);

%% update the inverse matrix here
X_eeg_clean.icawinv=pinv(X_eeg_clean.icaweights*X_eeg_clean.icasphere);

%% set up the channel locations here before applying ADJUST
X_eeg_clean.chanlocs = reshape(X_eeg_clean.chanlocs, 1, []);

%% readlocs is not working so use the customized one here
X_eeg_clean = assign_xyz_from_elc(X_eeg_clean, [pwd '/eeglab_current/eeglab2026.0.0/plugins/dipfit/standard_BEM/elec/standard_1020.elc']);
%% add the radius, theta and psi for the spherical representation
X_eeg_clean = pop_chanedit(X_eeg_clean, 'convert', {'cart2all'});
X_eeg_clean = eeg_checkset(X_eeg_clean);

%X_eeg_clean = pop_chanedit(X_eeg_clean, 'lookup', ...
%    [pwd '/eeglab_current/eeglab2026.0.0/plugins/dipfit/standard_BEM/elec/standard_1020.elc']);

rmpath(genpath('/home/jmm/Borjigin_Lab/eeglab_current/eeglab2026.0.0/plugins/Biosig3.8.5/NaN/inst/'))
%% apply Artifact recognition using heuristic classification
% [bad_ic, feat] = auto_reject_ica_no_topo(X_eeg_clean);
[bad_ic, feat] = auto_reject_ica_no_topo_minimal(X_eeg_clean, 3);


disp(feat)
disp('Bad ICs:')
disp(bad_ic.')

%% remove the bad ics here..
X_eeg_clean = pop_subcomp(X_eeg_clean, bad_ic, 0);

%% remove the ECG signal using a linear regression
ECG_signal_ds = resample(double(ECG_signal), X_eeg_clean.srate, old_srate);
[EEG_CLEAN_ECG, beta_all] = remove_ecg_regression_lagged(X_eeg_clean.data, ECG_signal_ds, X_eeg_clean.srate, 0.12, 1e-2);
X_eeg_clean.data = EEG_CLEAN_ECG;


%% Apply ADJUST here only if considered necessary - only if you want ADJUSt to be applicable
%%[art_channels_test]=ADJUST(X_eeg_clean, 'report_test_ADJUST.txt');
if process_status_separated == 0
   %% plot the signal here just to observe it
   stacked_eeg_plot(X_eeg_clean, 1:size(X_eeg_clean.data,1), 1:min(5120, size(X_eeg_clean.data,2)), 20, 20, path_patient);

   %% plot the frequency plots for each channels just checking
   plot_frequency_spectrum_all_channels(X_eeg_clean, [path_patient, suffix]);

   %% do the preliminary GC analysis 
   MVGC_application(X_eeg_clean, path_patient, 'custom');
else
   %% split the events in the new eeglab structure based on the event latencies
   split_EEGlab_large_events(X_eeg_clean, path_patient, suffix);
   %% save the preprocessed and splitted EEGlab structures
   disp(["Saving ", path_patient, " structures..."]);
   save([pwd '/preprocessed_save/JF_20250225_clean_and_splited', suffix, '.mat'], ...
            'X_eeg_clean', '-v7.3');
   %% continue here reading the structiures if necessary **UNCOMMENT THIS ONLY IF NECESSARY**
   %%for stages = 1:11
   %%   X_sub_eeg = load([path_patient, '_clean_and_splitted_S', num2str(stages), '.mat']);
   %%   MVGC_application(X_sub_eeg.EEG_k, [path_patient, '_S', num2str(stages)], 'custom');
   %%end
end