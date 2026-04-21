function EEG_k = split_EEGlab_large_events_update_corrected(X_eeg_clean, X_eeg_corr, path_patient, suffix)
%% Description:
%%   Splits EEG into segments based on events and saves each segment.
%%
%% Syntax:
%%   EEG_k = split_EEGlab_large_events(X_eeg_clean, path_patient, suffix)
%%
%% Inputs:
%%   X_eeg_clean - EEGLAB struct with events.
%%   X_eeg_corr - cells with corrected structures
%%   path_patient- Output prefix.
%%   suffix      - Filename suffix.
%%
%% Outputs:
%%   EEG_k - Last segment.
%%
%% Notes:
%%   Saves segments to disk, and use it for subsequent EEGlab structure definition
%% -------------------------------------------------------------------------
for k = 1:length(X_eeg_clean.event)

    evt = X_eeg_clean.event(k);

    t1 = (evt.latency - 1) / X_eeg_clean.srate;
    t2 = (evt.latency + evt.duration - 1) / X_eeg_clean.srate;
    disp(["Segmenting event ", string(evt.type), 'with duration ', num2str(t2-t1)])

    EEG_k = pop_select(X_eeg_clean, 'point', [evt.latency evt.latency + evt.duration - 1]);

    EEG_k.setname = sprintf('%s_%s', X_eeg_clean.setname, char(string(evt.type)));
    %% update the corrected data here
    EEG_k.data = X_eeg_corr{k}.data;

    %% transform string to char type
    evt_name = char(string(evt.type));
    evt_name = regexprep(evt_name, '[^a-zA-Z0-9]', '_'); % clean filename using regex rules

    %% create the filename with the adequate data types.
    filename = sprintf('%s_clean_and_splitted_%s%s.mat', path_patient, evt_name, suffix);
    %% get the preprocessed folder for saving interim preprocessing
    fullpath = fullfile(pwd, 'preprocessed_save', filename);

    save(fullpath, 'EEG_k', '-v7.3');
    %%pop_saveset(EEG_k, 'filename', [EEG_k.setname '.set']);
end