function EEG = reading_large_file(path_patient_file, varargin)
%% Convert EDF/EDF+ file into a minimal EEGLAB-style structure
%% using MATLAB edfread/edfinfo, and import EDF+ annotations as EEG.event.
%%
%% OUTPUT:
%%   EEG.data   -> [channels x time]
%%   EEG.event  -> imported from EDF+ annotations
%%   EEG.urevent
%%
%% USAGE:
%%   EEG = edfread_to_eeglab_with_annotations('JF_20250225.edf');
%%
%%   EEG = edfread_to_eeglab_with_annotations('JF_20250225.edf', ...
%%           'SelectedChannels', 1:20, ...
%%           'ForceSrate', []);
%%
%% OPTIONAL NAME-VALUE PAIRS:
%%   'SelectedChannels' : vector of channel indices to keep (default = all)
%%   'ForceSrate'       : manually force sampling rate (default = auto)
%%
%% NOTES:
%% - This assumes all kept channels have the same samples-per-record.
%% - EDF+ annotations are taken from edfinfo(...).Annotations.
%% - Requires MATLAB edfread / edfinfo.
%%
%% JMM-style version: practical, minimal, and robust in terms of Matlab and
%% local computing resources limitations
%% -------------------------------------------------------------------------
%% Description:
%%   Reads large EDF files, converts to EEGLAB-like struct, and imports events.
%%
%% Syntax:
%%   EEG = reading_large_file(path_patient_file)
%%
%% Inputs:
%%   path_patient_file - EDF file path.
%%
%% Outputs:
%%   EEG - EEGLAB-style struct.
%%
%% Notes:
%%   Supports channel selection and forced sampling rate. Use this for the large input EEG file**
%% -------------------------------------------------------------------------

% ----------------------------
% Parse inputs
% ----------------------------
p = inputParser;
addRequired(p, 'path_patient_file', @(x) ischar(x) || isstring(x));
addParameter(p, 'SelectedChannels', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'ForceSrate', [], @(x) isempty(x) || isnumeric(x));
parse(p, path_patient_file, varargin{:});

path_patient_file = char(p.Results.path_patient_file);
selected_channels = p.Results.SelectedChannels;
force_srate = p.Results.ForceSrate;

% ----------------------------
% Read EDF data and metadata
% ----------------------------
info = edfinfo(path_patient_file);
TT   = edfread(path_patient_file);   % timetable object
ann  = info.Annotations;             % timetable with Onset / Annotations / Duration use this to infer where to stop the stages annotation

% ----------------------------
% Channel selection
% ----------------------------
all_channel_names = TT.Properties.VariableNames;
nrec  = height(TT);
nchan_total = width(TT);

if isempty(selected_channels)
    selected_channels = 1:nchan_total;
end

channel_names = all_channel_names(selected_channels);
nchan = numel(channel_names);

% ----------------------------
% Infer samples per record and sampling rate
% ----------------------------
% Each timetable cell should contain a vector for one record.
first_cell = TT{1, selected_channels(1)};

% Some MATLAB versions return numeric vectors directly, some cell-wrapped.
if iscell(first_cell)
    first_vec = first_cell{1};
else
    first_vec = first_cell;
end

samples_per_record = numel(first_vec);

% Determine record duration in seconds.
% For EDF readouts like yours, each row is typically 1 second.
% We estimate it from row times if possible; otherwise fall back to 1 sec.
record_duration_sec = [];
try
    if nrec >= 2
        dt = seconds(TT.Properties.RowTimes(2) - TT.Properties.RowTimes(1));
        if ~isempty(dt) && isfinite(dt) && dt > 0
            record_duration_sec = dt;
        end
    end
catch
end

if isempty(record_duration_sec)
    % fallback for standard EDF+C cases like yours
    record_duration_sec = 1;
end

if isempty(force_srate)
    fs = samples_per_record / record_duration_sec;
else
    fs = force_srate;
end

% ----------------------------
% Flatten timetable to channels x time matrix
% ----------------------------
total_samples = nrec * samples_per_record;
X = zeros(nchan, total_samples, 'single');

for c = 1:nchan
    col = TT.(channel_names{c});

    % Convert timetable column into one long vector
    if iscell(col)
        tmp = cell2mat(col.');   % [samples_per_record x nrec] expected
    else
        % Sometimes timetable vars can already be matrix-like
        tmp = cell2mat(num2cell(col).');
    end

    X(c, :) = single(tmp(:)');
end

disp('--- The annotations obtained from reading the large edf file  ---');
disp(ann)

%% read the default structure here
X_default_structure = pop_biosig(path_patient_file, 'blockrange', [0 1000]);

 %% Use the default EEGLAB structure as template to substitute the data from edfread
EEG = X_default_structure;

% Replace signal data with full flattened matrix from edfread
EEG.data   = X;
EEG.nbchan = size(X, 1);
EEG.pnts   = size(X, 2);
EEG.trials = 1;
EEG.srate  = fs;
EEG.xmin   = 0;
EEG.xmax   = (EEG.pnts - 1) / EEG.srate;
EEG.times  = (0:EEG.pnts - 1) / EEG.srate * 1000;   % ms

% Reset fields that may no longer match the original small excerpt
EEG.icaact      = [];
EEG.icawinv     = [];
EEG.icasphere   = [];
EEG.icaweights  = [];
EEG.icachansind = [];
EEG.epoch       = [];
EEG.specdata    = [];
EEG.specicaact  = [];
EEG.stats       = [];
EEG.reject      = [];

% Update set/file naming
[filepath_, filename_, ext_] = fileparts(path_patient_file);
EEG.setname   = [filename_ '_edfread_full'];
EEG.filename  = [filename_ ext_];
EEG.filepath  = filepath_;

% Update channel labels if needed
if length(EEG.chanlocs) >= nchan
    for c = 1:nchan
        EEG.chanlocs(c).labels = channel_names{c};
    end
    if length(EEG.chanlocs) > nchan
        EEG.chanlocs = EEG.chanlocs(1:nchan);
    end
else
    EEG.chanlocs = struct('labels', '');
    for c = 1:nchan
        EEG.chanlocs(c).labels = channel_names{c};
    end
end

%% Rebuild EEG.event and EEG.urevent from EDF annotations
EEG.event   = struct('type', {}, 'latency', {}, 'duration', {}, 'urevent', {});
EEG.urevent = struct('type', {}, 'latency', {}, 'duration', {});

%% define herethe latencies of each stages S1-S11 in raw stage
if ~isempty(ann) && height(ann) > 0
    for k = 1:height(ann)
        onset_sec = seconds(ann.Onset(k));

        % annotation label
        evt_type = char(string(ann.Annotations(k)));

        % duration
        dur_sec = 0;
        try
            dur_sec = seconds(ann.Duration(k));
            if isnan(dur_sec)
                dur_sec = 0;
            end
        catch
            dur_sec = 0;
        end

        evt_latency  = round(onset_sec * fs) + 1;
        evt_duration = round(dur_sec * fs);

        % clamp inside bounds
        evt_latency  = max(1, min(evt_latency, EEG.pnts));
        evt_duration = max(0, evt_duration);

        EEG.event(k).type     = evt_type;
        EEG.event(k).latency  = evt_latency;
        EEG.event(k).duration = evt_duration;
        EEG.event(k).urevent  = k;

        EEG.urevent(k).type     = evt_type;
        EEG.urevent(k).latency  = evt_latency;
        EEG.urevent(k).duration = evt_duration;
    end
end

%% Optional consistency check
if exist('eeg_checkset', 'file') == 2
    EEG = eeg_checkset(EEG, 'eventconsistency');
    EEG = eeg_checkset(EEG);
end