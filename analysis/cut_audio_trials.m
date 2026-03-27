%%% cut the continuous audio file for a dbs-learn run into audio files for each trial

op.sub = 'sml001'; 
op.ses = 'subsyl'; 


% op.task = 'test1'; 
% op.run = 5; %%%% need to use 'run 56' ? 

op.task = 'test2'; 
op.run = 7; 

% op.channel_to_cut = 'mic';
op.channel_to_cut = 'headphone'; 

op.trial_event_start = 't_go_aud_on'; % field name in trial table
op.trial_event_end = 't_stim_aud_off'; % field name in trial table (of the next trial)
op.last_trial_duration = 6; % sec 


start_offset = 25; % optional - set the bounds for session audio outside of trials' boundaries
end_offset = 5; % optional - set the bounds for session audio outside of trials' boundaries

%% set paths
setpaths_dbs_learn()

paths.src_sub = [paths.data, filesep,'sourcedata', filesep, op.sub]; 
paths.src_ses = [paths.src_sub, filesep,'ses-',op.ses]; 
paths.audvid = [paths.src_ses, filesep, 'audio-video']; 
paths.beh = [paths.src_ses, filesep, 'beh']; 

paths.der_sub = [paths.data, filesep, 'derivatives', filesep, op.sub]; 
paths.trial_audio = [paths.der_sub, filesep, 'trial-audio']; 
paths.trial_audio_task = [paths.trial_audio, filesep, op.task]; % assume only 1 run per task
paths.annot = [paths.der_sub, filesep, 'annot']
paths.landmarks_file = [paths.annot filesep filestr,  '_sync-landmarks.tsv']; 


%%%% this string gets used in a variety of files associated with this run
% filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_' 'step-',op.step_id '_']; 
filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_']; 


audiofile_full_run = [paths.audvid, filesep, filestr,'recording-', op.channel_to_cut, '.wav']; 
audinfo = audioinfo(audiofile_full_run); 
[fpath, fname, fext] = fileparts(audiofile_full_run); 

%% load audio and timing info
% make trial audio directories
system(['mkdir ' paths.trial_audio]);
system(['mkdir ' paths.trial_audio_task]);

% load trial data
trials = readtable([paths.beh, filesep, filestr, 'trials.tsv'], 'FileType','text','Delimiter','tab');
landmarks = readtable([paths.beh, filesep, filestr, 'sync_landmarks.tsv'], 'FileType','text','Delimiter','tab');


% add vars to trial table
ntrials = height(trials);
cellcol = cell(ntrials,1);
nancol = NaN(ntrials,1); 
trials = [trials, table(cellcol,cellcol,nancol, 'VariableNames',...
                         {'dir','wavname','local_file_id'} )];

%%% get sample vs. audio time info
    roi = table;
    roi.name{1} = [fname, fext]; 
    roi.folder{1} = paths.audvid;
    % roi.chantype{1} = 'directionalmic';
    % roi.filetype{1} = 'audio.wav';
    roi.Fs(1) = audinfo.SampleRate; 
    roi.nSamples(1) = audinfo.TotalSamples; 
    % get approximate landmark sample using time and sample rate
    roi.s1(1) = 1; 
    roi.t1(1) = 1; 
    roi.s2(1) = landmarks.landmark_time(2) * audinfo.SampleRate; % instead of using index 2, we should actually look for rowmatch
    roi.t2(1) = landmarks.landmark_time(2);

% sync
audio_time_minus_trialtab_time = landmarks.landmark_time(2) - landmarks.landmark_time(1);  % instead of using index 1 and 2, we should actually look for rowmatch
trials.audfile_start = trials{:,op.trial_event_start} + audio_time_minus_trialtab_time; %%% not same time time coord system as other times in this table
trials.audfile_end(1:end-1) = trials{2:end,op.trial_event_end} + audio_time_minus_trialtab_time;
trials.audfile_end(end) = trials.audfile_end(1:end-1) + op.last_trial_duration; 

    % cut audio trials
%%%% load session audio
% this section can be run without having done syncing in protocol A04
cfg = [];
cfg.roi = roi; 
cfg.roi.starts(1) = min(trials.audfile_start) - start_offset; % set the bounds for session audio outside of trials' boundaries
cfg.roi.ends(1) = max(trials.audfile_end) + end_offset; % set the bounds for session audio outside of trials' boundaries
session_aud = bml_load_continuous(cfg);
% % % check_this_trial = 0; % initialize; if this var == 1, open the trial in praat