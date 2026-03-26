%%% cut the continuous audio file for a dbs-learn run into audio files for each trial

op.sub = 'sml001'; 
op.ses = 'subsyl'; 


op.task = 'test1'; 
op.run = 5; %%%% need to use 'run 56' ? 

% op.task = 'test2'; 
% op.run = 7; 


%% set paths
setpaths_dbs_learn()

paths.src_sub = [paths.data, filesep, op.sub]; 
paths.src_ses = [paths.src_sub, filesep, op.ses]; 
paths.audvid = [paths.src_ses, filesep, 'audio-video']; 
paths.beh = [paths.src_ses, filsep, 'beh']; 

paths.der_sub = [paths.data, filesep, 'derivatives', filesep, op.sub]; 
paths.trial_audio = [paths.der_sub, filesep, 'trial-audio']; 
paths.trial_audio_task = [paths.trial_audio, filesep, op.task]; % assume only 1 run per task

%%%% this string gets used in a variety of files associated with this run
% filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_' 'step-',op.step_id '_']; 
filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_']; 


audiofile_full_run = [paths.audvid, filesep, filestr,'recording-headphone.wav']; 

%% load audio and timing info
% make trial audio directories
system(['mkdir ' paths.trial_audio]);
system(['mkdir ' paths.trial_audio_task]);

% load trial data
trials = readtable([paths.beh, filesep, filestr, 'trials.tsv'], 'FileType','text','Delimiter','tab']);

% add vars to trial table
ntrials = height(trials);
cellcol = cell(ntrials,1);
nancol = NaN(ntrials,1); 
trials = [trials, table(cellcol,cellcol,nancol, 'VariableNames',...
                         {'dir','wavname','local_file_id'} )];