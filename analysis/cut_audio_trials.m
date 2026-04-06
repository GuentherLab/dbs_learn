%% NB: this script uses an 'adjusted sample rate' to account for different speeds of the audfile vs the 
%   clock time in trials table, using the differences between the durations between t1 and t2 in the two files
%%%%% this is very suboptimal - we should a 'synced' version of the trial table

%%% cut the continuous audio file for a dbs-learn run into audio files for each trial

function cut_audio_trials(op)

%% params
vardefault('op',struct);
field_default('op','sub','sml002');
field_default('op','ses','subsyl');

field_default('op','channels_to_cut',{'mic','headphone'});

field_default('op','trial_event_start','t_go_aud_on'); % field name in trial table
field_default('op','trial_event_end','t_stim_aud_off'); % field name in trial table (of the next trial)
field_default('op','last_trial_duration',6); % seconds

field_default('op','start_offset',25); % optional - set the bounds for session audio outside of trials' boundaries
field_default('op','end_offset',5); % optional - set the bounds for session audio outside of trials' boundaries

% this channel should always be listed in the landmarks sync table
landmark_audio_chan = 'headphone'; 

% manual editing of specific times
% use this table to modify the start/stop times of specific trials
%   ... for use on specific trials when sub starts early or trigger was sent at wrong time
%%%% indicate which trials to modify with 'inds'
%%%% for each ind, use 'sec' to specify how much to move start times [col 1] and stop times [col 2]
%%%%%%%% negative values move timepoints earlier, positive move later
field_default('op','trials_to_modify_inds',[]); 
    field_default('op','trials_to_modify_sec',[]); 

%% 
% set paths, load sync table 
paths = setpaths_dbs_learn(op); 
sync = readtable(paths.landmarks_file, 'FileType','text','Delimiter','tab');

switch op.ses
    case 'subsyl'
        tasks = {'famil';'pretest';'trainA';'trainB';'test1';'test2'}; 
        trialvars_to_copy = {'stim_group','is_native','name'}; % add these to audiofiles table
    case 'multisyl'
        tasks = {'fds';'famil';'assess';'pretest';'trainA';'trainB';'test1';'test2'}; 
        trialvars_to_copy = {'stim_group','name','n_syllables'}; % add these to audiofiles table
    otherwise 
        error('session not recognized')
end

n_syncrows = height(sync);

%% cut trial audio for each run with a trial-based task, for each audio channel to cut
for i_syncrow = 1:n_syncrows
    op.task = sync.task{i_syncrow}; 
    op.run = sync.run(i_syncrow);
    op.step = sync.step{i_syncrow}; 
    if ~any(string(op.task) == tasks) 
        fprintf(['   skipping task ''', thistask, ''' because it is not a trial-based task \n'])
        continue
    end
    paths = setpaths_dbs_learn(op); % add run-specific paths

    % syncing
    %%%% because headphone and mic have exactly the same timepoints (both collected simultaneously on focusrite)....
    %%%% .... we just load syncing for headphone channel and apply it to both
    landmark_audiofile = [paths.src_audvid, filesep, paths.filestr,'recording-', landmark_audio_chan, '.wav']; 
    [landmark_aud_path, landmark_aud_name, landmark_aud_ext] = fileparts(landmark_audiofile); 
    lnd_audiorow = sync.run == op.run & strcmp(sync.task,op.task) & cellfun(@(x)contains(x,getfname(landmark_audiofile)), sync.filename);
    lnd_trialsrow = sync.run == op.run & strcmp(sync.task,op.task) & cellfun(@(x)contains(x,'trials.tsv'), sync.filename);

    % skip runs for which landmarks haven't yet been marked
    if isnan(sync.t1(lnd_audiorow)) || isnan(sync.t2(lnd_audiorow))
        fprintf(['missing landmark time(s) for ', sync.filename{lnd_audiorow}, '.... skipping run \n'])
        continue 
    end

    audio_time_minus_trialtab_time = sync.t1(lnd_audiorow) - sync.t1(lnd_trialsrow);  
    
        % load trial data
    trials = readtable([paths.beh, filesep, paths.filestr_step, 'trials.tsv'], 'FileType','text','Delimiter','tab');
    ntrials = height(trials);
    cellcol = cell(ntrials,1);
    nancol = NaN(ntrials,1);

    % make table listing trialwise audio files
    audiofiles = [trials(:,trialvars_to_copy),...
                   table(cellcol,cellcol,nancol,nancol,nancol, 'VariableNames',...
                         {'filename','dir','starts','ends','duration'} )];
    audiofiles.starts = trials{:,op.trial_event_start} + audio_time_minus_trialtab_time; %%% starts is not in same time coord system as other times in trials table
    audiofiles.ends(1:end-1) = trials{2:end,op.trial_event_end} + audio_time_minus_trialtab_time;
    audiofiles.ends(end) = audiofiles.starts(end) + op.last_trial_duration; 
    audiofiles.duration = audiofiles.ends - audiofiles.starts; 
 
    for ichan = 1:length(op.channels_to_cut)
        
        
        
        % make trial audio directories
        cutchan = op.channels_to_cut{ichan}; 
        paths.trial_audio_task_chan = [paths.trial_audio, filesep, op.task,'_',cutchan]; % trial audio clips for this run and this recording chan
        if ~exist(paths.trial_audio, 'dir'); system(['mkdir ' paths.trial_audio]); end
        if ~exist(paths.trial_audio_task_chan, 'dir'); system(['mkdir ' paths.trial_audio_task_chan]); end
        
        % make channel-specific audiofile trial table
        audiofiles_this_channel = audiofiles; 
        audiofiles_this_channel.dir = repmat([paths.trial_audio_task_chan],ntrials,1); 
        
        % load full run audio
        audiofile_full_run = [paths.src_audvid, filesep, paths.filestr_step,'recording-', cutchan, '.wav']; 
        [run_aud, fs] = audioread(audiofile_full_run); 
        
        % adjust for differences in file/clock speeds
        fs_adj = fs * [(sync.t2(lnd_audiorow) - sync.t1(lnd_audiorow)) / (sync.t2(lnd_trialsrow) - sync.t1(lnd_trialsrow))];
        
        for itrial = 1:ntrials
            audiofiles_this_channel.filename{itrial} = ['trial', num2str(itrial), '_', cutchan '_', trials.name{itrial}, '.wav'];
        
            % apply manual edits to trialtimes
            if ~isempty(op.trials_to_modify_inds) && any(itrial == op.trials_to_modify_inds)
                matchrow = itrial == op.trials_to_modify_inds;
                start_time_edit = trials_to_modify_sec(matchrow, 1);
                end_time_edit = trials_to_modify_sec(matchrow, 2);
                audiofiles_this_channel.starts(itrial) = audiofiles_this_channel.starts(itrial) + start_time_edit; 
                audiofiles_this_channel.ends(itrial) = audiofiles_this_channel.ends(itrial) + end_time_edit; 
            end
        
            % Convert times to sample indices
            startSample = max(1, round(audiofiles_this_channel.starts(itrial) * fs_adj));
            endSample   = min(size(run_aud, 1), round(audiofiles_this_channel.ends(itrial) * fs_adj));
    
            trialaud = run_aud(startSample:endSample, :);
        
            audiowrite([paths.trial_audio_task_chan filesep audiofiles_this_channel.filename{itrial}], trialaud, fs)
            % % % if check_this_trial && open_problematic_trials_in_praat % load created audio file in praat
            % % %     cmd = ['praat --open ',  [audiofiles_this_channel.dir{itrial} filesep audiofiles_this_channel.filename{itrial}], ' &'];
            % % %     system(cmd);
            % % %     check_this_trial = 0; 
            % % % end
    
        end
    
        audiofiles_table_filename = [paths.annot, filesep, paths.filestr, 'audiofiles-',cutchan, '.tsv']; 
        audiofiles_this_channel = renamevars(audiofiles_this_channel,{'starts','ends','duration'},...
                            {'audfile_start','audfile_end','audfile_dur'}); 
        writetable(audiofiles_this_channel, audiofiles_table_filename, 'FileType','text','Delimiter','tab');
    end
end

