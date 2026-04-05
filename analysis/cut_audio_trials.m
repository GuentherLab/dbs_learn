%% NB: this script uses an 'adjusted sample rate' to account for different speeds of the audfile vs the 
%   clock time in trials table, using the differences between the durations between t1 and t2 in the two files
%%%%% this is very suboptimal - we should a 'synced' version of the trial table

%%% cut the continuous audio file for a dbs-learn run into audio files for each trial

clear

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

%% manual editing of specific times
% use this table to modify the start/stop times of specific trials
%   ... for use on specific trials when sub starts early or trigger was sent at wrong time
%%%% indicate which trials to modify with 'inds'
%%%% for each ind, use 'sec' to specify how much to move start times [col 1] and stop times [col 2]
%%%%%%%% negative values move timepoints earlier, positive move later
trials_to_modify_inds = []; 
    trials_to_modify_sec = []; 

%% set paths
paths = setpaths_dbs_learn(); 

audiofile_full_run = [paths.audvid, filesep, paths.filestr,'recording-', op.channel_to_cut, '.wav']; 
audinfo = audioinfo(audiofile_full_run); 
[audfull_path, audfull_name, audfull_ext] = fileparts(audiofile_full_run); 

%% load audio and timing info
% make trial audio directories
system(['mkdir ' paths.trial_audio]);
system(['mkdir ' paths.trial_audio_task]);

% load trial data
trials = readtable([paths.beh, filesep, paths.filestr, 'trials.tsv'], 'FileType','text','Delimiter','tab');
landmarks = readtable(paths.landmarks_file, 'FileType','text','Delimiter','tab');


% make audiofile trial table
ntrials = height(trials);
cellcol = cell(ntrials,1);
nancol = NaN(ntrials,1); 
audiofiles = [trials(:,{'stim_group','is_native','name'}),...
                   table(cellcol,cellcol,nancol,nancol,nancol, 'VariableNames',...
                         {'dir','wavname','starts','ends','duration'} )];

% % % % %%% get sample vs. audio time info
% % % % roi = table;
% % % % roi.name{1} = [audfull_name, fext]; 
% % % % roi.folder{1} = paths.audvid;
% % % % % roi.chantype{1} = 'directionalmic';
% % % % roi.filetype{1} = 'audio.wav';
% % % % roi.Fs(1) = audinfo.SampleRate; 
% % % % roi.nSamples(1) = audinfo.TotalSamples; 
% % % % % get approximate landmark sample using time and sample rate
% % % % roi.s1(1) = 1; 
% % % % roi.t1(1) = 1; 
% % % % roi.s2(1) = landmarks.t1(2) * audinfo.SampleRate; % instead of using index 2, we should actually look for rowmatch
% % % % roi.t2(1) = landmarks.t1(2);

% sync
lnd_audiorow = landmarks.run == op.run & strcmp(landmarks.task,op.task) & cellfun(@(x)contains(x,audfull_name), landmarks.name);
lnd_trialsrow = landmarks.run == op.run & strcmp(landmarks.task,op.task) & cellfun(@(x)contains(x,'trials.tsv'), landmarks.name);
audio_time_minus_trialtab_time = landmarks.t1(lnd_audiorow) - landmarks.t1(lnd_trialsrow);  
audiofiles.starts = trials{:,op.trial_event_start} + audio_time_minus_trialtab_time; %%% starts is not in same time coord system as other times in trials table
audiofiles.ends(1:end-1) = trials{2:end,op.trial_event_end} + audio_time_minus_trialtab_time;
audiofiles.ends(end) = audiofiles.starts(end) + op.last_trial_duration; 
audiofiles.duration = audiofiles.ends - audiofiles.starts; 

    % cut audio trials
%%%% load session audio
% % % cfg = [];
% % % cfg.roi = roi; 
% % % cfg.roi.starts(1) = min(audiofiles.starts) - start_offset; % set the bounds for session audio outside of trials' boundaries
% % % cfg.roi.ends(1) = max(audiofiles.ends) + end_offset; % set the bounds for session audio outside of trials' boundaries
% % % session_aud = bml_load_continuous(cfg);
% % % check_this_trial = 0; % initialize; if this var == 1, open the trial in praat
[run_aud, fs] = audioread(audiofile_full_run); 

% adjust for differnces in file/clock speeds
fs_adj = fs * [(landmarks.t2(lnd_audiorow) - landmarks.t1(lnd_audiorow)) / (landmarks.t2(lnd_trialsrow) - landmarks.t1(lnd_trialsrow))];

for itrial = 1:ntrials
    audiofiles.filename{itrial} = ['trial', num2str(itrial,'%03.f'), '_', trials.name{itrial}, '.wav'];
    trialdur = audiofiles.duration(itrial);
    
    % % % % if isnan(trialdur) || trialdur <= 0 % if erroneous keypress_time
    % % % %     audiofiles.keypress_time(itrial) = NaN; % erase erroneous time
    % % % %     % for trials with negative or NaN duration, set trial end as X many seconds before 
    % % % %     % the following trial's visual_onset
    % % % %     if itrial < ntrials % if not last trial
    % % % %         new_endtime = audiofiles.visual_onset(itrial+1) - corrected_trial_end_prestim; 
    % % % %     % if erroneous trial is the last trial, use X seconds as trial duration
    % % % %     elseif itrial == ntrials
    % % % %         new_endtime = audiofiles.starts(itrial) + last_trial_duration_correction;
    % % % %     end
    % % % %     audiofiles.ends(itrial) = new_endtime; clear new_endtime
    % % % %     audiofiles.duration(itrial) = audiofiles.ends(itrial) - audiofiles.starts(itrial); % correct the duration
    % % % %     check_this_trial = 1; 
    % % % %     warning([audiofiles.filename{itrial}, ' duration was originally erroneously labeled as ',...
    % % % %         num2str(trialdur), ' seconds. Replacing with ', num2str(audiofiles.duration(itrial)),...
    % % % %         'sec. See resulting audio file in Praat.'])
    % % % % elseif trialdur > max_trialdur
    % % % % audiofiles.ends(itrial) = audiofiles.starts(itrial) + max_trialdur; % cut trial duration
    % % % % audiofiles.duration(itrial) = max_trialdur; % correct the duration
    % % % % check_this_trial = 1; 
    % % % % warning([audiofiles.filename{itrial}, ' duration is unexpectedly long (',...
    % % % %     num2str(trialdur), ' seconds). Cutting down to ', num2str(max_trialdur), ' seconds. ', ...
    % % % %     'See resulting audio file in Praat.'])
    % % % % end

    % apply manual edits to trialtimes
    if ~isempty(trials_to_modify_inds) && any(itrial == trials_to_modify_inds)
        matchrow = itrial == trials_to_modify_inds;
        start_time_edit = trials_to_modify_sec(matchrow, 1);
        end_time_edit = trials_to_modify_sec(matchrow, 2);
        audiofiles.starts(itrial) = audiofiles.starts(itrial) + start_time_edit; 
        audiofiles.ends(itrial) = audiofiles.ends(itrial) + end_time_edit; 
    end

    % % % cfg=[];
    % % % cfg.epoch=audiofiles(itrial,:);
    % % % % % % cfg.epoch.ends(1) = cfg.epoch.ends(1) + post_keypress_extension_sec; % stoptrial buffer
    % % % % % % cfg.epoch.duration(1) = cfg.epoch.duration(1) + post_keypress_extension_sec; % stoptrial buffer
    % % % thistrial = bml_redefinetrial(cfg,session_aud);
    % % % fs = thistrial.fsample; 
    % % % trialaud = thistrial.trial{1};

    % Convert times to sample indices
    startSample = max(1, round(audiofiles.starts(itrial) * fs_adj));
    endSample   = min(size(run_aud, 1), round(audiofiles.ends(itrial) * fs_adj));


    % % % % make sure columns are channels, rows are data points
    % % % [~, longdim] = max(size(trialaud));         [~, shortdim] = min(size(trialaud)); 
    % % % trialaud = permute(trialaud, [longdim, shortdim]); 

    trialaud = run_aud(startSample:endSample, :);

    audiowrite([paths.trial_audio_task filesep audiofiles.filename{itrial}], trialaud, fs)
    % % % if check_this_trial && open_problematic_trials_in_praat % load created audio file in praat
    % % %     cmd = ['praat --open ',  [audiofiles.dir{itrial} filesep audiofiles.filename{itrial}], ' &'];
    % % %     system(cmd);
    % % %     check_this_trial = 0; 
    % % % end
end

% audiofiles = renamevars(audiofiles,{'starts','ends','duration'},...
%                         {'audfile_start','audfile_end','audfile_dur'}); 