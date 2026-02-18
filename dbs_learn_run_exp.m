% main script for running DBS-Learn multisyllabic and subsyllabic experiments
%   after starting this script, select appropriate config file for the experimental phase you intend to run
%
% capitalized times like "TIME_VIS_STIM_OFF" are time targets and are computed each trial based on previous timepoints; they do not get saved
% .... for each target time, we measure (approximately) the actual time of that event
% .......... and save it in the trial table with formatting like "trials.t_vis_stim_off(itrial)"
%
% what this script saves for each run: 
%   trials table .tsv - contains event times, stim info; saved every trial
%   op .mat - contains options structure and other info not specific to individual trials
%   beacon times and beacon event table file

function dbs_learn_run_exp(op)

starting_trial = 1; %%% if not ==1, currently causes a crash when compiling trial data at end of run
FLAG_SEND_EVENT_STIM_ONSET = 1;


% set priority for matlab to high for running experiments
system(sprintf('wmic process where processid=%d call setpriority "high priority"', feature('getpid')));

close all force

parpool_obj = gcp('nocreate');
if ~isempty(parpool_obj)
    cancel(parpool_obj.FevalQueue.RunningFutures);   % new-style futures queue (R2022b+)
end

paths = struct; setpaths_dbs_learn();
op.task_computer = getenv('COMPUTERNAME'); 

%% audio device setup
field_default('op','sub','qqq'); 
field_default('op','sestask', 'subsyl_famil'); 
% field_default('op','ses','subsyl'); % 'subsyl' or 'multisyl'
% field_default('op','task', 'famil'), 
field_default('op','max_repeated_trials', 3); 
field_default('op','record_audio', 1); 
field_default('op','require_keypress_every_trial',0); % if true, experimenter must press any key at end of trial to proceed to next trial
field_default('op','vis_offset_to_go',[0.25, 0.75]); % min and max of delay (jittered) between visual offset and GO cue presentation
field_default('op','ntrials_between_breaks',50); 
field_default('op','ortho_font_size',75); % 80 font size is max that will fit on 1920x1080 screen for 7syl
field_default('op','visual','orthography'); 
field_default('op','is_dbs_run',1); % if yes, will try to send beacon pulses for syncing
field_default('op','visual', 'orthography'), 
field_default('op','deviceHead','')      
field_default('op','beepoffset',0.1)   
field_default('op','subj_n_syls',7)   

% starting clock
CLOCKp = ManageTime('start');
TIME_PREPARE = 0.5; % Waiting period before experiment begin (sec)
runtimer = tic; % timer for elapsed time in this run

[trials, op] = generate_trial_table(op); 
paths.data_sub = [paths.data, filesep, op.sub]; 
paths.data_ses_beh = [paths.data_sub, filesep, op.ses, filesep, 'beh']; % behavioral data folder for the session (most outputs of this script)
paths.data_ses_audvid = [paths.data_sub, filesep, op.ses, filesep, 'audio-video']; 

% set session-specific timing params
%%% op.gobeep_to_next_trial is the speech window - time between GO cue onset and when the next trial's stim is presented
switch op.ses
    case 'subsyl'
        op.vis_stim_dur = 1.5; % 1.5s is similar to vis stim from intrasurgical experiment in richardson lab
        op.gobeep_to_next_trial = 4; 
    case 'multisyl'  
        op.vis_stim_dur = 3; % maximum audio stim length (at 7 syllables) is 2.8sec
        op.gobeep_to_next_trial = 6; 
end

% if forward digit span task, always require keypress on every trial 
%%%%%  (so we can score it online and end when an error occurs)
% do not show orthography during FDS (this task usually is audio-only) 
if op.task=="fds"
    op.require_keypress_every_trial = 1;
    op.visual = 'fixpoint'; 
end

%%%% parameters for sync pulses sent to CED computer [set by Zeyang]
op.pulse.interval = 0.4; 
op.pulse.count = 3; 
op.pulse.duration = 0.05; 

if ~exist(paths.data_ses_beh, 'dir')
    mkdir(paths.data_ses_beh)
end


% Check existing trial info files to increment run ID for this session
allEventFiles = dir([paths.data_ses_beh, filesep, '*_trials.tsv']);
if ~isempty(allEventFiles)
    prevRunIds = regexp({allEventFiles.name}, '_run-(\d+)_', 'tokens', 'ignorecase');
    prevRunIds = cellfun(@(x) str2double(x{1,1}), prevRunIds, 'UniformOutput', true);
    op.run = max(prevRunIds) + 1;
else
    op.run = 1;
end

filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_']; % this string gets used in a variety of files associated with this run

% specifying paths for this run
paths.run_exp_op_file = [paths.data_ses_beh, filesep, filestr,'run-exp-op.mat'];
paths.audio_stim_ses = [paths.code_dbs_learn, filesep, 'stimuli', filesep, 'audio-',op.ses]; 
paths.audio_stim_task = [paths.audio_stim_ses, filesep, op.task]; % contains the main audio stim files for this task
paths.trial_info_file = [paths.data_ses_beh, filesep, filestr,'trials.tsv'];
writetable(trials, paths.trial_info_file , 'Delimiter', '\t', 'FileType', 'text')


% visual setup
annoStr = setUpVisAnnot_HW([0 0 0]);
annoStr.Stim.FontSize = op.ortho_font_size; 
set(annoStr.Stim, 'String', 'Preparing...');
set(annoStr.Stim, 'Visible','on');


%% AUDIO SETUP
% modify the below function based on the task computer and audio devices you are using
aud = setup_audio_devices(); 

if op.record_audio
    if ~exist(paths.data_ses_audvid, 'dir')
        mkdir(paths.data_ses_audvid) %%%% recording function fails if the dir doesn't already exist
    end
    paths.aud_record_filename_1 = [paths.data_ses_audvid, filesep, filestr,'audrec_1.wav'];
    paths.aud_record_filename_2 = [paths.data_ses_audvid, filesep, filestr,'audrec_2.wav'];

    aud_record_obj_1 = parfeval(@recordMonoDevice, 0, ...
       aud.device_in_1, paths.aud_record_filename_1);
    aud_record_obj_2 = parfeval(@recordMonoDevice, 0, ...
       aud.device_in_2, paths.aud_record_filename_2);

end


 % load the stim audio in this table and play them, rather than loading them on every trial
stim_audio = table;
stim_audio.filename = cellfun(@(x)[paths.audio_stim_task, filesep, x,'.wav'], unique(trials.name, 'stable'), 'uni', 0); 
stim_audio.name = cellfun(@(x)getfname(x), stim_audio.filename, 'uni', 0); 

ok=cellfun(@(x)exist(x,'file'), stim_audio.filename); 
assert(all(ok), 'unable to find files %s', sprintf('%s ',stim_audio.filename{~ok}));

[stim_audio.audio,stim_audio.fs]=cellfun(@audioread, stim_audio.filename,'uni',0);
stim_audio.stimreads=cell(size(stim_audio.filename));
stim_audio.stimreads=cellfun(@(x)dsp.AudioFileReader(x, 'SamplesPerFrame', 2048),stim_audio.filename,'uni',0);
sileread = dsp.AudioFileReader(fullfile(paths.stim, 'silent.wav'), 'SamplesPerFrame', 2048);
stim_audio.dur=cellfun(@(a,b)numel(a)/b, stim_audio.audio, stim_audio.fs);

% set up sound output players
audiodevreset;
[beep_wav, beep_fs] = audioread([paths.stim, filesep,'flvoice_run_beep.wav']);
op.beep_dur = numel(beep_wav)/beep_fs;

beepread = dsp.AudioFileReader([paths.stim, filesep,'flvoice_run_beep.wav'], 'SamplesPerFrame', 2048);
headwrite = audioDeviceWriter('SampleRate',beepread.SampleRate,'Device',aud.device_out);

ok=ManageTime('wait', CLOCKp, TIME_PREPARE);
set(annoStr.Stim, 'Visible','off');     % Turn off preparation page
TIME_PREPARE_END=ManageTime('current', CLOCKp);

%% prep screen
set(annoStr.Stim, 'String', 'READY');
set(annoStr.Stim, 'Visible','on');
while ~isDone(sileread); sound=sileread();headwrite(sound);end;release(sileread);reset(headwrite);
ok=ManageTime('wait', CLOCKp, TIME_PREPARE_END+2);
set(annoStr.Stim, 'Visible','off');     % Turn off preparation page
CLOCK=[];                               % Main clock (not yet started)


%% BEACON SYNC SETUP......  can we make beacon times into a tsv file instead of mat file? 
if op.is_dbs_run
    %%% send signal to percept here
    %%%%% give option for experimenter to send more pulses to calibrate
    repeat_beacon = 1;
    beacon_times = []; 
    % >>> ZY mod: moved outside of while-loop, so all events saved correctly
    paths.beacon_times_fname = [paths.data_ses_beh, filesep, filestr,'beacon-times.mat'];
    %<<<
    while repeat_beacon
        beacon_times = [beacon_times, test_Beacon(op.pulse.interval,op.pulse.duration,op.pulse.count)];
        save(paths.beacon_times_fname, 'beacon_times');
        answer = questdlg('Repeat pulse or proceed to experiment?','','Repeat pulse','Proceed to experiment','Repeat pulse');
        if char(answer) == "Proceed to experiment"
            repeat_beacon = 0;
        end
    end
   
    % pause()
end
%>>> ZY addition
% # Initialize Taskcontrol
taskState = struct('task_isRunning',true,'pause_requested',false,'pause_isActive',false);
figTC=taskControlGUI_release(taskState,op.pulse,paths.beacon_times_fname,beacon_times);

% # Initialize EvtTime
if FLAG_SEND_EVENT_STIM_ONSET
    evt = []; evtCode = [];
    paths.trig_events_tab_fname = [paths.data_ses_beh, filesep, filestr,'trig-events.mat']
end
%<<<

% save paths and audio info into ops struct
op.paths_run_exp = paths; 
op.audiodev = aud; 
op % show options on command line
save(paths.run_exp_op_file, 'op'); % save ops structure, including paths

%% Main trial loop
for itrial = starting_trial:op.ntrials

    % set up trial (see subfunction at end of script)



    if isempty(CLOCK)
        CLOCK = ManageTime('start');                        % resets clock to t=0 (first-trial start-time)
        TIME_STIM_START = 0;
    end
    
    ok=ManageTime('wait', CLOCK, TIME_STIM_START);

    if (mod(itrial,op.ntrials_between_breaks) == 0) && (itrial ~= op.ntrials)  % Break after every X trials, but not on the last
        fprintf(['Scheduled break at every ', num2str(op.ntrials_between_breaks), ' trial\n'])

        if op.is_dbs_run
            %%% send signal to percept here
            fprintf([ 'Press any key to send pulses.... \n'])
            pause()

            beacon_times = [beacon_times, test_Beacon(op.pulse.interval,op.pulse.duration,op.pulse.count)];
            save(paths.beacon_times_fname,'beacon_times'); %%%% redundant - ok to delete

            writematrix(beacon_times,strrep(paths.beacon_times_fname,'.mat','.tsv'),'FileType','text','Delimiter','tab')
            
        end
        fprintf([ 'Press any key to resume next trial. \n'])
        pause()
    end
    %% >>>> ZY addition
    % check task control
    if ~exist('figTC','var') || ~ishandle(figTC)
        figTC=taskControlGUI_release(taskState,op.pulse,paths.beacon_times_fname,beacon_times);
    end
    ud = get(figTC,'UserData'); taskState = ud.taskState;
    if taskState.pause_requested
        % Create a timer that fires every second
        % pause requested => start pause
        taskState.pause_isActive = 1;
        ud.taskState = taskState; set(figTC,'UserData',ud);
        ud.updateGUIBasedOnTaskState();
        fprintf('Task Paused from TaskControl_Panel ...\n')
        while taskState.pause_requested 
            pause(0.2)
            ud = get(figTC,'UserData'); taskState = ud.taskState;
        end
        % resuming
        ud = get(figTC,'UserData'); taskState = ud.taskState;
        taskState.pause_isActive = 0; % no longer in pause
        ud.taskState = taskState; set(figTC,'UserData',ud); 
        ud.updateGUIBasedOnTaskState();
        % saving updated beaconTimes
        if length(ud.beacon_times)>=length(beacon_times)
            beacon_times = ud.beacon_times;
            save(paths.beacon_times_fname,'beacon_times'); %%%% redundant - ok to delete
        end
        fprintf('Task Resumed ...\n')
    end
    %<<<
    %%
    set(annoStr.Plus, 'Visible','on');

    audiorow = strcmp(trials.name{itrial}, stim_audio.name); % find the row of the audio file to play
    stimread = stim_audio.stimreads{audiorow};

    %% >>>> ZY addition
    % I put this before AM defines TIME_TRIAL_START, only so that I don't add additional delay to actual presentation
    % but we can figure out a better way to integrate

    % can tblEvt be saved as a .tsv table rather than mat? can it be combined with the beacon times table? 

    if FLAG_SEND_EVENT_STIM_ONSET
        % code 3 => sending on DI port 3 of Dev 2
        [evt_,evtCode_] = send_event([3],[],0.1,0.04,1,'Dev2'); 
        evt = cat(1,evt,evt_); evtCode = cat(1,evtCode,evtCode_);
        % ## Ideally, we should move the following to the end of trial ##
        tblEvt = table(evt,evtCode,'VariableNames',{'EventTime_dn','EventCode'});
        save(paths.trig_events_tab_fname, 'tblEvt');
    end
    %% <<<
    


   
    % show stim orthography
    if strcmp(op.visual, 'orthography')
        set(annoStr.Plus, 'Visible','off');
        set(annoStr.Stim, 'String', trials.name{itrial});
        set(annoStr.Stim, 'Visible','on');
    elseif strcmp(op.visual, 'fixpoint')
        set(annoStr.Plus, 'Visible','on'); % leave fixcross on screen
        set(annoStr.Plus, 'color', 'w'); % turn cross white while playing audio stim; indicate to not speak during stim
    end
    drawnow;
    trials.t_stim_vis_on(itrial) = ManageTime('current', CLOCK);


    % play stim sound - ASAP after visual onset
     %%% AM note 2026/2/15: if time is measured before starting the while loop or after the first iteration of the while loop...
     % %  ..... this is generally a 20ms difference (likely due to the size of the audio chunk being written)
    %  % ..... commnent in the tt and ii lines to test this, and compare to trials.t_stim_on(itrial)
     % % (running on stryx laptop)
    trials.t_stim_aud_on(itrial) = ManageTime('current', CLOCK);
    % ii = 0; 
    while ~isDone(stimread); 
        sound=stimread();
        headwrite(sound);
        % ii=ii+1;
        % if ii==1
            % tt(ii) = ManageTime('current', CLOCK); 
        % end
    end;
    release(stimread);
    reset(headwrite);

    % print progress to window
    fprintf(['Run ',num2str(op.run), ', trial ',num2str(itrial), '/',num2str(op.ntrials), ' ....... ', trials.name{itrial}, '\n'])

    if ~ok, fprintf('i am late for this trial TIME_STIM_START\n'); end

    % assume audio stim-off time based on audio stim-on time and audio stim duration
    % ... this stim-off time is exactly as accurate as the stim-on measurement
    audstimdur = stim_audio.dur(audiorow); 
    trials.t_stim_aud_off(itrial) = trials.t_stim_aud_on(itrial) + audstimdur; 

    % wait for a fixed time after visual onset, then switch back to fixation cross
    TIME_VIS_STIM_OFF = trials.t_stim_vis_on(itrial) + op.vis_stim_dur; 
    ManageTime('wait', CLOCK, TIME_VIS_STIM_OFF);
    if strcmp(op.visual, 'orthography')
        set(annoStr.Stim, 'Visible','off');
        set(annoStr.Plus, 'color',[1 1 1]);
        set(annoStr.Plus, 'Visible','on');
        drawnow;
        trials.t_stim_vis_off(itrial) = ManageTime('current', CLOCK);
    end

    % plan GO signal start time based on visual offset time
    % .... do not base off of audio offset, because this varies depending on the number of syllables in multisyl
    trials.go_delay(itrial) = min(op.vis_offset_to_go) + rand*diff(op.vis_offset_to_go); % get jittered time between stim offset and go cue
    TIME_GOSIGNAL_START = trials.t_stim_vis_off(itrial)  + trials.go_delay(itrial);          % GO signal target time
    ok=ManageTime('wait', CLOCK, TIME_GOSIGNAL_START);     % waits for GO signal time


    %%% AM note 2026/2/15: if time is measured before or after playing the go beep, there is 5-25ms difference (running on stryx laptop)
    %     this may be because it reads in the beep as a whole chunk 
    trials.t_go_aud_on(itrial) = ManageTime('current', CLOCK); 
    while ~isDone(beepread); sound=beepread();headwrite(sound);end;reset(beepread);reset(headwrite);

    % assume beep-off time based on beep-on time and beep duration
    % ... this beep-off time is exactly as accurate as the beep-on measurement
    trials.t_go_aud_off(itrial) = trials.t_go_aud_on(itrial) + op.beep_dur; 

    if strcmp(op.visual, 'fixpoint') || strcmp(op.visual, 'orthography')
        set(annoStr.Plus, 'color','g');
        drawnow;
    end
    if ~ok, fprintf('i am late for this trial TIME_GOSIGNAL_START\n'); end

    % show green cross ASAP after starting to play GO beep
    set(annoStr.Plus, 'color','g');
    trials.t_go_vis_on(itrial) = ManageTime('current', CLOCK); 

    % set the target start time for the start of the next trial's onset time
    TIME_STIM_START = trials.t_go_aud_on(itrial) + op.gobeep_to_next_trial; 

    % save timing info for this trial into the tsv 
    writetable(trials, paths.trial_info_file , 'Delimiter', '\t', 'FileType', 'text')


    if op.require_keypress_every_trial
        if itrial ~= op.ntrials
            fprintf('\n Press any key to proceed to next trial \n')
        end
        
        pause()
    end

end






%% end of experiment

% experiment time
op.elapsed_time = toc(runtimer)/60;    % elapsed time of the experiment
save(paths.run_exp_op_file, 'op'); % save ops structure, including paths
fprintf('\nElapsed Time: %f (min)\n', op.elapsed_time)

fprintf('Press any key to send final sync pulses and end this experimental phase')
pause()
beacon_times = [beacon_times, test_Beacon(op.pulse.interval,op.pulse.duration,op.pulse.count)];
save(paths.beacon_times_fname,'beacon_times');

% pause to make sure audio stim have finished playing, then close audio writer objects
pause(1)
release(headwrite);
release(beepread);

cancel(aud_record_obj_1)
cancel(aud_record_obj_2)

close all

   
end







function out = ManageTime(option, varargin)
% MANAGETIME time-management functions for real-time operations
%
% CLCK = ManageTime('start');             initializes clock to t=0
% T = ManageTime('current', CLCK);        returns current time T (measured in seconds after t=0)
% ok = ManageTime('wait', CLCK, T);       waits until time = T (measured in seconds after t=0)
%                                         returns ok=false if we were already passed T
%
% e.g.
%  CLCK = ManageTime('start');
%  ok = ManageTime('wait', CLCK, 10);
%  disp(ManageTime('current', CLCK));
%  disp(ManageTime('current', CLCK));
%  disp(ManageTime('current', CLCK));
%

DEBUG=false;
switch(lower(option))
    case 'start',
        out=clock;
    case 'wait'        
        out=true;
        t0=varargin{1};
        t1=varargin{2};
        t2=etime(clock,t0);
        if DEBUG, fprintf('had %f seconds to spare\n',t1-t2); end
        if t2>t1, out=false; return; end % i am already late
        if t1-t2>1, pause(t1-t2-1); end % wait until ~1s to go (do not waste CPU)
        while etime(clock,t0)<t1, end % wait until exactly the right time (~ms prec)
    case 'current'
        t0=varargin{1};
        out=etime(clock,t0);
end
end

function checkAndStopTimer(hFig, t)
    % Ensure the figure exists to avoid errors if the figure is closed
    if ishandle(hFig)
        % Check if UserData.A is 1
        if hFig.UserData.taskState.pause_requested == 0
            disp('Resuming ..... (by task control)');
            stop(t);  % Stop the timer
            delete(t);  % Clean up timer object
        end
    else
        disp('Figure handle is not valid. Stopping timer.');
        stop(t);  % Stop the timer if the figure is closed
        delete(t);  % Clean up timer object
    end
end