% main script for running DBS-Learn multisyllabic and subsyllabic experiments
%   after starting this script, select appropriate config file for the experimental phase you intend to run


function dbs_learn_run_exp(op)



starting_trial = 1; %%% if not ==1, currently causes a crash when compiling trial data at end of run

FLAG_SEND_EVENT_STIM_ONSET = 1;


% set priority for matlab to high for running experiments
system(sprintf('wmic process where processid=%d call setpriority "high priority"', feature('getpid')));

close all force

paths = struct; setpaths_dbs_learn();
op.task_computer = getenv('COMPUTERNAME'); 

%% audio device setup



field_default('op','ses','subsyl'); % 'subsyl' or 'multisyl'
field_default('op','require_keypress_every_trial',0); % if true, experimenter must press any key at end of trial to proceed to next trial
field_default('op','ntrials_between_breaks',50); 
field_default('op','ortho_font_size',50);
field_default('op','is_dbs_run',0); % if yes, will try to send beacon pulses for syncing
field_default('op','visual', 'orthography'), 
field_default('op','subject','TEST01'), 
field_default('op','task', 'famil'), 
field_default('op','deviceHead','')      
field_default('op','beepoffset',0.1)   

% starting clock
CLOCKp = ManageTime('start');
TIME_PREPARE = 0.5; % Waiting period before experiment begin (sec)
runtimer = tic; % timer for elapsed time in this run

% timing params inherited from FLVoice_run.m.... a lot of these descriptions won't be meaningful outside the fmri context :
%       timePostStim                : time (s) from end of the audio stimulus presentation to the GO signal (D1 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [.25 .75] 
%       timePostOnset               : time (s) from subject's voice onset to the scanner trigger (or to pre-stimulus segment, if scan=false) (D2 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [4.5] 
%       timeMax                     : maximum time (s) before GO signal and scanner trigger (or to pre-stimulus segment, if scan=false) (D3 in schematic above) (recording portion in a trial may end before this if necessary to start scanner) [5.5] 
%       timePreStim                 : time (s) from end of scan to start of next trial stimulus presentation (D5 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [.25] 
%       timePreSound                : time (s) from start of orthographic presentation to the start of sound stimulus (D6 in schematic above) [.5]
%       timePostSound               : time (s) from end of sound stimulus to the end of orthographic presentation
field_default('op','timePostStim', [.25 .75]), 
field_default('op','timePostOnset', 4.5), 
field_default('op','timePreStim', .25), 
field_default('op','timeMax', 5.5), 
field_default('op','timePreSound', .5), 
field_default('op','timePostSound', .47), 



paths.data_sub = [paths.data, filesep, op.sub]; 
paths.data_ses_beh = [paths.data_sub, filesep, op.ses, filesep, 'beh']; % behavioral data folder for the session (most outputs of this script)
[trials, op] = generate_trial_table(op); 

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
paths.run_exp_op_file = [paths.data_ses_beh, filesep, filestr,'run-exp-op.tsv'];
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
% % % % % % % % info=audiodevinfo;
% % % % % % % % strOUTPUT={info.output.Name};
% % % % % % % % if ~ismember(aud.device_out, strOUTPUT), 
% % % % % % % %     aud.device_out=strOUTPUT{find(strncmp(lower(aud.device_out),lower(strOUTPUT),numel(aud.device_out)),1)}; 
% % % % % % % % end
% % % % % % % % assert(ismember(aud.device_out, strOUTPUT), 'unable to find match to deviceHead name %s',aud.device_out);
% % % % % % % % % % % % % [ok,ID]=ismember(aud.device_out, strOUTPUT);
[beep_wav, beep_fs] = audioread([paths.stim, filesep,'flvoice_run_beep.wav']);
op.beep_dur = numel(beep_wav)/beep_fs;

beepread = dsp.AudioFileReader([paths.stim, filesep,'flvoice_run_beep.wav'], 'SamplesPerFrame', 2048);
headwrite = audioDeviceWriter('SampleRate',beepread.SampleRate,'Device',aud.device_out);


%% timing checks inherited from FLvoice_run.... AM not sure what these are actually for
assert(all(isfinite(op.timePostStim))&ismember(numel(op.timePostStim),[1,2]), 'timePostStim field must have one or two elements');
assert(all(isfinite(op.timePostOnset))&ismember(numel(op.timePostOnset),[1,2]), 'timePostOnset field must have one or two elements');
assert(all(isfinite(op.timePreStim))&ismember(numel(op.timePreStim),[1,2]), 'timePreStim field must have one or two elements');
assert(all(isfinite(op.timeMax))&ismember(numel(op.timeMax),[1,2]), 'timeMax field must have one or two elements');
if numel(op.timePostStim)==1, op.timePostStim=op.timePostStim+[0 0]; end
if numel(op.timePostOnset)==1, op.timePostOnset=op.timePostOnset+[0 0]; end
if numel(op.timePreStim)==1, op.timePreStim=op.timePreStim+[0 0]; end
if numel(op.timeMax)==1, op.timeMax=op.timeMax+[0 0]; end
op.timePostStim=sort(op.timePostStim);
op.timePostOnset=sort(op.timePostOnset);
op.timePreStim=sort(op.timePreStim);
op.timeMax=sort(op.timeMax);

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
        beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
        save(paths.beacon_times_fname, 'beacon_times');
        answer = questdlg('Repeat pulse or proceed to experiment?','','Repeat pulse','Proceed to experiment','Repeat pulse');
        if char(answer) == "Proceed to experiment"
            repeat_beacon = 0;
        end
    end
   
    pause()
end
%>>> ZY addition
% # Initialize Taskcontrol
taskState = struct('task_isRunning',true,'pause_requested',false,'pause_isActive',false);
figTC=taskControlGUI_release(taskState);
% # Initialize EvtTime
if FLAG_SEND_EVENT_STIM_ONSET
    evt = []; evtCode = [];
    paths.trig_events_tab_fname = [paths.data_ses_beh, filesep, filestr,'trig-events.mat']
end
%<<<

% save paths and audio info into ops struct
op.paths_run_exp = paths; 
op.audiodev = aud; 
save(paths.run_exp_op_file, 'op'); % save ops structure, including paths

%% Main trial loop
for itrial = starting_trial:op.ntrials

    % set up trial (see subfunction at end of script)
    % print progress to window
    fprintf(['Run ',num2str(op.run), ', trial ',num2str(itrial), '/',num2str(op.ntrials), ' ....... ', trials.name{itrial}, '\n'])

    if (mod(itrial,op.ntrials_between_breaks) == 0) && (itrial ~= op.ntrials)  % Break after every X trials  , but not on the last
        pause()

        if op.is_dbs_run
            %%% send signal to percept here
            beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
            save(paths.beacon_times_fname,'beacon_times');
            pause()
        end

    end
    % >>>> ZY addition
    % check task control
    if ~exist('figTC','var') || ~ishandle(figTC)
        figTC=taskControlGUI_release(taskState);
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
        fprintf('Task Resumed ...\n')
    end
    %<<<
    %
    set(annoStr.Plus, 'Visible','on');

  
    trials.time_post_stim(itrial) = op.timePostStim(1) + diff(op.timePostStim).*rand; 

    % these 4 values are exactly the same across trials
    TIME_PRE_STIM(itrial) = op.timePreStim(1) + diff(op.timePreStim).*rand; 
    TIME_MAX(itrial) = op.timeMax(1) + diff(op.timeMax).*rand; 
    TIME_POST_SOUND(itrial) = op.timePostSound;
    TIME_PRE_SOUND(itrial) = op.timePreSound;


    audiorow = strcmp(trials.name{itrial}, stim_audio.name); % find the row of the audio file to play
    stimread = stim_audio.stimreads{audiorow};
    trials.time_stim(itrial) = size(stim_audio.audio{audiorow},1)/stim_audio.fs{audiorow}; 

    if isempty(CLOCK)
        CLOCK = ManageTime('start');                        % resets clock to t=0 (first-trial start-time)
        TIME_TRIAL_START = 0;
        TIME_STIM_START = 0;
    else
        TIME_TRIAL_START = ManageTime('current', CLOCK);
    end
    
    % >>>> ZY addition
    % I put this before AM defines TIME_TRIAL_START, only so that I don't add additional delay to actual presentation
    % but we can figure out a better way to integrate

    % can tblEvt be saved as a .tsv table rather than mat? can it be combined with the beacon times table? 

    if FLAG_SEND_EVENT_STIM_ONSET
        try
            [evt_,evtCode_] = send_event([1],[],0.1,0.04,1,'Dev2');
            evt = cat(1,evt,evt_); evtCode = cat(1,evtCode,evtCode_);
            % ## Ideally, we should move the following to the end of trial ##
            tblEvt = table(evt,evtCode,'VariableNames',{'EventTime_dn','EventCode'});
            save(paths.trig_events_tab_fname, 'tblEvt');
        end
    end
    % <<<

    ok=ManageTime('wait', CLOCK, TIME_STIM_START);
    TIME_STIM_ACTUALLYSTART = ManageTime('current', CLOCK);
   
    set(annoStr.Plus, 'Visible','off');
    set(annoStr.Stim, 'String', trials.name{itrial});
    set(annoStr.Stim, 'Visible','on');
    drawnow;
    trials.time_stim_vis_on(itrial) = ManageTime('current', CLOCK);

    if ~ok, fprintf('i am late for this trial TIME_STIM_START\n'); end

    TIME_SOUND_START = TIME_STIM_ACTUALLYSTART + TIME_PRE_SOUND(itrial);
    %ok=ManageTime('wait', CLOCK, TIME_SOUND_START - stimoffset);
    ok=ManageTime('wait', CLOCK, TIME_SOUND_START);

     %%% AM note 2026/2/15: if time is measured before starting the while loop or after the first iteration of the while loop...
     % %  ..... this is generally a 20ms difference (likely due to the size of the audio chunk being written)
    %  % ..... commnent in the tt and ii lines to test this, and compare to trials.time_stim_on(itrial)
     % % (running on stryx laptop)
    trials.time_stim_aud_on(itrial) = ManageTime('current', CLOCK);
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

    % assume stim-off time based on stim-on time and stim duration
    % ... this stim-off time is exactly as accurate as the stim-on measurement
    stimdur = stim_audio.dur(audiorow); 
    trials.time_stim_aud_off(itrial) = trials.time_stim_aud_on(itrial) + stimdur; 


    TIME_SOUND_END = trials.time_stim_aud_on + trials.time_stim(itrial);           % stimulus ends
    if ~ok, fprintf('i am late for this trial TIME_SOUND_START\n'); end

    TIME_ALLSTIM_END = TIME_SOUND_END + TIME_POST_SOUND(itrial);
    ok=ManageTime('wait', CLOCK, TIME_ALLSTIM_END);
    if strcmp(op.visual, 'orthography')
        set(annoStr.Stim, 'Visible','off');
        set(annoStr.Plus, 'Visible','on');
        drawnow;
        trials.time_stim_vis_off(itrial) = ManageTime('current', CLOCK);
    end
    if ~ok, fprintf('i am late for this trial TIME_ALLSTIM_END\n'); end        

    TIME_GOSIGNAL_START = TIME_ALLSTIM_END + trials.time_post_stim(itrial);          % GO signal time


    ok=ManageTime('wait', CLOCK, TIME_GOSIGNAL_START - op.beepoffset);     % waits for recorder initialization time
    
    ok=ManageTime('wait', CLOCK, TIME_GOSIGNAL_START);     % waits for GO signal time


    %%% AM note 2026/2/15: if time is measured before or after playing the go beep, there is 5-25ms difference (running on stryx laptop)
    %     this may be because it reads in the beep as a whole chunk 
    trials.time_go_on(itrial) = ManageTime('current', CLOCK); 
    while ~isDone(beepread); sound=beepread();headwrite(sound);end;reset(beepread);reset(headwrite);

    % assume beep-off time based on beep-on time and beep duration
    % ... this beep-off time is exactly as accurate as the beep-on measurement
    trials.time_go_off(itrial) = trials.time_go_on(itrial) + op.beep_dur; 

    if strcmp(op.visual, 'fixpoint'),set(annoStr.Plus, 'color','g');drawnow;end
    if strcmp(op.visual, 'orthography'),set(annoStr.Plus, 'color','g');drawnow;end
    if ~ok, fprintf('i am late for this trial TIME_GOSIGNAL_START\n'); end


    TIME_SCAN_START = trials.time_go_on(itrial) + TIME_MAX(itrial);


   
    switch op.visual
        case 'fixpoint'
            set(annoStr.Plus, 'color','w');
        case 'orthography'
            set(annoStr.Plus, 'color','w');
    end

        NEXTTRIAL = TIME_SCAN_START + TIME_PRE_STIM(itrial);

    TIME_STIM_START = NEXTTRIAL; 

    % save timing info for this trial into the tsv 
    writetable(trials, paths.trial_info_file , 'Delimiter', '\t', 'FileType', 'text')


    if op.require_keypress_every_trial
        if itrial ~= op.ntrials
            fprintf('\n Press any key to proceed to next trial \n')
        end
        
        pause()
    end

end

release(headwrite);
release(beepread);



%% end of experiment
close all

% experiment time
op.elapsed_time = toc(runtimer)/60;    % elapsed time of the experiment
save(paths.run_exp_op_file, 'op'); % save ops structure, including paths
fprintf('\nElapsed Time: %f (min)\n', op.elapsed_time)

fprintf('Press any key to send final sync pulses and end this experimental phase')
pause()
beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
    save(paths.beacon_times_fname,'beacon_times');


   
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