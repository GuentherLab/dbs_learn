% main script for running DBS-Learn multisyllabic and subsyllabic experiments
%   after starting this script, select appropriate config file for the experimental phase you intend to run


function dbs_learn_run_exp(op)

show_mic_trace_figure = 0; % if false, make mic trace figure invisible

ntrials_between_breaks = 34; 

ortho_font_size = 50;

is_dbs_run = 0; 

starting_trial = 1; %%% if not ==1, currently causes a crash when compiling trial data at end of run

FLAG_SEND_EVENT_STIM_ONSET = 1;

% if true, experimenter must press any key at end of trial to proceed to next trial
require_keypress_every_trial = 0; 

CMRR = true;

% set priority for matlab to high for running experiments
system(sprintf('wmic process where processid=%d call setpriority "high priority"', feature('getpid')));

beepoffset = 0.100;

close all force

%
% timing params inherited from FLVoice_run.m :
%       timePostStim                : time (s) from end of the audio stimulus presentation to the GO signal (D1 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [.25 .75] 
%       timePostOnset               : time (s) from subject's voice onset to the scanner trigger (or to pre-stimulus segment, if scan=false) (D2 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [4.5] 
%       timeMax                     : maximum time (s) before GO signal and scanner trigger (or to pre-stimulus segment, if scan=false) (D3 in schematic above) (recording portion in a trial may end before this if necessary to start scanner) [5.5] 
%       timePreStim                 : time (s) from end of scan to start of next trial stimulus presentation (D5 in schematic above) (one value for fixed time, two values for minimum-maximum range of random times) [.25] 
%       timePreSound                : time (s) from start of orthographic presentation to the start of sound stimulus (D6 in schematic above) [.5]
%       timePostSound               : time (s) from end of sound stimulus to the end of orthographic presentation
%

paths = struct; setpaths_dbs_learn();

%% audio device setup
% modify the below function based on the task computer and audio devices you are using
aud = setup_audio_devices(); 

%%
ET = tic;
if ispc, [nill,host]=system('hostname');
else [nill,host]=system('hostname -f');
end
host=regexprep(host,'\n','');

if strcmp(host, '677-GUE-WL-0009')
    default_fontsize = 10;
else
    default_fontsize = 15;
end


       

expParams = op; 

field_default('expParams','ses','subsyl'), 
field_default('expParams','visual', 'orthography'), 
field_default('expParams','root', pwd), 
% field_default('expParams','audiopath', fullfile(pwd, 'stimuli', ['audio-',stimset])), 
field_default('expParams','subject','TEST01'), 
field_default('expParams','session', 1), 
% % % field_default('expParams','run', 1), 
field_default('expParams','task', 'test'), 
field_default('expParams','gender', 'unspecified'), 
field_default('expParams','timePostStim', [.25 .75]), 
field_default('expParams','timePostOnset', 4.5), 
field_default('expParams','timePreStim', .25), 
field_default('expParams','timeMax', 5.5), 
field_default('expParams','timePreSound', .5), 
field_default('expParams','timePostSound', .47), 
field_default('expParams','rmsThresh', .02),  %'rmsThresh', .05,
field_default('expParams','rmsBeepThresh', .1), 
field_default('expParams','rmsThreshTimeOnset', .02), % 'rmsThreshTimeOnset', .10,
field_default('expParams','rmsThreshTimeOffset', [.25 .25]), 
field_default('expParams','minVoiceOnsetTime', 0.4), 
field_default('expParams','deviceMic',''), 
field_default('expParams','deviceHead','')             

paths.data_sub = [paths.data, filesep, op.sub]; 
paths.data_ses = [paths.data_sub, filesep, op.ses]; 
[trials, op] = generate_trial_table(op); 
paths.data_task = [paths.data_ses, filesep, 'beh', filesep, op.task]; 

if ~exist(paths.data_task, 'dir')
    mkdir(paths.data_task)
end

% Check existing trial info files to increment run ID for this session
allEventFiles = dir([paths.data_task, filesep, '*_trials.tsv']);
if ~isempty(allEventFiles)
    prevRunIds = regexp({allEventFiles.name}, '_run-(\d+)_', 'tokens', 'ignorecase');
    prevRunIds = cellfun(@(x) str2double(x{1,1}), prevRunIds, 'UniformOutput', true);
    op.run = max(prevRunIds) + 1;
else
    op.run = 1;
end

filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run)]; % this string gets used in a variety of files associated with this run


paths.trial_info_file = [paths.data_task, filesep, filestr,'_trials.tsv'];

writetable(trials, paths.trial_info_file , 'Delimiter', '\t', 'FileType', 'text')

expParams.computer = host;
paths.audio_stim_ses = [paths.code_dbs_learn, filesep, 'stimuli', filesep, 'audio-',op.ses]; 
paths.audio_stim_task = [paths.audio_stim_ses, filesep, op.task]; % contains the main audio stim files for this task


audiodevreset;
info=audiodevinfo;
strOUTPUT={info.output.Name};
try, a=audioDeviceWriter('Device','asdf'); 
catch me; str=regexp(regexprep(me.message,'.*Valid values are:',''),'"([^"]*)"','tokens'); 
    strOUTPUT=[str{:}]; 
end;


    
strVisual={'figure', 'fixpoint', 'orthography'};


    annoStr = setUpVisAnnot_HW([0 0 0]);
% % % % % % % % % end
annoStr.Stim.FontSize = ortho_font_size; 

CLOCKp = ManageTime('start');
TIME_PREPARE = 0.5; % Waiting period before experiment begin (sec)
set(annoStr.Stim, 'String', 'Preparing...');
set(annoStr.Stim, 'Visible','on');



% locate files and generate trials 
taskpath = fullfile(expParams.root, sprintf('sub-%s',expParams.subject), sprintf('ses-%d',expParams.session),'beh', expParams.task);





Input_files = cellfun(@(x)[paths.audio_stim_task, filesep, x,'.wav'], unique(trials.name), 'UniformOutput', false); 
% % % % Input_files = cellfun(@(x)[x,'.wav'], unique(trials.name), 'UniformOutput', false);



ok=cellfun(@(x)exist(x,'file'), Input_files); 
assert(all(ok), 'unable to find files %s', sprintf('%s ',Input_files{~ok}));

                % % % % dirFiles=cellfun(@dir, Input_files(NoNull), 'uni', 0);
                % % % % NoNull=NoNull(cellfun(@(x)x.bytes>0, dirFiles));


Input_sound=cell(size(Input_files));
Input_fs=num2cell(ones(size(Input_files)));

 

[Input_sound,Input_fs]=cellfun(@audioread, Input_files,'uni',0);
[silent_sound,silent_fs]=audioread(fullfile(paths.stim, 'silent.wav'));
stimreads=cell(size(Input_files));
stimreads=cellfun(@(x)dsp.AudioFileReader(x, 'SamplesPerFrame', 2048),Input_files,'uni',0);
% % % % stimreads(setdiff(1:numel(stimreads), NoNull))=arrayfun(@(x)dsp.AudioFileReader(fullfile(paths.audio_stim_ses, 'silent.wav'), 'SamplesPerFrame', 2048),1:numel(Input_files(setdiff(1:numel(stimreads), NoNull))),'uni',0);
sileread = dsp.AudioFileReader(fullfile(paths.stim, 'silent.wav'), 'SamplesPerFrame', 2048);



Input_duration=cellfun(@(a,b)numel(a)/b, Input_sound, Input_fs);
%meanInput_duration=mean(Input_duration(Input_duration>0));
silence_dur=size(silent_sound,1)/silent_fs;
[Input_sound{Input_duration==0}]=deal(zeros(ceil(44100*silence_dur),1)); % fills empty audiofiles with average-duration silence ('NULL' CONDITIONS)
[Input_fs{Input_duration==0}]=deal(44100);
[Input_conditions{Input_duration==0}]=deal('NULL');

 

% set up device reader settings for accessing audio signal during recording
expParams.sr = 48000;            % sample frequency (Hz)
frameDur = .050;                 % frame duration in seconds
expParams.frameLength = expParams.sr*frameDur;      % framelength in samples
deviceReader = audioDeviceReader(...
    'Device', aud.device_in, ...
    'SamplesPerFrame', expParams.frameLength, ...
    'SampleRate', expParams.sr, ...
    'BitDepth', '24-bit integer');    
% set up sound output players
if ~ismember(expParams.deviceHead, strOUTPUT), expParams.deviceHead=strOUTPUT{find(strncmp(lower(expParams.deviceHead),lower(strOUTPUT),numel(expParams.deviceHead)),1)}; end
assert(ismember(expParams.deviceHead, strOUTPUT), 'unable to find match to deviceHead name %s',expParams.deviceHead);
[ok,ID]=ismember(expParams.deviceHead, strOUTPUT);
[twav, tfs] = audioread([paths.stim, filesep,'flvoice_run_beep.wav']);
beepdur = numel(twav)/tfs;
%stimID=info.output(ID).ID;
%beepPlayer = audioplayer(twav*0.2, tfs, 24, info.output(ID).ID);
beepread = dsp.AudioFileReader([paths.stim, filesep,'flvoice_run_beep.wav'], 'SamplesPerFrame', 2048);
%headwrite = audioDeviceWriter('SampleRate',beepread.SampleRate,'Device',expParams.deviceHead, 'SupportVariableSizeInput', true, 'BufferSize', 2048);
headwrite = audioDeviceWriter('SampleRate',beepread.SampleRate,'Device',expParams.deviceHead);

% checks values of timing variables
expParams.beepoffset = beepoffset;

assert(all(isfinite(expParams.timePostStim))&ismember(numel(expParams.timePostStim),[1,2]), 'timePostStim field must have one or two elements');
assert(all(isfinite(expParams.timePostOnset))&ismember(numel(expParams.timePostOnset),[1,2]), 'timePostOnset field must have one or two elements');
assert(all(isfinite(expParams.timePreStim))&ismember(numel(expParams.timePreStim),[1,2]), 'timePreStim field must have one or two elements');
assert(all(isfinite(expParams.timeMax))&ismember(numel(expParams.timeMax),[1,2]), 'timeMax field must have one or two elements');
if numel(expParams.timePostStim)==1, expParams.timePostStim=expParams.timePostStim+[0 0]; end
if numel(expParams.timePostOnset)==1, expParams.timePostOnset=expParams.timePostOnset+[0 0]; end
if numel(expParams.timePreStim)==1, expParams.timePreStim=expParams.timePreStim+[0 0]; end
if numel(expParams.timeMax)==1, expParams.timeMax=expParams.timeMax+[0 0]; end
expParams.timePostStim=sort(expParams.timePostStim);
expParams.timePostOnset=sort(expParams.timePostOnset);
expParams.timePreStim=sort(expParams.timePreStim);
expParams.timeMax=sort(expParams.timeMax);
rmsThresh = expParams.rmsThresh; % params for detecting voice onset %voiceCal.rmsThresh; % alternatively, run a few iterations of testThreshold and define rmsThreshd here with the resulting threshold value after convergence
rmsBeepThresh = expParams.rmsBeepThresh;
% nonSpeechDelay = .75; % initial estimate of time between go signal and voicing start
nonSpeechDelay = .5; % initial estimate of time between go signal and voicing start

%%%%% set up figure for real-time plotting of audio signal of next trial
if show_mic_trace_figure
    rtfig = figure('units','norm','position',[.1 .2 .4 .5],'menubar','none', 'Visible',show_mic_trace_figure);
    micSignal = plot(nan,nan,'-', 'Color', [0 0 0.5]);
    micLine = xline(0, 'Color', [0.984 0.352 0.313], 'LineWidth', 3);
    micLineB = xline(0, 'Color', [0.46 1 0.48], 'LineWidth', 3);
    micTitle = title('', 'Fontsize', default_fontsize-1, 'interpreter','none');
    xlabel('Time(s)');
    ylabel('Sound Pressure');
end

pause(1);
% save(Output_name, 'expParams');

%Initialize trialData structure
trialData = struct;

ok=ManageTime('wait', CLOCKp, TIME_PREPARE);
set(annoStr.Stim, 'Visible','off');     % Turn off preparation page
TIME_PREPARE_END=ManageTime('current', CLOCKp);

set(annoStr.Stim, 'String', 'READY');
set(annoStr.Stim, 'Visible','on');
while ~isDone(sileread); sound=sileread();headwrite(sound);end;release(sileread);reset(headwrite);
ok=ManageTime('wait', CLOCKp, TIME_PREPARE_END+2);
set(annoStr.Stim, 'Visible','off');     % Turn off preparation page
CLOCK=[];                               % Main clock (not yet started)
expParams.timeNULL = expParams.timeMax(1) + diff(expParams.timeMax).*rand;
intvs = [];



if is_dbs_run
    %%% send signal to percept here
    %%%%% give option for experimenter to send more pulses to calibrate
    repeat_beacon = 1;
    beacon_times = []; 
    % >>> ZY mod: moved outside of while-loop, so all events saved correctly
    beacon_times_fname = fullfile(taskpath,sprintf('sub-%s_ses-%d_run-%d_task-%s_beacon-times.mat',expParams.subject, expParams.session, expParams.run, expParams.task));
    %<<<
    while repeat_beacon
        beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
        save(beacon_times_fname, 'beacon_times');
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
    taskEvent_fname = fullfile(taskpath,sprintf('sub-%s_ses-%d_run-%d_task-%s_taskEvents.mat',...
        expParams.subject, expParams.session, expParams.run, expParams.task));
end
%<<<

%% Main trial loop
for itrial = starting_trial:op.ntrials

    % set up trial (see subfunction at end of script)
    %[trialData, annoStr] = setUpTrial(expParams, annoStr, stimName, condition, trialData, ii);
    % print progress to window
    fprintf('\nRun %d, trial %d/%d\n', expParams.run, itrial, op.ntrials);

    if (mod(itrial,ntrials_between_breaks) == 0) && (itrial ~= op.ntrials)  % Break after every X trials  , but not on the last
        pause()

        if is_dbs_run
            %%% send signal to percept here
            beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
            save(beacon_times_fname,'beacon_times');
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
        %assert(~taskState.pause_requested,'pause_requested should be 0, as task is resuming')
        taskState.pause_isActive = 0; % no longer in pause
        ud.taskState = taskState; set(figTC,'UserData',ud); 
        ud.updateGUIBasedOnTaskState();
        fprintf('Task Resumed ...\n')
    end
    %<<<
    %
    set(annoStr.Plus, 'Visible','on');
    trialData(itrial).stimName = trials.name{itrial};
    trialData(itrial).display = trials.name{itrial};
    
                        % trialData(itrial).condLabel = Input_conditions{itrial};
    
                        if strcmp(trialData(itrial).display, 'NULL'); trialData(itrial).display = 'yyy'; end
                        %     trialData(ii).timeStim = numel(Input_sound{ii})/Input_fs{ii}; 
    trialData(itrial).timeStim = size(Input_sound{itrial},1)/Input_fs{itrial}; 
    trialData(itrial).timePostStim = expParams.timePostStim(1) + diff(expParams.timePostStim).*rand; 
    trialData(itrial).timePostOnset = expParams.timePostOnset(1) + diff(expParams.timePostOnset).*rand; 
    % trialData(itrial).timeScan = expParams.timeScan(1) + diff(expParams.timeScan).*rand; 
    trialData(itrial).timePreStim = expParams.timePreStim(1) + diff(expParams.timePreStim).*rand; 
    trialData(itrial).timeMax = expParams.timeMax(1) + diff(expParams.timeMax).*rand; 
    trialData(itrial).timePostSound = expParams.timePostSound;
    trialData(itrial).timePreSound = expParams.timePreSound;
    %stimPlayer = audioplayer(Input_sound{ii},Input_fs{ii}, 24, stimID);
    stimread = stimreads{itrial};
                            % SpeechTrial=~strcmp(trialData(itrial).condLabel,'NULL');
                        %     SpeechTrial=~strcmp(trialData(ii).condLabel,'S');

    % set up variables for audio recording and voice detection
    % % % % % % % prepareScan=0.250*(expParams.scan~=0); % if scanning, end recording 250ms before scan trigger
    prepareScan = 0; 
    recordLen= trialData(itrial).timeMax-prepareScan; % max total recording time
    % % % % % % % % % % % % recordLenNULL = expParams.timeNULL-prepareScan;
    nSamples = ceil(recordLen*expParams.sr);
    % % % % % % % % % % % % % % % % % nSamplesNULL = ceil(recordLenNULL*expParams.sr);
    time = 0:1/expParams.sr:(nSamples-1)/expParams.sr;
    recAudio = zeros(nSamples,1);       % initialize variable to store audio
    % % % % % % % % % % % % nMissingSamples = 0;                % cumulative n missing samples between frames
    beepDetected = 0;
    voiceOnsetDetected = 0;             % voice onset not yet detected
    frameCount = 1;                     % counter for # of frames (starting at first frame)
    endIdx = 0;                         % initialize idx for end of frame
    voiceOnsetState = [];
    beepOnsetState = [];
        
    % set up figure for real-time plotting of audio signal of next trial
    if show_mic_trace_figure
        figure(rtfig)
        set(micTitle,'string',sprintf('%s %s run %d trial %d condition: %s', expParams.subject, expParams.task, expParams.run, itrial));
    end

    
    %t = timer;
    %t.StartDelay = 0.050;   % Delay between timer start and timer function
    %t.TimerFcn = @(myTimerObj, thisEvent)play(beepPlayer); % Timer function plays GO signal
    setup(deviceReader) % note: moved this here to avoid delays in time-sensitive portion
    
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
    if FLAG_SEND_EVENT_STIM_ONSET
        try
            [evt_,evtCode_] = send_event([1],[],0.1,0.04,1,'Dev2');
            evt = cat(1,evt,evt_); evtCode = cat(1,evtCode,evtCode_);
            % ## Ideally, we should move the following to the end of trial ##
            tblEvt = table(evt,evtCode,'VariableNames',{'EventTime_dn','EventCode'});
            save(taskEvent_fname, 'tblEvt');
        end
    end
    % <<<

    ok=ManageTime('wait', CLOCK, TIME_STIM_START);
    TIME_STIM_ACTUALLYSTART = ManageTime('current', CLOCK);
    

        set(annoStr.Plus, 'Visible','off');
        set(annoStr.Stim, 'String', {trialData(itrial).display});
        set(annoStr.Stim, 'Visible','on');
        drawnow;

    if ~ok, fprintf('i am late for this trial TIME_STIM_START\n'); end

    TIME_SOUND_START = TIME_STIM_ACTUALLYSTART + trialData(itrial).timePreSound;
    %ok=ManageTime('wait', CLOCK, TIME_SOUND_START - stimoffset);
    ok=ManageTime('wait', CLOCK, TIME_SOUND_START);

    TIME_SOUND_ACTUALLYSTART = ManageTime('current', CLOCK);
    while ~isDone(stimread); sound=stimread();headwrite(sound);end;release(stimread);reset(headwrite);
    TIME_SOUND_END = TIME_SOUND_ACTUALLYSTART + trialData(itrial).timeStim;           % stimulus ends
    if ~ok, fprintf('i am late for this trial TIME_SOUND_START\n'); end

    TIME_ALLSTIM_END = TIME_SOUND_END + trialData(itrial).timePostSound;
    %TIME_ALLSTIM_END = TIME_SOUND_RELEASED + trialData(ii).timePostSound;
    ok=ManageTime('wait', CLOCK, TIME_ALLSTIM_END);
    if strcmp(expParams.visual, 'orthography')
        set(annoStr.Stim, 'Visible','off');
        set(annoStr.Plus, 'Visible','on');
        drawnow;
    end
    if ~ok, fprintf('i am late for this trial TIME_ALLSTIM_END\n'); end        

    TIME_GOSIGNAL_START = TIME_ALLSTIM_END + trialData(itrial).timePostStim;          % GO signal time
    if  show_mic_trace_figure
        set(micLine,'visible','off');set(micLineB,'visible','off');
    end
    
    drawnow;


    ok=ManageTime('wait', CLOCK, TIME_GOSIGNAL_START - beepoffset);     % waits for recorder initialization time
    [nill, nill] = deviceReader(); % note: this line may take some random initialization time to run; audio signal start (t=0) will be synchronized to the time when this line finishes running
    if ~ok, fprintf('i am late for this trial TIME_GOSIGNAL_START - beepoffset\n'); end
    
    ok=ManageTime('wait', CLOCK, TIME_GOSIGNAL_START);     % waits for GO signal time
    while ~isDone(beepread); sound=beepread();headwrite(sound);end;reset(beepread);reset(headwrite);
    %TIME_GOSIGNAL_RELEASED = ManageTime('current', CLOCK);
    %TIME_GOSIGNAL_ACTUALLYSTART = TIME_GOSIGNAL_RELEASED - beepdur; % actual time for GO signal 
    TIME_GOSIGNAL_ACTUALLYSTART = ManageTime('current', CLOCK); % actual time for GO signal 
    if strcmp(expParams.visual, 'fixpoint'),set(annoStr.Plus, 'color','g');drawnow;end
    if strcmp(expParams.visual, 'orthography'),set(annoStr.Plus, 'color','g');drawnow;end
    if ~ok, fprintf('i am late for this trial TIME_GOSIGNAL_START\n'); end
    TIME_VOICE_START = TIME_GOSIGNAL_ACTUALLYSTART + nonSpeechDelay;                   % expected voice onset time


    TIME_SCAN_START = TIME_GOSIGNAL_ACTUALLYSTART + trialData(itrial).timeMax;



    while endIdx < nSamples
        % find beginning/end indices of frame
        begIdx = (frameCount*expParams.frameLength)-(expParams.frameLength-1);
        endIdx = (frameCount*expParams.frameLength);

        % read audio data
        [audioFromDevice, numOverrun] = deviceReader();     % read one frame of audio data % note: audio t=0 corresponds to first call to deviceReader, NOT to time of setup(...)
        numOverrun = double(numOverrun);    % convert from uint32 to type double
        % % % % % % if numOverrun > 0, recAudio(begIdx:begIdx+numOverrun-1) = 0; end      % set missing samples to 0
        % % % % % % % % % recAudio(begIdx+numOverrun:endIdx+numOverrun) = audioFromDevice;    % save frame to audio vector
        % % % % % % % % % % % % % nMissingSamples = nMissingSamples + numOverrun;     % keep count of cumulative missng samples between frames

        % plot audio data
        drawnow()

        % voice onset exclusion
        minVoiceOnsetTime = max(0, expParams.minVoiceOnsetTime-(begIdx+numOverrun)/expParams.sr);
        
        % detect beep onset
        if beepDetected == 0 && expParams.minVoiceOnsetTime > (begIdx+numOverrun)/expParams.sr
            % look for beep onset
            [beepDetected, bTime, beepOnsetState]  = detectVoiceOnset(recAudio(begIdx+numOverrun:endIdx+numOverrun), expParams.sr, expParams.rmsThreshTimeOnset, rmsBeepThresh, 0, beepOnsetState);
            if beepDetected
                beepTime = bTime + (begIdx+numOverrun)/expParams.sr; 
                 if show_mic_trace_figure
                    set(micLineB,'value',beepTime,'visible','on');
                 end
            end
        elseif voiceOnsetDetected == 0,% && frameCount > onsetWindow/frameDur
            if ~beepDetected; 
                beepTime = 0; 
                % disp('Beep not detected. Assign beepTime = 0.'); % this warning needs to be moved - it prints endlessly until speech onset
            end
            trialData(itrial).beepTime = beepTime;

            % look for voice onset in previous onsetWindow
            [voiceOnsetDetected, voiceOnsetTime, voiceOnsetState]  = detectVoiceOnset(recAudio(begIdx+numOverrun:endIdx+numOverrun), expParams.sr, expParams.rmsThreshTimeOnset, rmsThresh, minVoiceOnsetTime, voiceOnsetState);


        end

        frameCount = frameCount+1;

    end
    if voiceOnsetDetected == 0, fprintf('warning: voice was expected but not detected (rmsThresh = %f)\n',rmsThresh); end
    release(deviceReader); % end recording
    
    switch expParams.visual
        case 'fixpoint'
            set(annoStr.Plus, 'color','w');
        case 'orthography'
            set(annoStr.Plus, 'color','w');
    end


    %% save voice onset time and determine how much time left before sending trigger to scanner
    if voiceOnsetDetected == 0 %if voice onset wasn't detected
        trialData(itrial).onsetDetected = 0;
        trialData(itrial).voiceOnsetTime = NaN;
        trialData(itrial).nonSpeechDelay = nonSpeechDelay;
    else
        trialData(itrial).onsetDetected = 1;
        trialData(itrial).voiceOnsetTime = voiceOnsetTime;
        trialData(itrial).nonSpeechDelay = NaN;
    end


    % % % % % % % % % % else
        TIME_SCAN_ACTUALLYSTART=nan;
        %TIME_TRIG_RELEASED = nan;
        TIME_SCAN_END = nan;
        NEXTTRIAL = TIME_SCAN_START + trialData(itrial).timePreStim;
    % % % % % % % % % % % % % % % % end
        
    trialData(itrial).timingTrial =    [TIME_TRIAL_START;TIME_STIM_START;TIME_STIM_ACTUALLYSTART;TIME_SOUND_START;TIME_SOUND_ACTUALLYSTART;TIME_SOUND_END;TIME_ALLSTIM_END;TIME_GOSIGNAL_START;TIME_GOSIGNAL_ACTUALLYSTART;TIME_VOICE_START];
    expParams.timingTrialNames = split('TIME_TRIAL_START;TIME_STIM_START;TIME_STIM_ACTUALLYSTART;TIME_SOUND_START;TIME_SOUND_ACTUALLYSTART;TIME_SOUND_END;TIME_ALLSTIM_END;TIME_GOSIGNAL_START;TIME_GOSIGNAL_ACTUALLYSTART;TIME_VOICE_START', ';');

    TIME_STIM_START = NEXTTRIAL; 



    %% save for each trial
    trialData(itrial).s = recAudio(1:nSamples);
    trialData(itrial).fs = expParams.sr;
    if voiceOnsetDetected, trialData(itrial).reference_time = voiceOnsetTime;
    else trialData(itrial).reference_time = nonSpeechDelay;
    end
    % % % % % % % % trialData(itrial).percMissingSamples = (nMissingSamples/(recordLen*expParams.sr))*100;

    %JT save update test 8/10/21
    % save only data from current trial
    tData = trialData(itrial);

    % fName_trial will be used for individual trial files (which will
    % live in the run folder)
    fName_trial = fullfile(taskpath,sprintf('sub-%s_ses-%d_run-%d_task-%s_trial-%d.mat',expParams.subject, expParams.session, expParams.run, expParams.task,itrial));
    save(fName_trial,'tData');

    if require_keypress_every_trial
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
expParams.elapsed_time = toc(ET)/60;    % elapsed time of the experiment
fprintf('\nElapsed Time: %f (min)\n', expParams.elapsed_time)
save(Output_name, 'expParams', 'trialData');

% number of trials with voice onset detected
onsetCount = nan(op.ntrials,1);
for j = 1: op.ntrials
    onsetCount(j) = trialData(j).onsetDetected;
end
numOnsetDetected = sum(onsetCount);    

fprintf('Voice onset detected on %d/%d trials\n', numOnsetDetected, op.ntrials);

fprintf('Press any key to send final sync pulses and end this experimental phase')
pause()
beacon_times = [beacon_times, test_Beacon(0.4,0.05,5)];
    save(beacon_times_fname,'beacon_times');

end




function [voiceOnsetDetected, voiceOnsetTime, state]  = detectVoiceOnset(samples, Fs, onDur, onThresh, minVoiceOnsetTime, state)
% function [voiceOnsetDetected, voiceOnsetTime]  = detectVoiceOnset(samples, Fs, onDur, onThresh, minVoiceOnsetTime)
% 
% Function to detect onset of speech production in an audio recording.
% Input samples can be from a whole recording or from individual frames.
%
% INPUTS        samples                 vector of recorded samples
%               Fs                      sampling frequency of samples
%               onDur                   how long the intensity must exceed the
%                                       threshold to be considered an onset (s)
%               onThresh                onset threshold
%               minVoiceOnsetTime       time (s) before which voice onset
%                                       cannot be detected (due to
%                                       anticipation errors, throat
%                                       clearing etc) - often set to .09 at
%                                       beginning of recording/first frame
%
% OUTPUTS       voiceOnsetDetected      0 = not detected, 1 = detected
%               voiceOnsetTime          time (s) when voice onset occurred
%                                       (with respect of first sample)
%
% Adapted from ACE study in Jan 2021 by Elaine Kearney (elaine-kearney.com)
% Matlab 2019a 
%
%%

% set up parameters
winSize = ceil(Fs * .002);                  % analysis window = # samples per 2 ms (matched to Audapter frameLen=32 samples @ 16KHz)
Incr = ceil(Fs * .002);                     % # samples to increment by (2 ms)
rmsFF = 0.90;                               % forgetting factor (Audapter's ma_rms computation)
BegWin = 1;                                 % first sample in analysis window
EndWin = BegWin + winSize -1;               % last sample in analysis window
voiceOnsetDetected = false;                 % voice onset (not yet detected)
voiceOnsetTime = [];                        % variable for storing voice onset time
if nargin<5||isempty(minVoiceOnsetTime), minVoiceOnsetTime = 0; end % minimum voice onset time
if nargin<6||isempty(state), state=struct('rms',0,'count',0); end   % rms: last-call rms value; count: last-call number of supra-threshold rms values

% main loop
while EndWin <= length(samples)
    
    dat = samples(BegWin:EndWin);
    %dat = detrend(dat, 0);                 % legacy step: removes mean value from data in analysis window
    %dat = convn(dat,[1;-.95]);             % legacy step: applies a high pass filter to the data
                                            % and may reduce sensitivity to production onset,
                                            % especially if stimulus starts with a voiceless consonant
    int = mean(dat.^2);                      % mean of squares
    state.rms =  rmsFF*state.rms+(1-rmsFF)*sqrt(int);            % calculate moving-average rms (calcRMS1 function)
    if state.rms > onThresh && BegWin > minVoiceOnsetTime*Fs
        state.count = state.count + 1;
    else
        state.count = 0;
    end
    
    % criteria for voice onset:
    % (1) onThresh must have been continuously exceeded for the duration of onDur
    % (2) time of voice onset is greater than minVoiceOnsetTime
    
    if state.count >= onDur*Fs/Incr
        voiceOnsetDetected = true;                              % onset detected
        voiceOnsetTime = (BegWin-(state.count-1)*Incr-1)/Fs;    % time when onThresh was first reached
        break
    end
    
    % increment analysis window and iteration by 1 (until voice onset detected)
    BegWin = BegWin + Incr;          
    EndWin = EndWin + Incr;                
end

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