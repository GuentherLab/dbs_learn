 % generate an annot table which translates stim computer trial times into Percept file times
 %%% first need to manually find and save landmarks

 %%%% to create landmarks file, we need to find the CED and stim-computer times of an audible event
 % easiest will usually be stim audio onset of trial 1, which in stim computer data is the 'TIME_SOUND_ACTUALLYSTART' value in tData
 % .... this will be the tData.timeingTrial value in the slot labeled 'TIME_VOICE_START' in exParams.timingTrialNames
 % to find this time CED, pick a start sample and run: 
 % ...... startsamp = $;  player = audioplayer(physdat.CED.voice.trial{1}(startsamp:startsamp+1e5),physdat.CED.voice.fsample); play(player); physdat.CED.voice.time{1}(startsamp)
 % plot the audio waveform in the sample that was just played:
 % ...... figure; plot(physdat.CED.voice.time{1}(startsamp:startsamp+1e5), physdat.CED.voice.trial{1}(startsamp:startsamp+1e5))
 % or save as WAV and open in Audacity
 %%%%% this method is somewhat suboptimal because we are assuming that TIME_VOICE_START is accurate (that the stim computer got a good estimate of voice onset)

clear

 op.sub = 'dbs001'; 
 op.ses = 1; 
 op.run = 1; 

 op.task = 'trainA'; op.task_dir = 'train-a'; 
  % op.task = 'trainB'; op.task_dir = 'train-b'; 

physio_preproc_file = 'train-a_preproc_240417-131548.mat';
% physio_preproc_file = 'train-b_preproc_240417-132735.mat';

op.ylimits = [];
op.ylimits = [3 11] * 1e4; %%% limits for DBS-ON mode - dbs001 train A
% op.ylimits = [0.5 3] * 1e5; %%% limits dbs001 train B for DBS-OFF mode

op.dbs_response_freq = 17; % hz
op.wavtrans_out_freq = 100; 
op.wavtrans_wav_width = 7; % number of wavelets

op.chan_ind = 5; % chan ind to plot - in future, analyze all chans

%%% remove datapoints greater than this value as artifacts..... suboptimal method - replace with standard deviation method
% op.dbs_hard_thresh = 5e5; % reasonable thresh for dbs001 trainA (DBS ON)
op.dbs_hard_thresh = 1.1e6; % reasonable thresh for dbs001 trainB (DBS OFF)


op.do_high_pass_filter = 'yes'; % should a high pass filter be applied
    op.high_pass_filter_freq = 1; %cutoff frequency of high pass filter, if one is to be used
op.do_bsfilter = 'no';  % should a band stop filter be applied
    op.line_noise_harm_freqs = [60 120 180 240]; % for notch filters for 60hz harmonics, if band stop filter is to be used


trial_window_peri_speech = [-7.5 8]; 


op.smooth_method = 'gaussian'; 
op.smooth_windowsize = 3; 

 op.normalize_by_pretest = 1; 
 op.exclude_no_response_trials = 1; 




setpaths_dbsmulti(); 

paths.sub = [paths.data, filesep, 'sub-', op.sub];
paths.ses = [paths.sub, filesep, 'ses-', num2str(op.ses)];
paths.physio = [paths.ses, filesep, 'physio']; 
paths.beh = [paths.ses, filesep, 'beh']; 
paths.beh_task = [paths.beh, filesep, op.task_dir]; 

landmarks_file = [paths.ses, filesep, ['sub-',op.sub,'_ses-',num2str(op.ses),'_task-',op.task,'_landmarks.tsv']]; 

physdat = load([paths.physio, filesep, physio_preproc_file]);



%%%%%%%%%%% translate Percept times into CED times
physdat = ft_checkconfig(physdat,'renamed',{'preprocEMG','CED'});
physdat = ft_checkconfig(physdat,'renamed',{'preprocLFP','PCP'});
ppCED = physdat.CED;
ppBHV = physdat.bhv;
ppPCP = physdat.PCP;
m_Align    = physdat.m_Align;
pcp = ppPCP.data;
t_PCP = pcp.time{1};
t_CED = ppCED.stimArt.time{1};
t_BHV = [];
% aligned times
fcn_translateTime = m_Align.fcn_translateTime; % function to apply the translation
tTrans_PCP_to_CED = m_Align.TimeTrans_PCP_to_CED_byEvent; % translating PCP time into CED time
tTrans_BHV_to_CED = m_Align.TimeTrans_BHV_to_CED; % translating BHV time into CED time
%
if ~isempty(tTrans_PCP_to_CED)
    t_PCP_inCED = fcn_translateTime(t_PCP,tTrans_PCP_to_CED);
else
    t_PCP_inCED = [];
end
pcp.oldtime = pcp.time;
pcp.time{1} = t_PCP_inCED;




%%%% load behavioral data
load_accuracy_results()


%%% dbs responses
trials_resp = subs{op.sub, ['tr_',op.task]}{1};
ntrials = height(trials_resp); 
resp.timecourse = cell(ntrials,1); 


%%%% find stimcomp-CED adjustment
ldmrk_tab = readtable(landmarks_file,'FileType','text','Delimiter','tab');
ldmrk_tab.Properties.RowNames = ldmrk_tab.computer;
ced_minus_stimcomp_time = ldmrk_tab{'CED','time'} - ldmrk_tab{'stim','time'}; 


%% get dbs responses at freq of interest

% Applying High Pass Filter and line noise filter
cfg=[];
cfg.hpfilter=op.do_high_pass_filter;
    cfg.hpfreq=op.high_pass_filter_freq;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.bsfilter = op.do_bsfilter;
    cfg.bsfreq= [op.line_noise_harm_freqs-1; op.line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
D_hpf = ft_preprocessing(cfg,pcp);

% get wave power
cfg = [];
cfg.out_freq = op.wavtrans_out_freq;
cfg.wav_freq = op.dbs_response_freq; 
cfg.wav_width = op.wavtrans_wav_width;
D_wavpow = bml_envelope_wavpow(cfg,D_hpf);


% artifact rejection.... need to use method with more systematically determined thresh
artf_inds = D_wavpow.trial{1} > op.dbs_hard_thresh;
D_wavpow.trial{1}(artf_inds) = nan; 

%% get dbs response at each trial
for itrial = 1:ntrials
    trialnum = trials_resp.trial(itrial); %%% may be different from itrial if some trials are excluded
    load([paths.beh_task, filesep, 'sub-',op.sub,'_ses-',num2str(op.ses),'_run-',num2str(op.run),'_task-',op.task_dir,'_trial-',num2str(trialnum),'.mat'])

    % get key trial times in the CED time coordinates
    %%%% see expParams for timepoint labels
    trials_resp.audio_stim_on(itrial) = tData.timingTrial(5) + ced_minus_stimcomp_time; 
    trials_resp.audio_stim_off(itrial) = tData.timingTrial(6) + ced_minus_stimcomp_time; 
    trials_resp.gobeep_start(itrial) = tData.timingTrial(9) + ced_minus_stimcomp_time; 
    trials_resp.sp_on(itrial) = tData.timingTrial(10) + ced_minus_stimcomp_time; % speech on time

    % speech-aligned responses
    time_window_this_trial = trial_window_peri_speech + trials_resp.sp_on(itrial); %%% CED time coordinates
    time_inds_this_trial = D_wavpow.time{1} > time_window_this_trial(1) & D_wavpow.time{1} < time_window_this_trial(2); 
    trials_resp.timecourse{itrial} = D_wavpow.trial{1}(op.chan_ind,time_inds_this_trial); 


end

%% plotting

%%%%% need to baseline normalize

op.smooth_windowsize = 50;


% trials_to_plot = 1:33;
% trials_to_plot = 34:66;
% trials_to_plot = 67:ntrials;
trials_to_plot = 1:50;
% trials_to_plot = 51:ntrials;

% trials_to_plot = 1:ntrials; 


% close all

% hfig = figure;  

mean_stimon_to_speech = mean(trials_resp.sp_on - trials_resp.audio_stim_on); 
mean_stimoff_to_speech = mean(trials_resp.sp_on - trials_resp.audio_stim_off); 
mean_rt = mean(trials_resp.sp_on - trials_resp.gobeep_start); 

timecourses = trials_resp.timecourse;
min_trial_length_samples = min(cellfun(@length,timecourses)); 
tc_truncated = cellfun(@(x)x(1:min_trial_length_samples),timecourses,'UniformOutput',false); 
timecourses = cell2mat(tc_truncated);

timecourses_to_plot = timecourses(trials_to_plot,:); 
timecourses_to_plot = smoothdata(timecourses_to_plot, 2, op.smooth_method, op.smooth_windowsize); 

tc_mean = nanmean(timecourses_to_plot,1);
tc_sem = std(timecourses_to_plot,1) / sqrt(size(timecourses_to_plot,1)); 
tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem]; 

xtime = linspace(trial_window_peri_speech(1), trial_window_peri_speech(2), min_trial_length_samples); 


   % hold off
    hfill = fill([xtime, fliplr(xtime)], [tc_sem_lims(1,:), fliplr(tc_sem_lims(2,:))], [0.8 0.8 0.8],...
        'HandleVisibility','off'); % standard error
%         hfill.LineStyle = 'none'; % no border
        hfill.EdgeColor = [0.8 0.8 0.8]; 
    hold on 
    hplot = plot(xtime, tc_mean);
        hplot.LineWidth = 1;

 xline_stimon = xline(-mean(mean_stimon_to_speech),'k','HandleVisibility','off');
 xline_stimoff = xline(-mean(mean_stimoff_to_speech),'k','HandleVisibility','off');
 xline_gobeep_on = xline(-mean(mean_rt),'g','HandleVisibility','off');
 xline_sp_on = xline(0,'LineStyle','--','Color',[0.5 0.5 0.5],'HandleVisibility','off');

ylabel('beta power')
xlabel('time after speech onset (sec)')

if ~isempty(op.ylimits)
     ylim(op.ylimits)
end

