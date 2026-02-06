%%%% analyze DBS responses in a single phase for one subject at a time

 %%%% to create landmarks file, we need to find the CED and stim-computer times of an audible event
 % easiest will usually be stim audio onset of trial 1, which in stim computer data is the 'TIME_SOUND_ACTUALLYSTART' value in tData
 % .... this will be the tData.timeingTrial value in the slot labeled 'TIME_VOICE_START' in exParams.timingTrialNames [probably slot 10]
 % to find this time CED, pick a start sample and run: 
 % ...... startsamp = $;  player = audioplayer(physdat.CED.voice.trial{1}(startsamp:startsamp+1e5),physdat.CED.voice.fsample); play(player); physdat.CED.voice.time{1}(startsamp)
 % plot the audio waveform in the sample that was just played:
 % ...... figure; plot(physdat.CED.voice.time{1}(startsamp:startsamp+1e5), physdat.CED.voice.trial{1}(startsamp:startsamp+1e5))
 % or save as WAV and open in Audacity
 %%%%% this method is somewhat suboptimal because we are assuming that TIME_VOICE_START is accurate (that the stim computer got a good estimate of voice onset)

clear

% op.sub = 'dbs001'; 
% op.sub = 'dbs002'; 
op.sub = 'dbs003'; 
% op.sub = 'dbs004'; 

 op.ses = 1; 
 op.run = 1; 

 % op.task = 'trainA'; 
  % op.task = 'trainB';
op.task = 'test'; 

op.ylimits = [];
% op.ylimits = [3 11] * 1e4; %%% limits for DBS-ON mode - dbs001 train A
% op.ylimits = [0.5 3] * 1e5; %%% limits dbs001 train B for DBS-OFF mode

op.dbs_response_freq = 17; % hz
op.wavtrans_out_freq = 125; % must be integral factor of sample rate (250hz)
op.wavtrans_wav_width = 7; % number of wavelets

op.chan_ind = 1; % chan ind to plot - in future, analyze all chans

%%% remove datapoints with wavpow [not raw] greater than this value as artifacts..... suboptimal method - replace with standard deviation method
% op.dbs_hard_thresh = 5e5; % reasonable thresh for dbs001 trainA (DBS ON)
% op.dbs_hard_thresh = 1.1e6; % reasonable thresh for dbs001 trainB (DBS OFF)
op.dbs_hard_thresh = 6e5; % reasonable thresh for dbs004 trainB test (DBS ON)


op.do_high_pass_filter = 'yes'; % should a high pass filter be applied
    op.high_pass_filter_freq = 1; %cutoff frequency of high pass filter, if one is to be used
op.do_bsfilter = 'no';  % should a band stop filter be applied
    op.line_noise_harm_freqs = [60 120 180 240]; % for notch filters for 60hz harmonics, if band stop filter is to be used


trial_window_peri_speech = [-7.5 8]; 

op.smooth_method = 'gaussian'; 

op.smooth_windowsize = 200; 
% op.smooth_windowsize = 1; 

 op.normalize_by_pretest = 1; 
 op.exclude_no_response_trials = 1; 


%% load data

setpaths_dbsmulti(); 

paths.sub = [paths.data, filesep, 'sub-', op.sub];
paths.ses = [paths.sub, filesep, 'ses-', num2str(op.ses)];
paths.physio = [paths.ses, filesep, 'physio']; 
paths.beh = [paths.ses, filesep, 'beh']; 
paths.beh_task = [paths.beh, filesep, task2dir(op.task)]; 

landmarks_file = [paths.annot, filesep, ['sub-',op.sub,'_ses-',num2str(op.ses),'_task-',op.task,'_landmarks.tsv']]; 


%%%% load behavioral data
load_accuracy_results()

%%% dbs responses
[physdat, pcp] = load_ced_and_percept_data(op)

trials_resp = subs{op.sub, ['tr_',op.task]}{1};
ntrials = height(trials_resp); 
trials_resp.timecourse = cell(ntrials,1); 


%%%% find stimcomp-CED adjustment
ldmrk_tab = readtable(landmarks_file,'FileType','text','Delimiter','tab');
ldmrk_tab.Properties.RowNames = ldmrk_tab.computer;
ced_minus_stimcomp_time = ldmrk_tab{'CED','time'} - ldmrk_tab{'stim','time'}; 


%% get dbs responses at freq of interest

%%%% replace nans with zero and note locations
% count timepoints in which data for any channel is missing as missing for all channels
pcp.missingdat = any(isnan(pcp.trial{1}),1);  % record locations so we can exclude later
pcp.trial{1}(:,pcp.missingdat) = 0; % replace w/ zero so we can run wavelet transform

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

% sometimes (if there is only 1 channel), bml_envelope_wavpow transposes 'trial' field in output
%%%% if so, transpose it back
for iblock = 1:length(D_wavpow.trial)
    if size(D_wavpow.trial{iblock},2) ~= size(D_wavpow.time{iblock},2) &&...
       size(D_wavpow.trial{iblock},2) == size(D_wavpow.time{iblock},1) 
            D_wavpow.trial{iblock} = D_wavpow.trial{iblock}'; 
    end
end

% mask missing data with nans
downsample_points = round(linspace(1, size(pcp.missingdat,2), size(D_wavpow.trial{1},2)));
D_wavpow.missingdat = pcp.missingdat(downsample_points); 
D_wavpow.trial{1}(:,D_wavpow.missingdat) = nan; % restore nans

% artifact rejection.... need to use method with more systematically determined thresh
artf_inds = D_wavpow.trial{1} > op.dbs_hard_thresh;
D_wavpow.trial{1}(artf_inds) = nan; 

%% get dbs response at each trial
for itrial = 1:ntrials
    trialnum = trials_resp.trial(itrial); %%% may be different from itrial if some trials are excluded
    load([paths.beh_task, filesep, 'sub-',op.sub,'_ses-',num2str(op.ses),'_run-',num2str(op.run),'_task-',task2dir(op.task),'_trial-',num2str(trialnum),'.mat'])

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


% trials_to_plot = 1:33;
% trials_to_plot = 34:66;
% trials_to_plot = 67:ntrials;
% trials_to_plot = 1:50;
% trials_to_plot = 51:ntrials;
trials_to_plot = 1:ntrials; 


% close all

hfig = figure('Color',[1 1 1]);  

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

box off

title([op.sub, ' ', D_wavpow.label{op.chan_ind}, ' ', op.task, ' (DBS-', upper(get_dbs_state(op.sub, op.task, paths)), ')'])

if ~isempty(op.ylimits)
     ylim(op.ylimits)
end

