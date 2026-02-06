 %%% generate trialwise table of response power at a specific freq for a specific channel 
 % 'pcp' is the fieldtrip struct containing Percept data

function trials_resp = trialwise_resp(pcp, subs, op)



setpaths_dbsmulti()
paths.sub = [paths.data, filesep, 'sub-', op.sub];
paths.ses = [paths.sub, filesep, 'ses-', num2str(op.ses)];
paths.physio = [paths.ses, filesep, 'physio']; 
paths.beh = [paths.ses, filesep, 'beh']; 
paths.beh_task = [paths.beh, filesep, task2dir(op.task)]; 

landmarks_file = [paths.annot, filesep, ['sub-',op.sub,'_ses-',num2str(op.ses),'_task-',op.task,'_landmarks.tsv']]; 


ldmrk_tab = readtable(landmarks_file,'FileType','text','Delimiter','tab');
ldmrk_tab.Properties.RowNames = ldmrk_tab.computer;
ced_minus_stimcomp_time = ldmrk_tab{'CED','time'} - ldmrk_tab{'stim','time'}; 

trials_resp = subs{op.sub, ['tr_',op.task]}{1};
ntrials = height(trials_resp); 
trials_resp.timecourse = cell(ntrials,1); 
    

%%%% replace nans with zero and note locations
% count timepoints in which data for any channel is missing as missing for all channels
pcp.missingdat = any(isnan(pcp.trial{1}),1);  % record locations so we can exclude later
pcp.trial{1}(:,pcp.missingdat) = 0; % replace w/ zero so we can run wavelet transform

cfg = []; 
cfg.channel = op.chan
D_sel = ft_selectdata(cfg,pcp); 

% Applying High Pass Filter and line noise filter
cfg=[];
cfg.hpfilter=op.do_high_pass_filter;
    cfg.hpfreq=op.high_pass_filter_freq;
cfg.hpfilttype='but';
cfg.hpfiltord=5;
cfg.hpfiltdir='twopass';
cfg.bsfilter = op.do_bsfilter;
    cfg.bsfreq= [op.line_noise_harm_freqs-1; op.line_noise_harm_freqs+1]'; % notch filters for 60hz harmonics
D_hpf = ft_preprocessing(cfg,D_sel);

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
    time_window_this_trial = op.trial_window_peri_speech + trials_resp.sp_on(itrial); %%% CED time coordinates
    time_inds_this_trial = D_wavpow.time{1} > time_window_this_trial(1) & D_wavpow.time{1} < time_window_this_trial(2); 
    trials_resp.timecourse{itrial} = D_wavpow.trial{1}(1,time_inds_this_trial); 

    % baseline normalization
    baseline_inds_this_trial = D_wavpow.time{1} > time_window_this_trial(1) & D_wavpow.time{1} < trials_resp.audio_stim_on(itrial); 
    trials_resp.base(itrial) = mean(D_wavpow.trial{1}(1,baseline_inds_this_trial));

    if op.norm_subtract_baseline
        trials_resp.timecourse{itrial} = trials_resp.timecourse{itrial} - trials_resp.base(itrial);
    end

    if op.norm_divide_abs_baseline 
        trials_resp.timecourse{itrial} = trials_resp.timecourse{itrial} ./ abs(trials_resp.base(itrial)); 
    end

end

% trial-relative times
trials_resp.rt = trials_resp.sp_on - trials_resp.gobeep_start; 
trials_resp.stimon_to_speech = trials_resp.sp_on - trials_resp.audio_stim_on; 
trials_resp.stimoff_to_speech = trials_resp.sp_on - trials_resp.audio_stim_off; 

% % truncate the end of each trial so that all trials are the same length in samples
min_trial_length = min(cellfun(@length,trials_resp.timecourse)); 
trials_resp.timecourse = cellfun(@(x)x(1:min_trial_length),trials_resp.timecourse,'UniformOutput',false); 
