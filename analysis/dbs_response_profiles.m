 %%% for each channel in each subject, get response profiles during each phase

clear
setpaths_dbsmulti(); 

 %% params

 op.ses = 1; 
 op.run = 1; 

 op.tasks_to_analyze = {'trainA','trainB','test'}; 
op.savepath = [paths.groupanalysis, filesep, 'resp_all']; 

%%%%%%%%%%% behavioral analysis params
op.exclude_no_response_trials = 1; % if true, exclude trials in which subject did not produce spoken response; otherwise these trials will be counted as having zero accuracy
op.normalize_by_pretest = 1; % subtract pretest accuracy from later accuracy... only matter for beh. analysis, not response profiles
op.n_early_trials_per_stimid = 5; % trials when the majority of learning within a phase is expected to occur


%%%%%%%%%%%%% dbs analysis params
% wavelet transform params
 op.dbs_response_freq = 17; % hz
op.wavtrans_out_freq = 125; % must be integral factor of sample rate (250hz)
op.wavtrans_wav_width = 7; % number of wavelets

% filtering params
op.do_high_pass_filter = 'yes'; % should a high pass filter be applied
    op.high_pass_filter_freq = 1; %cutoff frequency of high pass filter, if one is to be used
op.do_bsfilter = 'no';  % should a band stop filter be applied
    op.line_noise_harm_freqs = [60 120 180 240]; % for notch filters for 60hz harmonics, if band stop filter is to be used

% remove datapoints with wavpow [not raw] greater than this value as artifacts..... suboptimal method - replace with standard deviation method
op.dbs_hard_thresh = 8e5; 

op.trial_window_peri_speech = [-7.5 8];  % window in seconds around speech onset

op.norm_subtract_baseline = 1; % subtract trialwise baseline data from response values
op.norm_divide_abs_baseline = 1; % divide each trial timecourse by abs val of baseline; do so after subtract_baseline if that op is true

%% get response timecourses 

%%%% load behavioral data
load_accuracy_results()

ntasks = length(op.tasks_to_analyze); 

resp = table({},{},'VariableNames',{'sub', 'chan'});

warning('off','MATLAB:table:RowsAddedExistingVars') % we will be building the table 1 row at a time

for isub = 1:nsubs
    thissub = subs.sub{isub}

    paths.sub = [paths.data, filesep, 'sub-', thissub];
    paths.ses = [thissub, filesep, 'ses-', num2str(op.ses)];
    paths.physio = [paths.ses, filesep, 'physio']; 
    paths.beh = [paths.ses, filesep, 'beh']; 

    for itask = 1:ntasks
        thistask = op.tasks_to_analyze{itask}; 
        task_dir = task2dir(thistask);
        paths.beh_task = [paths.beh, filesep, task_dir]; 
        

        % load dbs data for this task in this sub
        cfg = op; 
        cfg.sub = thissub; 
        cfg.task = thistask; 
        [physdat, pcp] = load_ced_and_percept_data(cfg);

        thesechans = pcp.label;
        nchans = length(thesechans);
        for ichan = 1:nchans
            thischan = thesechans{ichan}; 

            resprow = find(resp.sub==string(thissub) & resp.chan==string(thischan));
            if isempty(resprow) % if we don't already have a row for this chan
                resprow = height(resp)+1; 
                resp.sub{resprow} = thissub; % add row
                resp.chan{resprow} = thischan; 
            end

            cfg = op;
            cfg.sub = thissub; 
            cfg.task = thistask;             
            cfg.chan = thischan; 
            trialtab = trialwise_resp(pcp, subs, cfg); 

            taskresp = struct;
            taskresp.trials = trialtab; 

            resp{resprow,thistask} = {taskresp}; 
        end
    end
end

%% sort responses by stimulus
for isub = 1:nsubs
    thissub = subs.sub{isub};
    for itask = 1:ntasks
        thistask = op.tasks_to_analyze{itask}; 
        thesechans = unique(resp.chan(resp.sub==string(thissub))); 
        nchans = length(thesechans);
        for ichan = 1:nchans
            thischan = thesechans{ichan}; 
            resprow = find(resp.sub==string(thissub) & resp.chan==string(thischan));
            if ~isempty(resp{resprow,thistask}{1}) % if we have recordings from this chan in this task
                trials_resp = resp{resprow,thistask}{1}.trials; 
                stimlist = unique(trials_resp.name);
                nstim = length(stimlist); 
                resp_stim = table(stimlist,cell(nstim,1),cell(nstim,1),cell(nstim,1),'VariableNames',{'name','con','stim_group','trials'},'RowNames',stimlist);
                for iname = 1:nstim % for each unique stim in this task
                    thisstim = stimlist{iname};
                    matchtrialrows = find(trials_resp.name == string(thisstim));
                    resp_stim.stim_group{iname} = trials_resp.stim_group{matchtrialrows(1)}; 
                    if string(resp_stim.stim_group{iname}) == "novel";
                        resp_stim.con{iname} = 'novel';
                    else
                        resp_stim.con{iname} = get_dbs_state(thissub, resp_stim.stim_group{iname}, paths); 
                    end
                    resp_stim.trials{iname} = trials_resp(matchtrialrows,:); 
                end

                % get early test phase stim sorted by condition
                %%%% to be more accurate about trial timing averages, we could exclude trials when this channel doesn't have data [see fieldtrip struct]
                if string(thistask)=="test";
                    conlist = unique(trials_resp.traincon);
                    ncons = length(conlist);
                    resp_con = table(conlist,cell(ncons,1),cell(ncons,1),cell(ncons,1),cell(ncons,1),cell(ncons,1),cell(ncons,1),'VariableNames',{'name','stim_group','trials','resp_mean','resp_sem','stimresp_mean','stimresp_sem'},'RowNames',conlist);
                    for icon = 1:ncons
                        thiscon = conlist{icon};
                        matchrows = string(resp_stim.con) == thiscon;
                        early_stim_tables_this_con = cellfun(@(x)x(1:op.n_early_trials_per_stimid,:),resp_stim.trials(matchrows),'UniformOutput',0); 
                        resp_con.trials{icon} = cat(1,early_stim_tables_this_con{:}); % stack trials from early instance of stim in this training condition
                        respcontrials = cell2mat(resp_con.trials{icon}.timecourse);
                        resp_con.resp_mean{icon} = nanmean(respcontrials,1); 
                        resp_con.resp_sem{icon} = nanstd(respcontrials,1) / sqrt(size(respcontrials,1)); 

                        %%% mean trial timing
                        resp_con.rt_mean(icon) = mean(resp_con.trials{icon}.rt);
                        resp_con.stimon_to_speech_mean(icon) = mean(resp_con.trials{icon}.stimon_to_speech);
                        resp_con.stimoff_to_speech_mean(icon) = mean(resp_con.trials{icon}.stimoff_to_speech);

                        %  response during stim
                        xtime = linspace(op.trial_window_peri_speech(1),op.trial_window_peri_speech(2), size(resp_con.trials{icon}.timecourse{1},2));  
                        stiminds = xtime > -resp_con{icon,'stimon_to_speech_mean'} & xtime < -resp_con{icon,'stimoff_to_speech_mean'};
                        % % % % % % % % % % % % % % % stimresp_not_normed = mean(respcontrials(:,stiminds),2); 
                        % % % % % % % % % % % % % % % % resp_con{icon,'trials'}{1}.stimresp = stimresp_not_normed ./ abs(resp_con.trials{icon}.base);
                        resp_con{icon,'trials'}{1}.stimresp = mean(respcontrials(:,stiminds),2); 
                        resp_con{icon,'stimresp'} = nanmean(resp_con{icon,'trials'}{1}.stimresp); 
                        

                    end
                    resp{resprow,thistask}{1}.con = resp_con; 
   
                end
                resp{resprow,thistask}{1}.rt_mean = mean(trials_resp.rt);
                resp{resprow,thistask}{1}.stimon_to_speech_mean = mean(trials_resp.stimon_to_speech);
                resp{resprow,thistask}{1}.stimoff_to_speech_mean = mean(trials_resp.stimoff_to_speech);
                resp{resprow,thistask}{1}.stim = resp_stim;
            end
            
        end
    end
end

%

% save results
mkdir(paths.groupanalysis)
save(op.savepath, 'resp', 'subs', 'op')






