  %%%% plot responses from early trials in test phase sorted by training condition

 setpaths_dbsmulti()

 % load([paths.groupanalysis, filesep, 'resp_all']) %%%% data created by dbs_response_profiles.m
   load([paths.groupanalysis, filesep, 'resp_all_5early']) %%%% data created by dbs_response_profiles.m
  % load([paths.groupanalysis, filesep, 'resp_all_10early']) %%%% data created by dbs_response_profiles.m

close all


%%
op.smooth_method = 'gaussian'; 

op.smooth_windowsize = 250; 
% op.smooth_windowsize = 1; 

op.plot_test_individual_chans = 0; 
op.plot_test_avg_all_chans = 1; 
op.plot_trainphase_first_vs_last_half = 1; 

op.dbs_on_off_colors = [0.6 0.6 0;... 
                        0 0 0.7];

conlist = {'on','off'};
% conlist = {'on','off','novel'};


op.sem_from_all_trials = 1; 

% approximate average speech duration in seconds to use
op.sp_dur_aprx_train_fh_on = 2.63; % first half dbs on
op.sp_dur_aprx_train_fh_off = 2.74; % first half dbs off
op.sp_dur_aprx_train_lh_on = 2.72; % last half dbs on
op.sp_dur_aprx_train_lh_off = 2.81; % last half dbs off
op.sp_dur_aprx_test_early = 2.575; % early test phase trials

set(0,'defaulttextinterpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare data
 % keep only channels with recordings during test phase
 resp = resp(~cellfun(@isempty,resp.test),:);

 ncons = length(conlist);
nchans = height(resp); 
nsubs = height(subs);

allchans_stimon_to_speech = [];
allchans_stimoff_to_speech = [];
allchans_rt = [];

onoff_test = table(repmat({nan(nchans,1)},2,1),nan(2,1),nan(2,1),'VariableNames',{'stimresp','stimresp_mean','stimresp_sem'},'RowNames',{'on','off'}); % for compare on vs. off in training

dbs_state_list = {'on','off'}

for ichan = 1:nchans
    allchans_stimon_to_speech = [allchans_stimon_to_speech; mean(resp.test{ichan}.trials.stimon_to_speech)]; 
    allchans_stimoff_to_speech = [allchans_stimoff_to_speech; mean(resp.test{ichan}.trials.stimoff_to_speech)]; 
    allchans_rt = [allchans_rt; mean(resp.test{ichan}.trials.rt)]; 

   for istate = 1:2
       dbs_state = dbs_state_list{istate};
       onoff_test{dbs_state,'stimresp'}{1}(ichan) = resp{ichan,'test'}{1}.con{dbs_state,'stimresp'};
   end
    
end

for istate = 1:2
    dbs_state = dbs_state_list{istate};
    onoff_test.stimresp_mean(istate) = mean(onoff_test{dbs_state,'stimresp'}{1});
    onoff_test.stimresp_sem(istate) = std(onoff_test{dbs_state,'stimresp'}{1}) / sqrt(nchans);
    onoff_test.stimresp_std(istate) = std(onoff_test{dbs_state,'stimresp'}{1});
end

plotcolors = table(op.dbs_on_off_colors,'VariableNames',{'color'},'RowNames',{'on','off'});

%% test phase, early trials - plot individual channels

for ichan = 1:nchans
      thischan = resp.chan{ichan}; 
    thissub = resp.sub{ichan}; 
    contab = resp.test{ichan}.con;       
    for icon = 1:ncons
        thiscon = conlist{icon}; 
     tc_mean = contab{thiscon,'resp_mean'}{1}; 
                resp{ichan,['mean_',thiscon]} = tc_mean;

    end
end

if op.plot_test_individual_chans

    hfig = figure('Color','w'); 
    sgtitle ('Beta responses by training condition (Test phase)')
    
    for ichan = 1:nchans
        thischan = resp.chan{ichan}; 
        thissub = resp.sub{ichan}; 
        contab = resp.test{ichan}.con; 
        subplot(2,4,ichan)
        hold on
        for icon = 1:ncons
            thiscon = conlist{icon}; 
           
    
            plot( smoothdata(tc_mean, 2, op.smooth_method, op.smooth_windowsize) )
            
        end
        legend(conlist)
        title([thissub, ' ', thischan])
        ylabel('beta power')
    
    
        box off
    end

end

%% test phase - average over subjects
if op.plot_test_avg_all_chans 
    hfig = figure('Color','w'); 
    
    for icon = 1:ncons
            thiscon = conlist{icon}; 
    
            % NB - might need to truncate some of these timecourses so they're all the same length
            tc_allchans = resp.(['mean_',thiscon]);
            tc_mean = nanmean(tc_allchans);
    
            if op.sem_from_all_trials
                tc_sem = nanstd(tc_allchans) / sqrt(nchans*2*op.n_early_trials_per_stimid);
            else
                tc_sem = nanstd(tc_allchans) / sqrt(nchans);
            end
            tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem]; 
            
    
            xtime = [1:size(tc_mean,2)] ./ op.wavtrans_out_freq + op.trial_window_peri_speech(1); 
    
            tc_sem_lims_to_plot = smoothdata(tc_sem_lims, 2, op.smooth_method, op.smooth_windowsize);
    
            hold on
            hfill = fill([xtime, fliplr(xtime)], [tc_sem_lims_to_plot(1,:), fliplr(tc_sem_lims_to_plot(2,:))], [0.8 0.8 0.8],...
                'HandleVisibility','off'); % standard error
    %         hfill.LineStyle = 'none'; % no border
            hfill.EdgeColor = [0.8 0.8 0.8]; 
    
            ylabel('beta power')
            
    end
    
    for icon = 1:ncons
            thiscon = conlist{icon}; 
    
            % NB - might need to truncate some of these timecourses so they're all the same length
            tc_allchans = resp.(['mean_',thiscon]);
            tc_mean = nanmean(tc_allchans);
    
            xtime = [1:size(tc_mean,2)] ./ op.wavtrans_out_freq + op.trial_window_peri_speech(1); 
    
    
            hplot = plot(xtime, smoothdata(tc_mean, 2, op.smooth_method, op.smooth_windowsize) );
            hplot.LineWidth = 2; 
            
    end
    
     % xline_stimon = xline(-mean(allchans_stimon_to_speech),'k','HandleVisibility','off');
    xline_stimon = xline(-allchans_stimon_to_speech,'k','LineStyle', ":",'HandleVisibility','off');
    
     xline_stimoff = xline(-mean(allchans_stimoff_to_speech),'k','HandleVisibility','off');
      xline_gobeep_on = xline(-mean(allchans_rt),'g','HandleVisibility','off');
     xline_sp_on = xline(0,'LineStyle','--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    
    legend(conlist)
    title(['All channels average - test phase'])
    
    
    box off
end

%% early vs. late  training




timecourse_linecolor = [0 0 0]; 
linestyles_first_last = {'-','--'}; 

op.trainphase_halfway_trial = 50; % last half starts after this trial

nan2 = nan(2,1);
trialtiming = table(nan2,nan2,nan2,'VariableNames',{'stimon_to_speech','stimoff_to_speech','rt'},   'RowNames',{'firsthalf','lasthalf'});

onoff_train = table('RowNames',{'on','off'}); % for compare on vs. off in training

xtime = linspace(op.trial_window_peri_speech(1),op.trial_window_peri_speech(2),size(resp.trainA{1}.trials.timecourse{1},2)); 



    if op.plot_trainphase_first_vs_last_half
        % close all
        hfig = figure('Color','w'); 
    end

for istate = 1:2
    dbs_state = dbs_state_list{istate};

    % in the following arrays, first value = timing in first half of trials, second value = last half
    allchans_stimon_to_speech = nan(nsubs,2); 
    allchans_stimoff_to_speech = nan(nsubs,2); 
    allchans_rt = nan(nsubs,2); 


   resptab = resp(:,{'sub','chan'});
    resptab.tc = cell(nchans,1); 
    fhlh = table(repmat({resptab},2,1),'VariableNames',{'resp'},'RowNames',{'fh','lh'}); 


    for ichan = 1:nchans
        thischan = resp.chan{ichan}; 
        thissub = resp.sub{ichan}; 
        this_task = dbs_state_to_trainphase(thissub, dbs_state, paths); 
    
        fhidx = 1:op.trainphase_halfway_trial;
        lhidx = op.trainphase_halfway_trial+1: length(resp{ichan,this_task}{1}.trials.timecourse); 
    
        fhlh{'fh','resp'}{1}{ichan,'tc'} = {cell2mat(resp{ichan,this_task}{1}.trials.timecourse(fhidx))};
        fhlh{'lh','resp'}{1}{ichan,'tc'} = {cell2mat(resp{ichan,this_task}{1}.trials.timecourse(lhidx))};

        fhlh{'fh','resp'}{1}{ichan,'base'} = {(resp{ichan,this_task}{1}.trials.base(fhidx))};
        fhlh{'lh','resp'}{1}{ichan,'base'} = {(resp{ichan,this_task}{1}.trials.base(lhidx))};

        % behavioral timing
        fh_trials = resp{ichan,this_task}{1}.trials(fhidx,:);
        lh_trials = resp{ichan,this_task}{1}.trials(lhidx,:);
        fhlh{'fh','resp'}{1}{ichan,'stimon_to_speech'} = mean(fh_trials.stimon_to_speech); 
        fhlh{'lh','resp'}{1}{ichan,'stimon_to_speech'} = mean(lh_trials.stimon_to_speech); 
        fhlh{'fh','resp'}{1}{ichan,'stimoff_to_speech'} = mean(fh_trials.stimoff_to_speech); 
        fhlh{'lh','resp'}{1}{ichan,'stimoff_to_speech'} = mean(lh_trials.stimoff_to_speech); 
        fhlh{'fh','resp'}{1}{ichan,'rt'} = mean(fh_trials.rt); 
        fhlh{'lh','resp'}{1}{ichan,'rt'} = mean(lh_trials.rt); 

        % response during stim
        for ihalf = 1:2
            stiminds = xtime > -fhlh{ihalf,'resp'}{1}.stimon_to_speech(ichan) & xtime < -fhlh{ihalf,'resp'}{1}.stimoff_to_speech(ichan);
             % mean resp of this chan during stim in this half
             %%%% normalize by baseline magnitude 
            fhlh{ihalf,'resp'}{1}.stimresp(ichan) = nanmean(nanmean(fhlh{ihalf,'resp'}{1}.tc{ichan}(:,stiminds),2) ./ abs(fhlh{ihalf,'resp'}{1}{ichan,'base'}{1}));
        end
    
    end

    for ihalf = 1:2
        fhlh{ihalf,'mean_tc'} = {nanmean(cell2mat(fhlh{ihalf,'resp'}{1}.tc))};
    end
   
    
    mean_firsthalf_plot = smoothdata(fhlh{'fh','mean_tc'}{1}, 2, op.smooth_method, op.smooth_windowsize);
    mean_lastthalf_plot = smoothdata(fhlh{'lh','mean_tc'}{1}, 2, op.smooth_method, op.smooth_windowsize);
    
    if op.plot_trainphase_first_vs_last_half
        subplot(1,2,istate)

        %%%%% hfill error bars here


        hold on
        hplot = plot(xtime,mean_firsthalf_plot,'Color',timecourse_linecolor,'LineWidth',2,'LineStyle',linestyles_first_last{1});
        hplot2 = plot(xtime,mean_lastthalf_plot,'Color',timecourse_linecolor,'LineWidth',2,'LineStyle',linestyles_first_last{2});




        ylabel('Beta power')
        
         xline_stimoff = xline(-mean(fhlh{'fh','resp'}{1}.stimoff_to_speech,1),'k','LineStyle',linestyles_first_last{1},'HandleVisibility','off');
         xline_stimoff = xline(-mean(fhlh{'lh','resp'}{1}.stimoff_to_speech,1),'k','LineStyle',linestyles_first_last{2},'HandleVisibility','off');

          xline_stimon = xline(-mean(fhlh{'fh','resp'}{1}.stimon_to_speech,1),'k','LineStyle',linestyles_first_last{1},'HandleVisibility','off');
         xline_stimonn = xline(-mean(fhlh{'lh','resp'}{1}.stimon_to_speech,1),'k','LineStyle',linestyles_first_last{2},'HandleVisibility','off');

          xline_gobeep_on = xline(-mean(fhlh{'fh','resp'}{1}.rt,1),'g','LineStyle',linestyles_first_last{1},'HandleVisibility','off');
         xline_gobeep_on = xline(-mean(fhlh{'lh','resp'}{1}.rt,1),'g','LineStyle',linestyles_first_last{2},'HandleVisibility','off');

         xline_sp_on = xline(0,'LineStyle','--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        
        
        title(['Train DBS-',char(upper(dbs_state))])
        
        legend({'First half trials','Last half trials'})
    end

    onoff_train{dbs_state,'fhlh'} = {fhlh};
    onoff_train{dbs_state,'stimresp_fh_minus_lh'} = {fhlh{'fh','resp'}{1}.stimresp - fhlh{'lh','resp'}{1}.stimresp};
    onoff_train{dbs_state,'mean_stimresp_fh_minus_lh'} = mean(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 
    onoff_train{dbs_state,'sem_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}) / sqrt(nchans); 
    onoff_train{dbs_state,'std_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 

end
    
    


%% plot mean stimresp firsthalf vs lasthalf for dbs on vs. dbs off
%%%% .... and test phase stimresp dbs on vs. off
% close all
fig = figure('Color','w');


subplot(1,2,1)
    hbar = bar(-onoff_train.mean_stimresp_fh_minus_lh,'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(-onoff_train.mean_stimresp_fh_minus_lh, -onoff_train.sem_stimresp_fh_minus_lh, 'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    title('Training phases','Second minus first half responses')
    ylabel({'Beta power','during stim (normed)'})
    % ylim(ylimits)
    hax = gca;
    hax.XTickLabels = {'DBS-on','DBS-off'};    set(gca,'TickLabelInterpreter','none')
    box off

subplot(1,2,2)
    hbar = bar(onoff_test.stimresp_mean,'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(onoff_test.stimresp_mean,onoff_test.stimresp_sem, 'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    title('Test phase')
    ylabel({'Beta power','during stimnormed)'})
    % ylim(ylimits)
    hax = gca;
    hax.XTickLabels = {'DBS-on','DBS-off'};    set(gca,'TickLabelInterpreter','none')
    box off    


%% power analysis of during-stim beta effects

% number of subs needed to make test-phase during-stim beta train-on vs. train-off sgnf
op.planned_nsubs = 10; 

chans_per_sub = 2; %%% only consider channels available during DBS ON for purposes of power analysis
planned_nchans = 2*op.planned_nsubs; 

pwran.test_stimresp_std = std(cell2mat(onoff_test.stimresp));

pwran.test_stimresp_onvsoff_nsubs = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran.test_stimresp_std], onoff_test{'off','stimresp_mean'}, 0.8, [], 'Tail','right') / chans_per_sub;

pwran.test_stimresp_onvsoff_power = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran.test_stimresp_std], onoff_test{'off','stimresp_mean'}, [],planned_nchans, 'Tail','right');

pwran.train_stimresp_fhlh_onvsoff_power = sampsizepwr('t', [onoff_train{'on','mean_stimresp_fh_minus_lh'}, onoff_train{'on','std_stimresp_fh_minus_lh'}], ...
    onoff_train{'off','mean_stimresp_fh_minus_lh'}, [], planned_nchans, 'Tail','left');

pwran
