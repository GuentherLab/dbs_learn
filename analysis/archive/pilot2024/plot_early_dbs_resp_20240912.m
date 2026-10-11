 %%%% plot responses from early trials in test phase sorted by training condition

 setpaths_dbsmulti()

 % load([paths.groupanalysis, filesep, 'resp_all']) %%%% data created by dbs_response_profiles.m
   load([paths.groupanalysis, filesep, 'resp_all_5early']) %%%% data created by dbs_response_profiles.m
  % load([paths.groupanalysis, filesep, 'resp_all_10early']) %%%% data created by dbs_response_profiles.m

% close all


%% params
%%%% best stimresp train on vs. off quantitative result is with divide_by_trial_baseline off, op.trainphase_half_ntrials = 50

%%%% does not yet affect bars, only timecourse
op.divide_by_trial_baseline  = 0; 
op.divide_by_chan_baseline = 0; %%% need to reprogram to affect timecourses

op.sem_from_all_trials = 1; 

op.trainphase_half_ntrials = 50; % number of trials to use in each 'half'


op.smooth_method = 'gaussian'; 

op.smooth_windowsize = 250; 
% op.smooth_windowsize = 1; 

op.plot_test_individual_chans = 0; 
op.plot_test_avg_all_chans = 0; 
op.plot_trainphase_first_vs_last_half = 1; 
op.bar_stimresp = 1; 

op.dbs_on_off_colors = [0.6 0.6 0;... 
                        0 0 0.7];
op.ebar_facealpha = 0.5; % error bar opacity
op.xlabel = 'Time after speech onset (s)'; 
op.xlimits = [-6 6.5]; 

op.annot_height = 0.9; % relative vertical position of xline labels
op.annot_bkg_color = 1*[1 1 1]; 
op.annot_edge_col = 0*[1 1 1]; 
op.annot_fontsize = 10; 

op.annot_gobeep_xoffset_sec = 0.25; 

op.annot_text_stim = 'stim';
op.annot_text_gobeep = {'GO', 'cue'}; 
op.annot_text_sp = 'speech'; 

op.sp_off_color = [1 0.5 0.5]; 

op.leg_height = 0.2; 

op.timecourse_linecolor = [0 0 0]; 

conlist = {'on','off'};
% conlist = {'on','off','novel'}; 



op.linestyles_fhlh = {'-','--'}; % line styles for first half, last half



% approximate average speech duration in seconds to use
%%% derived from random sample of durations across subs/tasks
op.sp_dur_aprx_train_fh_on = 2.63; % first half dbs on
op.sp_dur_aprx_train_fh_off = 2.74; % first half dbs off
op.sp_dur_aprx_train_lh_on = 2.72; % last half dbs on
op.sp_dur_aprx_train_lh_off = 2.86; % last half dbs off
op.sp_dur_aprx_test_early = 2.575; % early test phase trials


% number of subs needed to make test-phase during-stim beta train-on vs. train-off sgnf
pwran.planned_nsubs = 10; 

set(0,'defaulttextinterpreter','none')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare data
 % keep only channels with recordings during test phase
 resp = resp(~cellfun(@isempty,resp.test),:);

 ncons = length(conlist);
nchans = height(resp); 
nsubs = height(subs);

allchans_stimon_to_speech = [];
allchans_stimoff_to_speech = [];
allchans_rt = [];

onoff_test = table(repmat({nan(nchans,1)},2,1),nan(2,1),nan(2,1),'VariableNames',{'stimresp','stimresp_mean','stimresp_sem'},'RowNames',{'on','off'}); % for compare on vs. off in training

dbs_state_list = {'on','off'};

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
     tc_mean = contab{thiscon,'resp_mean'}{1}; 
                resp{ichan,['mean_',thiscon]} = tc_mean;    
            plot( smoothdata(tc_mean, 2, op.smooth_method, op.smooth_windowsize) )
            
        end
        legend(conlist)
        title([thissub, ' ', thischan])
        ylabel('Beta power')
    
    
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
                error('not implemented')
                % % tc_sem = nanstd(tc_allchans) / sqrt(nchans);           % should average over channel averages, not all trials
            end
            tc_sem_lims = [tc_mean-tc_sem; tc_mean+tc_sem]; 
            
    
            xtime = [1:size(tc_mean,2)] ./ op.wavtrans_out_freq + op.trial_window_peri_speech(1); 
    
            tc_sem_lims_to_plot = smoothdata(tc_sem_lims, 2, op.smooth_method, op.smooth_windowsize);
    
            hold on
            hfill = fill([xtime, fliplr(xtime)], [tc_sem_lims_to_plot(1,:), fliplr(tc_sem_lims_to_plot(2,:))], [0.8 0.8 0.8],...
                'FaceAlpha', op.ebar_facealpha, 'HandleVisibility','off'); % standard error
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
    
    stimon_to_speech = -mean(allchans_stimon_to_speech);
    stimoff_to_speech = -mean(allchans_stimoff_to_speech); 

     xline_stimon = xline(stimon_to_speech,'k','HandleVisibility','off');
    % xline_stimon = xline(-allchans_stimon_to_speech,'k','LineStyle', ":",'HandleVisibility','off');
    
    gobeep = -mean(allchans_rt);

     xline_stimoff = xline(stimoff_to_speech,'k','HandleVisibility','off');
      xline_gobeep_on = xline(gobeep,'g','HandleVisibility','off');
     xline_sp_on = xline(0,'LineStyle','--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    
    xline_sp_off = xline(op.sp_dur_aprx_test_early, 'LineStyle','-','Color',op.sp_off_color,'HandleVisibility','off');

    ylimits = ylim; 
    annot_yloc = ylimits(1) + op.annot_height*diff(ylimits);
    annot_stim_xloc = mean([stimon_to_speech, stimoff_to_speech]); 
    annot_gobeep_xloc = gobeep + op.annot_gobeep_xoffset_sec; 
    annot_speech_xloc = op.sp_dur_aprx_test_early / 2; 

    hannot_stim = text(annot_stim_xloc, annot_yloc, op.annot_text_stim, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
    hannot_gobeep = text(annot_gobeep_xloc, annot_yloc, op.annot_text_gobeep, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
    hannot_sp = text(annot_speech_xloc, annot_yloc, op.annot_text_sp, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');


    xlabel(op.xlabel)
    xlim(op.xlimits)

    hleg = legend(conlist);
            hleg.Position(2) = op.leg_height; 
    title(['All channels average - test phase'])
    
    
    box off
end

%% early vs. late  training


op.linestyles_fhlh = {'-','--'}; % line styles for first half, last half

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
    
        fhidx = 1:op.trainphase_half_ntrials;
        lhidx = length(resp{ichan,this_task}{1}.trials.timecourse) - op.trainphase_half_ntrials : length(resp{ichan,this_task}{1}.trials.timecourse); 
    
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
             %% normalize by baseline magnitude 




               if op.divide_by_trial_baseline 
                    

                    % normalize (divide) magnitude by baseline in each trial
                    fhlh{ihalf,'resp'}{1}.tc{ichan} = fhlh{ihalf,'resp'}{1}.tc{ichan} ./ abs(fhlh{ihalf,'resp'}{1}{ichan,'base'}{1});
                                % tc_all_chans_trials = tc_all_chans_trials ./ cell2mat(fhlh{ihalf,'resp'}{1}.base);
        
                elseif op.divide_by_chan_baseline 
                    tc_all_chans_trials = [];
                    error('need to reprogram to affect timecourses')
                    for ichan = 1:nchans
                        tc_all_chans_trials = [tc_all_chans_trials; fhlh{ihalf,'resp'}{1}.tc{ichan} / fhlh{ihalf,'resp'}{1}.base_mean(ichan)]; 
                    end
                end



            fhlh{ihalf,'resp'}{1}.stimresp(ichan) = nanmean(nanmean(fhlh{ihalf,'resp'}{1}.tc{ichan}(:,stiminds),2));

        end
    
    end

    for ihalf = 1:2
        tc_all_chans_trials = cell2mat(fhlh{ihalf,'resp'}{1}.tc); 

        fhlh{ihalf,'resp'}{1}.base_mean = cellfun(@nanmean,fhlh{ihalf,'resp'}{1}.base); % channel baseline means across trials
        ntrials = size(fhlh{ihalf,'resp'}{1}.tc{1},1);
        
        
        fhlh{ihalf,'mean_tc'} = {nanmean(tc_all_chans_trials)};
        if op.sem_from_all_trials 
            
            fhlh{ihalf,'sem_tc'} = {nanstd(tc_all_chans_trials,1) / sqrt(size(tc_all_chans_trials,1))} ;
        else
            error('not implemented')
        end
    end

    if op.plot_trainphase_first_vs_last_half
        subplot(1,2,istate)
        hold on

         % draw error bars first
        for ihalf = 1:2
            tc_sem_lims = [fhlh{ihalf,'mean_tc'}{1} - fhlh{ihalf,'sem_tc'}{1}; fhlh{ihalf,'mean_tc'}{1} + fhlh{ihalf,'sem_tc'}{1}];
            tc_sem_lims_plot = smoothdata(tc_sem_lims, 2, op.smooth_method, op.smooth_windowsize);
            hfill = fill([xtime, fliplr(xtime)], [tc_sem_lims_plot(1,:), fliplr(tc_sem_lims_plot(2,:))], [0.8 0.8 0.8],...
                'FaceAlpha',op.ebar_facealpha,'HandleVisibility','off'); % standard error
    %         hfill.LineStyle = 'none'; % no border
            hfill.EdgeColor = [0.8 0.8 0.8]; 
        
        end

         % mean timecourses
        mean_firsthalf_plot = smoothdata(fhlh{'fh','mean_tc'}{1}, 2, op.smooth_method, op.smooth_windowsize);
        mean_lastthalf_plot = smoothdata(fhlh{'lh','mean_tc'}{1}, 2, op.smooth_method, op.smooth_windowsize);
        hplot = plot(xtime,mean_firsthalf_plot,'Color',op.timecourse_linecolor,'LineWidth',2,'LineStyle',op.linestyles_fhlh{1});
        hplot2 = plot(xtime,mean_lastthalf_plot,'Color',op.timecourse_linecolor,'LineWidth',2,'LineStyle',op.linestyles_fhlh{2});

        ylabel('Beta power')
        
        stimoff_to_sp_fh = -mean(fhlh{'fh','resp'}{1}.stimoff_to_speech,1);
        stimoff_to_sp_lh = -mean(fhlh{'lh','resp'}{1}.stimoff_to_speech,1);

        stimon_to_sp_fh = -mean(fhlh{'fh','resp'}{1}.stimon_to_speech,1);
        stimon_to_sp_lh = -mean(fhlh{'lh','resp'}{1}.stimon_to_speech,1);

        gobeep_fh = -mean(fhlh{'fh','resp'}{1}.rt,1); 
        gobeep_lh = -mean(fhlh{'lh','resp'}{1}.rt,1);

         xline_stimoff(1) = xline(stimoff_to_sp_fh,'k','LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         xline_stimoff(2) = xline(stimoff_to_sp_lh,'k','LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');

          xline_stimon(1) = xline(stimon_to_sp_fh,'k','LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         xline_stimon(2) = xline(stimon_to_sp_lh,'k','LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');

          xline_gobeep_on(2) = xline(gobeep_fh,'g','LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         xline_gobeep_on(2) = xline(gobeep_lh,'g','LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');

         xline_sp_on = xline(0,'LineStyle','-','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        
         % get speech off time
        if strcmp(dbs_state,'on')
            xline_sp_off_fh = xline(op.sp_dur_aprx_train_fh_on, 'LineStyle',op.linestyles_fhlh{1},'Color',[1 0.1 0.1],'HandleVisibility','off');
            xline_sp_off_lh = xline(op.sp_dur_aprx_train_lh_on, 'LineStyle',op.linestyles_fhlh{2},'Color',op.sp_off_color,'HandleVisibility','off');
            annot_speech_xloc = mean([op.sp_dur_aprx_train_fh_on, op.sp_dur_aprx_train_lh_on]) / 2; 
        elseif strcmp(dbs_state,'off')
            xline_sp_off_fh = xline(op.sp_dur_aprx_train_fh_off, 'LineStyle',op.linestyles_fhlh{1},'Color',[1 0.1 0.1],'HandleVisibility','off');
            xline_sp_off_lh = xline(op.sp_dur_aprx_train_lh_off, 'LineStyle',op.linestyles_fhlh{2},'Color',op.sp_off_color,'HandleVisibility','off');
            annot_speech_xloc = mean([op.sp_dur_aprx_train_fh_off, op.sp_dur_aprx_train_lh_off]) / 2; 
        end

        ylimits = ylim; 
        annot_yloc = ylimits(1) + op.annot_height*diff(ylimits);
        annot_stim_xloc = mean([stimoff_to_sp_fh, stimoff_to_sp_lh, stimon_to_sp_fh, stimon_to_sp_lh]);
        annot_gobeep_xloc = mean([gobeep_fh, gobeep_lh]) + op.annot_gobeep_xoffset_sec;
        
        hannot_stim = text(annot_stim_xloc, annot_yloc, op.annot_text_stim, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
        hannot_gobeep = text(annot_gobeep_xloc, annot_yloc, op.annot_text_gobeep, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
        hannot_sp = text(annot_speech_xloc, annot_yloc, op.annot_text_sp, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');

        title(['Train DBS-',char(upper(dbs_state))])
        
        xlabel(op.xlabel)
        xlim(op.xlimits)

        hleg = legend({'First half trials','Last half trials'});
            hleg.Position(2) = op.leg_height; 
            % hleg.Position(1) = 0.8; 

    end

    onoff_train{dbs_state,'fhlh'} = {fhlh};
    onoff_train{dbs_state,'stimresp_fh_minus_lh'} = {fhlh{'fh','resp'}{1}.stimresp - fhlh{'lh','resp'}{1}.stimresp};
    onoff_train{dbs_state,'mean_stimresp_fh_minus_lh'} = mean(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 
    onoff_train{dbs_state,'sem_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}) / sqrt(nchans); 
    onoff_train{dbs_state,'std_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 

end
    
    


%% bargraph - mean stimresp firsthalf vs lasthalf for dbs on vs. dbs off
%%%% .... and test phase stimresp dbs on vs. off

% close all

if op.bar_stimresp

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
    ylabel({'Beta power','during stim (normed)'})
    % ylim(ylimits)
    hax = gca;
    hax.XTickLabels = {'DBS-on','DBS-off'};    set(gca,'TickLabelInterpreter','none')
    box off   

end


%% power analysis of during-stim beta effects

pwran.chans_per_sub = 2; 

chans_per_sub = 2; %%% only consider channels available during DBS ON for purposes of power analysis
planned_nchans = pwran.chans_per_sub*pwran.planned_nsubs; 

pwran.test_stimresp_std = std(cell2mat(onoff_test.stimresp));

pwran.test_stimresp_onvsoff_nsubs = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran.test_stimresp_std], onoff_test{'off','stimresp_mean'}, 0.8, [], 'Tail','right') / chans_per_sub;

pwran.test_stimresp_onvsoff_power = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran.test_stimresp_std], onoff_test{'off','stimresp_mean'}, [],planned_nchans, 'Tail','right');

pwran.train_stimresp_fhlh_onvsoff_power = sampsizepwr('t', [onoff_train{'on','mean_stimresp_fh_minus_lh'}, onoff_train{'on','std_stimresp_fh_minus_lh'}], ...
    onoff_train{'off','mean_stimresp_fh_minus_lh'}, [], planned_nchans, 'Tail','left');

[~, pwran.train_stimresp_fhlh_onvsoff_p] = ttest(onoff_train.stimresp_fh_minus_lh{1},onoff_train.stimresp_fh_minus_lh{2}); 

pwran
