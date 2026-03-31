 %%%% make plots for 2024 ro1 PD grant application
 % also output power analysis and t tests

clear

 setpaths_dbsmulti()

 % load([paths.groupanalysis, filesep, 'resp_all']) %%%% data created by dbs_response_profiles.m
   % load([paths.groupanalysis, filesep, 'resp_all_5early_base_subtract']) %%%% data created by dbs_response_profiles.m
   load([paths.groupanalysis, filesep, 'resp_all_5early_base_subtract_divide']) %%%% data created by dbs_response_profiles.m
  % load([paths.groupanalysis, filesep, 'resp_all_10early_base_subtract']) %%%% data created by dbs_response_profiles.m



close all


%% params
%%%% best stimresp train on vs. off quantitative result is with divide_by_trial_baseline off, op.trainphase_half_ntrials = 50

op.sem_from_all_trials = 1; 

op.trainphase_half_ntrials = 50; % number of trials to use in each 'half'


op.smooth_method = 'gaussian'; 

op.smooth_windowsize = 250; 
% op.smooth_windowsize = 1; 

op.plot_test_individual_chans = 0; % separate figure; not for grant


op.dbs_on_off_colors = [0.6 0.6 0;... 
                        0 0 0.7];
op.ebar_facealpha = 0.5; % error bar opacity

% accuracy bargraphs params
op.normalize_by_pretest = 1; % subtract pretest accuracy from later accuracy
    ops.ylims_if_acc_normalized = [0. 0.35]; % if normalizing, use these ylims (not implented for all plots yet)
op.plot_early_acc_on_vs_off_detailed_title = 0; 



op.tc_xlabel = 'Time after speech onset (s)'; % beta timecourse xlabel
op.tc_xlimits = [-6 6.5]; 

op.annot_height = 0.93; % relative vertical position of xline labels
op.annot_bkg_color = 1*[1 1 1]; 
op.annot_edge_col = 0*[1 1 1]; 
op.annot_fontsize = 10; 

op.annot_gobeep_xoffset_sec = 0.25; 

op.annot_text_stim = 'stim';
op.annot_text_gobeep = {'GO', 'cue'}; 
op.annot_text_sp = 'speech'; 

op.gobeep_color = [0 0 0]; 

% op.sp_off_color = [1 0.5 0.5]; 
op.sp_off_color = [0 0 0]; 

op.leg_height = 0.2; 

op.timecourse_linecolor = [0 0 0]; 

op.train_tc_ylims = [-0.25, 0.3]; 

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

% power analysis params; also applies to t test
op.planned_nsubs = 10; % number of subs needed to make test-phase during-stim beta train-on vs. train-off sgnf
op.chans_per_sub = 2; 

pwran_test.stimresp_onvsoff_tail = 'both';

pwran_train.stimresp_fhlh_onvsoff = 'left'; 

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

% indices of the first channel for each  subject
[~,unqsubinds] = unique(resp.sub);

for istate = 1:2
    dbs_state = dbs_state_list{istate};
    onoff_test.stimresp_mean(istate) = mean(onoff_test{dbs_state,'stimresp'}{1});
    onoff_test.stimresp_sem(istate) = std(onoff_test{dbs_state,'stimresp'}{1}) / sqrt(nchans);
    onoff_test.stimresp_std(istate) = std(onoff_test{dbs_state,'stimresp'}{1});

    onoff_test.rt_early(istate) = {cell2mat(cellfun(@(x)x.con{dbs_state,'rt_mean'},resp.test(unqsubinds),UniformOutput=false))};
    onoff_test.rt_early_mean(istate) = mean(onoff_test{dbs_state,'rt_early'}{1});
    onoff_test.rt_early_sem(istate) = std(onoff_test{dbs_state,'rt_early'}{1}) / sqrt(nsubs);
    onoff_test.rt_early_std(istate) = std(onoff_test{dbs_state,'rt_early'}{1});
end

plotcolors = table(op.dbs_on_off_colors,'VariableNames',{'color'},'RowNames',{'on','off'});

% individual channel responses sorted by condition
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

load_accuracy_results()

hfig = figure('Color',[1 1 1]); 
htl = tiledlayout(2,3);

%% test phase acc bargraph
  if op.normalize_by_pretest
        test_rows_to_plot = {'trainon','trainoff'};
        test_xticklabels = {'DBS ON','DBS OFF'}; 
        ylab = {'proportion syllables accurate','(normalized)'}; 
        ylimits = ops.ylims_if_acc_normalized; 
    else
        test_rows_to_plot = {'trainon','trainoff', 'novel'};
        test_xticklabels = {'DBS ON','DBS OFF','Novel'}; 
        ylab = {'proportion syllables accurate'}; 
        ylimits = [0 1]; 
    end

    nexttile(tilenum(htl,1,1))

    hax = gca;      
    hbar = bar(con.test.acc_early_mean(test_rows_to_plot,:),'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(con.test.acc_early_mean(test_rows_to_plot,:), con.test.acc_early_sem(test_rows_to_plot,:), 'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    if     op.plot_early_acc_on_vs_off_detailed_title
        title({['Test phase, first ' num2str(op.n_early_trials_per_stimid) ' of each stim'], ['n = ', num2str(nsubs)]})
    else 
        % title({'Early Testing phase',''})
        title({'Testing phase',''})
    end
    ylabel(ylab)
    ylim(ylimits)
    hax.XTickLabels = test_xticklabels;    set(gca,'TickLabelInterpreter','none')
    box off

%% test phase - average over subjects
    nexttile(tilenum(htl,1,2))
    
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

     % xline_stimoff = xline(stimoff_to_speech,'k','HandleVisibility','off');
      xline_gobeep_on = xline(gobeep,'g','HandleVisibility','off');
     xline_sp_on = xline(0,'LineStyle','--','Color',[0.5 0.5 0.5],'HandleVisibility','off');
    
    xline_sp_off = xline(op.sp_dur_aprx_test_early, 'LineStyle','-','Color',op.sp_off_color,'HandleVisibility','off');

    ylimits = ylim; 
    annot_yloc = ylimits(1) + op.annot_height*diff(ylimits);
    annot_stim_xloc = mean([stimon_to_speech, stimoff_to_speech]); 
    annot_gobeep_xloc = gobeep + op.annot_gobeep_xoffset_sec; 
    annot_speech_xloc = op.sp_dur_aprx_test_early / 2; 

    hannot_stim = text(annot_stim_xloc, annot_yloc, op.annot_text_stim, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
    % hannot_gobeep = text(annot_gobeep_xloc, annot_yloc, op.annot_text_gobeep, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
    hannot_sp = text(annot_speech_xloc, annot_yloc, op.annot_text_sp, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');


    xlabel(op.tc_xlabel)
    xlim(op.tc_xlimits)

    hleg = legend(cellfun(@(x)['trained-',upper(x)],conlist,'UniformOutput',0));
            hleg.Position(2) = op.leg_height + 0.4; % adjust height up for tiledlayout 
            hleg.Position(1) = hleg.Position(1) + .12;
    title(['Test phase'])
    
    

%% early vs. late  training




nan2 = nan(2,1);
trialtiming = table(nan2,nan2,nan2,'VariableNames',{'stimon_to_speech','stimoff_to_speech','rt'},   'RowNames',{'firsthalf','lasthalf'});

onoff_train = table('RowNames',{'on','off'}); % for compare on vs. off in training

xtime = linspace(op.trial_window_peri_speech(1),op.trial_window_peri_speech(2),size(resp.trainA{1}.trials.timecourse{1},2)); 


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

        nexttile(tilenum(htl,2,istate)) % second row, first and second columns
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

        if ~isempty(op.train_tc_ylims)
            ylim(op.train_tc_ylims) 
        end

        stimoff_to_sp_fh = -mean(fhlh{'fh','resp'}{1}.stimoff_to_speech,1);
        stimoff_to_sp_lh = -mean(fhlh{'lh','resp'}{1}.stimoff_to_speech,1);

        stimon_to_sp_fh = -mean(fhlh{'fh','resp'}{1}.stimon_to_speech,1);
        stimon_to_sp_lh = -mean(fhlh{'lh','resp'}{1}.stimon_to_speech,1);

        gobeep_fh = -mean(fhlh{'fh','resp'}{1}.rt,1); 
        gobeep_lh = -mean(fhlh{'lh','resp'}{1}.rt,1);

         % xline_stimoff(1) = xline(stimoff_to_sp_fh,'k','LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         % xline_stimoff(2) = xline(stimoff_to_sp_lh,'k','LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');

         %  xline_stimon(1) = xline(stimon_to_sp_fh,'k','LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         % xline_stimon(2) = xline(stimon_to_sp_lh,'k','LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');
         xline_stimon(2) = xline(stimon_to_sp_lh,'k','LineStyle','-','HandleVisibility','off');

         %  xline_gobeep_on(2) = xline(gobeep_fh,op.gobeep_color,'LineStyle',op.linestyles_fhlh{1},'HandleVisibility','off');
         % xline_gobeep_on(2) = xline(gobeep_lh,op.gobeep_color,'LineStyle',op.linestyles_fhlh{2},'HandleVisibility','off');
         xline_gobeep_on(2) = xline(gobeep_lh,'Color',op.gobeep_color,'LineStyle','-','HandleVisibility','off');

         xline_sp_on = xline(0,'LineStyle','-','Color',[0.5 0.5 0.5],'HandleVisibility','off');
        
         % get speech off time
        if strcmp(dbs_state,'on')
            % xline_sp_off_fh = xline(op.sp_dur_aprx_train_fh_on, 'LineStyle',op.linestyles_fhlh{1},'Color',op.sp_off_color,'HandleVisibility','off');
            xline_sp_off_lh = xline(op.sp_dur_aprx_train_lh_on, 'LineStyle',op.linestyles_fhlh{2},'Color',op.sp_off_color,'HandleVisibility','off');
            annot_speech_xloc = mean([op.sp_dur_aprx_train_fh_on, op.sp_dur_aprx_train_lh_on]) / 2; 
        elseif strcmp(dbs_state,'off')
            % xline_sp_off_fh = xline(op.sp_dur_aprx_train_fh_off, 'LineStyle',op.linestyles_fhlh{1},'Color',op.sp_off_color,'HandleVisibility','off');
            xline_sp_off_lh = xline(op.sp_dur_aprx_train_lh_off, 'LineStyle',op.linestyles_fhlh{2},'Color',op.sp_off_color,'HandleVisibility','off');
            annot_speech_xloc = mean([op.sp_dur_aprx_train_fh_off, op.sp_dur_aprx_train_lh_off]) / 2; 
        end

        ylimits = ylim
        annot_yloc = ylimits(1) + op.annot_height*diff(ylimits);
        annot_stim_xloc = mean([stimoff_to_sp_fh, stimoff_to_sp_lh, stimon_to_sp_fh, stimon_to_sp_lh]);
        annot_gobeep_xloc = mean([gobeep_fh, gobeep_lh]) + op.annot_gobeep_xoffset_sec;
        
        hannot_stim = text(annot_stim_xloc, annot_yloc, op.annot_text_stim, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
        % hannot_gobeep = text(annot_gobeep_xloc, annot_yloc, op.annot_text_gobeep, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');
        hannot_sp = text(annot_speech_xloc, annot_yloc, op.annot_text_sp, 'BackgroundColor',op.annot_bkg_color , 'EdgeColor',op.annot_edge_col, 'FontSize',op.annot_fontsize, 'HorizontalAlignment','Center');

        title(['Train DBS-',char(upper(dbs_state))])
        
        xlabel(op.tc_xlabel)
        xlim(op.tc_xlimits)

        hleg = legend({'First half trials','Last half trials'});
            hleg.Position(2) = op.leg_height - 0.08; 
            hleg.Position(1) = hleg.Position(1) + 0.01;



    onoff_train{dbs_state,'fhlh'} = {fhlh};
    onoff_train{dbs_state,'stimresp_fh_minus_lh'} = {fhlh{'fh','resp'}{1}.stimresp - fhlh{'lh','resp'}{1}.stimresp};
    onoff_train{dbs_state,'mean_stimresp_fh_minus_lh'} = mean(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 
    onoff_train{dbs_state,'sem_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}) / sqrt(nchans); 
    onoff_train{dbs_state,'std_stimresp_fh_minus_lh'} = std(onoff_train{dbs_state,'stimresp_fh_minus_lh'}{1}); 

end
    
    


 %% bargraph - mean stimresp firsthalf vs lasthalf for dbs on vs. dbs off


nexttile(tilenum(htl,1,3))

    hbar = bar(onoff_test.stimresp_mean,'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(onoff_test.stimresp_mean,onoff_test.stimresp_sem, 'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    title('Test phase')
    ylabel({'Beta power','during stim (normed)'})
    % ylim(ylimits)
    hax = gca;
    hax.XTickLabels = {'DBS-on','DBS-off'};    set(gca,'TickLabelInterpreter','none')
     
box off


%% bargraph - training phase stimresp dbs on vs. off




nexttile(tilenum(htl,2,3))

    hbar = bar(-onoff_train.mean_stimresp_fh_minus_lh,'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(-onoff_train.mean_stimresp_fh_minus_lh, -onoff_train.sem_stimresp_fh_minus_lh, 'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    title('Training phases','Second minus first half responses')
    ylabel({'Beta power','during stim (normed)'})
    % ylim(ylimits)
    hax = gca;
    hax.XTickLabels = {'DBS-on','DBS-off'};    set(gca,'TickLabelInterpreter','none')
    box off





%% power analysis of during-stim beta effects
%%% only consider channels available during DBS ON for purposes of power analysis


op.planned_nchans = op.chans_per_sub*op.planned_nsubs; 

pwran_test.stimresp_std = std(cell2mat(onoff_test.stimresp));
[~, pwran_test.stimresp_onvsoff_p] = ttest(onoff_test{'on','stimresp'}{1}, onoff_test{'off','stimresp'}{1}, 'Tail',pwran_test.stimresp_onvsoff_tail);
pwran_test.stimresp_onvsoff_nsubs = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran_test.stimresp_std], onoff_test{'off','stimresp_mean'}, 0.8, [], 'Tail',pwran_test.stimresp_onvsoff_tail) / op.chans_per_sub;
pwran_test.stimresp_onvsoff_power = sampsizepwr('t', [onoff_test{'on','stimresp_mean'}, pwran_test.stimresp_std], onoff_test{'off','stimresp_mean'}, [],op.planned_nchans, 'Tail',pwran_test.stimresp_onvsoff_tail);
[~, pwran_test.rt_onvsoff_p] = ttest(onoff_test{'on','rt_early'}{1}, onoff_test{'off','rt_early'}{1}); % could do more sophisticated analysis to analyze single-trial or single-stim data.... currently just using subject means

pwran_train.stimresp_fhlh_onvsoff_power = sampsizepwr('t', [onoff_train{'on','mean_stimresp_fh_minus_lh'}, onoff_train{'on','std_stimresp_fh_minus_lh'}], ...
    onoff_train{'off','mean_stimresp_fh_minus_lh'}, [], op.planned_nchans, 'Tail',pwran_train.stimresp_fhlh_onvsoff);
[~, pwran_train.stimresp_fhlh_onvsoff_p] = ttest(onoff_train.stimresp_fh_minus_lh{2},onoff_train.stimresp_fh_minus_lh{1}, 'Tail',pwran_train.stimresp_fhlh_onvsoff ); 

pwran_test
pwran_train




%% test phase, early trials - plot individual channels


if op.plot_test_individual_chans

    hfig_individual_chan = figure('Color','w'); 
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

%% timecourse plot formatting function


function format_timecourse(cfg)

% add shading for stim period and speech period
color code dbs on off tc lines

end