%%%% compile accuracies from dbsmulti subjects, look for training condition effects on learning

clear
close all

warning('off','MATLAB:table:RowsAddedNewVars')

%% params
% number of trials from test phase to analyze from each stimulus ID
op.n_early_trials_per_stimid = 5;

% if true, exclude trials in which subject did not produce spoken response
%%%% these trials will be counted as having zero accuracy
op.exclude_no_response_trials = 1; 

% subtract pretest accuracy from later accuracy
op.normalize_by_pretest = 0; 

%%%%% analyses and plots to generate
op.plot_bar_mean_sem = 0; 
op.power_analysis_plot = 0; 
    % op.power_test_tail = 'right';
    op.power_test_tail = 'both'; 
op.plot_sub_training_timecourse = 0; 
op.plot_all_learning_curve_fits = 0;
op.early_train_acc_on_vs_off = 1; 
op.plot_learn_params_on_vs_off = 0; 
op.plot_learn_params_A_vs_B = 0; 


op.dbs_order_colors = [0.6 0.6 0;... % for plotting: ON-OFF, OFF-ON
                        0 0 0.7]; 
%%

%%% load data
setpaths_dbsmulti() 

subtable_file = [paths.code_analysis, filesep, 'dbsmulti_subs_master.csv']; 

subs = readtable(subtable_file);

%%%%%%%% only analyze certain subjects
subs = subs(1:4,:);
% subs = subs(2:4,:);

nsubs = height(subs); 
stim_master = readtable([paths.code_dbsmulti, filesep, 'stim_master_dbsmulti.xlsx']); 

phases_chron = table({'assess','pretest','trainA','trainB','test'}','VariableNames',{'phase'});
nphases_chron = height(phases_chron);
phasetypes = {'train','test'};
nphasestypes = length(phasetypes); 

subcel = cell(nsubs,1); 
subnan = nan(nsubs,1); 
subnan2 = nan(nsubs,2);
subnan3 = nan(nsubs,3); 
sublog2 = false(nsubs,2); 
subs = [subs, table(subnan3,subcel,subcel,subcel,subcel,subcel,'VariableNames',...
                {'dbs_order_color','tr_assess','tr_pretest','tr_trainA','tr_trainB','tr_test'})];

con = struct; % conditions accuracy results
            
for isub = 1:nsubs
   thissub = subs.sub{isub};
   
   %%% import trialtables for this subject
   thissub_acc_filename = [paths.data, filesep, 'sub-', thissub, '_acc.xlsx']; 
   for iphase = 1:nphases_chron
       thisphase = phases_chron.phase{iphase};
      thisphasetab = readtable(thissub_acc_filename,'Sheet',thisphase);  % get trials for this phase for this sub 
      thisphasetab = thisphasetab(thisphasetab.unusable_trial ~= 1, :);  % remove unusable trials
      
      if op.exclude_no_response_trials && contains('no_response', thisphasetab.Properties.VariableNames)
          thisphasetab = thisphasetab(thisphasetab.no_response ~= 1, :);  % remove unusable trials      
      end
      
      [~,stimmaster_rows] = ismember(thisphasetab.name,stim_master.name); 
            thisphasetab.stim_group = stim_master.stim_group(stimmaster_rows); % label stim group, in case it's missing
      thisphasetab.propacc = thisphasetab.n_correct_syls ./ subs.n_syls(isub); % get proportion of correct syls

        if op.normalize_by_pretest
            if any(strcmp(thisphase,{'trainA','trainB','test'}))
                pretest_grpstats = grpstats(subs.tr_pretest{isub}(:,{'name','propacc'}), 'name', {'mean'});
                [~,pretest_grpstats_rows] = ismember(thisphasetab.name,pretest_grpstats.name); 
                pretest_acc_to_subtract = nan(height(thisphasetab),1); % initialize
                    pretest_acc_to_subtract(pretest_grpstats_rows~=0) = pretest_grpstats.mean_propacc(pretest_grpstats_rows(pretest_grpstats_rows~=0)); 
                thisphasetab.propacc = thisphasetab.propacc - pretest_acc_to_subtract; % normalize by pretest vals on each stim
            end
        end

       subs{isub,['tr_',thisphase]} = {thisphasetab};
       subs{isub,['stim_',thisphase]} = {grpstats(thisphasetab(:,{'name','propacc'}), 'name', {'mean','sem','std'}, ...% stim stats in this phase, this sub
                                                    'VarNames',{'name','GroupCount','mean_acc','sem_acc','std_acc'})}; 
            [~,stimmaster_rows] = ismember(subs{isub,['stim_',thisphase]}{1}.name,stim_master.name); 
            subs{isub,['stim_',thisphase]}{1}.con = stim_master.stim_group(stimmaster_rows); % label condition for each stim
            subs{isub,['stim_',thisphase]}{1} = sortrows(subs{isub,['stim_',thisphase]}{1},'RowNames'); % alphabetize

        %%%%% get early trials and add to sub.stim_ tables for each stim
        if any(strcmp(thisphase,{'trainA','trainB','test'}))
            phasestim = subs{isub,['stim_',thisphase]}{1}.name; 
            nphasestim = length(phasestim);
            for istim = 1:nphasestim
                thisstim = phasestim{istim};
                subs{isub,['stim_',thisphase]}{1}.early_inds{istim} = find(string(thisphasetab.name) == thisstim, op.n_early_trials_per_stimid); 
            end
            all_early_inds = sort(cell2mat(subs{isub,['stim_',thisphase]}{1}.early_inds)); % first X trials of each unique stim in this phase
    
            grpstats_early = grpstats(thisphasetab(all_early_inds,{'name','propacc'}), 'name', {'mean','sem','std'},...
                'VarNames',{'name','GroupCount','mean_acc_early','sem_acc_early','std_acc_early'});
                grpstats_early.GroupCount = [];
            subs{isub,['stim_',thisphase]} = {join(subs{isub,['stim_',thisphase]}{:}, grpstats_early)}; 
        end

       if any(strcmp(thisphase,{'test'})) % if this phase has multiple stim conditions... add pretest later
            subs{isub,['con_',thisphase]} = {grpstats(thisphasetab(:,{'stim_group','propacc'}), 'stim_group', {'mean','sem','std'})};  % condition in stats this phase, this sub
       end
   end

   if string(subs.train_order{isub})=='on_off'
        subs.dbs_order_color(isub,:) = op.dbs_order_colors(1,:);
   elseif string(subs.train_order{isub})=='off_on'
        subs.dbs_order_color(isub,:) = op.dbs_order_colors(2,:);
   end
end

%% get accuracy for each condition for training (if applicable) and testing
for iphasetype = 1:nphasestypes % train vs. test
   thisphasetype = phasetypes{iphasetype}; 
   switch thisphasetype 
       case 'train'
           cons_to_analyze = {'trainon','trainoff'};
       case 'test'
           cons_to_analyze = {'trainon','trainoff','novel'};
   end
   ncons_to_analyze = length(cons_to_analyze); 
   concel = cell(ncons_to_analyze,1); 
   connan = nan(ncons_to_analyze,1); 
   fit_tab_cells = repmat({table(subs.sub,subnan2,subnan2,subnan2,sublog2,'VariableNames',{'sub','max_acc','acc_increase','exp_learn_rate','starts_at_ceil'},'RowNames',subs.sub)},ncons_to_analyze,1); 
   con.(thisphasetype) = table(cons_to_analyze',concel,connan,connan,connan,fit_tab_cells,...
              'VariableNames',{'con','acc_early','acc_early_mean','acc_early_std','acc_early_sem','learnfit'}, 'RowNames',cons_to_analyze');
   
   for icon = 1:ncons_to_analyze
      thiscon = cons_to_analyze{icon}; 

      con.(thisphasetype).acc_early{icon} = subnan;
      for isub = 1:nsubs
         thissub = subs.sub{isub};

         %%%% determine which train phase was dbs on vs. off for this sub
           if string(subs.train_order{isub}) == 'on_off'
               train_on_label = 'A'; train_off_label = 'B'; 
           elseif string(subs.train_order{isub}) == 'off_on'
              train_on_label = 'B'; train_off_label = 'A'; 
           else 
               error('train_order not recognized')
           end

          switch thiscon
              case 'trainon'
                  conlab = ['train',train_on_label];
              case 'trainoff'
                  conlab = ['train',train_off_label];
              case 'novel'
                  conlab = 'novel';
          end

          % get name of the phase we are looking at
           if strcmp(thisphasetype,'train')
               phaselab = conlab;
           elseif strcmp(thisphasetype,'test')
                phaselab = 'test';
           end

        sub_phase_stim_tab = subs{isub,['stim_',phaselab]}{1};
        sub_phase_stim_rows = strcmp(sub_phase_stim_tab.con, conlab); 
         con.(thisphasetype).acc_early{icon}(isub) = mean(sub_phase_stim_tab.mean_acc_early(sub_phase_stim_rows)); % mean acc for this sub, this condition, this phase
       con.(thisphasetype).acc_early_mean(icon) = mean(con.(thisphasetype).acc_early{icon}); % group mean across subjects
       con.(thisphasetype).acc_early_std(icon) = std(con.(thisphasetype).acc_early{icon}); % group mean across subjects
       con.(thisphasetype).acc_early_sem(icon) = std(con.(thisphasetype).acc_early{icon}) ./ sqrt(nsubs); % group mean across subjects
      end

   end

end

% get phase info for later
for iphase = 1:nphases_chron
    thisphase = phases_chron.phase{iphase};
    phases_chron.n_unq_stim(iphase) = max(cellfun(@height,subs{:,['stim_',thisphase]})); 
end
trainphases = phases_chron(contains(phases_chron.phase,'train'),:); 
n_train_phases = height(trainphases); % train A and train B
fit_tab_cells = repmat({table(subs.sub,subnan2,subnan2,subnan2,sublog2,'VariableNames',{'sub','max_acc','acc_increase','exp_learn_rate','starts_at_ceil'},'RowNames',subs.sub)},n_train_phases,1); 
fit_tab_cells_to_plot = repmat({table(subs.sub,subnan,subnan,subnan,'VariableNames',{'sub','max_acc','acc_increase','exp_learn_rate'},'RowNames',subs.sub)},n_train_phases,1); 
trainphases = [trainphases, table(fit_tab_cells, fit_tab_cells_to_plot, 'VariableNames',{'learnfit','learnfit_to_plot'})]; 
trainphases.Properties.RowNames = trainphases.phase; 

%%%%%%%%%%%%%%%%%%% stats tests and plotting %%%%%%%%%%%%%%%%%%%%

%% - dbs-on vs. dbs-off, group level


acc_for_anova = cell2mat(con.test.acc_early);
con_plot_labels = strrep( con.test.con, 'train', 'train_dbs_');
cons_for_anova = repelem(con_plot_labels,  nsubs); 
[anova_p,anova_tbl,anova_stats] = anova1(acc_for_anova, cons_for_anova,'off');

[~, p_dbs_on_vs_off] = ttest(con.test.acc_early{1}, con.test.acc_early{2});

%% make bar plot of accuracy
if op.plot_bar_mean_sem
    hfig = figure('Color','w');
    hax = gca;      
    hbar = bar(con.test.acc_early_mean,'FaceColor','flat','CData',[0.4 0.4 0.4]);
    hold on
    hebar = errorbar(con.test.acc_early_mean,con.test.acc_early_sem,'LineWidth',2,'LineStyle','none','Color',[0 0 0 ]);
    title({['Test phase, first ' num2str(op.n_early_trials_per_stimid) ' of each stim'], ['n = ', num2str(nsubs)]})
    ylabel('proportion syllables accurate')
    hax.XTickLabels = con_plot_labels;    set(gca,'TickLabelInterpreter','none')
    box off
end



%% power analysis - detecting difference of early trials vs train-on vs. train-off

if op.power_analysis_plot
    numsubs_vals = 1:100;
    acc_early_all = [con.test.acc_early{'trainon'};  con.test.acc_early{'trainoff'}]
    % % % acc_early_mean = mean(acc_early_all);
    % % % acc_early_effect_size = con.test.acc_early_mean('trainon') - con.test.acc_early_mean('trainoff');
    acc_early_std = std(acc_early_all); 
    
    nout = sampsizepwr('t',[con.test.acc_early_mean('trainoff') acc_early_std], con.test.acc_early_mean('trainon') ,0.80, [], 'Tail',op.power_test_tail)
    
    pwrout = sampsizepwr('t',[con.test.acc_early_mean('trainoff') acc_early_std], con.test.acc_early_mean('trainon') ,[],numsubs_vals, 'Tail',op.power_test_tail);
    
    hfig = figure('Color','w'); 
    plot(numsubs_vals,pwrout,'b-',nout,0.8,'ro')
    title({'Power versus Sample Size', ['using n = ', num2str(nsubs) ' DBS pilot subjects'],['Tail = ',op.power_test_tail]})
    xlabel('Sample Size')
    ylabel('Power')
    box off
end

%% plot performance accuracy timecourse for each subject/phase in individual panes
%     used for getting an impression of what function would be good to fit to training timecourses
if op.plot_sub_training_timecourse
    hfig = figure('Color','w'); 
    htl = tiledlayout(n_train_phases,nsubs); 
    for isub = 1:nsubs
        thissub = subs.sub{isub};
        for iphase = 1:n_train_phases
            thisphase = trainphases.phase{iphase};
            nexttile(tilenum(htl,iphase,isub),[1,1])

            acc_trials_plot_op.intable = subs{isub,['tr_',thisphase]}{1}; % accuracy table
            acc_trials_plot_op.sortvar = 'name'; % sort by stim name
            acc_trials_plot_op.plotvar = 'propacc'; % plot proportion accuracy
            acc_trials_plot_op.windowsize = 1; % plot individual trials, don't do windowing
            acc_trials_plot_op.endpoint_method = 'discard';
            % acc_trials_plot_op.line_style = '-'; % 
                acc_trials_plot_op.line_style = 'none'; % do not connect points 
            acc_trials_plot_op.marker_style = '.';
            acc_trials_plot_op.marker_size = 20; 
            acc_trials_plot_op.ylimits = [0 1.05]; 
            acc_trials_plot_op.jitter_prop_of_yrange = 0.005; 
            acc_trials_plot_op.legend_location = 'southeast';
            acc_trials_plot_op.newfig = 0; 
            plot_windowed_timecourse(acc_trials_plot_op); 

            htitle = title([thissub '...' thisphase '..... trial-window size = ' num2str(op.windowsize)],'Interpreter','none');
        end
    end
end

%% fit learning curves, make plots if specified
if op.plot_all_learning_curve_fits
    curve_fit_plot_op.sortvar = 'name'; % sort by stim name
    curve_fit_plot_op.plotvar = 'propacc'; % plot proportion accuracy
    curve_fit_plot_op.windowsize = 1; % plot individual trials, don't do windowing
    curve_fit_plot_op.endpoint_method = 'discard';
    curve_fit_plot_op.line_style = 'none'; % do not connect points 
    curve_fit_plot_op.marker_style = '.';
    curve_fit_plot_op.marker_size = 20; 
    curve_fit_plot_op.ylimits = [0 1.05]; 
    curve_fit_plot_op.jitter_prop_of_yrange = 0.005; 
    curve_fit_plot_op.legend_location = 'southeast';
    curve_fit_plot_op.newfig = 0; 

    hfig = figure('Color','w'); 
    htl = tiledlayout(sum(trainphases.n_unq_stim),nsubs); 
    sgtitle(['train phase curve fitting...trial-window size = ' num2str(curve_fit_plot_op.windowsize)])
end

for isub = 1:nsubs
    thissub = subs.sub{isub};
    rowcounter = 0; 
    for iphase = 1:n_train_phases
        thisphase = trainphases.phase{iphase};
        sub_phase_stim = subs{isub,['stim_',thisphase]}{1}.name;

       if string(subs.train_order{isub}) == 'on_off'
           if string(thisphase) == 'trainA';
                dbs_con = 'on';
           else string(thisphase) == 'trainB';
                dbs_con = 'off';
           end
       elseif string(subs.train_order{isub}) == 'off_on'
           if string(thisphase) == 'trainA';
                dbs_con = 'off';
           else string(thisphase) == 'trainB';
                dbs_con = 'on';
           end
       end

        for istim = 1:trainphases.n_unq_stim(iphase)
            rowcounter = rowcounter + 1; 
            thisstim = sub_phase_stim{istim}; 
            thisphasetab = subs{isub,['tr_',thisphase]}{1};
            phasestimtab = thisphasetab(string(thisphasetab.name)==thisstim, :); 

            % example paper where this exponential model was used + explained: https://doi.org/10.1177/1545968316662526
            %%% in paper above.... A = max_acc, B = acc_increase (reversed sign because we're using accuracy instead of error)
            %%% .................. alpha = exp_learn_rate
            expModel = fittype('max_acc - acc_increase * exp(-exp_learn_rate * [trialnum])', 'independent', 'trialnum', 'coefficients', {'max_acc','acc_increase', 'exp_learn_rate'});
            fit_ops =  fitoptions('Method', 'NonlinearLeastSquares',...
                'StartPoint', [0.93, 0.7, 0.12],... % start at reasonably interpretable coef values
                'Lower', [0, -1, -inf],...  % restrict accuracy between 0 and 1
                'Upper', [1, 1, inf]... % restrict accuracy between 0 and 1
                );
            [fit_result, gof] = fit([1:height(phasestimtab)]', phasestimtab.propacc, expModel, fit_ops);
            
            % mark whether this phase/stim started at ceiling accuracy
            %%%% use the op.n_early_trials_per_stimid to decide how many trials should be tested
            %%%% ... this requirement might not catch all stim/phases which could be reasonably described as starting at ceiling
                starts_at_ceiling = all( phasestimtab.n_correct_syls(1:op.n_early_trials_per_stimid ) == subs.n_syls(isub) ); 

            %%%%%%%%%%% add data to con.train table, where phases are organized by DBS ON vs. OFF
            con.train.learnfit{['train',dbs_con],1}.max_acc(isub,istim) = fit_result.max_acc; 
            con.train.learnfit{['train',dbs_con],1}.acc_increase(isub,istim) = fit_result.acc_increase; 
            con.train.learnfit{['train',dbs_con],1}.exp_learn_rate(isub,istim) = fit_result.exp_learn_rate;
            con.train.learnfit{['train',dbs_con],1}.starts_at_ceil(isub,istim) = starts_at_ceiling;
                

            
            %%%%%%%%%%%%%%%%% add data to trainphases table, where phases are organized by chronological order
            trainphases.learnfit{iphase}.max_acc(isub,istim) = fit_result.max_acc; 
            trainphases.learnfit{iphase}.acc_increase(isub,istim) = fit_result.acc_increase; 
            trainphases.learnfit{iphase}.exp_learn_rate(isub,istim) = fit_result.exp_learn_rate;
            trainphases.learnfit{iphase}.starts_at_ceil(isub,istim) = starts_at_ceiling;

            if op.plot_all_learning_curve_fits
                nexttile(tilenum(htl,rowcounter,isub),[1,1])
                curve_fit_plot_op.intable = phasestimtab; % accuracy table
                plot_windowed_timecourse(curve_fit_plot_op); 
                hold on
                hplot = plot(fit_result);
                htitle = title([thissub '...' thisphase], 'Interpreter','none');
            end
        end
    end
end

%% organize fit results by dbs-on vs. dbs-off
% max_acc
max_acc_on_all = con.train.learnfit{'trainon',1}.max_acc;  
      max_acc_on_all(con.train.learnfit{'trainon',1}.starts_at_ceil) = nan;  % exclude stim/phases which start at ceiling
      max_acc_on = nanmean(max_acc_on_all,2);  % average within subjects
max_acc_off_all = con.train.learnfit{'trainoff',1}.max_acc;   
      max_acc_off_all(con.train.learnfit{'trainoff',1}.starts_at_ceil) = nan; % exclude stim/phases which start at ceiling
      max_acc_off = nanmean(max_acc_off_all,2);  % average within subjects

% acc_increase
acc_increase_on_all = con.train.learnfit{'trainon',1}.acc_increase;  
      acc_increase_on_all(con.train.learnfit{'trainon',1}.starts_at_ceil) = nan;  % exclude stim/phases which start at ceiling
      acc_increase_on = nanmean(acc_increase_on_all,2);  % average within subjects
acc_increase_off_all = con.train.learnfit{'trainoff',1}.acc_increase;   
      acc_increase_off_all(con.train.learnfit{'trainoff',1}.starts_at_ceil) = nan; % exclude stim/phases which start at ceiling
      acc_increase_off = nanmean(acc_increase_off_all,2);  % average within subjects

% exp_learn_rate
%%%% when comparing learning rates, exclude negative/flat learning slope... 
%%%% .... because our interpretation of alpha assumes a positive learning slope (i.e. higher alpha indicates faster learning gains)
positive_acc_increase_on = acc_increase_on_all > 0; % determine which phase/stim had positive accuracy changes
exp_learn_rate_on_all = con.train.learnfit{'trainon',1}.exp_learn_rate;  
      exp_learn_rate_on_all(con.train.learnfit{'trainon',1}.starts_at_ceil) = nan;  % exclude stim/phases which start at ceiling
      exp_learn_rate_on_all(~positive_acc_increase_on) = nan; 
      exp_learn_rate_on = nanmean(exp_learn_rate_on_all,2);  % average within subjects
positive_acc_increase_off = acc_increase_off_all > 0; % determine which phase/stim had positive accuracy changes
exp_learn_rate_off_all = con.train.learnfit{'trainoff',1}.exp_learn_rate;   
      exp_learn_rate_off_all(con.train.learnfit{'trainoff',1}.starts_at_ceil) = nan; % exclude stim/phases which start at ceiling
      exp_learn_rate_off_all(~positive_acc_increase_off) = nan; 
      exp_learn_rate_off = nanmean(exp_learn_rate_off_all,2);  % average within subjects

[~, p_acc_increase_on_vs_off] = ttest(acc_increase_on, acc_increase_off);
[~, p_max_acc_on_vs_off] = ttest(max_acc_on, max_acc_off);
[~, p_exp_learn_rat_on_vs_off] = ttest(exp_learn_rate_on, exp_learn_rate_off);

%% organize fit results by chron. order
for iphase = 1:n_train_phases
    thisphase = trainphases.phase{iphase};

    max_acc_all = trainphases{thisphase,'learnfit'}{1}.max_acc;  
    max_acc_all(trainphases{thisphase,'learnfit'}{1}.starts_at_ceil) = nan;    % exclude stim/phases which start at ceiling
    trainphases{thisphase,'learnfit_to_plot'}{1}.max_acc = nanmean(max_acc_all,2);

    acc_increase_all = trainphases{thisphase,'learnfit'}{1}.acc_increase;  
    acc_increase_all(trainphases{thisphase,'learnfit'}{1}.starts_at_ceil) = nan;    % exclude stim/phases which start at ceiling
    trainphases{thisphase,'learnfit_to_plot'}{1}.acc_increase = nanmean(acc_increase_all,2);

    %%%% when comparing learning rates, exclude negative/flat learning slope... 
    %%%% .... because our interpretation of alpha assumes a positive learning slope (i.e. higher alpha indicates faster learning gains)
    positive_acc_increase_on = acc_increase_all > 0; % determine which phase/stim had positive accuracy changes
    exp_learn_rate_all = trainphases{thisphase,'learnfit'}{1}.exp_learn_rate;  
    exp_learn_rate_all(trainphases{thisphase,'learnfit'}{1}.starts_at_ceil) = nan;    % exclude stim/phases which start at ceiling
    exp_learn_rate_all(~positive_acc_increase_on) = nan; 
    trainphases{thisphase,'learnfit_to_plot'}{1}.exp_learn_rate = nanmean(exp_learn_rate_all,2);
end


if op.early_train_acc_on_vs_off
    hfig = figure('Color','w'); 
    hplot = plot([con.train{'trainon','acc_early'}{1}, con.train{'trainoff','acc_early'}{1}]');
        ylabel('proportion syllables accurate')
        title(['Training phases, first ' num2str(op.n_early_trials_per_stimid) ' of each stim'])
        format_dbs_on_off_plot(subs)
end

%% plot learning fit parameters in dbs-on vs. dbs-off
if op.plot_learn_params_on_vs_off
    hfig = figure('Color','w'); 
    
    subplot(1,3,1)
    plot([max_acc_on, max_acc_off]')
        title('Max accuracy')
        ylabel('Max proportion syl acc. (fitted)')    
         format_dbs_on_off_plot(subs)
    
    subplot(1,3,2)
    plot([acc_increase_on, acc_increase_off]')
        title('Training phase accuracy increases')
        ylabel('Proportion syl acc. increase (fitted)')
        format_dbs_on_off_plot(subs)
    
    subplot(1,3,3)
    plot([exp_learn_rate_on, exp_learn_rate_off]')
        title('Training learning rates')
        ylabel('Exponential learning rate')    
        format_dbs_on_off_plot(subs)

end



%% plot learning fit parameters in train-A vs. train-B while color coding dbs order
if op.plot_learn_params_A_vs_B
        hfig = figure('Color','w'); 
    
    subplot(1,3,1)
    hplot = plot(cell2mat(cellfun(@(x)x.max_acc,trainphases.learnfit_to_plot,'UniformOutput',false)')');
        title('Max accuracy')        
        ylabel('Max proportion syl acc. (fitted)')    
        format_train_AB_plot(subs,hplot)

    subplot(1,3,2)
    hplot = plot(cell2mat(cellfun(@(x)x.acc_increase,trainphases.learnfit_to_plot,'UniformOutput',false)')');
        title('Training phase accuracy increases')
        ylabel('Proportion syl acc. increase (fitted)')
        format_train_AB_plot(subs,hplot)

    subplot(1,3,3)
    hplot = plot(cell2mat(cellfun(@(x)x.exp_learn_rate,trainphases.learnfit_to_plot,'UniformOutput',false)')');
        title('Training learning rates')
        ylabel('Exponential learning rate')   
         format_train_AB_plot(subs,hplot)

end





%% subfunctions
function format_dbs_on_off_plot(subs)
        hax = gca;
        hax.XTick = [1 2]; 
        hax.XTickLabel = {'DBS On','DBS Off'};
        legend(cell2mat([subs.sub, repmat({'...'},height(subs),1), subs.train_order]),Interpreter="none")
        xlim([0.9 2.1])
        box off
end

function format_train_AB_plot(subs,hplot)
        hax = gca;
        hax.XTick = [1 2]; 
        hax.XTickLabel = {'Train A','Train B'};
        arrayfun(@(k) set(hplot(k), 'Color', subs.dbs_order_color(k,:)), 1:numel(hplot));
        legend(cell2mat([subs.sub, repmat({'...'},height(subs),1), subs.train_order]),Interpreter="none")
        xlim([0.9 2.1])
        box off
end