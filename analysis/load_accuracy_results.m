 %% load and process accuracy data

warning('off','MATLAB:table:RowsAddedNewVars')

 setpaths_dbsmulti() 


vardefault('op',struct);
field_default('op','exclude_no_response_trials',1);
field_default('op','normalize_by_pretest',1);
field_default('op','dbs_order_colors',[0.6 0.6 0;... % for plotting: ON-OFF, OFF-ON
                                         0 0 0.7]   ); 
field_default('op','n_early_trials_per_stimid', 5); 

subtable_file = [paths.annot, filesep, 'dbsmulti_subs_master.csv']; 

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
subs = [subs, table(subnan3,subcel,subcel,subcel,subcel,subcel,'RowNames',subs.sub,'VariableNames',...
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

      % if test phase, label trials by whether the stim is trained-on, trained-off, or novel
        if string(thisphase) == "test";
            thisphasetab.traincon = cell(height(thisphasetab),1);
            novel_rows = string(thisphasetab.stim_group) == "novel"; 
            thisphasetab.traincon(novel_rows,1) = {'novel'};
            thisphasetab.traincon(~novel_rows,1) = cellfun(@(x)get_dbs_state(thissub,x,paths), thisphasetab.stim_group(~novel_rows,1), 'UniformOutput', false); % DBS state of trained stim
        end

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
