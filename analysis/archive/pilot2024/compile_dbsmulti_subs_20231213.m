%%% compile accuracy data from pilot subjects 

clear

datapath = 'C:\dbsmulti';

performance_var = 'syl_acc_prop';  %% proportion of correct syllables within each stim
% performance_var = 'correct'; %%% logical - no errors vs. any errors

cond_windowsize = 10; 
stimname_windowsize = 4; 


subtab = table({'pilot_001';'pilot_003';'pilot_004';'pilot_005'},[1;1;1;1],[1;1;1;1],{'ab';'ba';'ab';'ba'},...
    'VariableNames',{'sub',                         'sess','run',     'trainorder'});

nsubs = height(subtab);
subtab.acc_trained = nan(nsubs,1);
subtab.acc_novel = nan(nsubs,1);
subtab.acc_trained_1 = nan(nsubs,1); % accuracy of the first trained stim group
subtab.acc_trained_2 = nan(nsubs,1); % accuracy of the second trained stim group
subtab.condplot = cell(nsubs,1);
subtab.nameplot = cell(nsubs,1);



%% import and compile subjects
for isub = 1:nsubs
    thissub = subtab.sub{isub};
    thissess = num2str(subtab.sess(isub)); 
    thisrun = num2str(subtab.sess(isub)); 
    trials_test = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-test_acc.csv']);

    % look for recency effects; sort subs by which stim set they trained on first/second
    if 'ab' == string(subtab.trainorder{isub})
        subtab.acc_trained_1(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_a'}),performance_var});
        subtab.acc_trained_2(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_b'}),performance_var});
        trials_train_1 = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-train-a_acc.csv']);
        trials_train_2 = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-train-b_acc.csv']);
    elseif 'ba' == string(subtab.trainorder{isub})
        subtab.acc_trained_1(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_b'}),performance_var});
        subtab.acc_trained_2(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_a'}),performance_var});
        trials_train_1 = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-train-b_acc.csv']);
        trials_train_2 = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-train-a_acc.csv']);
    end
    trials_train = [trials_train_1; trials_train_2]; 
  
    % accuracy and learning effects
    subtab.acc_trained(isub) =  mean(trials_test{contains(trials_test.stim_group,{'test_a','test_b'}),performance_var});
    subtab.acc_novel(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_novel'}),performance_var});

    %% organize data into windowed accuracy timecourses, sorted by condition or stimulus ID
    % condition
    unq_cons = {'test_a';'test_b';'test_novel'};
    ncons = length(unq_cons);
    conacc = cell(ncons,1); 
    for icond = 1:ncons
        thisname = unq_cons{icond};
        temptrials = trials_test(strcmp(trials_test.stim_group,thisname),:);
        ncontrials = height(temptrials);
    
        xvals = cond_windowsize:ncontrials; % index of trial relative to all trials within this condition
        ntrials_to_plot = length(xvals);
        conacc{icond} = nan(ntrials_to_plot,1); 
        for itrial = 1:ntrials_to_plot
            trial_con_relative = xvals(itrial); 
            conacc{icond}(itrial) = mean(temptrials{trial_con_relative-cond_windowsize+1:trial_con_relative, performance_var});
        end
    end
    subtab.condplot{isub} = table(unq_cons, conacc,'VariableNames',{'cond','trialacc'});
        subtab.condplot{isub}.ntrials = cellfun(@length,subtab.condplot{isub}.trialacc);
    
    % stim name
    unq_names = unique(trials_test.name); 
    nnames = length(unq_names);
    nameacc = cell(nnames,1); 
    subtab.nameplot{isub} = table(unq_names,cell(nnames,1),nan(nnames,3),'VariableNames',{'name','stim_group','color'});
        clear unq_names
    %%%%%%% sort by condition
    for iname = 1:nnames
        thisname = subtab.nameplot{isub}.name{iname};
        subtab.nameplot{isub}.stim_group{iname} = trials_test.stim_group{find(strcmp(trials_test.name,thisname),1)};
    end
    subtab.nameplot{isub} = sortrows(subtab.nameplot{isub},'stim_group');
    
    for iname = 1:nnames
        thisname = subtab.nameplot{isub}.name{iname};
    
        temptrials = trials_test(strcmp(trials_test.name,thisname),:);
        nnametrials = height(temptrials);
    
        xvals = stimname_windowsize:nnametrials; % index of trial relative to all trials of this name
        ntrials_to_plot = length(xvals);
        nameacc{iname} = nan(ntrials_to_plot,1); 
        for itrial = 1:ntrials_to_plot
            trial_con_relative = xvals(itrial); 
            nameacc{iname}(itrial) = mean(temptrials{(trial_con_relative-stimname_windowsize+1:trial_con_relative), performance_var});
        end
    end
    subtab.nameplot{isub}.trialacc = nameacc; 
        subtab.nameplot{isub}.ntrials = cellfun(@length,subtab.nameplot{isub}.trialacc);
end

allsubs_cond = table(subtab.condplot{1}.cond,cell(3,1),cell(3,1),cell(3,1),'VariableNames',{'label','trialacc','acc_mean','acc_sem'});
allsubs_name = table(subtab.nameplot{1}.name,cell(3,1),cell(3,1),cell(3,1),'VariableNames',{'label','trialacc','acc_mean','acc_sem'});


switch performance_var
    case 'correct'
        ylab = 'Accuracy (full stimulus)';
    case 'syl_acc_prop'
        ylab = 'Proportion of accurate syllables';
end
set(0, 'DefaultTextInterpreter', 'none')

%% look for recency effect
% close all
% linewidth = 2; 
% 
% hfig = figure('Color','w');
% yvals_to_plot = [subtab.acc_trained_1, subtab.acc_trained_2]';
% hplot = plot(repmat([1;2],1,nsubs),yvals_to_plot,'LineWidth',linewidth);
% box('off')
% legend(subtab.sub)
% ylim([min(yvals_to_plot(:))-0.05,1.05])
% xlim([0.9 2.1])
% hax = gca;
%     hax.XTick = [1 2];
%     hax.XTickLabel = {'First trained set','Second traineed set'};
%     hax.FontWeight = 'bold';
% ylabel(ylab)
% title({'Accuracy of 1st vs. 2nd trained stim sets', 'during Testing phase'})
% 
% [~, p_recency] = ttest(subtab.acc_trained_1, subtab.acc_trained_2); 
% hanot = annotation('textbox','String',['p=',num2str(p_recency)],'FitBoxToText','on');
%     hanot.Position(1:2) = [0.5 0.15];


%% plot timecourse accuracies by stimname and syllable, compiled from all subjects



%% subfunction for processing trialtable files
function tabout = import_trialtable(filename)
    tabout = readtable(filename);
    tabout = tabout(~(tabout.unusableTrial==1 & ~isnan(tabout.unusableTrial==1)),:);
    
    tabout.syl_acc_prop = tabout.n_correct_syls ./ tabout.n_syllables; 
    
    ntrials = height(tabout);
    if ~iscell(tabout.notes(1)) % if notes are empty, they might be numeric; convert to cell
        tabout.notes = cell(ntrials,1);
    end
end

%% subfunction for compiling accuracies across subjects into plot-able matrices


