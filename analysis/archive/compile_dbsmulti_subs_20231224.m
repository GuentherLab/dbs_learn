%%% compile accuracy data from pilot subjects 
% includes options to plot results across subjects

% clear

datapath = 'C:\dbsmulti';

performance_var = 'syl_acc_prop';  %% proportion of correct syllables within each stim
% performance_var = 'correct'; %%% logical - no errors vs. any errors

vardefault('cond_windowsize',10); 
vardefault('stimname_windowsize',4); 

% when restricting analysis to only the early trials of an experimental phase, use this many trials of each label
%%%% note that this is the number of early trials of each LABEL (condition or stimname), not the total number of early trials across labels
vardefault('n_early_trials_cond',7);
vardefault('n_early_trials_stimname',4);

vardefault('recency_plot_individual_subs',false);
vardefault('plot_allsubs_timecourse_condition',false);
vardefault('plot_allsubs_timecourse_stimname',false);

subtab = table({'pilot_001';'pilot_003';'pilot_004';'pilot_005'},[1;1;1;1],[1;1;1;1],{'ab';'ba';'ab';'ba'},...
    'VariableNames',{'sub',                         'sess',         'run',                'trainorder',   });

nsubs = height(subtab);
subtab.acc_trained = nan(nsubs,1);
subtab.acc_novel = nan(nsubs,1);
subtab.acc_trained_1 = nan(nsubs,1); % accuracy of the first trained stim group
subtab.acc_trained_2 = nan(nsubs,1); % accuracy of the second trained stim group
subtab.condplot = cell(nsubs,1);
subtab.nameplot = cell(nsubs,1);
subtab.trials_trainphase_a = cell(nsubs,1); 
subtab.trials_trainphase_b = cell(nsubs,1); 



%% import and compile subjects
for isub = 1:nsubs
    thissub = subtab.sub{isub};
    thissess = num2str(subtab.sess(isub)); 
    thisrun = num2str(subtab.sess(isub)); 
    trials_test = import_trialtable([datapath filesep 'sub-' thissub filesep...
                    'sub-' thissub '_ses-' thissess '_run-' thisrun '_task-test_acc.csv']);

    % look for recency effects; sort subs by which stim set they trained on first/second
    %%% also extract training-phase trials for each stim group, including "early" training trials
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
    traingroups = {'a','b'}; ntraingroups = length(traingroups);
    for igroup = 1:ntraingroups
        thisgroup = traingroups{igroup};
        subtab{isub,['trials_trainphase_',thisgroup]} = {trials_train(strcmp(trials_train.stim_group,['test_',thisgroup]), :)};
    end
  
    % accuracy and learning effects
    subtab.acc_trained(isub) =  mean(trials_test{contains(trials_test.stim_group,{'test_a','test_b'}),performance_var});
    subtab.acc_novel(isub) = mean(trials_test{contains(trials_test.stim_group,{'test_novel'}),performance_var});


    %% organize data into windowed accuracy timecourses, sorted by condition or stimulus ID
    % condition
    unq_cons = {'test_a';'test_b';'test_novel'};
    nconds = length(unq_cons);
    conacc = cell(nconds,1); 
    for icond = 1:nconds
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
    subtab.condplot{isub} = table(repmat({thissub},nconds,1),unq_cons, conacc,'VariableNames',{'sub','label','trialacc'});
        subtab.condplot{isub}.ntrials = cellfun(@length,subtab.condplot{isub}.trialacc);
    
    % stim name
    unq_names = unique(trials_test.name); 
    nnames = length(unq_names);
    nameacc = cell(nnames,1); 
    subtab.nameplot{isub} = table(repmat({thissub},nnames,1),unq_names,cell(nnames,1),nan(nnames,3),'VariableNames',{'sub','label','stim_group','color'});

        clear unq_names
    %%%%%%% sort by condition
    for iname = 1:nnames
        thisname = subtab.nameplot{isub}.label{iname};
        subtab.nameplot{isub}.stim_group{iname} = trials_test.stim_group{find(strcmp(trials_test.name,thisname),1)};
    end
    subtab.nameplot{isub} = sortrows(subtab.nameplot{isub},'stim_group');
    
    for iname = 1:nnames
        thisname = subtab.nameplot{isub}.label{iname};
    
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

switch performance_var
    case 'correct'
        ylab = 'Accuracy (full stimulus)';
    case 'syl_acc_prop'
        ylab = 'Proportion of accurate syllables';
end
set(0, 'DefaultTextInterpreter', 'none')

% condition - average across subs
ops = []; 
ops.varname = 'condplot';
ops.n_early_trials = n_early_trials_cond; 

if plot_allsubs_timecourse_condition
    ops.linecolor = [  0.6 0.6 1;...       % test a
                                    0 0.7 0.4;...   % test b
                                    1 0.1 0.6];      % test novel 
    ops.windowsize = cond_windowsize; 
    ops.xjittermax = 0.5;
    ops.yjittermax = 0; 
    ops.show_figure = true; 
    allsubs_cond = compile_timecourses(subtab,ops);
        hylab = ylabel (ylab);
        hxlab = xlabel('trial (within-condition)');
else
    ops.show_figure = false; 
    allsubs_cond = compile_timecourses(subtab,ops);
end

% stimname - average across subs
ops = []; 
ops.varname = 'nameplot';
ops.n_early_trials = n_early_trials_stimname; 

    % early trainphase trials
    ops.trialtab_varname = 'name';
    ops.performance_var = performance_var; 
    subtab = early_trainphase_acc(subtab, ops); 

if plot_allsubs_timecourse_stimname
    ops.linecolor = [  0.6 0.6 1;...       % test a
                                    0.1 0.1 0.7;...   % test a
                                    0.6 1 0;...       % test b
                                    0 0.7 0.4;...   % test b
                                    1 0 0;...       % test novel
                                    1 0.1 0.6];      % test novel    
    ops.windowsize = stimname_windowsize; 
    ops.xjittermax = 0.5;
    ops.yjittermax = 0; 
    ops.show_figure = true; 
    allsubs_name = compile_timecourses(subtab,ops);
        hylab = ylabel (ylab);
        hxlab = xlabel('trial (within-stim-name)');

else
    ops.show_figure = false; 
    allsubs_name = compile_timecourses(subtab,ops);
end
allsubs_name.stim_group = {'test_a';'test_a';'test_b';'test_b';'test_novel';'test_novel'};


%% look for recency effect
if recency_plot_individual_subs

    % close all
    linewidth = 2; 
    
    hfig = figure('Color','w');
    yvals_to_plot = [subtab.acc_trained_1, subtab.acc_trained_2]';
    hplot = plot(repmat([1;2],1,nsubs),yvals_to_plot,'LineWidth',linewidth);
    box('off')
    legend(subtab.sub)
    ylim([min(yvals_to_plot(:))-0.05,1.05])
    xlim([0.9 2.1])
    hax = gca;
        hax.XTick = [1 2];
        hax.XTickLabel = {'First trained set','Second traineed set'};
        hax.FontWeight = 'bold';
    ylabel(ylab)
    title({'Accuracy of 1st vs. 2nd trained stim sets', 'during Testing phase'})
    
    [~, p_recency] = ttest(subtab.acc_trained_1, subtab.acc_trained_2); 
    hanot = annotation('textbox','String',['p=',num2str(p_recency)],'FitBoxToText','on');
        hanot.Position(1:2) = [0.5 0.15];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction for plotting timecourse accuracies by stimname and syllable, compiled from all subjects
% compile accuracies across subjects into plot-able matrices
%
% also compute mean 'early trial' accuracy

function tabout = compile_timecourses(subtab,ops)

field_default('ops','linecolor',[0 0 0]);
field_default('ops','show_figure',true);
field_default('ops','linewidth',2); 

% close all


nsubs = height(subtab);
labels = subtab{1,ops.varname}{1}.label;
nlabels = length(labels); 
tabout = table(labels,cell(nlabels,1),cell(nlabels,1),cell(nlabels,1),cell(nlabels,1),cell(nlabels,1),...
 'VariableNames',{'label', 'trialacc',   'acc_mean', 'acc_sem',   'early_acc_submean', 'ntrials',       });


for ilabel = 1:nlabels
    for isub = 1:nsubs
        tabout.ntrials{ilabel}(isub) = subtab{isub,ops.varname}{1}.ntrials(ilabel);
    end
    
    ntrials_to_use = min(tabout.ntrials{ilabel});

    for isub = 1:nsubs
        tabout.trialacc{ilabel}(isub,1:ntrials_to_use) = subtab{isub,ops.varname}{1}.trialacc{ilabel}(1:ntrials_to_use);
    end
    
    tabout.acc_mean{ilabel} = mean(tabout.trialacc{ilabel});
    tabout.acc_sem{ilabel} = std(tabout.trialacc{ilabel}) ./ sqrt(nsubs);
    tabout.early_acc_submean{ilabel} = mean(tabout.trialacc{ilabel}(:,1:ops.n_early_trials),2); 
end

if ops.show_figure
    hfig = figure('Color','w');
   hold on
    for ilabel = 1:nlabels
        ntrials_to_use = min(tabout.ntrials{ilabel});
        % plot timecourse for this label
        xvals = [1:ntrials_to_use] + ops.windowsize; % index of trial relative to all trials within this condition
            xvals = xvals + [rand(size(xvals))-0.5]*ops.xjittermax;        
        yvals = tabout.acc_mean{ilabel};
            yvals = yvals + [rand(size(yvals))-0.5]*ops.yjittermax;
        hplot = plot(xvals,yvals,'LineWidth',ops.linewidth,'Color',ops.linecolor(ilabel,:));
        herb = errorbar(xvals,yvals,tabout.acc_sem{ilabel},'HandleVisibility','off','Color',ops.linecolor(ilabel,:));
    
    end
    
    ileg = legend(labels,'Interpreter','none');
    % title(['Test phase accuracy.... trial-window size = ' num2str(ops.windowsize)])
    title(['Test phase accuracy'])
end


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% subfunction for computing early-trial accuracy during training phases
function subtab_out = early_trainphase_acc(subtab_in, ops)
    subtab_out = subtab_in;
    nsubs = height(subtab_out);
    
    % get indices of the first $ trials for each unique label
    traingroups = {'a','b'}; ntraingroups = length(traingroups);
    for isub = 1:nsubs
        for igroup = 1:ntraingroups % for each training-phase stimulus group
            thisgroup = traingroups{igroup};
            trialtab =  subtab_out{isub,['trials_trainphase_',thisgroup]}{1}; 
            labels = unique(trialtab{:,ops.trialtab_varname});
            % % % % % nlabels = length(labels); 
            early_inds = cell2mat(arrayfun(@(x) find(string(trialtab{:,ops.trialtab_varname})==x,ops.n_early_trials),labels,'un',0))';
            early_trials_acc = trialtab{early_inds,ops.performance_var}; 
            subtab_out{isub,['acc_trainphase_early_',thisgroup]} = mean(early_trials_acc); 
        end
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
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




