
%% plot train and test timecourses
% sml002 subsyl
op.sub = 'sml002';
op.ses = 'subsyl'; 


tasks = {'trainA','trainB','test1','test2'};
ntasks = length(tasks);

op.frac_trials_to_plot = 1; 

paths = setpaths_dbs_learn(op);
runs = readtable(paths.src_runs_table,'FileType','text','Delimiter','tab');

close all
hfig = figure('Color','w'); 

for itask = 1:ntasks
    op.task = tasks{itask};
    op.run = runs.run(strcmp(runs.task,op.task));
    paths = setpaths_dbs_learn(op); % get run-specific info

    if any(strcmp(op.task, {'trainA','trainB'}))
        
        op.sortvar = 'name'; 
    elseif any(strcmp(op.task, {'test1','test2'}))
        op.sortvar = 'stim_group';
    end

    op.newfig = 0; 
    subplot(2,2,itask)
    dbs_learn_timecourse(op,paths)



end


%% multisyl - sml001 and sml002 

frac_trials_to_plot = .25;

close all hidden
hfig  = figure('Color','w')

subplot(2,2,1)
t = readtable('sub-sml001_ses-multisyl_task-test1_run-12_trials_accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = t.n_correct_syls ./ t.n_syllables;
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'train-ON','train-OFF','novel'};
ylabel({'proportion correct syllables'})
title({'multisyllabic accuracy during Testing-DBS-ON - subject 1',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off



%% subsyl
frac_trials_to_plot = .25;

close all
hfig = figure;

subplot(2,2,1)
t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-subsyl_task-test2_run-7_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.rime_error);
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);
g = sortrows(g, 'stim_group')

% close all
% hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-ON - subject 1',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off





subplot(2,2,2)
t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-subsyl_task-test1_run-56_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = [isnan(t.onset_error) | t.onset_error==7] & strcmp(t.rime_error,'');
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);
g = sortrows(g, 'stim_group')

% hfig = figure; 
% close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 1',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0.5 1])
box off



t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot\sub-sml002_ses-subsyl_task-test1_run-9_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = sum([isnan(t.onset_error) | t.onset_error==7, isnan(t.coda_error) | t.coda_error==7], 2) / 2; 
t.acc_binary = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.coda_error) 
g = grpstats(t,"stim_group","mean","DataVars",["acc","acc_binary"]);
g = sortrows(g, 'stim_group')

subplot(2,2,3)
hbar = bar(g.mean_acc);
hax = gca;
xticklabels = {,'train-ON-nn','novel-nn','train-OFF-nn', 'novel-nat'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-ON - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off


subplot(2,2,4)
t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot\sub-sml002_ses-subsyl_task-test2_run-11_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = sum([isnan(t.onset_error) | t.onset_error==7, isnan(t.coda_error) | t.coda_error==7], 2) / 2; 
t.acc_binary = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.coda_error) 
g = grpstats(t,"stim_group","mean","DataVars",["acc","acc_binary"]);
g = sortrows(g, 'stim_group')

hbar = bar(g.mean_acc);
hax = gca;
xticklabels = {,'train-OFF-nn','novel-nat','novel-nn', 'train-ON-nn'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off













function dbs_learn_timecourse(op,paths)

    trials_acc = readtable([paths.annot, filesep, paths.filestr, 'trials-accuracy.tsv'],'FileType','text','Delimiter','tab');

    dbs_state = get_dbs_state(op); 

    field_default('op','sortvar','stim_group'); 
    if strcmp(op.ses,'multisyl')
        trials_acc.prop_correct = trials_acc.n_correct_syls ./ trials_acc.n_syllables;

         
    elseif strcmp(op.ses,'subsyl')
        trials_acc.prop_correct = sum([isnan(trials_acc.onset_error) | trials_acc.onset_error==7, isnan(trials_acc.coda_error) | trials_acc.coda_error==7], 2) / 2; 

    end

    op.intable = trials_acc;
    op.plotvar = 'prop_correct';
    outtab = plot_windowed_timecourse(op); 
    
    hlgd = findobj(gcf, 'Type', 'Legend');
    % hlgd.String = {'novel','train on', 'train off'};
    
    title([op.task,'_dbs-',dbs_state],'Color','k')


end