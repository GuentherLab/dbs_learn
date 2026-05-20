
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
trials = readtable('sub-sml001_ses-multisyl_task-test1_run-12_trials_accuracy.tsv','FileType','text','Delimiter','tab')
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = trials.n_correct_syls ./ trials.n_syllables;
g = grpstats(trials,"stim_group","mean","DataVars",["acc"]);

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
trials = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-subsyl_task-test2_run-7_trials-accuracy.tsv','FileType','text','Delimiter','tab')
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = [isnan(trials.onset_error) | trials.onset_error==7] & isnan(trials.rime_error);
g = grpstats(trials,"stim_group","mean","DataVars",["acc"]);
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
trials = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-subsyl_task-test1_run-56_trials-accuracy.tsv','FileType','text','Delimiter','tab')
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = [isnan(trials.onset_error) | trials.onset_error==7] & strcmp(trials.rime_error,'');
g = grpstats(trials,"stim_group","mean","DataVars",["acc"]);
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



trials = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot\sub-sml002_ses-subsyl_task-test1_run-9_trials-accuracy.tsv','FileType','text','Delimiter','tab')
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = sum([isnan(trials.onset_error) | trials.onset_error==7, isnan(trials.coda_error) | trials.coda_error==7], 2) / 2; 
trials.acc_binary = [isnan(trials.onset_error) | trials.onset_error==7] & isnan(trials.coda_error) 
g = grpstats(trials,"stim_group","mean","DataVars",["acc","acc_binary"]);
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
trials = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot\sub-sml002_ses-subsyl_task-test2_run-11_trials-accuracy.tsv','FileType','text','Delimiter','tab')
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = sum([isnan(trials.onset_error) | trials.onset_error==7, isnan(trials.coda_error) | trials.coda_error==7], 2) / 2; 
trials.acc_binary = [isnan(trials.onset_error) | trials.onset_error==7] & isnan(trials.coda_error) 
g = grpstats(trials,"stim_group","mean","DataVars",["acc","acc_binary"]);
g = sortrows(g, 'stim_group')

hbar = bar(g.mean_acc);
hax = gca;
xticklabels = {,'train-OFF-nn','novel-nat','novel-nn', 'train-ON-nn'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off




%% sml003 subsyl
% scored by Logan
op.sub = 'sml003';
op.ses = 'subsyl';
op.task = 'test1';

frac_trials_to_plot = .5;

paths = setpaths_dbs_learn(op); 

subs = readtable(paths.subject_list_master,'VariableNamingRule', 'preserve'); subinfo = subs(op.sub==string(subs.sub), :); 

close all hidden
hfig  = figure('Color','w'); 

close all
hfig = figure;

subplot(2,2,1)
trials = readtable(paths.beh_annot_table); 
trials.stim_group = strrep(strrep(trials.stim_group,'novel1','novel'),'novel2','novel');
trials = trials(1:round(frac_trials_to_plot*height(trials)),:);
trials.acc = [trials.cons1_accuracy + trials.cons2_accuracy] / 4; 
cond_order = {'novel_nat','trainA','trainB','novel_nn'}; 
g = grpstats(trials,"stim_group","mean","DataVars",["acc"]);
g = g(cond_order,:); 
nconds = height(g);

% look up the DBS state that each condition was learned in [not the state it's being tested in]
for icond = 1:nconds
    thiscond = g.stim_group{icond};
    cfg = [];    cfg.sub = op.sub;    cfg.ses = op.ses;     cfg.task = thiscond; 
    g.dbs_state_learn{icond} = get_dbs_state(cfg);
    if ~isempty(g.dbs_state_learn{icond})
        g.xlab{icond} = [g.stim_group{icond}, '_DBS-', upper(g.dbs_state_learn{icond})];
    else
        g.xlab{icond} = g.stim_group{icond}; 
    end
end

% close all
% hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;


ylabel({'proportion correct', '(ignoring epenthesis)'})
title({[op.ses, ' accuracy during ', op.task, ' DBS-', upper(get_dbs_state(op)), ' - subject ', op.sub, ' (',upper(subinfo.dx{1}),')'],...
    ['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = g.xlab;
ylim([0 1])
box off



%%


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