%%%%% sml 002

% cd('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot')
cd('C:\Users\amsmeier\Dropbox\R01-SML_data_shared\derivatives\sml001\annot')

%% test-on subsyl - cohort c

frac_trials_to_plot = 1;

t = readtable('sub-sml001_ses-subsyl_task-test1_run-56_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = [isnan(t.onset_error) | t.onset_error==7] & strcmp(t.rime_error,'');
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

hfig = figure; 
close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 1',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0.5 1])
box off

%% test-off subsyl - cohort c

frac_trials_to_plot = 0.25;

t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-subsyl_task-test2_run-7_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.rime_error);
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

% close all
hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-ON - subject 1',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

%% sml001 test-on multisyl - cohort a
% close all

clear

frac_trials_to_plot = 0.25;

t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot\sub-sml002_ses-multisyl_task-test1_run-14_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = t.n_correct_syls ./ t.n_syllables;
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'train-ON','train-OFF','novel'};
ylabel({'proportion correct syllables'})
title({'multisyllabic accuracy during Testing-DBS-ON - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

%% sml002 test-on multisyl - cohort a
% close all

clear

frac_trials_to_plot = 0.25;

t = readtable('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot\sub-sml001_ses-multisyl_task-test2_run-12_trials_accuracy.tsv','FileType','text','Delimiter','tab')
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



%% sml002
cd('C:\Dropbox\R01-SML_data_shared\derivatives\sml002\annot')
set(0, 'DefaultTextInterpreter', 'none')

%% test-on subsyl 

frac_trials_to_plot = .25;

t = readtable('sub-sml002_ses-subsyl_task-test1_run-9_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = sum([isnan(t.onset_error) | t.onset_error==7, isnan(t.coda_error) | t.coda_error==7], 2) / 2; 
t.acc_binary = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.coda_error) 
g = grpstats(t,"stim_group","mean","DataVars",["acc","acc_binary"]);

% close all
hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;
xticklabels = {,'train-ON-nn','novel-nn','train-OFF-nn', 'novel-nat'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-ON - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

hfig = figure; 
hbar = bar(g.mean_acc_binary);
hax = gca;
xticklabels = {,'train-ON-nn','novel-nn','train-OFF-nn', 'novel-nat'};
ylabel({'proportion of trials correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-ON - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

%% test-off subsyl 

frac_trials_to_plot = .25;

t = readtable('sub-sml002_ses-_task-test2_run-11_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:round(frac_trials_to_plot*height(t)),:);
t.acc = sum([isnan(t.onset_error) | t.onset_error==7, isnan(t.coda_error) | t.coda_error==7], 2) / 2; 
t.acc_binary = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.coda_error) 
g = grpstats(t,"stim_group","mean","DataVars",["acc","acc_binary"]);

close all
hfig = figure; 
hbar = bar(g.mean_acc);
hax = gca;
xticklabels = {,'train-OFF-nn','novel-nat','novel-nn', 'train-ON-nn'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

hfig = figure; 
hbar = bar(g.mean_acc_binary);
hax = gca;
xticklabels = {,'train-OFF-nn','novel-nat','novel-nn', 'train-ON-nn'};
ylabel({'proportion of clusters correct', '(ignoring epenthesis)'})
title({'subsyllabic accuracy during Testing-DBS-OFF - subject 2',['first ', num2str(100*frac_trials_to_plot), 'pct of trials']})
hax.XTickLabels = xticklabels;
ylim([0 1])
box off

