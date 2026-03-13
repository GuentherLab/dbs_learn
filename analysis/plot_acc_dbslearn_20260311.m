cd('C:\Dropbox\R01-SML_data_shared\derivatives\sml001\annot')

%% test-off subsyl - cohort c

close all

t = readtable('sub-sml001_ses-subsyl_task-test1_run-56_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t.acc = [isnan(t.onset_error) | t.onset_error==7] & strcmp(t.rime_error,'');
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title('subsyllabic accuracy during Testing-DBS-ON - subject 1')
hax.XTickLabels = xticklabels;
ylim([0.5 1])
box off

%% test-on subsyl - cohort c
close all

t = readtable('sub-sml001_ses-subsyl_task-test2_run-7_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t.acc = [isnan(t.onset_error) | t.onset_error==7] & isnan(t.rime_error);
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-ON-nn','train-OFF-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title('subsyllabic accuracy during Testing-DBS-OFF - subject 1')
hax.XTickLabels = xticklabels;
ylim([0.5 1])
box off

%% test-on multisyl - cohort a
close all
clear

t = readtable('sub-sml001_ses-multisyl_task-test1_run-12_trials_accuracy.tsv','FileType','text','Delimiter','tab')
t = t(1:19,:);
t.acc = t.n_correct_syls ./ t.n_syllables;
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'train-ON','train-OFF','novel'};
ylabel({'proportion correct syllables'})
title('multisyllabic accuracy during Testing-DBS-ON - subject 1')
hax.XTickLabels = xticklabels;
ylim([0 1])
box off