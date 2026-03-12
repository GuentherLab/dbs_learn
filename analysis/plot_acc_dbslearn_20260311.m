 t = readtable('sub-sml001_ses-subsyl_task-test1_run-56_trials-accuracy.tsv','FileType','text','Delimiter','tab')
t.acc = [isnan(t.onset_error) | t.onset_error==7] & strcmp(t.rime_error,'');
g = grpstats(t,"stim_group","mean","DataVars",["acc"]);

close all
hbar = bar(g.mean_acc);
hax = gca;

set(0, 'DefaultTextInterpreter', 'none')
xticklabels = {'novel-nn','novel-nat','train-OFF-nn','train-ON-nn'};
ylabel({'proportion correct', '(ignoring epenthesis)'})
title('subsyllabic accuracy during Testing-DBS-ON - subject 1')
hax.XTickLabels = xticklabels;
ylim([0.5 1])
box off