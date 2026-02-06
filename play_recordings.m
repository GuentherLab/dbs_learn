% Script to play recording from one single trial
% Note that you should not load the combined result file (*audio.mat)

% Load the trial result you want to listen

%     'trial-70.mat']);
% load(['C:\dbsmulti\sub-qqq\ses-1\beh\famil\sub-qqq_ses-1_run-1_task-famil_',...
%     'trial-1.mat']);
% load(['C:\dbsmulti\sub-qqq\ses-1\beh\famil\sub-qqq_ses-1_run-1_task-famil_',...
%     'trial-2.mat']);
% load(['C:\dbsmulti\sub-qqq\ses-1\beh\test\\sub-qqq_ses-1_run-1_task-test_',...
%     'trial-1.mat']);
% load(['C:\dbsmulti\sub-dbs999\ses-1\beh\famil\sub-dbs999_ses-1_run-1_task-famil_',...
%     'trial-1.mat']);
% load(['C:\dbsmulti\sub-dbs003\ses-1\beh\famil\sub-dbs003_ses-1_run-1_task-famil_',...
%     'trial-1.mat']);
% load(['C:\dbsmulti\sub-dbs003\ses-1\beh\pretest\sub-dbs003_ses-1_run-1_task-pretest_',...
%     'trial-1.mat']);
% load(['C:\dbsmulti\sub-dbs003\ses-1\beh\train-a\sub-dbs003_ses-1_run-1_task-train-a_',...
%     'trial-99.mat']);
load(['C:\dbsmulti\sub-dbs003\ses-1\beh\train-b\sub-dbs003_ses-1_run-1_task-train-b_',...
    'trial-99.mat']);


info=audiodevinfo; 
player=audioplayer(tData.s, tData.fs, 24, info.output(1).ID);
play(player)