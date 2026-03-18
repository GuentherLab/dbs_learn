

op.sub = 'qqq';
op.ses = 'multisyl'; 
op.is_dbs_run = 1; 
op.record_audio = 1; 

%% famil
op.task='famil'; 
op.dbs_state = 'on';
op.step_id = 'A101';
dbs_learn_run_exp(op)

%% asses
op.task='asses'; 
op.dbs_state = 'on';
op.step_id = 'A201';
dbs_learn_run_exp(op)

%% set number of syls
op.subj_n_syls = 0;

%% pretest
op.task='pretest'; 
op.dbs_state = 'on';
op.step_id = 'A301';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step_id = 'A401';
dbs_learn_run_exp(op)

%====================================== 
%% Washout-1
%% trainA
op.task='trainA'; 
op.dbs_state = 'on';
op.step_id = 'B201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step_id = 'B301';
dbs_learn_run_exp(op)

%====================================== washout starts

%% free base
op.task='base-free'; 
op.dbs_state = 'off';
op.step_id = 'C001';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step_id = 'C401';
dbs_learn_run_exp(op)

%% trainB
op.task='trainB'; 
op.dbs_state = 'off';
op.step_id = 'C201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step_id = 'C301';
dbs_learn_run_exp(op)



%====================================== washout starts

%% free base
op.task='base-free'; 
op.dbs_state = 'on';
op.step_id = 'D001';
dbs_learn_run_exp(op)

%% test1
op.task='test1'; 
op.dbs_state = 'on';
op.step_id = 'D201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step_id = 'D301';
dbs_learn_run_exp(op)

%======================================
%% test2
op.task='test2'; 
op.dbs_state = 'off';
op.step_id = 'E201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step_id = 'E301';
dbs_learn_run_exp(op)


%% TELL
op.task='reading'; 
op.dbs_state = 'off';
op.step_id = 'K501';
dbs_learn_run_exp(op)