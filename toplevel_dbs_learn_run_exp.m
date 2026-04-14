

op.sub = 'qqq';
op.ses = 'multisyl'; 
op.is_dbs_run = 0; 
op.record_audio = 0; 

%% fds
op.task='fds'; 
op.dbs_state = 'on';
op.step = 'K401';
dbs_learn_run_exp(op)

%% famil
op.task='famil'; 
op.dbs_state = 'on';
op.step = 'A101';
dbs_learn_run_exp(op)

%% asses
op.task='asses'; 
op.dbs_state = 'on';
op.step = 'A201';
dbs_learn_run_exp(op)

%% set number of syls
op.subj_n_syls = 0;

%% pretest
op.task='pretest'; 
op.dbs_state = 'on';
op.step = 'A301';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step = 'A401';
dbs_learn_run_exp(op)

%====================================== Washout_01
%% trainA
op.task='trainA'; 
op.dbs_state = 'on';
op.step = 'B201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step = 'B301';
dbs_learn_run_exp(op)

%====================================== washout_02 starts

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'C401';
dbs_learn_run_exp(op)

%% trainB
op.task='trainB'; 
op.dbs_state = 'off';
op.step = 'C201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'C301';
dbs_learn_run_exp(op)



%====================================== washout_03 starts

%% test1
op.task='test1'; 
op.dbs_state = 'on';
op.step = 'D201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step = 'D301';
dbs_learn_run_exp(op)

%====================================== washout_04 starts
%% test2
op.task='test2'; 
op.dbs_state = 'off';
op.step = 'E201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'E301';
dbs_learn_run_exp(op)
