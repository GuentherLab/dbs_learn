

op.sub = 'qqq';
op.ses = 'subsyl'; 
op.is_dbs_run = 1; 
op.record_audio = 1; 
op.syl_struct = 'ccvcc';

%% famil
op.task='famil'; 
op.dbs_state = 'on';
op.step = 'A101';
dbs_learn_run_exp(op)

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

%====================================== 
%% Washout-1
%% trainA
op.task='trainA'; 
op.dbs_state = 'off';
op.step = 'B201';
dbs_learn_run_exp(op)

%% ctr base - Percept survey mode
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'B301';
dbs_learn_run_exp(op)

%% ctr base - Percept sense mode
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'B401';
dbs_learn_run_exp(op)
%====================================== washout_02 starts

%% trainB
op.task='trainB'; 
op.dbs_state = 'on';
op.step = 'C201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'on';
op.step = 'C301';
dbs_learn_run_exp(op)

%====================================== washout03 starts



%% test1
op.task='test1'; 
op.dbs_state = 'off';
op.step = 'D201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'D301';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'D401';
dbs_learn_run_exp(op)

%====================================== washout04 starts
%% test2
op.task='test2'; 
op.dbs_state = 'on';
op.step = 'E201';
dbs_learn_run_exp(op)

%% ctr base
op.task='base-controlled'; 
op.dbs_state = 'off';
op.step = 'E301';
dbs_learn_run_exp(op)


