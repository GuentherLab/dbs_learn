function setup_subj_DBSMulti_pretest_train_test(subjID, session, task, run, subj_n_syls)

% Write experiment desc text files for Pretest, Train, and Test phases of DBS-Learn (multisyllabic and subsyllabic experiments)
% By Andrew Meier @ Guenther Lab, 2026

%% initial setup and parameters

vardefault('session',1);
vardefault('run',1);

%%% specify the number of repetitions for each stim group within each task
% no entries for famil or assses because we just present the lists once in order
stim_reps_by_task = {...

    %subsyl tasks
    'subsyl', 'trainA', 'trainA', 75;...

    'subsyl', 'trainB', 'trainB', 75;...

    'subsyl', 'test1',   'trainA', 20;...
    'subsyl', 'test1',   'trainB', 20;...
    'subsyl', 'test1',  'novel1_nat', 10;...
    'subsyl', 'test1',  'novel1_nn', 10;...

    'subsyl', 'test2',   'trainA' 20;...
    'subsyl', 'test2',   'trainB' 20;...
    'subsyl', 'test2',  'novel2_nat', 10;...
    'subsyl', 'test2',  'novel2_nn', 10;...


    % multisyl tasks
    'multisyl', 'trainA', 'trainA', 75;...

    'multisyl', 'trainB', 'trainB', 75;...

    'multisyl', 'test1', 'trainA', 20;... 
    'multisyl', 'test1', 'trainB', 20;... 
    'multisyl', 'test1', 'novel1', 20;... 

    'multisyl', 'test2', 'trainA', 20;... 
    'multisyl', 'test2', 'trainB', 20;... 
    'multisyl', 'test2', 'novel2', 20;... 

    }



% ntrials_pretest_per_unq_stim = 2; % number of times to repeat each unique training-phase stim during the pretest
% ntrials_train = 100; % total trials in each training phase (divided between trained stim)
% ntrials_test = 150; % total trials in each testing phase (divided between trained stim and test stim)

% switch sestask
%     case 'subsyl_famil'
% 
%     case 'sub_pretest'
% 
%     case 'sub_trainA'
% 
%     case 'sub_trainB'
% 
%     case 'sub_test1'
% 
%     case 'sub_test2'
% 
%     case 'multi_famil'
% 
%     case 'multi_assess'
% 
%     case 'multisyl_pretest'
% 
%     case 'multi_trainA'
% 
%     case 'multi_trainB'
% 
%     case 'multi_test1'
% 
%     case 'multi_test2'
% 
% end

op.ses = regexp(task,'^[^_]*','match'); 
op.task = regexp(sestas,'(?<=_).*', 'match'); 

op.stim_master_file = [paths.code_dbs_multi, filesep, 'stim_master_', op.ses, '.tsv']; 
stim_master = readtable('op.stim_master_file','FileType','text'); 


if op.task == "famil" || op.task == "assess" % if familiarization or assessment phase, just take the exact listed stim in order
    stimtab = stim_master(stim_master.stim_group == op.task, : )
else % for training/testing, stim list needs to be sorted / multiplied / shuffled 
    
    % for multisyl, keep only stim with this subject's selected syllable count (does not apply to famil/assess)
    if op.ses == "multisyl"
        stim_master = stim_master(stim_master.n_syllables == op.subj_n_syls, :); 
    end

    



end



%% pretest
taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','pretest');
if ~isfolder(taskpath); mkdir(taskpath); end
stim_thistask = stim_master(contains(stim_master.stim_group,{'test_a';'test_b'}),:); % pretest the train-phase stim, not the novel stim

% select the stim set with the number of syllables selected for this subject
stim_thistask = stim_thistask(stim_thistask.n_syllables==subj_n_syls,:); 

pretest_tab = repmat(stim_thistask, ntrials_pretest_per_unq_stim, 1); 
    ntrials_pretest = height(pretest_tab);
pretest_tab = pretest_tab(randperm(ntrials_pretest),:); % shuffle
spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-pretest_trials.csv',subjID,session,run)); 

%%%% check for overwriting
if exist(spath,'file')
    overwrite_ok = input('Stimuli for this subject/session/run already exists; overwrite? (Enter y if OK to delete)','s');
    if strcmp(overwrite_ok,'y')
        fprintf(['Overwriting pretest/train/test stim for subject ', subjID, ', session ', num2str(session) , ', run ', num2str(run), '\n'])
    else 
        error('Not overwriting stim files; quitting setup_subj script')
    end
end

writetable(pretest_tab,spath) % save table with condition labels

pretest_words = pretest_tab.name;

spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-pretest_desc-stimulus.txt',subjID,session,run));
writetable(table(pretest_words),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions

%% train a
taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','train-a');
if ~isfolder(taskpath); mkdir(taskpath); end

stim_thistask = stim_master(strcmp(stim_master.stim_group,'trainA'),:); 

% select the stim set with the number of syllables selected for this subject
stim_thistask = stim_thistask(stim_thistask.n_syllables==subj_n_syls,:); 

trainwords_1_tab = repmat(stim_thistask,ntrials_train,1);
trainwords_1_tab = trainwords_1_tab(1:ntrials_train,:); 
trainwords_1_tab = trainwords_1_tab(randperm(ntrials_train),:); % shuffle
spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-train-a_trials.csv',subjID,session,run)); 
writetable(trainwords_1_tab,spath) % save table with condition labels

trainwords_1 = trainwords_1_tab.name; 

spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-train-a_desc-stimulus.txt',subjID,session,run));
writetable(table(trainwords_1),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions




%% train b
taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','train-b');
if ~isfolder(taskpath); mkdir(taskpath); end

stim_thistask = stim_master(strcmp(stim_master.stim_group,'trainB'),:); 


% select the stim set with the number of syllables selected for this subject
stim_thistask = stim_thistask(stim_thistask.n_syllables==subj_n_syls,:); 

trainwords_2_tab = repmat(stim_thistask,ntrials_train,1);
trainwords_2_tab = trainwords_2_tab(1:ntrials_train,:); 
trainwords_2_tab = trainwords_2_tab(randperm(ntrials_train),:); % shuffle
spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-train-b_trials.csv',subjID,session,run)); 
writetable(trainwords_2_tab,spath) % save table with condition labels

trainwords_2 = trainwords_2_tab.name; 

spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-train-b_desc-stimulus.txt',subjID,session,run));
writetable(table(trainwords_2),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions




%% test
 
%% update to include 2 test sessions
%% update to work differently for multi and subsyl

taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','test');
if ~isfolder(taskpath); mkdir(taskpath); end

stimrowmatch = strcmp(stim_master.stim_group,'trainA') | strcmp(stim_master.stim_group,'trainB') | strcmp(stim_master.stim_group,'novel'); 
stim_thistask = stim_master(stimrowmatch,:); 

% select the stim set with the number of syllables selected for this subject

stim_thistask = stim_thistask(stim_thistask.n_syllables==subj_n_syls,:); 

testwords_tab = repmat(stim_thistask,ntrials_test,1);
testwords_tab = testwords_tab(1:ntrials_test,:); 
testwords_tab = testwords_tab(randperm(ntrials_test),:); % shuffle

testwords = testwords_tab.name; 

spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-test_trials.csv',subjID,session,run)); 
writetable(table(testwords),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions

testwords = testwords_tab.name; 
spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-test_desc-stimulus.txt',subjID,session,run));
writecell(testwords,spath)



