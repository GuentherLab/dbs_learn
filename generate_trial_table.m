function [trials, op] = generate_trial_table(op)

% Write experiment desc text files for all phases of DBS-Learn (multisyllabic and subsyllabic experiments)
% By Andrew Meier @ Guenther Lab, 2026

% NB: for training phases, where there are only 2 unique stim, it can take a few seconds to find a stim order when op.max_repeated_trials == 3
% ..... op.max_repeated_trials = 2 will likely lag infinitely
% ..... further optimization could probably speed up this randomization for low values of op.max_repeated_trials

%% initial setup and parameters

setpaths_dbs_learn()

%%% specify the number of repetitions for each stim group within each task
% no entries for famil or assses because we just present the lists once in order
stim_reps_by_task = {...

    %subsyl tasks
    'subsyl','pretest','trainA', 2;...
    'subsyl','pretest','trainB', 2;...
    'subsyl','pretest','novel1_nat', 2;...
    'subsyl','pretest','novel1_nn', 2;...
    'subsyl','pretest','novel2_nat', 2;...
    'subsyl','pretest','novel2_nn', 2;...
    
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
    'multisyl','pretest','trainA', 2;...
    'multisyl','pretest','trainB', 2;...
    'multisyl','pretest','novel1', 2;...
    'multisyl','pretest','novel2', 2;...

    'multisyl', 'trainA', 'trainA', 75;...

    'multisyl', 'trainB', 'trainB', 75;...

    'multisyl', 'test1', 'trainA', 25;... 
    'multisyl', 'test1', 'trainB', 25;... 
    'multisyl', 'test1', 'novel1', 10;... 

    'multisyl', 'test2', 'trainA', 25;... 
    'multisyl', 'test2', 'trainB', 25;... 
    'multisyl', 'test2', 'novel2', 10;... 

    };


stim_reps_by_task = cell2table(stim_reps_by_task, 'VariableNames',{'ses','task','stim_group','reps_per_unq_stim'});

op.ses = regexp(op.sestask,'^[^_]*','match'); op.ses = op.ses{1}; 
op.task = regexp(op.sestask,'(?<=_).*', 'match'); op.task = op.task{1}; 
op = rmfield(op,'sestask'); 

op.stim_master_file = [paths.code_dbs_learn, filesep, 'stim_master_', op.ses, '.tsv']; 
stim_master = readtable(op.stim_master_file,'FileType','text'); 


if op.task == "famil" || op.task == "assess" % if familiarization or assessment phase, just take the exact listed stim in order
    trials = stim_master(stim_master.stim_group == string(op.task), : );
else % for training/testing, stim list needs to be sorted / multiplied / shuffled 
    trials = table; 

    % for multisyl, keep only stim with this subject's selected syllable count (does not apply to famil/assess)
    if op.ses == "multisyl"
        stim_master = stim_master(stim_master.n_syllables == op.subj_n_syls, :); 
    end

    stim_reps_tab = stim_reps_by_task(stim_reps_by_task.ses == string(op.ses)  & ...
                                      stim_reps_by_task.task == string(op.task), :); 

    n_stim_groups = height(stim_reps_tab);
    for igroup = 1:n_stim_groups
        thisgroup = stim_reps_tab.stim_group{igroup};
        stim_in_group = stim_master(stim_master.stim_group==string(thisgroup), :); 

        % make copies of stim in this group and add it to the trials table; 
        trials = [trials; repmat(stim_in_group, stim_reps_tab.reps_per_unq_stim(igroup), 1)];
    end
    op.ntrials = height(trials); 

    % shuffle stim and make sure that we don't exceed the max allowed number of sequential repetitions of a stim
    max_reps_in_a_row = inf; 
    while max_reps_in_a_row > op.max_repeated_trials
        trials = trials(randperm(op.ntrials), :); % shuffle trials
        [~,~,stiminds] = unique(trials.name); 

        max_reps_in_a_row = 1; 
        streak = 1; 
        for itrial = 2:op.ntrials
            if stiminds(itrial)==stiminds(itrial-1)
                streak = streak+1;
                max_reps_in_a_row = max(max_reps_in_a_row, streak); 
            else
                streak = 1; 
            end
        end
    end


end

nancol = nan(op.ntrials,1);
trials = [trials, table(nancol,nancol,nancol,nancol,nancol,nancol,'VariableNames',...
    {'time_stim_vis_on','time_stim_vis_off','time_stim_aud_on','time_stim_aud_off','time_beep_on','time_beep_off'})]; 