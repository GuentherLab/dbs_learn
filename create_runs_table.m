
%% NB - we need to change this process.... instead of creating this table, we should ....
.... get the full 'steps' table from Zeyang at the start of the session - have zeyang add a 'run' column to it

%%%%% this script could then be run afterward to merge what AM has added during the run to ZY's notes to create a finalized table






% create table listing the runs containing usable experimental data (1 run per ses/task)
%%% this script just creates a table with a blank column for run which must be manually filled out afteward
%
%%% this script can be run at the beginning of an acquisition session, in which case you will fill out step ID
%%% .... and load it at each run to apply the correct step ID
%%% .... in this case, also fill in run ID during the session

%%% alternatively, if runs table was not created before/during acquisition, run afterward
%%% ... in this case, step ID will already have been assigned by some other means
%%% ... and it can be read from a trials table file and automatically added to the table

function create_runs_table(op)

field_default('op','sub','qqq');

paths = setpaths_dbs_learn(op);
paths.data = paths.data_local;

% either run on the specified session or (default) on both sessions
if isfield(op,'ses')
    seslist = {op.ses}; 
elseif ~isfield(op,'ses')
    seslist = {'subsyl','multisyl'};
end

for thisses = seslist
    op.ses = thisses{1};
    [paths] = setpaths_dbs_learn(op,paths);
    paths.src_runs_table = [paths.src_ses, filesep, 'sub-',op.sub,'_ses-',op.ses, '_runs.tsv']; 

    if ~exist(paths.src_ses,'dir')
        mkdir(paths.src_ses)
    end

    if exist(paths.src_runs_table, 'file')
        fprintf(['\nRuns table file for ses==',op.ses ,' already exists - skipping \n..........%s'],paths.src_runs_table)
    elseif ~exist(paths.landmarks_file,'file')

        switch op.ses
            case 'subsyl'
                tasks = {'famil';'pretest';'trainA';'trainB';'test1';'test2'}; 
            case 'multisyl'
                tasks = {'fds';'famil';'assess';'pretest';'trainA';'trainB';'test1';'test2'}; 
        end
    
        ntasks = length(tasks); 
        nancol = nan(ntasks,1);
        celcol = cell(ntasks,1);
        runs = table( tasks,nancol,celcol,celcol, 'VariableNames',...
                    {'task','run','step','notes'}); 
    
        for itask = 1:ntasks
            op.task = tasks{itask};

            % if data has already been acquired, get the step label by finding a trials file which matches this task
            %%% we don't need to know run number for this,...
            %%% ...  because this task always matches a specific step label within a given subject
            if exist(paths.beh,'dir')
                beh_files = {dir(paths.beh).name}';
                matchfile = beh_files{find(matches(beh_files, ...
                    ['sub-'  op.sub  '_ses-'  op.ses  '_task-'  op.task '_run-'] +...
                    wildcardPattern  + '_step-' + wildcardPattern + 'trials.tsv'),      1)}; 
                steplabel = extractBetween(matchfile,'_step-','_trials.tsv'); 
                runs.step{itask} = steplabel{1}; 
            end
    
            fprintf(['Creating empty runs table to fill in here:\n    %s,\n'],paths.src_runs_table)
            writetable(runs, paths.src_runs_table, "FileType","text", "Delimiter","tab"); 
        end
    end

end