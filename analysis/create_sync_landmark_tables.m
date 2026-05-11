%%% create table for filling in manual landmarks to be used for syncing for a subject
% need to create and fill out runs table for this session first - create_runs_table.m
% 
% also make a table for manually altering where audio/video trials will be cut

function create_sync_landmark_tables(op)

vardefault('op',struct);
field_default('op','sub','sml003');

field_default('op','ses','subsyl');
    % field_default('op','ses','multisyl');

[paths, compname] = setpaths_dbs_learn(op);

runs = readtable(paths.src_runs_table,'Filetype','text', 'Delimiter','tab'); 

switch op.ses
    case 'subsyl'
        tasks = {'famil';'pretest';'trainA';'trainB';'test1';'test2'}; 
    case 'multisyl'
        tasks = {'fds';'famil';'assess';'pretest';'trainA';'trainB';'test1';'test2'}; 
end

filetypes_suffix = {'trials.tsv','trials.tsv';...
    'audio.wav','recording-headphone.wav';...
    'preproc.mat','';.... % need to get expected filenames - probably include step ID 
    'audio.wav','recording-gopro.mp4';... %%%%% list the audio stripped from gopro video, because video won't be stored in shared dropbox
    }; 
nfiletypes = size(filetypes_suffix, 1); 
nruns = length(tasks);
nrows = nruns * nfiletypes;
nancol = nan(nrows,1);
celcol = cell(nrows,1); 
runs_col =  nan(nrows,1);
tasks_col = repelem(tasks,nfiletypes,1);
step_col =  nan(nrows,1);
filetype_col = repmat(filetypes_suffix(:,1),nruns,1);
suffix_col = repmat(filetypes_suffix(:,2),nruns,1);

sync = table(tasks_col, runs_col, step_col,  filetype_col,suffix_col,celcol,celcol,nancol,nancol,celcol,celcol,  'VariableNames',...
            {'task',	'run', 'step',	'filetype','suffix',	'dir','filename'	't1'	't2'	't1_description','t2_description'});
sync.t1_description = repmat({'first go beep onset'},nrows,1);
sync.t2_description = repmat({'last go beep onset'},nrows,1);


for isynctabrow = 1:nrows
    op.task = sync.task{isynctabrow};
    % % % % % op.run = sync.run(isynctabrow); 
    runrowmatch = string(op.task) == runs.task %%%%% & op.run == runs.run;

    if  ismember('step', sync.Properties.VariableNames) && iscell(sync.step) % in subject sml001, we didn't add 'step' to sourcedata filenames
        op.step = runs.step{runrowmatch};
    end

    op.run = runs.run(runrowmatch); 
    paths = setpaths_dbs_learn(op); % get paths.filestr

    switch sync.filetype{isynctabrow} % set dir 
        case 'trials.tsv'
            sync.dir{isynctabrow} = paths.beh;

            if  ismember('step', sync.Properties.VariableNames) && iscell(sync.step) % in subject sml001, we didn't add 'step' to sourcedata filenames            
                    sync.filename{isynctabrow} = [paths.filestr_step, 'trials.tsv'];
            else
                sync.filename{isynctabrow} = [paths.filestr, 'trials.tsv'];
            end

            sync.t1_description{isynctabrow} = 'trials.t_go_aud_on(1)..... first trial go beep onset'; 

            % open trial table to add GO beep datapoints
            trials = readtable([sync.dir{isynctabrow}, filesep, sync.filename{isynctabrow}], 'Filetype','text', 'Delimiter','tab');
            sync.t1(isynctabrow) = trials.t_go_aud_on(1);

            %% trial boundary adjustments
            % %    create boundary adjustments trial file to fill in if it doesn't exist yet
            % 
            %%% note that these adjustments will not change any subsequent statistical analyses....
            %%% ... they change where the trial-wise audio/video files will be cut to make it easier to do annotations
            %%% generally this will be useful for when the sub answers early and you want the file to start earlier
            if ~exist(paths.trialfile_boundary_adjustments, 'file')
                adjust_start = zeros(size(trials,1),1); 
                adjust_end = zeros(size(trials,1),1); 
                trials_bound_adj = [table(adjust_start,adjust_end), trials]; 
                writetable(trials_bound_adj, paths.trialfile_boundary_adjustments, 'FileType','text','Delimiter','tab'); 
            end


            %%
            % for FDS task we do not usually expect to complete all trials, but it has to get to trial 3, because we wait for 3 error trials
            %%% other tasks are expected to complete the final trial in the trial table
            switch op.task
                case 'fds'
                    mintrials_fds = 3; 
                    sync.t2_description{isynctabrow} = ['trials.t_go_aud_on(',num2str(mintrials_fds),')..... trial ',num2str(mintrials_fds),' go beep onset']; 
                    sync.t2(isynctabrow) = trials.t_go_aud_on(mintrials_fds); 
                otherwise
                    sync.t2_description{isynctabrow} = 'trials.t_go_aud_on(end)..... last trial go beep onset'; 
                    sync.t2(isynctabrow) = trials.t_go_aud_on(end); 
            end

        case 'audio.wav'
            sync.dir{isynctabrow} = paths.src_audvid; 
            sync.filename{isynctabrow} = [paths.filestr,'recording-headphone.wav'];
        case 'preproc.mat'
            % need to specify name 

        case 'video'
            sync.dir{irow} = paths.src_gopro; 
            sync.filename{irow} = [paths.filestr,'recording-video.wav'];
    end
end

% save new landmarks sync file - check to avoid overwriting
if ~exist(paths.landmarks_file,'file')
    if ~exist(paths.der_sub,'dir')
        mkdir(paths.der_sub)
    end
    writetable(sync,paths.landmarks_file,'FileType','text','Delimiter','tab')
    fprintf(['Landmarks file already exists - aborting \n..........%s'],paths.landmarks_file)
elseif exist(paths.landmarks_file,'file')
    fprintf(['Landmarks file already exists - not overwriting \n..........%s \n'],paths.landmarks_file)
end



