%%% create table for filling in manual landmarks to be used for syncing for a subject
% need to create and fill out runs table for this session first - create_runs_table.m

function create_sync_landmark_tables(op)

vardefault('op',struct);
field_default('op','sub','sml002');
field_default('op','ses','subsyl');

[paths, compname] = setpaths_dbs_learn(op);

if exist(paths.landmarks_file,'file')
    fprintf(['Landmarks file already exists - aborting \n..........',paths.landmarks_file])
    return
end

runs = readtable(paths.src_runs_table,'Filetype','text', 'Delimiter','tab'); 
tasks = runs.task; 

filetypes = {'trials.tsv','trials.tsv';...
    'audio.wav','recording-headphone.wav';...
    'preproc.mat','';.... % need to get expected filenames - probably include step ID 
    'video','video';... %%%%% need to get correct filenames and types - gopro
    }; 
nfiletypes = size(filetypes, 1); 
nruns = length(tasks);
nrows = nruns * nfiletypes;
nancol = nan(nrows,1);
celcol = cell(nrows,1); 
runs_col = repelem(runs.run,nfiletypes,1);
tasks_col = repelem(tasks,nfiletypes,1);
step_col = repelem(runs.step,nfiletypes,1);
filetype_col = repmat(filetypes,nruns,1);

sync = table(tasks_col, runs_col, step_col,  filetype_col,celcol,celcol,nancol,nancol,celcol,celcol,  'VariableNames',...
            {'task',	'run', 'step',	'filetype',	'dir','filename'	't1'	't2'	't1_description','t2_description'});
sync.t1_description = repmat({'first go beep onset'},nrows,1);
sync.t2_description = repmat({'last go beep onset'},nrows,1);


for irow = 1:nrows
    op.task = sync.task{irow};
    op.run = sync.run(irow); 
    runrowmatch = string(op.task) == runs.task & op.run == runs.run;
    op.step = runs.step{runrowmatch};
    paths = setpaths_dbs_learn(op); % get paths.filestr

    switch sync.filetype{irow} % set dir 
        case 'trials.tsv'
            sync.dir{irow} = paths.beh;
            sync.filename{irow} = [paths.filestr_step, 'trials.tsv'];
            sync.t1_description{irow} = 'trials.t_go_aud_on(1)..... first trial go beep onset'; 
            sync.t2_description{irow} = 'trials.t_go_aud_on(end)..... last trial go beep onset'; 

            % open trial table to add go beep datapoints
            trials = readtable([sync.dir{irow}, filesep, sync.filename{irow}], 'Filetype','text', 'Delimiter','tab');
            sync.t1(irow) = trials.t_go_aud_on(1);
            sync.t2(irow) = trials.t_go_aud_on(end); 

        case 'audio.wav'
            sync.dir{irow} = paths.src_audvid; 
            sync.filename{irow} = [paths.filestr,'recording-headphone.wav'];
        case 'preproc.mat'
            % need to specify name 

        case 'video'

    end
end

% save new file
if ~exist(paths.landmarks_file,'file')
    if ~exist(paths.der_sub,'dir')
        mkdir(paths.der_sub)
    end
    writetable(sync,paths.landmarks_file,'FileType','text','Delimiter','tab')
    fprintf(['writing sync table to manually complete: \n     ',paths.landmarks_file])
end

