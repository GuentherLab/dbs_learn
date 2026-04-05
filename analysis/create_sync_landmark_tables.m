%%% create table for filling in manual landmarks to be used for syncing for a subject
% need to create and fill out runs table for this session first - create_runs_table.m

function create_sync_landmark_tables(op)

op.sub = 'sml001';
op.ses = 'subsyl'; 


[paths, compname] = setpaths_dbs_learn(op);



if exist(paths.landmarks_file,'file')
    fprintf(['Landmarks file already exists - aborting \n..........',paths.landmarks_file])
    return
end


replace cmted out part with loading paths.src_runs_table
% % % % % switch op.ses
% % % % %     case 'subsyl'
% % % % %         tasks = {'famil','pretest','trainA','trainB','test1','test2'}; 
% % % % %     case 'multisyl'
% % % % %         tasks = {'fds','famil','assess','pretest','trainA','trainB','test2'}; 
% % % % % end
% % % % % runs = readtable([paths.annot, filesep, 'sub-',op.sub,'_ses-',op.ses, '_runs.tsv']); 





filetypes = {'trials.tsv','trials.tsv';...
    'audio.wav','recording-headphone.wav';...
    'preproc.mat','';....
    'video','video';... %%%%% need to get correct filenames and types - gopro
    }; 
nfiletypes = size(filetypes, 1); 
% % % % % ntasks = length(tasks);
nrows = ntasks * nfiletypes;
nancol = nan(nrows,1);
celcol = cell(nrows,1); 
tasks_col = repelem(tasks',nfiletypes,1);
filetype_col = repmat(filetypes,ntasks,1);

sync = table(tasks_col, nancol, celcol,  filetype_col,celcol,celcol,nancol,nancol,celcol,celcol,  'VariableNames',...
            {'task',	'run', 'step',	'filetype',	'dir','filename'	't1'	't2'	't1_description','t2_description'});
sync.t1_description = repmat({'first go beep onset'},nrows,1);
sync.t2_description = repmat({'last go beep onset'},nrows,1);


for irow = 1:nrows
    op.task = sync.task{irow};
    op.run = sync.run{irow}; 
    paths = setpaths_dbs_learn(op); % get paths.filestr

    switch sync.filetype{irow} % set dir 
        case 'trials.tsv'
            sync.dir{irow} = paths.beh;
            sync.filename{irow} = [filestr_step, 'trials.tsv']

            % open trial file table to add 
            trials = readtable([sync.dir{irow}, filesep, sync.filename{irow}, 'Filetype','text', 'Delimiter','tab']);
            % sync.

        case 'audio.wav'
            sync.dir{irow} = [paths.src_audvid, filesep, paths.filestr,'recording-headphone.wav'];
            sync.filename{irow} = paths.data_ses_audvid; 
        case 'preproc.mat'
            % need to specify name 

        case 'video'

    end
end

if ~exist(exist(paths.landmarks_file,'file'))
% save new file
end

