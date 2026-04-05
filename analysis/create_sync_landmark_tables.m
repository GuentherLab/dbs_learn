%%% create table for filling in manual landmarks to be used for syncing for a subject

function create_sync_landmark_tables(op)

op.sub = 'sml001';
op.ses = 'subsyl'; 


[paths, compname] = setpaths_dbs_learn;

paths.aud_record_filename_headphone = [paths.src_audvid, filesep, filestr,'recording-headphone.wav']
paths.aud_record_filename_mic = [paths.data_audvid, filesep, filestr,'recording-mic.wav'];

if exist(paths.landmarks_file,'file')
    fprintf(['Landmarks file already exist - aborting \n..........',paths.landmarks_file])
end

switch op.ses
    case 'subsyl'
        tasks = {'famil','pretest','trainA','trainB','test1','test2'}; 
    case 'multisyl'
        tasks = {'fds','famil','assess','pretest','trainA','trainB','test2'}; 
end
filetypes = {'trials.tsv','trials.tsv';...
    'audio.wav','recording-headphone.wav';...
    'preproc.mat','';....
    'video','video';... %%%%% need to get correct filenames and types - gopro
    }; 
nfiletypes = size(filetypes, 1); 
ntasks = length(tasks);
nrows = ntasks * nfiletypes;
nancol = nan(nrows,1);
celcol = cell(nrows,1); 
tasks_col = repelem(tasks',nfiletypes,1);
filetype_col = repmat(filetypes,ntasks,1);

sync = table(tasks_col, nancol,filetype_col,celcol,celcol,nancol,nancol,celcol,celcol,  'VariableNames',...
    {'task',	'run',	'filetype',	'dir','filename'	't1'	't2'	't1_description','t2_description'});
sync.t1_description = repmat({'first go beep onset'},nrows,1);
sync.t2_description = repmat({'last go beep onset'},nrows,1);


for irow = 1:nrows
    switch sync.filetype{irow} % set dir 
        case 'trials.tsv'
            sync.dir{irow} = paths.beh;
            sync.filename{irow} = [filestr_step, 'trials.tsv']

            % open trial file table to add 
            trials = readtable([sync.dir{irow}, filesep, sync.filename{irow}, 'Filetype','text', 'Delimiter','tab']);
            % sync.

        case 'audio.wav'
            sync.dir{irow} = paths.aud_record_filename_headphone;
            sync.filename{irow} = paths.data_ses_audvid; 
        case 'preproc.mat'
            % need to specify name 

        case 'video'

    end
end

if ~exist(exist(paths.landmarks_file,'file'))
% save new file
end

