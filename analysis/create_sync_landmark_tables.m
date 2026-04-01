%%% create table for filling in manual landmarks to be used for syncing for a subject


op.sub = 'sml001'


%% subsyl 
tasks = {'famil','pretest','trainA','trainB','test1','test2'}; 
filetypes = {'trials.tsv','trials.tsv';...
    'audio.wav','recording-headphone.wav';...
    'preproc.mat','';....
    'video','video';... %%%%% need to get correct filenames and types - gopro
    }; 
nfiletypes = size(filetypes, 1); 
ntasks = length(tasks);
nrows = ntasks *
nancol = nan(nrows,1);
celcol = cell(nrows,1); 

% specify that landmark_1 should be 'first go beep onset', landmark_2 'last go beep onset'

sync = table(tasks', nancol,celcol,celcol,celcol,nancol,nancol,celcol,  'VariableNames',...
    {'task',	'run',	'filetype',	'dir','name'	't1'	't2'	'landmark_description'});

for irow = 1:nrows
    switch sync.filetype % set dir 
        case 'trials.tsv'

        case 'audio.wav'

        case 'preproc.mat'
            % need to specify name 

        case 'video'

    end
end

