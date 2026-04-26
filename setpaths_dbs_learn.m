%%% set matlab paths for DBS-Multisyllabic experiment 

function [paths, compname] = setpaths_dbs_learn(op,paths)

% if a paths variable wasn't provided, create a new one
if ~exist('paths','var')
    paths = struct; 
end

compname = getenv('COMPUTERNAME');

switch compname 
    case '677-GUE-WL-0012' %%% thinkpad to stay in dbs-learn experiment room
        paths.code = 'C:\code';
    case {'677-GUE-WL-0010', 'AMSMEIER'} % AM thinkpad, strix 
        paths.code = 'C:\docs\code'; 
    case 'ORDI'
        paths.code = 'C:\Users\lwolf\github_repos'; 
end

paths.data_local  = 'C:\dbs_learn_data'; % local location to save when running experiment
paths.data_remote = 'C:\Dropbox\R01-SML_data_shared'; % remote data - upload here after a session and analyze from here
    paths.groupanalysis = [paths.data_remote, filesep, 'groupanalysis']; 


if ~isfield(paths,'data')
    paths.data = paths.data_remote;
end
    

paths.code_dbs_learn = [paths.code, filesep, 'dbs_learn']; 
paths.code_analysis = [paths.code_dbs_learn, filesep, 'analysis'];
paths.annot = [paths.code_analysis, filesep, 'annot']; % small annotation files
paths.stim = [paths.code_dbs_learn, filesep, 'stimuli'];

% external toolboxes
paths.spm = [paths.code, filesep, 'spm12']; %%%% 
paths.fieldtrip_toolbox = [paths.code, filesep, 'fieldtrip']; % only used in analysis, not running experiment
paths.bml = [paths.code, filesep, 'bml']; % RM Richardson lab toolbox

% paths to specific reference files
paths.subject_list_master = [paths.data_remote, filesep, 'subject_list_master.xlsx'];
paths.cohort_dbs_states = [paths.data_remote, filesep, 'cohort_dbs_states.tsv']; 

 paths_to_add =  {paths.code_dbs_learn;...
     [paths.code_analysis];...
     [paths.annot];...
     [paths.code_dbs_learn, filesep, 'util'];...
    [paths.stim];...
    paths.spm;... %%%% use version of spm in fieldtrip
    paths.fieldtrip_toolbox;...
};

 addpath(paths_to_add{:}); clear paths_to_add

addpath(genpath([paths.code_dbs_learn, filesep, 'NIMH_daqtoolbox_(Apr-7-2016)'])) % needed for sending beacon

set(0,'defaultTextInterpreter','none') 

%% if more details are provided, output relevant paths

if exist('op','var')
    if isfield(op,'sub') % if sub specified
        paths.src_sub = [paths.data, filesep,'sourcedata', filesep, op.sub]; 
        paths.der_sub = [paths.data, filesep, 'derivatives', filesep, op.sub]; 
        paths.annot = [paths.der_sub, filesep, 'annot']; 
        if isfield(op,'ses') % if session specified
            paths.src_ses = [paths.src_sub, filesep,'ses-',op.ses]; 
            paths.beh = [paths.src_ses, filesep, 'beh'];
            paths.src_audvid = [paths.src_ses, filesep, 'audio-video'];
            paths.landmarks_file = [paths.annot filesep, 'sub-',op.sub, '_ses-',op.ses,  '_sync-landmarks.tsv']; 
            paths.trial_audio = [paths.der_sub, filesep, 'trial-audio']; 
            paths.src_runs_table = [paths.src_ses, filesep, 'sub-',op.sub,'_ses-',op.ses, '_runs.tsv']; 
            if isfield(op,'task')
                if isfield(op,'run')

                    %%%% the following string gets used in a variety of files associated with this run
                    % filestr_step includes step number - use only for sourcedata - useful for Zeyang preprocessing
                    % filestr omits step number for all further processing as it's redundant with session+task
                    paths.filestr = ['sub-',op.sub, '_ses-',op.ses, '_task-',op.task, '_run-',num2str(op.run), '_']; 
                    if isfield(op,'step')
                        paths.filestr_step = [paths.filestr 'step-',op.step '_']; 
                    elseif ~isfield(op,'step') && isfield(paths,'filestr_step')
                        paths = rmfield(paths,'filestr_step'); % remove unless it's clearly been added to op
                    end
                    
                    % trial boundary adjustments - gets created during create_sync_landmark_tables.m
                    % .... gets used during audio/video trial cutting
                    %%% note that these adjustments will not change any subsequent statistical analyses....
                    %%% ... they change where the trial-wise audio/video files will be cut to make it easier to do annotations
                    %%% generally this will be useful for when the sub answers early and you want the file to start earlier
                    paths.trialfile_boundary_adjustments = [paths.trial_audio, filesep, paths.filestr, 'boundary_adjustments.tsv'];

                end
            end

               
            
        end
    end
end

