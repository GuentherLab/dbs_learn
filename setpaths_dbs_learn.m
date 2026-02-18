%%% set matlab paths for DBS-Multisyllabic experiment 


%% ADD MONKEYLOGIC TO DBSLEARN REPO AND ADD TO THIS SCRIPT

paths.data  = 'C:\dbs_learn_data'; 
    paths.groupanalysis = [paths.data, filesep, 'groupanalysis']; 
paths.code = 'C:\docs\code'; 
    paths.code_dbs_learn = [paths.code, filesep, 'dbs_learn']; 
    paths.code_analysis = [paths.code_dbs_learn, filesep, 'analysis'];
    paths.annot = [paths.code_analysis, filesep, 'annot']; % small annotation files
    paths.stim = [paths.code_dbs_learn, filesep, 'stimuli'];

    % external toolboxes
    paths.spm = [paths.code, filesep, 'spm12']; %%%% use version of spm in fieldtrip
    paths.fieldtrip_toolbox = [paths.code, filesep, 'fieldtrip']; % only used in analysis, not running experiment
    paths.bml = [paths.code, filesep, 'bml']; % RM Richardson lab toolbox

if ~isfolder(paths.data)
    mkdir(paths.data)
end

 paths_to_add =  {paths.code_dbs_learn;...
     [paths.code_analysis];...
     [paths.annot];...
     [paths.code_dbs_learn, filesep, 'util'];...
    [paths.stim];...
    paths.spm;... %%%% use version of spm in fieldtrip
    % paths.bml;...
    paths.fieldtrip_toolbox;...
};

 addpath(paths_to_add{:}); clear paths_to_add

addpath(genpath([paths.code_dbs_learn, filesep, 'NIMH_daqtoolbox_(Apr-7-2016)'])) % needed for sending beacon

% bml_defaults()
% ft_defaults()


set(0,'defaulttextInterpreter','none') 

 clear paths_to_add