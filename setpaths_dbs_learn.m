%%% set matlab paths for DBS-Multisyllabic experiment 

paths.data  = 'C:\dbsmulti'; 
    paths.groupanalysis = [paths.data, filesep, 'groupanalysis']; 
paths.code = 'C:\docs\code'; 
    paths.code_dbsmulti = [paths.code, filesep, 'dbs_learn']; 
    paths.code_analysis = [paths.code_dbsmulti, filesep, 'analysis'];
    paths.annot = [paths.code_analysis, filesep, 'annot']; % small annotation files
    paths.config = [paths.code_dbsmulti, filesep, 'config']; 
    paths.spm = [paths.code, filesep, 'spm12']; %%%% use version of spm in fieldtrip
    paths.fieldtrip_toolbox = [paths.code, filesep, 'fieldtrip']; % only used in analysis, not running experiment
    paths.bml = [paths.code, filesep, 'bml']; % RM Richardson lab toolbox

if ~isfolder(paths.data)
    mkdir(paths.data)
end

 paths_to_add =  {paths.code_dbsmulti;...
     [paths.code_analysis];...
     [paths.annot];...
     [paths.code_dbsmulti, filesep, 'util'];...
    [paths.code_dbsmulti, filesep, 'stimuli'];...
    [paths.code_dbsmulti, filesep, 'config'];...
    [paths.code_dbsmulti, filesep, 'NIMH_daqtoolbox_(Apr-7-2016)'];... % needed for sending beacon
    paths.spm;... %%%% use version of spm in fieldtrip
    paths.bml;...
    paths.fieldtrip_toolbox;...
};

 addpath(paths_to_add{:})

bml_defaults()
% ft_defaults()

addpath('C:\Users\amsmeier\Downloads\NIMH_MonkeyLogic_2.2')


set(0,'defaulttextInterpreter','none') 

 clear paths_to_add