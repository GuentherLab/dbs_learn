%% 2024 pilot version

function setup_subj_DBSMulti_fam_assess(subjID, session, stimset, run)

% Write experiment desc text files for Familiarization and Assessment phases of SEQ experiment - DBS-Multisyllabic branch.
% By Andrew Meier @ Guenther Lab, 2023

vardefault('session',1);
vardefault('run',1);

projpath = 'C:\dbs_learn';

stim_master_file = ['stim_master_', stimset, '.tsv'; 
stim_master = readtable('stim_master_multisyl.tsv','FileType','text'); 



%% familiarization
taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','famil');
if ~isfolder(taskpath); mkdir(taskpath); end

famil_words = stim_master.name(strcmp(stim_master.stim_group,'familiarization'),:); % use all assessment stim; no repetitions or randomization

spath = fullfile(taskpath, sprintf('sub-%s_ses-%d_run-%d_task-famil_desc-stimulus.txt',subjID,session, run));
writetable(table(famil_words),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions
    
if strcmp(stimset, 'multisyl') % subsyl has no assessment phase

    %% assessment
    taskpath = fullfile(projpath, sprintf('sub-%s', subjID),sprintf('ses-%d', 1),'beh','assess');
    if ~isfolder(taskpath); mkdir(taskpath); end
    
    
    assess_words = stim_master.name(strcmp(stim_master.stim_group,'assessment'),:); % use all assessment stim; no repetitions or randomization
    
    spath = fullfile(taskpath, sprintf('sub-%s_ses-1_run-1_task-assess_desc-stimulus.txt',subjID));
    writetable(table(assess_words),spath,'WriteVariableNames',false) % writetable rather than writecell for use w/ older matlab versions
end









end