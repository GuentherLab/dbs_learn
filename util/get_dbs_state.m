% output dbs state op struct based subject cohort
%%% op struct must include fields 'sub','ses','task'
% do not add to op structure - this risks that it remains there when it is no longer correct

function dbs_state = get_dbs_state(op)

paths = setpaths_dbs_learn(op);
subs = readtable(paths.subject_list_master,'VariableNamingRule', 'preserve');
subs.Properties.RowNames = subs.sub; 
cohort = subs{op.sub,['cohort_',op.ses]}{1};

cohort_dbs_states = readtable(paths.cohort_dbs_states, 'FileType','text','Delimiter','tab'); 
cohort_dbs_states.Properties.RowNames = cohort_dbs_states.task; 

% if the input isn't a listed task [e.g. if it's a trial condition like 'novel1_nn'], give empty output
if ~any(op.task == string(cohort_dbs_states.task))
    dbs_state = ''; 
elseif any(op.task == string(cohort_dbs_states.task))
    dbs_state = cohort_dbs_states{op.task,['cohort_',cohort]}{1}; 
end
