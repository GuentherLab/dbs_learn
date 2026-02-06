%       dbs_state = get_dbs_state(thissub, thistask, paths)
% 
%%% output the string 'off' or 'on' indicating the DBS state in the queried subject and task
% first run setpaths_dbsmulti, input the 'paths' variable)
% 
% assumes that DBS is ON for all tasks except for either trainA or trainB

function dbs_state = get_dbs_state(thissub, thistask, paths)

subtable_file = [paths.annot, filesep, 'dbsmulti_subs_master.csv']; 
subs = readtable(subtable_file);
subrow = strcmp(subs.sub, thissub);
train_order = subs.train_order{subrow};

switch thistask
    case  'assess',
        dbs_state = 'on';
    case 'pretest',
        dbs_state = 'on';
    case 'trainA',
        if train_order == "on_off"
            dbs_state = 'on';
        elseif train_order == "off_on"
            dbs_state = 'off';
        end
    case'trainB',
        if train_order == "on_off"
            dbs_state = 'off';
        elseif train_order == "off_on"
            dbs_state = 'on';
        end
    case'test'
        dbs_state = 'on';
    otherwise
        error('experimental phase name not recognized')
end