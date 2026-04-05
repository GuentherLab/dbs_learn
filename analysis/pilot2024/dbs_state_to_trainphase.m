%       dbs_state = get_dbs_state(thissub, thistask, paths)
% 
%%% output the string 'trainA' or 'trainB' indicating the DBS state in the queried subject and dbs state
% dbsstate = 'on' or 'off'
%
% first run setpaths_dbsmulti, input the 'paths' variable

function trainphase = dbs_state_to_trainphase(thissub, dbs_state, paths)

subtable_file = [paths.annot, filesep, 'dbsmulti_subs_master.csv']; 
subs = readtable(subtable_file);
subrow = strcmp(subs.sub, thissub);
train_order = subs.train_order{subrow};

switch dbs_state
    case 'on',
        if train_order == "on_off"
            trainphase = 'trainA';
        elseif train_order == "off_on"
            trainphase = 'trainB';
        end
    case'off',
        if train_order == "on_off"
            trainphase = 'trainB';
        elseif train_order == "off_on"
            trainphase = 'trainA';
        end
    otherwise
        error('DBS state name not recognized (must be either ''on'' or ''off'')')
end