%%% add dash after 'trial' in trial audio files

cd('C:\docs\code\dbs_learn') 
paths = setpaths_dbs_learn()

target_folder = 'C:\Dropbox\R01-SML_data_shared\derivatives\sml001\trial-audio'; 
% target_folder = 'C:\Dropbox\R01-SML_data_shared\derivatives\sml002\trial-audio'; 

cd(target_folder)

% 2. Get all contents and filter for subfolders
items = dir(target_folder);
dirFlags = [items.isdir] & ~ismember({items.name}, {'.', '..'});
subFolders = items(dirFlags);

% 3. Loop through and execute
for k = 1:length(subFolders)
    thisFolderName = [target_folder, filesep, subFolders(k).name];
    
    fprintf('Renaming in: %s\n', thisFolderName);
    
    % % Call your function
    batch_rename(thisFolderName, 'trial', 'trial-');
end