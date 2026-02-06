% convert all trial audio files for a subject/session to wavs

setpaths_dbsmulti()

%%% INPUTs ===============================================================
subjID = 'dbs004';
sesID = 1;

basePath = paths.data;
outputPath = paths.data;

start_run_num = 1;
end_run_num = 1;


%% convert files loop
dd = struct2table(dir(sprintf('%s/sub-%s/ses-%d/beh',basePath, subjID, sesID))); 
    dd = dd(dd.isdir,:);
tasklist = dd.name(3:end);
ntasks = length(tasklist); 

for itask = 1: ntasks
    task = tasklist{itask};
    
    % % % % % mkdir(sprintf('%s/sub-%s/ses-%d/beh/%s/',outputPath, subjID, sesID, task));
    for irun = start_run_num : end_run_num
        data = load(sprintf('%s/sub-%s/ses-%d/beh/%s/sub-%s_ses-%d_run-%d_task-%s_desc-audio.mat',basePath, subjID, sesID, task, subjID, sesID, irun, task));
        desc_file = sprintf('%s/sub-%s/ses-%d/beh/%s/sub-%s_ses-%d_run-%d_task-%s_desc-stimulus.txt',basePath, subjID, sesID, task, subjID, sesID, irun, task);
        desc = regexp(fileread(desc_file),'[\n\r]+','split');
        nTrial = size(data.trialData,2);
        f = waitbar(0, sprintf('Run %d', irun));
        for j = 1:nTrial
            sti = strrep(strrep(desc{j}, '/', '-'), '\', '-');
            fname = sprintf('%s/sub-%s/ses-%d/beh/%s/sub-%s_ses-%d_run-%d_task-%s_trial-%s_%s_response.wav',outputPath, subjID, sesID, task, subjID, sesID, irun, task, pad(num2str(j), 2, 'left', '0'),sti);
            audiowrite(fname, data.trialData(j).s, data.trialData(j).fs);
            waitbar(j/nTrial, f, sprintf('Progress: %d %%', floor(j/nTrial*100)));
        end
        close(f)
    end

end