%%% save audio recordings from CED computer as .wav files

clear
setpaths_dbsmulti()

phase_names = {'famil','assess','pretest','train-a','train-b','test'}; 
nphases = length(phase_names);

subtable_file = [paths.annot, filesep, 'dbsmulti_subs_master.csv']; 
 subs = readtable(subtable_file);
nsubs = height(subs); 



for isub = 1:nsubs
    cfg = [];
    cfg.sub = subs.sub{isub};
    cfg.ses = 1; 
    outfiledir = [paths.data, filesep, ['sub-',cfg.sub], filesep, ['ses-',num2str(cfg.ses)], filesep, 'audio-video']; 
    for iphase = 1:nphases
        cfg.task = phase_names{iphase};
        physiopath = [paths.data, filesep, ['sub-',cfg.sub], filesep, ['ses-',num2str(cfg.ses)], filesep, 'physio'];
        physiodirlist = struct2table(dir(physiopath)); 
        last_preproc_file_idx = find(contains(physiodirlist.name,[cfg.task,'_preproc']),1,'last'); % most updated preproc file
        
        if ~isempty(last_preproc_file_idx) % if CED and aligned Percept file (preproc file) exists
            [physdat, pcp] = load_ced_and_percept_data(cfg);
            outfilename = [outfiledir, filesep,  'sub-',cfg.sub, '_ses-',num2str(cfg.ses), '_task-',cfg.task, '_ced-audio.wav']; 
            audiowrite(outfilename, physdat.CED.voice.trial{1}, physdat.CED.voice.fsample)
        end
    end
end