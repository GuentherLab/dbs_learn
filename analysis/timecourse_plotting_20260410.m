%%% plot timecourse of accuracy

op.sub = 'sml002';
op.ses = 'multisyl'; 
op.run = 14;
op.task = 'test1';

paths = setpaths_dbs_learn(op);

trials_acc = readtable([paths.annot, filesep, paths.filestr, 'trials-accuracy.tsv'],'FileType','text','Delimiter','tab');
trials_acc.prop_correct = trials_acc.n_correct_syls ./ trials_acc.n_syllables;

op.intable = trials_acc;
op.plotvar = 'prop_correct';
op.sortvar = 'stim_group'; 
outtab = plot_windowed_timecourse(op); 

hlgd = findobj(gcf, 'Type', 'Legend');
hlgd.String = {'novel','train on', 'train off'};

title('test off','Color','k')


%%

op.sub = 'sml002';
op.ses = 'multisyl'; 
op.run = 15;
op.task = 'test2';

paths = setpaths_dbs_learn(op);

trials_acc = readtable([paths.annot, filesep, paths.filestr, 'trials-accuracy.tsv'],'FileType','text','Delimiter','tab');
trials_acc.prop_correct = trials_acc.n_correct_syls ./ trials_acc.n_syllables;

op.intable = trials_acc;
op.plotvar = 'prop_correct';
op.sortvar = 'stim_group'; 
outtab = plot_windowed_timecourse(op); 

hlgd = findobj(gcf, 'Type', 'Legend');
hlgd.String = {'novel','train on', 'train off'};
title('test on','Color','k')


%% plot timecourse of duration


op.sub = 'sml003';
op.ses = 'subsyl'; 
% op.plotvar = 'speech_dur';
    op.plotvar = 'vowel_dur';
op.newfig = 0; 

    close all
    hfig = figure('Color','w')

subplot(1,2,1)
    op.run = 10;
    op.task = 'test1';
    
    paths = setpaths_dbs_learn(op);
    
    trials = readtable(paths.beh_annot_table); 
        trials.speech_dur = trials.speech_offset - trials.speech_onset; 
        trials.vowel_dur = trials.vowel_offset - trials.vowel_onset;
    trials = trials(~isnan(trials.speech_onset),:);
    trials = trials(cellfun(@isempty,trials.disfluency),:);
    
    op.intable = trials;
    op.sortvar = 'stim_group'; 

    outtab = plot_windowed_timecourse(op); 
    
    hlgd = findobj(gcf, 'Type', 'Legend');
    hlgd.String = {'novel_nat','novel_nn', 'train A off','train B on'};
    subtitle(['Subject ', op.sub, ' - test off'],'Color','k')  


subplot(1,2,2)
    op.run = 12;
    op.task = 'test2';
    
    paths = setpaths_dbs_learn(op);
    
    trials = readtable(paths.beh_annot_table); 
    trials.speech_dur = trials.speech_offset - trials.speech_onset; 
    trials.vowel_dur = trials.vowel_offset - trials.vowel_onset;
    trials = trials(~isnan(trials.speech_onset),:);
    trials = trials(cellfun(@isempty,trials.disfluency),:);
    
    op.intable = trials;
    op.sortvar = 'stim_group'; 

    outtab = plot_windowed_timecourse(op); 
    
    hlgd = findobj(gcf, 'Type', 'Legend');
    % hlgd.String = {'novel_nat','novel_nn', 'train A off','train B on'};
    subtitle(['Subject ', op.sub, ' - test on'],'Color','k')    

    