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