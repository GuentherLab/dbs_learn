%%% plot accuracy of a single pilot subject

close all 

datapath = 'C:\dbsmulti';

% sub = 'pilot_001';
% sub = 'pilot_003';
% sub = 'pilot_004';
% sub = 'pilot_005';
% sub = 'dbs001';
% sub = 'dbs002';
sub = 'dbs003';


ses = '1'; 
run = '1';

% task = 'train-a';
% task = 'train-b';
task = 'test';

% performance_var = 'correct';
performance_var = 'syl_acc_prop'; 

add_jitter = 1; % add noise to traces so that they don't perfectly overlap
    jittermax = 0.005; 
linedwidth = 2; 

plot_by_condition = 1; 
    cond_windowsize = 10; 
plot_by_stimname = 1; 
    stimname_windowsize = 4; 


set(0, 'DefaultTextInterpreter', 'none')

%% load and prepare data
trials = readtable([datapath filesep 'sub-' sub filesep...
                    'sub-' sub '_ses-' ses '_run-' run '_task-' task '_acc.csv']);
trials = trials(~(trials.unusableTrial==1 & ~isnan(trials.unusableTrial==1)),:);
trials.syl_acc_prop = trials.n_correct_syls ./ trials.n_syllables; 

trials.stim_group = cell(height(trials), 1); 

for itrial = 1:height(trials)
    switch trials.name{itrial}
        case 'DOY-TOW-VAB-FO-VIL-AN'
            trials.stim_group{itrial} = 'test_a';
        case 'ILL-BER-AM-TAY-VOI-CHYG'
            trials.stim_group{itrial} = 'test_a';
        case 'SIK-AN-THOD-HES-TO-MEE'
            trials.stim_group{itrial} = 'test_b';
        case 'GINS-TA-BUL-JELL-AN-TIF'
            trials.stim_group{itrial} = 'test_b';
        case 'KAF-TOG-ROP-INK-IS-TA'
            trials.stim_group{itrial} = 'test_novel';
        case 'NOI-CHER-TYSS-HEEL-KO-SUM'
            trials.stim_group{itrial} = 'test_novel';
    end
end

%%%%%%%%%%%%%%% correct for pretest performance
normalize_by_pretest = 0; 
if normalize_by_pretest
    if strcmp(sub,'dbs001')
        for itrial = 1:height(trials)
           switch trials.name{itrial}
            case 'DOY-TOW-VAB-FO-VIL-AN-SHA'
                trials.n_correct_syls(itrial) = trials.n_correct_syls(itrial) - 3.25;
            case 'ILL-BER-AM-TAY-VOI-CHYG-BO'
                trials.n_correct_syls(itrial) = trials.n_correct_syls(itrial) - 3.75;
            case 'SIK-AN-THOD-HES-TO-MEE-ZIL'
                trials.n_correct_syls(itrial) = trials.n_correct_syls(itrial) - 6.5;
            case 'GINS-TA-BUL-JELL-AN-TIF-CHOO'
                trials.n_correct_syls(itrial) = trials.n_correct_syls(itrial) - 4.5;
            case 'KAF-TOG-ROP-INK-IS-TA-NAY'
                trials.stim_group{itrial} = 'test_novel';
            case 'NOI-CHER-TYSS-HEEL-KO-SUM-FEE'
                trials.stim_group{itrial} = 'test_novel';
            end 
        end
    end
end





% sort by condition
unq_cons = {'test_a';'test_b';'test_novel'};
ncons = length(unq_cons);
conacc = cell(ncons,1); 
for icon = 1:ncons
    thisname = unq_cons{icon};
    temptrials = trials(strcmp(trials.stim_group,thisname),:);
    ncontrials = height(temptrials);

    xvals = cond_windowsize:ncontrials; % index of trial relative to all trials within this condition
    ntrials_to_plot = length(xvals);
    conacc{icon} = nan(ntrials_to_plot,1); 
    for itrial = 1:ntrials_to_plot
        trial_con_relative = xvals(itrial); 
        conacc{icon}(itrial) = mean(temptrials{trial_con_relative-cond_windowsize+1:trial_con_relative, performance_var});
    end
end
contab = table(unq_cons, cellfun(@mean,conacc), conacc,'VariableNames',{'con','acc_mean','trialacc'});

% sort by stim name
unq_names = unique(trials.name); 
nnames = length(unq_names);
nameacc = cell(nnames,1); 
nametab = table(unq_names,cell(nnames,1),nan(nnames,3),'VariableNames',{'name','stim_group','color'});
    clear unq_names
%%%%%%% sort by stim name
for iname = 1:nnames
    thisname = nametab.name{iname};
    nametab.stim_group{iname} = trials.stim_group{find(strcmp(trials.name,thisname),1)};
end
nametab = sortrows(nametab,'stim_group');

for iname = 1:nnames
    thisname = nametab.name{iname};

    temptrials = trials(strcmp(trials.name,thisname),:);
    nnametrials = height(temptrials);

    xvals = stimname_windowsize:nnametrials; % index of trial relative to all trials of this name
    ntrials_to_plot = length(xvals);
    nameacc{iname} = nan(ntrials_to_plot,1); 
    for itrial = 1:ntrials_to_plot
        trial_con_relative = xvals(itrial); 
        nameacc{iname}(itrial) = mean(temptrials{(trial_con_relative-stimname_windowsize+1:trial_con_relative), performance_var});
    end
end
nametab.acc_mean = cellfun(@mean,nameacc);
nametab.trialacc = nameacc; 

switch performance_var
    case 'correct'
        ylab = 'Accuracy (full stimulus)';
    case 'syl_acc_prop'
        ylab = 'Proportion of accurate syllables';
end



% get accuracy from all trials, not windowed data
learned_accuracy =  mean(trials{contains(trials.stim_group,{'test_a','test_b'}),performance_var})
novel_accuracy = mean(trials{contains(trials.stim_group,{'test_novel'}),performance_var})

%% plot timecourses of each condition


if plot_by_condition 
    hfig = figure;
    hold on
    for icon = 1:ncons
        xvals = [1:length(contab.trialacc{icon})]' + cond_windowsize - 1; % index of trial relative to all trials within this condition
        yvals = contab.trialacc{icon};
        if add_jitter
            yvals = yvals + [rand(size(yvals))-0.5]*jittermax;
        end
        hplot = plot(xvals,yvals,'LineWidth',linedwidth);
    end
    
    ileg = legend(unq_cons,'Interpreter','none');
    hylab = ylabel (ylab);
    hxlab = xlabel('trial (within-condition)');
    title([sub '...phase = ' task '..... trial-window size = ' num2str(cond_windowsize)])
end



%% sort trials by stim name and plot timecourses
if plot_by_stimname
switch task
    case 'train-a'
        nametab.linecolor = [  0.6 0.6 1;...       % test a
                                0.1 0.1 0.7];   % test a
    case 'train-b'
        nametab.linecolor = [  0.6 1 0;...       % test b
                                0 0.7 0.4];   % test b
    case 'test'
        nametab.linecolor = [  0.6 0.6 1;...       % test a
                                0.1 0.1 0.7;...   % test a
                                0.6 1 0;...       % test b
                                0 0.7 0.4;...   % test b
                                1 0 0;...       % test novel
                                1 0.1 0.6];      % test novel
        nametab.linestyle = {   '-';...     % test a
                        '-';...     % test a
                        '-';...     % test b
                        '-';...     % test b
                        '-';...     % test novel
                        '-'};       % test novel
end


% make plot
hfig = figure;
hold on
for iname = 1:nnames
    xvals = [1:length(nametab.trialacc{iname})]' + stimname_windowsize - 1; % index of trial relative to all trials of this name
    yvals = nametab.trialacc{iname};
    if add_jitter
        yvals = yvals + [rand(size(yvals))-0.5]*jittermax;
    end
    hplot(iname) = plot(xvals,yvals,'Color',nametab.linecolor(iname,1:3), 'LineWidth',linedwidth);
    % hplot(iname).LineStyle = nametab.linestyle{iname};
end

leglabels = strcat(nametab.stim_group,repmat({'...'},nnames,1),nametab.name); 
ileg = legend(leglabels,'Interpreter','none');
hylab = ylabel (ylab);
hxlab = xlabel('trial (within-stim-name)');
title([sub '...phase = ' task '..... trial-window size = ' num2str(stimname_windowsize)])

end