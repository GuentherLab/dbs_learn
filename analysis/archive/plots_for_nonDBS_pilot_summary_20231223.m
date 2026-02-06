%%%% create plots for Frank Guenther to send to Mark Richardson and Todd Herrington

clear
close all

n_early_trials_cond = 6; 
n_early_trials_stimname = 3; 

recency_plot_individual_subs = 0;
plot_allsubs_timecourse_condition = 0; 
plot_allsubs_timecourse_stimname = 0; 
if ~exist('subtab','var')
    compile_dbsmulti_subs_20231224();
end



plotops.ylab = ylab; 

%% plot recency effects across all subs
% close all
hfig = figure('Color','w'); box off
    hfig.Position = [680   443   288   435]; 

plotdat = [subtab.acc_trained_1, subtab.acc_trained_2];
[herb, hax, hanot] = create_2datapoint_plot(plotdat, plotops);
title({'Accuracy of 1st vs. 2nd trained stim sets', 'during Testing phase'})
hax.XTickLabel = {'1st trained set','2nd traineed set'};
hax.XTickLabelRotation = 30;
    hanot.Position(1:2) = [0.23 0.23];
ylim([0.9 1]); 
herb.LineStyle = 'none';

[~, p_recency] = ttest(subtab.acc_trained_1, subtab.acc_trained_2); 

%% plot test-set accuracy in trained set A vs. trained set B across all subs

% close all
hfig = figure('Color','w'); box off
    hfig.Position = [680   443   288   435]; 


plotdat = [mean(allsubs_cond.trialacc{string(allsubs_cond.label)=="test_a"},2),...
           mean(allsubs_cond.trialacc{string(allsubs_cond.label)=="test_b"},2)];
[herb, hax, hanot] = create_2datapoint_plot(plotdat, plotops);
title({'Accuracy of trained stim sets A vs. B', 'during Testing phase'})
hax.XTickLabel = {'Trained set A','Traineed set B'};
hax.XTickLabelRotation = 30;
hanot.Position(1:2) = [0.5 0.13];
ylim([0.9 1]);
herb.LineStyle = 'none';




%% plot within-set training effect
% accuracy of training sets in the first X trials of their training phase compared to their first X trials in testing phase

% close all
hfig = figure('Color','w'); box off
    hfig.Position = [680   443   288   435]; 

early_acc_trainphase = [subtab.acc_trainphase_early_a, subtab.acc_trainphase_early_b];
early_acc_testphase = horzcat(allsubs_name.early_acc_submean{string(allsubs_name.stim_group)=="test_a" | string(allsubs_name.stim_group)=="test_b"});
plotdat = [mean(early_acc_trainphase,2), mean(early_acc_testphase,2)];
[herb, hax, hanot] = create_2datapoint_plot(plotdat, plotops);
title({'Accuracy at start of Training phase vs.', 'start of Testing phase'})
hax.XTickLabel = {'Training','Testing'};
hanot.Position(1:2) = [0.6 0.1];
ylim([0.8 1]);
herb.LineStyle = ':';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting function
function [herb, hax, hanot] = create_2datapoint_plot(plotdat, plotops)
    % hbar = bar(mean([plotdat])); hold on
nsubs = size(plotdat,1);
herb = errorbar(mean(plotdat), std(plotdat) ./ nsubs, 'LineStyle','none', 'Marker','o', 'Color',[0 0 0]);
    herb.LineWidth = 2;
    herb.MarkerSize = 10;
    herb.MarkerFaceColor = [0.4 0.4 0.4];
xlim([0.2 2.8])
hax = gca;
    hax.XTick = [1 2];
    hax.FontWeight = 'bold';
    hax.Box = 'off';
ylabel(plotops.ylab)

[~, pval] = ttest(plotdat(:,1), plotdat(:,2));
hanot = annotation('textbox','String',['p=',num2str(pval)],'FitBoxToText','on');



end