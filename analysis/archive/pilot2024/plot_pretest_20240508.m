hfig = figure;


hbar = bar([3.75, 4] / 6);

hax = gca;

% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)','novel (not pretested)'}
% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}
hax.XTickLabels = {'test_a (dbs-off)','test_b (dbs-on)'}



% ylabel('n syllables correct during pretest')
ylabel('prop syllables correct during pretest')

%% dbs 003

hfig = figure;


hbar = bar([mean([4.5 4.5 4.5 3]), mean([5 5.5 3 2.5])] / 6);

hax = gca;

% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)','novel (not pretested)'}
% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}
hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}



% ylabel('n syllables correct during pretest')
ylabel('prop syllables correct during pretest')

%% dbs 004

hfig = figure;


hbar = bar([mean([3 3.5 3.5 4]), mean([3.5 1.5 2.5 1.5])] / 7);

hax = gca;

% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)','novel (not pretested)'}
% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}
hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}



% ylabel('n syllables correct during pretest')
ylabel('prop syllables correct during pretest')