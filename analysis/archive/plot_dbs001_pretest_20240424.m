hfig = figure;

% hbar = bar([3.5, 5.5, 0] / 7);
hbar = bar([3.5, 5.5] / 7);

hax = gca;

% hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)','novel (not pretested)'}
hax.XTickLabels = {'test_a (dbs-on)','test_b (dbs-off)'}


% ylabel('n syllables correct during pretest')
ylabel('prop syllables correct during pretest')