 %%% run pd grant plotting function while varying parameters

close all

 for ii = [5 10 25 33 50]
op.trainphase_half_ntrials = ii;
pd_r01_figures_20240916;
nexttile(tilenum(htl,2,1))
title(['train ON, ntrials = ', num2str(ii)])
end