%%% import accuracy data from SEQ-Multisyllabic dataset
% this dataset (.xlsx) is expected to contain separate spreadsheet tabs for each subject
% ... this script will compile all subjects and organize accuracy by stim identity

seqm_results_filename = 'C:\Users\amsmeier\Downloads\SEQM_Results.xlsx'; 
subsheet_name_start_string = 'SEQM'; % all subject sheets must start with this string
    spreadsheet_range_to_import = "B1:T31"; % accuracy data must be contained within exactly this range across all subjects
    fmri_run_var_names = {'Run1','Run2','Run3','Run4','Run5','Run1_1','Run2_1','Run3_1','Run4_1','Run5_1'}; 
    subject_vars_that_must_not_be_nans = {'Pre','Post','Pre_1','Post_1'}; 
stimlist_sheet_name = 'ABC'; % stimuli divided into groups A, B, C should be contained in this sheet
    stimA_group_range = 'E2:E31';
    stimB_group_range = 'B2:B31';
    stimC_group_range = 'N2:N31';
output_savefile = [fileparts(seqm_results_filename) filesep 'seqm_stim_accuracy.xlsx']; 
plot_results = 1; 
    thresh_count_to_plot = 4; % do not plot stim which have fewer trials than this value

%% concatenate data across subjects
sheettabs = sheetnames(seqm_results_filename);
sub_tab_names = sort(sheettabs(contains(sheettabs, subsheet_name_start_string)));
nsubs = size(sub_tab_names,1);

rawtable = table; 
for isub = 1:nsubs
    thissheetname = sub_tab_names{isub};
    temptab = readtable(seqm_results_filename,'Sheet',thissheetname, Range=spreadsheet_range_to_import);
    nstimrows = height(temptab);
    temptab.sub = repmat(thissheetname,nstimrows,1);
    temptab = movevars(temptab,'sub','Before',1);

    % check whether there are any empty values in pretest/posttest scorings; if so, this subj has not been fully scored for pre/post and should be skipped
    prepost_var_types = cellfun(@ (x) class(temptab.(x)), subject_vars_that_must_not_be_nans, 'UniformOutput', false); 
    if all(strcmp('double',prepost_var_types) ) &&...
            ~any(isnan(temptab{:,subject_vars_that_must_not_be_nans}(:)))

        % remove subject-specific stim set names so that we can concat all subjs
        tabvars = temptab.Properties.VariableNames;
        newvars = cellfun(@(x)regexprep(x,'Set.','Set'),tabvars,'UniformOutput',false);
        temptab.Properties.VariableNames = newvars;
        
        % if this subject's fmri runs have not yet been scored, fill in with nans
        if iscell(temptab{:,fmri_run_var_names{1}})
            for ivar = 1:length(fmri_run_var_names)
                temptab.(fmri_run_var_names{ivar}) = nan(nstimrows,1);
            end
        end
    
        % add this sub's table to the overall table
        rawtable = [rawtable; temptab];

    end
end

%% sort by stim name and compute accuracies
% load stim lists
stimlist_A = readtable(seqm_results_filename,'Sheet',stimlist_sheet_name, Range=stimA_group_range, ReadVariableNames=false);
stimlist_B = readtable(seqm_results_filename,'Sheet',stimlist_sheet_name, Range=stimB_group_range, ReadVariableNames=false);
stimlist_C = readtable(seqm_results_filename,'Sheet',stimlist_sheet_name, Range=stimC_group_range, ReadVariableNames=false);
stim_grouped = table([repmat('a',height(stimlist_A),1); repmat('b',height(stimlist_B),1); repmat('c',height(stimlist_C),1)],...
    [stimlist_A{:,:}; stimlist_B{:,:}; stimlist_C{:,:}], 'VariableNames', {'stimset', 'stim'} ); 


% compute accuracies
stimnames = unique([rawtable.Set_Novel_Beh_; rawtable.Set_Trained_]);
nstim = length(stimnames); 
nancol = nan(nstim,1);
celcol = cell(nstim,1);
stimacc = table(stimnames, celcol,     nancol,    nancol,   nancol,    nancol, 'VariableNames', ...
                {'stim', 'stimset',  'pre_mean','pre_std','pre_sem','pre_count'}); 
for istim = 1:nstim
    thisstim = stimnames{istim};
    stimacc.stimset{istim} = stim_grouped.stimset(strcmp(thisstim,stim_grouped.stim)); 
    matchnovel = strcmp(thisstim,rawtable.Set_Novel_Beh_);
    matchtrained = strcmp(thisstim,rawtable.Set_Trained_);
    pre_vals = [rawtable.Pre(matchnovel); rawtable.Pre_1(matchtrained)]; % concat pretest trials from trained and novel
    stimacc.pre_mean(istim) = mean(pre_vals);
    stimacc.pre_std(istim) = std(pre_vals); 
    stimacc.pre_count(istim) = nnz(~isnan(pre_vals)); 
    stimacc.pre_sem(istim) = std(pre_vals) / sqrt(stimacc.pre_count(istim)); 
end

writetable(stimacc, output_savefile)

%% plotting
if plot_results
    stimacc_to_plot = stimacc(stimacc.pre_count >= thresh_count_to_plot, :); 
    nstim_plot = height(stimacc_to_plot); 
    hfig = figure;
    hplot = plot(stimacc_to_plot.pre_mean); 
        hplot.LineStyle = 'none';
    hax = gca;
        hax.XTick = 1:nstim_plot;
        hax.XTickLabel = stimacc_to_plot.stim;
        hax.XTickLabelRotation = 10;
        hax.XTickLabelRotation = 90;
    box off
    hold on
    hebar = errorbar(stimacc_to_plot.pre_mean,stimacc_to_plot.pre_sem);
        hebar.LineStyle = 'none';
end








