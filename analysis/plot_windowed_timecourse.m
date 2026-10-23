%%%%% sort table by specific variable (e.g. stim or condition)...
% ..... then plot timecourse of another variable (e.g. accuracy)  for each value of the sorting variable

function [outtab] = plot_windowed_timecourse(op)

field_default('op','intable',table); 
field_default('op','sortvar','name');
field_default('op','plotvar','propacc'); 
field_default('op','windowsize',1); 
field_default('op','endpoint_method','discard'); 
field_default('op','newfig',1); 
field_default('op','line_style','-'); 
field_default('op','line_width',2); 
field_default('op','marker_style','none'); 
field_default('op','marker_size',6); 
field_default('op','ylimits',[]); 
field_default('op','xlabel','trial'); 
field_default('op','jitter_prop_of_yrange',0.005); % jitter traces by this proportion of y axis range; set to zero to not jitter
field_default('op','legend_location', 'southeast'); % [left bottom width height]

unq_sortvals = unique(op.intable{:,op.sortvar},'sorted');
n_sortvals = length(unq_sortvals); 
celcol = cell(n_sortvals,1); 
outtab = table(unq_sortvals, celcol, celcol, celcol, 'RowNames', unq_sortvals, 'VariableNames', {op.sortvar, 'tab', op.plotvar, [op.plotvar, '_win']});

all_plotvals = op.intable{:,op.plotvar}; 
plotvalrange = range(all_plotvals); 

if op.newfig
    hfig = figure('Color','w'); 
end

task_order = {'trainA','trainB','novel1','novel2','novel1_nat','novel2_nat','novel1_nn','novel2_nn'}; 

%%
% 1. Get the current row names of your table
current_names = outtab.Properties.RowNames;

% 2. Find row names that are in task_order, preserving task_order's alignment
% 'stable' ensures the order of task_order is strictly maintained
ordered_matches = intersect(task_order, current_names, 'stable');

% 3. Find any remaining row names that weren't included in task_order
% 'stable' preserves their original relative order from the table
remaining_names = setdiff(current_names, ordered_matches, 'stable');

% 4. Combine them into the final desired row name order
final_order = [ordered_matches; remaining_names];

% 5. Index back into the table to reorder the rows
outtab = outtab(final_order, :);



%%
for isortval = 1:n_sortvals
    this_sortval = unq_sortvals{isortval};
    matchrows = string(op.intable{:,op.sortvar})==this_sortval; 
    outtab.tab{isortval} = op.intable(matchrows,:); % full table for this sortval
    outtab{isortval,op.plotvar} = {outtab.tab{isortval}{:,op.plotvar}}; % array of only the plotvar for this sortvar
    outtab{isortval, [op.plotvar, '_win']} = {movmean(outtab{isortval,op.plotvar}{1}, op.windowsize, 'Endpoints', op.endpoint_method)}; 
    xvals = op.windowsize : length(outtab{isortval,op.plotvar}{1}); % label x values by the latest value within the window

    if isempty(op.ylimits)
        yax_buffer_prop = 0.1; 
        op.ylimits = [min(all_plotvals)-yax_buffer_prop*plotvalrange, max(all_plotvals)+yax_buffer_prop*plotvalrange];
    end
    jittermax = op.jitter_prop_of_yrange * diff(op.ylimits);
    yvals = outtab{isortval, [op.plotvar, '_win']}{1}; 
        yvals = yvals + [rand(size(yvals))-0.5]*jittermax; % add jitter to avoid overlapping traces
    
    hplot = plot(xvals,yvals,'LineWidth',op.line_width, 'LineStyle',op.line_style, 'Marker',op.marker_style, 'MarkerSize',op.marker_size);
    ylim(op.ylimits)
    hold on

    if any(strcmp(this_sortval,{'trainA','trainB'}))
       cfg = [];
        cfg.sub = op.sub;
        cfg.ses = op.ses;
        cfg.task = this_sortval;
        outtab.dbs_train_state{isortval} = get_dbs_state(cfg); 
    else
        outtab.dbs_train_state{isortval} = '';
    end
    outtab.leglabel{isortval} = [this_sortval, ' ', outtab.dbs_train_state{isortval}]; 

end

hax = gca;
hax.Color = 'w';
hax.XColor = 'k';
hax.YColor = 'k';

hylab = ylabel(op.plotvar,'Interpreter','none');
hxlab = xlabel(op.xlabel,'Interpreter','none');
hleg = legend(outtab.leglabel, 'Location',op.legend_location, 'Interpreter','none');
box off

