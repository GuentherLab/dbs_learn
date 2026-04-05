%% timecourse plot formatting function
% called by pd_r01_figures script


function h_tfill = format_timecourse(cfg)
    hax = gca;
    ylimits = hax.YLim; 

% add shading for specified time periods (e.g. stim and speech)
        hold on

    n_fillblocks = size(cfg.tfill_xvals,1);
    for iblock = 1:n_fillblocks
        xpair = cfg.tfill_xvals(iblock,:);
        h_tfill(iblock) = fill([xpair(1),xpair(2),xpair(2),xpair(1)],[ylimits(1),ylimits(1),ylimits(2),ylimits(2)],  cfg.tfill_color,...
            'FaceAlpha', cfg.tfill_alpha, 'HandleVisibility','off'); % standard error
        h_tfill(iblock).LineStyle = cfg.tfill_LineStyle; % border
            h_tfill(iblock).EdgeColor = cfg.tfill_edge_color; 

    end

    box off

end