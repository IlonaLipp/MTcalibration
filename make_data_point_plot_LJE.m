function make_data_point_plot_LJE(confound, mtmaps, fT, fa_variation_pulse_strengths_in_rad, nominal_mt_angle_in_rad, caudate_roi_coords, splenium_roi_coords, figname, textfileprefix)

    roinames={'splenium (high myelin)','caudate (low myelin)'};
    shortroinames={'splenium','caudate'};
    ROIs{1}=splenium_roi_coords;
    ROIs{2}=caudate_roi_coords;

    arange=[min(fT(fT>0)),max(fT(fT>0))]*rad2deg(nominal_mt_angle_in_rad);
    ex_angles=linspace(0,deg2rad(1000));
    ylimits=[0,10];
    
    figure('Units','inches','Position',[4,5,8.5,2.3])
    colors=get(gca,'colororder');
    
    tiledlayout(1,2,'TileSpacing','compact','Padding','tight')
    
    nexttile
    plotschematic
    plot(xlim,[0,0],'k') % add black line for x axis as rectangles overlap it...
    title('A')
    
    nexttile
    hold on
    % show range of flip angles in sample
    rectangle('Position',[arange(1),ylimits(1),diff(arange),diff(ylimits)],'FaceColor',[0.9*[1 1 1],0.5],'LineStyle','none')
    
    for ROIidx=1:length(ROIs)
        ROI=ROIs{ROIidx};
        aobs=rad2deg(fa_variation_pulse_strengths_in_rad)*fT(ROI(1),ROI(2),ROI(3));
        %MTobs=squeeze(fT(ROI(1),ROI(2),ROI(3)).^2.*mtmaps(ROI(1),ROI(2),ROI(3),:));
        MTobs=squeeze(mtmaps(ROI(1),ROI(2),ROI(3),:));
        plot(aobs,MTobs,'x','Color',colors(3+ROIidx,:));
    end
    plot(rad2deg(nominal_mt_angle_in_rad)*[1,1],ylimits,'k--')
    xlabel('\beta_{loc} / degrees')
    ylabel('MTsat / p.u.')
    ylim(ylimits)
    xlim([0,1000])
    plot(xlim,[0,0],'k') % add black line for x axis as rectangles overlap it...
    legend(roinames,'Location','northwest')
    title('B')
    
    exportgraphics(gcf,figname,'Resolution',300)

    for r = 1:length(ROIs)

        coords_to_use = ROIs{r};
        mtsatloc = squeeze(mtmaps(coords_to_use(1), coords_to_use(2), coords_to_use(3), :));
        local_ft_brain = fT(coords_to_use(1), coords_to_use(2), coords_to_use(3));
        
        writemat = [confound', mtsatloc, fa_variation_pulse_strengths_in_rad, repmat(local_ft_brain, length(fa_variation_pulse_strengths_in_rad), 1)];
        writemat_header = {'sequencenumber' 'mtapp' 'fa' 'ft'};

        %%% write out mtsatloc
        textfilename = [textfileprefix,'_',shortroinames{r},'.csv'];
        delete(textfilename);
        fid = fopen(textfilename, 'w');
        commaheader = [writemat_header; repmat({', '}, 1, numel(writemat_header))];
        commaheader = commaheader(:)';
        fprintf(fid,'%s\n',cell2mat(commaheader));
        fclose(fid)
        dlmwrite(textfilename, writemat, '-append');
    end

end


