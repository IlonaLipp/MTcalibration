%% group analysis
%%% makes figures and tables for MTsat calibration analysis
close all

brain_idx = find(sum(all.MTsat_reference) ~= 0);
colors_to_use = parula(length(brain_names));

%% histograms
for_histogram = {'n','R2','A','C','MTsat_reference'}; 
for hist_to_make = 1:length(for_histogram)
   f = figure(hist_to_make);
    set(f,'Position',[100 100 1000 300]);
       varname = for_histogram{hist_to_make}
       clear lower_lim upper_lim fi xi
       for brain = 1:length(brain_names) %brain_idx
           try
               values = all.(varname)(:,brain);
               idx = find(all.voxelmask(:,brain) ~= 0); %%%
               values_to_use = values(idx);
               %lower_lim(brain) = prctile(values_to_use,2); %%% limit range of x-axis, otherwise difficult to compare distributions
               %upper_lim(brain) = prctile(values_to_use,98);
               lower_lim(brain) = 0; %%% limit range of x-axis, otherwise difficult to compare distributions
               upper_lim(brain) = max(values_to_use);
                mask = find(values_to_use >= lower_lim(brain) & values_to_use <= upper_lim(brain));
                test = values_to_use(mask);
                [fi, xi] = ksdensity(test); %%% make  Gaussians, so values can occur outside the range
                plot(xi,fi,'color',colors_to_use(brain,:),'LineWidth',2);
                hold on
                
                set(gca,'FontSize',15)
           end
       end
        xlim([min(lower_lim), max(upper_lim)]);
        if strcmp(varname,'C')
          aah = xline(1.2);
          aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
          xticks([0:0.2:2])
        end
        legend(brain_names(brain_idx),'interpreter','None','Location','NorthWest');
        xlabel(varname,'interpreter','None');
        ylabel('ProbabilityDensity')
        saveas(f,[outpath,'/Figure_',variation,'_',varname,'_histo.png'])
end

%% table
curr_c = 'C';
display('value table')
table_filename = [outpath,'/Table_values_',curr_c,'_',variation,'.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
%% get values for each brain
clear med mean_R2 std_R2
for brain = brain_idx
    values = all.(curr_c)(:,brain);
    %idx = find(all.ft_reference(:,brain) > 0 ); %%% define brain mask
    idx = find(all.voxelmask(:,brain) > 0 ); %%% define brain mask
    idx = find(all.voxelmask(:,brain) > 0 & all.C(:,brain) > 0 & all.C(:,brain) < 1.4 ); %%% constrain statistics to realistic values!
    %med(brain) = median(values(idx)); %%% median C across brain
    med(brain) = mean(values(idx));
    var(brain) = std(values(idx));
    var_perc(brain) = 100 * var(brain) / med(brain);
    pc = partialcorr([all.(curr_c)(idx,brain),all.ft_reference(idx,brain),all.MTsat_reference_const_corr_indopt(idx,brain)],'Type','Spearman','rows','complete');
    pc_CvsB1 = pc(1,2);
    pc_CvsMT = pc(1,3);
    %%% fit median
    %mean_R2(brain) = median(all.R2(idx,brain));
    %std_R2(brain) = iqr(all.R2(idx,brain));
    mean_R2(brain) = mean(all.R2(idx,brain));
    std_R2(brain) = std(all.R2(idx,brain)); 
    within_vox_mean(brain) = mean(all.C_SE(idx,brain));
    within_vox_std(brain) = std(all.C_SE(idx,brain));
    within_vox_mean_perc(brain) = 100 * nanmean(all.C_SE(idx,brain)./ all.C(idx,brain));
    within_vox_std_perc(brain) = 100 * nanstd(all.C_SE(idx,brain) ./ all.C(idx,brain));
    fprintf(fid,['%s & ',repmat('%.3f%s%.3f & ',1,5),'%.3f & %.3f','%s'],brain_names{brain},med(brain),'$\pm$',var(brain),100,'$\pm$',var_perc(brain),within_vox_mean(brain),'$\pm$',within_vox_std(brain), within_vox_mean_perc(brain),'$\pm$', within_vox_std_perc(brain),mean_R2(brain),'$\pm$',std_R2(brain),pc_CvsB1,pc_CvsMT,'\\');
    fprintf(fid,'\n');
end
table_filename = [outpath,'/Table_group_median_iqr_',curr_c,'_',variation,'.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
fprintf(fid,['%.3f %.3f'], median(med(med~=0)), iqr(med(med~=0)));
fprintf(fid,'\n');
fclose(fid)
table_filename = [outpath,'/Table_group_mean_std_',curr_c,'_',variation,'.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
fprintf(fid,['%.3f %.3f'], mean(med(med~=0)), std(med(med~=0)));
fprintf(fid,'\n');
fclose(fid)


%% table for comparing corrected maps
display('correction approach comparison table')
correction_approaches = {'uncorrected', 'voxel', 'individual', 'global'};
correction_name_in_structure = {'MTsat_reference','MTsat_reference_corr','MTsat_reference_const_corr_indopt','MTsat_reference_const_corr_groupopt'};
table_filename = [outpath,'/Table_correction_comparison_',curr_c,'_',variation,'.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
fprintf(fid,['%s'],'brain & uncorr & voxel & ind & group & uncorr vs voxel & uncorr vs ind & voxel vs ind & ind vs group'); 
fprintf(fid,'\n');
for brain = brain_idx
    %figure()
    idx = find(all.ft_reference(:,brain) > 0 ); %%% define brain mask
    fprintf(fid,['%s & '],brain_names{brain});
    %% calculate partial spearman correlation between mt map and b1, regressing out r2s
    for c = 1:length(correction_approaches)
        correction_approach = correction_approaches{c};
        curr_corrected_map_vec = all.(correction_name_in_structure{c})(:,brain); %%% get current map
        b1partcorr.(correction_approach)(brain) = partialcorr(curr_corrected_map_vec(idx),all.ft_reference(idx, brain), all.R2s(idx, brain),'type','Spearman');
        fprintf(fid,['%.3f & '], b1partcorr.(correction_approach)(brain));
    end
    %% calculate difference between the different maps
    combinations = {[1 2], [1 3], [2 3], [3 4]};
    for c2 = 1:length(combinations)
       vec1 = all.(correction_name_in_structure{combinations{c2}(1)})(idx,brain); %%% get first map
       vec2 = all.(correction_name_in_structure{combinations{c2}(2)})(idx,brain); %%% get second map
       percent_change = 100 * abs((vec1 - vec2) ./ vec1);
       comparison_name{c2} = [correction_approaches{combinations{c2}(1)},'_vs_',correction_approaches{combinations{c2}(2)}];
       absolute_diff_between_maps.(comparison_name{c2})(brain) = nanmedian(percent_change); %nanmean(percent_change); 
       absolute_diff_between_maps_iqr.(comparison_name{c2})(brain) = iqr(percent_change); %nanstd(percent_change); 
       absolute_diff_between_maps_max.(comparison_name{c2})(brain) = max(percent_change); 
       fprintf(fid,['%.2f%s%.2f%s%.2f%s & '], absolute_diff_between_maps.(comparison_name{c2})(brain), '$\pm$', absolute_diff_between_maps_iqr.(comparison_name{c2})(brain), ' \newline (max = ', round(absolute_diff_between_maps_max.(comparison_name{c2})(brain),2),')');
    end
    fprintf(fid,'\n');
end
fclose(fid)


%% correlation with MT and b1        
corrvars = {'C', 'ft_reference', 'MTsat_reference_const_corr_indopt'}
for var1 = 1:length(corrvars)
    var1name = corrvars{var1};
    lower_lim1 = min(prctile(all.(var1name)(:,brain),5));
    upper_lim1 = min(prctile(all.(var1name)(:,brain),95));
    for var2 = 1:length(corrvars)
        f5 = figure()
        set(f5,'position',[1 1 1500 1000]); %, 'units', 'normalized')
        count = 1;
        for brain = brain_idx
            subplot(2,3,count); count = count + 1; 
            var2name = corrvars{var2};
            lower_lim2 = min(prctile(all.(var2name)(:,brain),2));
            upper_lim2 = min(prctile(all.(var2name)(:,brain),98));
            idx = find(all.(var1name)(:,brain) > lower_lim1 & all.(var1name)(:,brain) < upper_lim1 & all.(var2name)(:,brain) > lower_lim2 & all.(var2name)(:,brain) < upper_lim2);
            binscatter(all.(var1name)(idx, brain),all.(var2name)(idx,brain))
            xlabel(var1name)
            ylabel(var2name)
        end
        set(f5,'position',[1 1 1000 800]); %, 'units', 'normalized')
        saveas(f5,[outpath,'/Figure_',variation,'_',var1name,'_v2_',var2name,'_scatters.png']);
    end
end


%% temperature table
display('temperature')
table_filename = [outpath,'/Table_temperature_',variation,'.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
%% get values for each person
for brain = brain_idx
   fprintf(fid,['%s & ',repmat('%.0f & ',1,2),' %s'],brain_names{brain},all.min_temp(brain),all.max_temp(brain),'\\');
   fprintf(fid,'\n');
end
fclose(fid)