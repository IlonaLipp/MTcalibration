clear all
close all

addpath(genpath('~/Documents/scripts/postmortembrain-mpm/'))
outpath = '/data/pt_02101/results/mt_calibration_7t';
mkdir('/data/pt_02101/results/mt_calibration_7t/hirescors/')

brain_id_slash_scan_ids = {'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB'...
    '018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB'...
    '025_C_W_RAVEL_TAI_S/mr/200923_Terra_7T_32Ch_WB'...
    '026_C_W_OSCAR_TAI_S/mr/200921_Terra_7T_32Ch_WB'...
    '032_C_C_SAMSON_SL_TAC/mr/210206_Terra_7T_32Ch_WB'};
brain_names = {'brain 1','brain 2','brain 3', 'brain 4', 'brain 5'}; %, 'brain 6', 'brain 7'};

%%% load global (median across individual)
global_C = importdata('/data/pt_02101/results/mt_calibration_7t/Table_group_median_iqr_C_MTvariation.txt');
global_C = 1.2

for scan_to_analyse = 1:length(brain_id_slash_scan_ids)
    
    brain_id_slash_scan_id = brain_id_slash_scan_ids{scan_to_analyse};
    brain_id_underscore_scan_id = strrep(brain_id_slash_scan_id,'/','_');
    processed_data_folder = ['/data/pt_02101/preprocessed/',brain_id_slash_scan_id,'_V2'];
    hires_mt_map = [processed_data_folder,'/MPM_calc/MTsat_real_TE0_WLS_rescaled_to_FA700_processed_s_rbc_u_0p3_run01_brain_masked.nii'];
    hires_r2s_map = [processed_data_folder,'/MPM_calc/R2s_WLS_processed_s_rbc_u_0p3_run01_brain_masked.nii'];
    hires_r1_map = [processed_data_folder,'/MPM_calc/R1_TE0_WLS_processed_s_rbc_u_0p3_run01_brain_masked.nii'];
    b1_map = [processed_data_folder,'/B1mapping/B1map_AL_processed_0p3_mean.nii'];
    
    if exist(hires_mt_map) == 2 
        
        %%% read in stuff
        mtvol = spm_vol(hires_mt_map);
        mtuncorr = spm_read_vols(spm_vol(hires_mt_map));
        mtuncorrvec = mtuncorr(:);
        
        b1 = spm_read_vols(spm_vol(b1_map));
        b1vec = b1(:);
        mask = find(b1vec > 0);
        
        r2s = spm_read_vols(spm_vol(hires_r2s_map));
        r2svec = r2s(:);
        medrel.r2smedian(scan_to_analyse) = median(r2svec(mask));
        medrel.r2siqr(scan_to_analyse) = iqr(r2svec(mask));
        
        r1 = spm_read_vols(spm_vol(hires_r1_map));
        r1vec = r1(:);
        medrel.r1median(scan_to_analyse) = median(r1vec(mask));
        medrel.r1iqr(scan_to_analyse) = iqr(r1vec(mask));
        
        %%% load individual constants
        clear ind_global_C
        
        load(['/data/pt_02101/results/mt_calibration_7t/',brain_id_underscore_scan_id,'_MTvariation_individual_calibration_constants.mat'])

        correction_approaches = {'individual', 'global', 'uncorrected'};
        
        corr_map.uncorrected = mtuncorr;
        corr_map.individual = b1_correct_mt_map_linear(mtuncorr, b1, round(ind_global_C,2));
        corr_map.global = b1_correct_mt_map_linear(mtuncorr, b1, round(global_C,2));

        %% make histogram
        f = figure();
        cs_to_use = [2 3];
        colors_to_use = {[115 20 122]./255,[79 177 84]./255,[39 34 138]./255}
        for cc = 1:length(cs_to_use)
            c = cs_to_use(cc);
           varname = correction_approaches{c};
           values = corr_map.(varname)(:);
           values_to_use = values(values ~= 0);
           lower_lim = min(prctile(values_to_use,2));
           upper_lim = min(prctile(values_to_use,98));
            mask = find(values > lower_lim & values < upper_lim & values ~= 0);
            [fi, xi] = ksdensity(values(mask));
            plot(xi,fi,'color',colors_to_use{cc},'LineWidth',3);
            hold on
            xlim([lower_lim, upper_lim]);
            xlabel('MTsat');
            ylabel('frequency');
            set(gca,'FontSize',15)
        end
        legend(correction_approaches(cs_to_use),'interpreter','None','Location','NorthWestOutside');
        title(brain_names{scan_to_analyse},'interpreter','None');
        saveas(f,[outpath,'/Figure_MTvar_',brain_id_underscore_scan_id,'_hires_correction_histo.png'])
        
        %%
        for c = 1:length(correction_approaches)
            correction_approach = correction_approaches{c};
            curr_corrected_map = corr_map.(correction_approach);
            curr_corrected_map_vec = curr_corrected_map(:);
            %%% save map
            mtvol.fname =  ['/data/pt_02101/results/mt_calibration_7t/hirescors/', brain_id_underscore_scan_id,'_hires_corr_', correction_approach,'.nii'];
            spm_write_vol(mtvol, curr_corrected_map);
            %%% analyse
            curr_corrected_map_vec = curr_corrected_map(:);
            %%% what is the difference to uncorrected? 
            idx = b1vec > 0; %mtuncorr(:) > 0;
            percent_change = 100 * abs((mtuncorr(idx) - curr_corrected_map_vec(idx)) ./ mtuncorr(idx) );
            absolute_diff_between_maps.(correction_approach)(scan_to_analyse) = nanmean(percent_change);
            absolute_diff_between_maps_std.(correction_approach)(scan_to_analyse) = nanstd(percent_change);
            %%% correlation with b1
            rho = corr(curr_corrected_map_vec(idx),b1vec(idx),'Type','Spearman');
            b1corr.(correction_approach)(scan_to_analyse) = rho; % r(1,2);
            b1partcorr.(correction_approach)(scan_to_analyse) = partialcorr(curr_corrected_map_vec(idx),b1vec(idx),r2svec(idx),'type','Spearman');
 
            %%% correlation with r2s
           % [r p] = corrcoef(curr_corrected_map_vec(idx),r2svec(idx));
            rho = corr(curr_corrected_map_vec(idx),r2svec(idx),'Type','Spearman');
            r2scorr.(correction_approach)(scan_to_analyse) = rho;% r(1,2);
        end
          %%% diff between individual and global
          indvec = corr_map.individual(:); 
          globvec = corr_map.global(:);
          percent_change = 100 * abs((indvec(idx) - globvec(idx)) ./ indvec(idx) );
          absolute_diff_between_maps.individual_vs_global(scan_to_analyse) = nanmean(percent_change); % nanmean(percent_change); 
          absolute_diff_between_maps_std.individual_vs_global(scan_to_analyse) = iqr(percent_change); %nanstd(percent_change); 
    end
end

%%% write out table
table_filename = [outpath,'/Table_hires_correction_evaluation_MTvar.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
%% get values for each person
for brain = 1:length(brain_id_slash_scan_ids)
   fprintf(fid,['%s & ',repmat('%.3f & ',1,3),repmat('%.2f%s%.2f & ',1,3), ' %s'],brain_names{brain},round(b1partcorr.uncorrected(brain),3),round(b1partcorr.individual(brain),3),round(b1partcorr.global(brain),3),round(absolute_diff_between_maps.individual(brain),2), '$\pm$', round(absolute_diff_between_maps_std.individual(brain),2), round(absolute_diff_between_maps.global(brain),2),'$\pm$', round(absolute_diff_between_maps_std.global(brain),2), round(absolute_diff_between_maps.individual_vs_global(brain),2),'$\pm$',round(absolute_diff_between_maps_std.individual_vs_global(brain),2),'\\');
   fprintf(fid,'\n');
end
fclose(fid)

%%% median table
table_filename = [outpath,'/Table_median_iqr_R2s_and_R1.txt'];
delete(table_filename)
fid = fopen(table_filename,'w');
%% get values for each person
for brain = 1:length(brain_id_slash_scan_ids)
   fprintf(fid,['%s & ',repmat('%.2f %s %.2f & ',1,2), ' %s'],brain_names{brain},round(medrel.r1median(brain),3),'$\pm$',round(medrel.r1iqr(brain),3),round(medrel.r2smedian(brain),3),'$\pm$',round(medrel.r2siqr(brain),3),'\\');
   fprintf(fid,'\n');
end
fclose(fid)

