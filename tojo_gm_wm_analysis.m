close all
clear all
clc
%%% makes GM and WM histograms for segmented brain


addpath(genpath('~/Documents/scripts/postmortembrain-mpm/'))
axes_to_flip = [];

mt_folder = '/data/pt_02101/results/mt_calibration_7t/018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB';
fs_folder = '/data/pt_02101/preprocessed/018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB_V2/freesurfer/MTsat_0p3mm_downsampled_to_0p7mm_FSV7_hires/';

%% load B1
b1file = [fs_folder,'/../../MPMs_to_use/B1_0p3_resampled_to_0p7_run01_brain_masked_reoriented.nii']
b1mat = spm_read_vols(spm_vol(b1file));
b1_vec = b1mat(:);

%% load cortex mask
lh_cortex_file = [fs_folder,'/aseg_volume_masks/lh_cortex.nii'];
cormat = spm_read_vols(spm_vol(lh_cortex_file));
rh_cortex_file = [fs_folder,'/aseg_volume_masks/rh_cortex.nii'];
cormat = cormat + spm_read_vols(spm_vol(rh_cortex_file));
cor_vec = cormat(:);
cor_idx = find(cor_vec > 0);

%% load WM mask
lh_wm_file = [fs_folder,'/aseg_volume_masks/lh_wm.nii'];
wmmat = spm_read_vols(spm_vol(lh_wm_file));
rh_wm_file = [fs_folder,'/aseg_volume_masks/rh_wm.nii'];
wmmat = wmmat + spm_read_vols(spm_vol(rh_wm_file));
wm_vec = wmmat(:);
wm_idx = find(wm_vec > 0);

%% load C 
highres_file = [mt_folder,'/../hirescors/018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_hires_corr_global.nii'];
cfile = [mt_folder,'/C_MTvariation_variation_018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_2p1_run01_brain_masked.nii'];
%%% to get it into the same space as masks, theoretically need to be
%%% upsampled to 0p3, reoriented the same way as mpms, and then downsample
%%% again
cfile_resampled = [mt_folder,'/C_MTvariation_variation_018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_2p1_run01_brain_masked_resampled.nii'];
lowres_mt = '/data/pt_02101/results/mt_calibration_7t/018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB/MTsat_reference_const_corr_groupopt_MTvariation_variation_018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_2p1_run01_brain_masked.nii';
system(['flirt -in ', lowres_mt, ' -ref ', highres_file, ' -dof 6 -omat ', mt_folder, '/reg2p1to0p3'])
%%% apply now to C
system(['flirt -in ', cfile, ' -ref ', highres_file, ' -out ', cfile_resampled, ' -interp nearestneighbour -applyxfm -init ', mt_folder, '/reg2p1to0p3'])
%%% reorient
cfile_reoriented = [mt_folder,'/C_MTvariation_variation_018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_2p1_run01_brain_masked_reoriented.nii'];
reorientation_vector = [2 3 1];
correct_PMB_orientation(cfile_resampled,cfile_reoriented,reorientation_vector,axes_to_flip);
%%% resample again
cfile_reoriented_to0p7 = [mt_folder,'/C_MTvariation_variation_018_C_C_TOJO_mr_200922_Terra_7T_32Ch_WB_2p1_run01_brain_masked_reoriented_resampled0p7.nii'];
system(['flirt -in ', cfile_reoriented, ' -ref ', lh_cortex_file, ' -out ', cfile_reoriented_to0p7, ' -interp nearestneighbour -nosearch -applyisoxfm 0.7']);
%%% finally load
cmat = spm_read_vols(spm_vol(cfile_reoriented_to0p7));
c_vec = cmat(:);

%% analyse
boxplot(make_matrix_for_boxplot({c_vec(cor_idx), c_vec(wm_idx)}));
ylim([0.8 1.6])
%%% get rid of outliers in c_vec
lb = prctile(c_vec(c_vec>0),1);
ub = prctile(c_vec(c_vec>0),99);
idx = find(c_vec > lb & c_vec < ub);

%% figures
f = figure()
 set(f,'Position',[100 100 1000 500]);
subplot(1,3,2)
    hold on
    
    v1 = c_vec(intersect(cor_idx,idx));
    [fi, xi] = ksdensity(v1); %%% make  Gaussians, so values can occur outside the range
    plot(xi,fi,'color','r','LineWidth',2);
    aah = xline(median(v1),'color','r');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    v2 = c_vec(intersect(wm_idx, idx));
    [fi, xi] = ksdensity(v2);
    plot(xi,fi,'color','b','LineWidth',2);
    aah = xline(median(v2));
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %legend({'cortex','wm'})
    legend({['cortex: ', num2str(round(median(v1),2)), '\pm', num2str(round(iqr(v1),2))],['WM: ', num2str(round(median(v2),2)), '\pm', num2str(round(iqr(v2),2))]}, 'Location', 'NorthOutside')

    xlim([1.1 1.3])
    xlabel(['C'])
    ylabel('Probability density')

%figure()
subplot(1,3,1)
    hold on
    
    v1 = b1_vec(intersect(cor_idx,idx));
    [fi, xi] = ksdensity(v1); %%% make  Gaussians, so values can occur outside the range
    plot(xi,fi,'color','r','LineWidth',2);
    aah = xline(median(v1),'color','r');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    v2 = b1_vec(intersect(wm_idx,idx));
    [fi, xi] = ksdensity(v2);
    plot(xi,fi,'color','b','LineWidth',2);
    aah = xline(median(v2),'color','b');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    legend({['cortex: ', num2str(round(median(v1))), '\pm', num2str(round(iqr(v1)))],['WM: ', num2str(round(median(v2))), '\pm', num2str(round(iqr(v2)))]}, 'Location', 'NorthOutside')
    xlim([40 120])
    xlabel(['B_1^+'])
    ylabel('Probability density')

subplot(1,3,3)
    hold on
    low_b1_idx = find(b1_vec > 0 & b1_vec < 80);
    mid_b1_idx = find(b1_vec > 80 & b1_vec < 100);
    high_b1_idx = find(b1_vec > 100);
    
    v1 = c_vec(intersect(idx,low_b1_idx));
    [fi, xi] = ksdensity(v1); %%% make  Gaussians, so values can occur outside the range
    plot(xi,fi,'color','c','LineWidth',2);
    aah = xline(median(v1),'color','c');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    v2 = c_vec(intersect(idx,mid_b1_idx));
    [fi, xi] = ksdensity(v2);
    plot(xi,fi,'color','m','LineWidth',2);
    aah = xline(median(v2),'color','m');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    v3 = c_vec(intersect(idx,high_b1_idx));
    [fi, xi] = ksdensity(v3);
    plot(xi,fi,'color','g','LineWidth',2);
    aah = xline(median(v3),'color','g');
    aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    legend({['B_1^+ < 80: ', num2str(round(median(v1),2)), '\pm', num2str(round(iqr(v1),2))],['B_1^+ 80-100: ', num2str(round(median(v2),2)), '\pm', num2str(round(iqr(v2),2))], ['B_1^+ 100-120:',  num2str(round(median(v3),2)), '\pm', num2str(round(iqr(v3),2))]}, 'Location', 'NorthOutside')
    xlim([1.1 1.3])
    xlabel(['C'])
    ylabel('Probability density')

    set(gca, 'Color', 'white');
    saveas(f,[mt_folder,'/Figure_Tojo_Histograms.png']);
    %%
    
    
    f = figure()
    set(f,'Position',[100 100 300 500]);
    cols_to_plot = {'c','m','g'};
    %subplot(1,3,3)
    hold on
    count = 1;
    for b = 1:3
        if b == 1
            b1_idx = find(b1_vec > 0 & b1_vec < 80); b1_str = 'B_1^+ < 80';
        elseif b == 2
            b1_idx = find(b1_vec > 80 & b1_vec < 100); b1_str = 'B_1^+ 80-100'; 
        elseif b == 3
            b1_idx = find(b1_vec > 100); b1_str = 'B_1^+ 100-120';
        end
        for t = 1:2
            if t == 1
                t_idx = wm_idx; t_str = 'WM';
            elseif t == 2
                t_idx = cor_idx; t_str = 'cortex';
            end
            b1_t_idx = intersect(b1_idx, t_idx);
            v = c_vec(intersect(idx,b1_t_idx));
            [fi, xi] = ksdensity(v); %%% make  Gaussians, so values can occur outside the range
            if t == 1
                plot(xi,fi,'color',cols_to_plot{b},'LineWidth',2);
            elseif t == 2
                aah = plot(xi,fi,[cols_to_plot{b},'--'],'LineWidth',2); 
                %aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            aah = xline(median(v),'color',cols_to_plot{b});
            aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
            legentr{count} = [t_str,' ',b1_str,': ', num2str(round(median(v),2)), '\pm', num2str(round(iqr(v),2))]
            count = count + 1;
        end
    end

    % legend({['B1^+ < 80: ', num2str(round(median(v1),2)), '\pm', num2str(round(iqr(v1),2))],['B1^+ 80-100: ', num2str(round(median(v2),2)), '\pm', num2str(round(iqr(v2),2))], ['B1^+ 100-120:',  num2str(round(median(v3),2)), '\pm', num2str(round(iqr(v3),2))]}, 'Location', 'NorthOutside')
    legend(legentr, 'Location', 'NorthOutside')
    xlim([1.1 1.3])
    xlabel(['C'])
    ylabel('Probability density')
     saveas(f,[mt_folder,'/Figure_Tojo_interaction_histograms.png']);
    

f = figure()
 set(f,'Position',[100 100 300 500]);
    subplot(2,1,1)  
        binscatter(b1_vec(cor_idx), c_vec(cor_idx), 'XLimits', [60 110], 'YLimits', [1.1 1.25])
        xlabel('B_1^+')
        ylabel('C')
        title('cortex')
    subplot(2,1,2)
        binscatter(b1_vec(wm_idx), c_vec(wm_idx), 'XLimits', [60 110], 'YLimits', [1.1 1.25])
        xlabel('B_1^+')
        ylabel('C')
        title('WM')
       set(gca, 'Color', 'white');
       saveas(f,[mt_folder,'/Figure_Tojo_HistScatters.png']);
