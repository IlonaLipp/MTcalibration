clear all
close all

addpath(genpath('~/Documents/scripts/postmortembrain-mpm/'))

outpath = '/data/pt_02101/results/mt_calibration_7t';
system(['mkdir ',outpath]);

warning('off'); %%% will otherwise take longer

brain_id_slash_scan_ids = {'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB'...
    '018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB'...
    '025_C_W_RAVEL_TAI_S/mr/200923_Terra_7T_32Ch_WB'...
    '026_C_W_OSCAR_TAI_S/mr/200921_Terra_7T_32Ch_WB'...
    '032_C_C_SAMSON_SL_TAC/mr/210206_Terra_7T_32Ch_WB'};

%%% careful, coords are from fsleyes, with 1 added because fsleyes starts
%%% to count at 0, so when checking in eyes, subtract a 1
caudate_roi_coords = {[24 28 21], [21 20 22], [20 24 21], [27 16 18], [24 32 20]}; %%% hoima tojo ravel oscar samson
splenium_roi_coords = {[21 26 22], [20 24 23], [17 27 20], [25 17 20], [22 27 25]}; %%% hoima tojo ravel oscar samson
lowb1_roi_coords = {[45 45 24],[32 7 23],[32 8 20],[25 4 15],[39 48 21]} %%% hoima tojo ravel oscar samson

brain_names = {'brain 1','brain 2','brain 3', 'brain 4', 'brain 5'}; %, 'brain 6', 'brain 7'};
global_C = 1.2;
nominal_mt_pulse_to_calibrate_to_in_rad = deg2rad(700);

clear all_A all_B  all_C all_r2 all_b1 all_MTsat_reference all_MTsat_reference_corr all_MTsat_reference_const_corr all_R2s all_R1 all_bubbles 
clear CV_vectors_globalC CV_vectors_localC CV_mtpulse mt_avg_curve mt_pulse_used CV_fromTV170 CV_fromTVreferences

variation = 'MTvariation';

workspacefile = [outpath,'/',variation,'_workspace.mat'];
%delete(workspacefile);

if exist(workspacefile) == 0

    for scan_to_analyse = 1:length(brain_id_slash_scan_ids) 

        clear mtmaps r1maps b1maps ft_mt_pulse b1map_simulated mtavg CV_localC CV_globalC 

        brain_id_slash_scan_id = brain_id_slash_scan_ids{scan_to_analyse};
        brain_id_underscore_scan_id = strrep(brain_id_slash_scan_id,'/','_');
        ind_outpath = [outpath,'/',brain_id_underscore_scan_id];

        %% TEMPERATURE AND TEMPORAL CONFOUND
        %%% get temperature for scanning session
        processed_data_folder = ['/data/pt_02101/preprocessed/',brain_id_slash_scan_id,'_',variation,'_V2'];
        if exist(processed_data_folder) ~= 0
            param = load([processed_data_folder,'/param_2p1.mat']);
            confound = param.param.folder_series_numbers.mtw; %%% when series was acquired respective to others
            tempfile = dir([processed_data_folder,'/*_temperature_trace.mat']);
            if length(tempfile) > 0
                tempinfo = load([tempfile(1).folder,'/',tempfile(1).name]);
                %%% filter series names by calibration series
                if strcmp(variation,'TVvariation')
                   filter = '_TV'; 
                elseif strcmp(variation,'MTvariation')
                   filter = 'pulse';
                end
                temp_vec = [];
                for s = 1:length(tempinfo.temp_per_seq)
                    if contains(tempinfo.temp_per_seq(s).seq, filter)
                        temp_vec = [temp_vec,tempinfo.temp_per_seq(s).temp];
                    end
                end
                all.min_temp(scan_to_analyse) = min(temp_vec);
                all.max_temp(scan_to_analyse) = max(temp_vec);
            end
        else
            confound = NaN;
        end

        %% BRAIN-SPECIFIC settings
        %%% set mt pulse strengths that were used
        if strcmp(brain_id_slash_scan_id,'008_C_W_HOIMA-1/mr/200624_Terra_7T_32Ch_WB')
            fa_variation_pulse_strengths_in_rad = deg2rad([200:100:700]'); %%% these are the MT amplitudes used 
        elseif strcmp(brain_id_slash_scan_id,'007_C_C_NEGRA_ID/mr/200708_Terra_7T_32Ch_WB')...
                || strcmp(brain_id_slash_scan_id,'016_C_C_ROSIE/mr/200715_Terra_7T_32Ch_WB')
            fa_variation_pulse_strengths_in_rad = deg2rad([210:40:730]'); %%% these are the MT amplitudes used 
        elseif strcmp(brain_id_slash_scan_id,'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB')
            fa_variation_pulse_strengths_in_rad = deg2rad([200:20:700]');
        else
            fa_variation_pulse_strengths_in_rad = deg2rad([200:20:760]'); %%% these are the MT amplitudes used 
        end

        no_runs = length(fa_variation_pulse_strengths_in_rad);
        refmaps = find(fa_variation_pulse_strengths_in_rad == nominal_mt_pulse_to_calibrate_to_in_rad);

        %% read data
        mergedfile = [outpath,'/',brain_id_underscore_scan_id,'_',variation,'_data.nii'];
        for run = 1:no_runs
           filename = [processed_data_folder,'/MPM_calc/MTsat_real_TE0_WLS_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'];
           filename_r1 = [processed_data_folder,'/MPM_calc/R1_TE0_WLS_b1uncorr_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'];

           if exist(filename) ~= 0
               dataexist = 1;
               %% read mtsat
               mtvol = spm_vol(filename);
               mtmap = spm_read_vols(mtvol);
               mtmaps(:,:,:,run) = mtmap; %%% concatenated MTsat maps calculated from these pulse strengths
               r1maps(:,:,:,run) = spm_read_vols(spm_vol(filename_r1));
               %% create simulated b1+ based on pulse strength
               ft_mt_pulse(:,:,:,run) = (fa_variation_pulse_strengths_in_rad(run) / nominal_mt_pulse_to_calibrate_to_in_rad) * ones(size(mtmap)); %%% one value per volume
               mtvals = mtmap(:);
               mtavg(run) = mean(mtvals(mtvals~=0));
               %% bubble map, r2* and r1
               if run == refmaps
                  TEMPMAPS.R2s = spm_read_vols(spm_vol([processed_data_folder,'/MPM_calc/R2s_WLS_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'])); 
                  r2s_ols_r2 = spm_read_vols(spm_vol([processed_data_folder,'/MPM_calc/OLSfit_R2_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'])); 
                  TEMPMAPS.bubbles = r2s_ols_r2 < 0.99;
                  TEMPMAPS.R1 = spm_read_vols(spm_vol([processed_data_folder,'/MPM_calc/R1_TE0_WLS_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'])); 
                  TEMPMAPS.voxelmask = (TEMPMAPS.R2s > 0 & TEMPMAPS.bubbles == 0);
               end
           else
               dataexist = 0;
           end

        end

        if dataexist == 1

            %%% save mtmaps as merged
            save_4d_nifti(mtmaps, mergedfile, mtvol, outpath);
            save_4d_nifti(r1maps, [mergedfile(1:end-4),'_r1tocheck.nii'], mtvol, outpath);

            %%% get b1+ map
            b1file = [processed_data_folder,'/B1mapping/B1map_AL_processed_2p1_mean.nii'];
            b1map = spm_read_vols(spm_vol(b1file));
            ft_brain = b1map / 100;
            result_map_file = [outpath,'/',brain_id_underscore_scan_id,'_',variation,'_resultmaps.mat'];
            delete(result_map_file)
            if exist(result_map_file) == 2
                load(result_map_file);
            else    
                MAPS = calculate_calibration_constants_MT_variation_linear_attempt2(mtmaps, ft_brain, fa_variation_pulse_strengths_in_rad, nominal_mt_pulse_to_calibrate_to_in_rad);
                save(result_map_file,'MAPS');
            end

            %% plot (also writes out individual datapoints)
            if dataexist == 1
               figname = [outpath,'/',brain_id_underscore_scan_id,'_',variation,'_datapointplot.jpg'];
               textfileprefix = [outpath,'/',brain_id_underscore_scan_id,'_',variation,'_datapoints'];
               if exist(figname) == 0
                make_data_point_plot_LJE(confound, mtmaps, ft_brain, fa_variation_pulse_strengths_in_rad, nominal_mt_pulse_to_calibrate_to_in_rad, caudate_roi_coords{scan_to_analyse}, splenium_roi_coords{scan_to_analyse}, figname, textfileprefix)
               end
            end


           MAPS.ft_reference = b1map/100;
           MAPS.MTsat_reference = mtmaps(:,:,:,refmaps);

           %%% get individually global values
           C_vec = MAPS.C(:);
           ind_global_C = round(mean(C_vec(TEMPMAPS.voxelmask(:)~=0)), 2); 
           save([outpath,'/',brain_id_underscore_scan_id,'_',variation,'_individual_calibration_constants'], 'ind_global_C');
           %% convert c range to percent
           MAPS.C_range_percent = 100 * MAPS.C_range ./ MAPS.C;
           %% correct map different ways and save everything
           [newmap, corrfac] = b1_correct_mt_map_linear(MAPS.MTsat_reference, b1map, MAPS.C); %%% voxelwise correction
           MAPS.MTsat_reference_corr = newmap; %%% save corrected map
           MAPS.MTsat_reference_corr_factor = corrfac; %%% also save correction factor
           [newmap, corrfac] = b1_correct_mt_map_linear(MAPS.MTsat_reference, b1map, global_C); %%% correction with 1.2
           MAPS.MTsat_reference_const_corr_groupopt = newmap;
           MAPS.MTsat_reference_const_corr_groupopt_factor = corrfac;
           [newmap, corrfac] = b1_correct_mt_map_linear(MAPS.MTsat_reference, b1map, ind_global_C); %%% correction with individual value
           MAPS.MTsat_reference_const_corr_indopt = newmap;
           MAPS.MTsat_reference_const_corr_indopt_factor = corrfac;
           MAPS.C_in_perc_from_mean = 100 * (MAPS.C - ind_global_C) / ind_global_C;
           MAPS.C_SE_in_perc = 100 * MAPS.C_SE ./ MAPS.C;

           var_names = fieldnames(TEMPMAPS);
            for var = 1:length(var_names);
                var_name = var_names{var};
               MAPS.(var_name) = TEMPMAPS.(var_name);
            end
            %% save all maps as vector for later stats
            var_names = fieldnames(MAPS);
            for var = 1:length(var_names)
                varname = var_names{var};
                %%% when doing so 
                all.(varname)(:,scan_to_analyse) = MAPS.(varname)(:);
            end

            ind_outpath = [outpath,'/',brain_id_underscore_scan_id];
            system(['mkdir ', ind_outpath]);
            save_and_brain_mask_files_from_MAPS_structure(MAPS, mtvol, [processed_data_folder,'/MPM_calc/brain_mask_2p1_run01.nii'], ind_outpath, [variation,'_variation_',brain_id_underscore_scan_id], 1, '2p1')

        end
    end
    save(workspacefile);
else
    load(workspacefile);
end

%% group analysis
MT_calibration_experiment_figures
