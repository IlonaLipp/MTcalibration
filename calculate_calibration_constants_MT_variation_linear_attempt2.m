function MAPS = calculate_calibration_constants_MT_variation_linear_attempt2(mtmaps, ft_brain, fa_variation_pulse_strengths_in_rad, nominal_mt_pulse_to_calibrate_to_in_rad, confound, useconfound)

    mtmap = mtmaps(:,:,:,1);
    %%% initialise output maps
    MAPS.n = zeros(size(mtmap));
    MAPS.R2 = zeros(size(mtmap));
    MAPS.A = zeros(size(mtmap));
    MAPS.C = zeros(size(mtmap));
    MAPS.C_SE = zeros(size(mtmap));
    MAPS.C_CI_low = zeros(size(mtmap));
    MAPS.C_CI_high = zeros(size(mtmap));
    
    for x = 1:size(mtmap,1)
       for y = 1:size(mtmap,2)
           for z = 1:size(mtmap,3)
               mt_data = squeeze(mtmaps(x,y,z,:));

               
               %%% identify data points to include:
               alpha_sat_all = fa_variation_pulse_strengths_in_rad;
               alpha_local_all = alpha_sat_all * ft_brain(x,y,z);
               idx = find(alpha_local_all > deg2rad(220) & mt_data > 0);
               
               %length(idx)
               if length(idx) > 2 %%% only consider voxels in which both MT and b1 are realistic values
                   %%% regress out confound from data
                   alpha_sat = fa_variation_pulse_strengths_in_rad(idx); %%% in radians
                   local_ft_brain = ft_brain(x,y,z);
                   alpha_local = alpha_sat * local_ft_brain; %%% in radians
                   
                   dv = mt_data(idx); 
                   iv = (alpha_local - nominal_mt_pulse_to_calibrate_to_in_rad);
          
                   linear_fit_to_insert
                   
                   
               end
               MAPS.n(x,y,z) = length(idx);
           end
       end
    end
   MAPS.C_range = MAPS.C_CI_high - MAPS.C_CI_low;

end
