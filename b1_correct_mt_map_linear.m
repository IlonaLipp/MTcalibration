function [newmap, corrfac] = b1_correct_mt_map_linear(mtmap, b1map, C)
%%% C = calibration constant C from Luke and Evgeniya's approach
%%% mt map is map to be calibrated
%%% b1 map in %
%%% newmap = b1 corrected mtmap
    map_of_ones = ones(size(mtmap));
    ft = b1map/100;
    bottom = 1 + (ft - map_of_ones).*C;
    corrfac = map_of_ones ./ bottom; %%% this is for not using the small angle approx
    newmap = mtmap .* corrfac; 
end
