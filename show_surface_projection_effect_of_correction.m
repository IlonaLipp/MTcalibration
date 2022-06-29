clear all
close all
clc

addpath(genpath('/data/tu_lippi/Software/BrewerMap'))

parameters = {'B1', 'MTsat_uncorr', 'MTsat'};
display_ranges = {[60 140],[3.5 5], [3.5 5]};

hemispheres = {'lh','rh'};
resultdir = '/data/pt_02101/results/mt_calibration_7t'

for h = 1:length(hemispheres)
    
    hemisphere = hemispheres{h}
    freesurfer_folder = '/data/pt_02101/preprocessed/018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB_V2/freesurfer/MTsat_0p3mm_downsampled_to_0p7mm_FSV7_hires';
    surface_file = [freesurfer_folder, '/surf/', hemisphere, '.inflated']

    for p = 1:length(parameters)

        parameter = parameters{p};
        display_range = display_ranges{p};

        projection_file = [freesurfer_folder, '/SurfaceProjections/', parameter, '_0p3_average_',hemisphere,'.mgh'];
        [vol, M, mr_parms, volsz] = load_mgh(projection_file);

        gg.cdata = single(vol); % + 0.01*rand(length(val_vec)));
        gg_newmat = gifti(surface_file); 

        rotation_angles = [90 -90];
        for ra = 1:length(rotation_angles)
            rotation_angle = rotation_angles(ra);

            colormap(brewermap(256, '*Spectral'));
            figure, plot(gg_newmat,gg); 
            caxis([display_range(1), display_range(end)]);
            %title([effects{e},' ',dstr])
            view([rotation_angle 0])

            hfig = light; set(hfig, 'position', [1 1 0.2]); lighting gouraud; material dull
            colormap(brewermap(256, '*Spectral'));
            hfig = colorbar();

            set(hfig,'XTick',display_range,'FontSize',10);
            ylabel(hfig,parameter,'FontSize',11,'interpreter','none');
            set(gcf,'color','w');

           figure_file_to_save_prefix = [resultdir,'/Figure_Tojo_',hemisphere,'_',parameter,'_brewermap'];
            figure_file_to_save = [figure_file_to_save_prefix,'_',num2str(ra),'.jpg'];
           saveas(gcf,figure_file_to_save);
        end
    end
end