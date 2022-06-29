close all
clear all
clc

variation = 'MTvariation'
nominal_mt_pulse_to_calibrate_to_in_rad = deg2rad(700);
figc = -2;
useconfound = 0
outpath = ['/data/pt_02101/results/mt_calibration_7t']
plotcols = {'r', 'g', 'b'}

warning('off')
regions = {'caudate','splenium'}; %,'low b1+'};

brain_id_slash_scan_ids = {'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB'...
    '018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB'...
    '025_C_W_RAVEL_TAI_S/mr/200923_Terra_7T_32Ch_WB'...
    '026_C_W_OSCAR_TAI_S/mr/200921_Terra_7T_32Ch_WB'...
    '032_C_C_SAMSON_SL_TAC/mr/210206_Terra_7T_32Ch_WB'};
nominal_mt_pulse_to_calibrate_to_in_rad = deg2rad(700);
colors=colororder;
colors_to_use = parula(length(brain_id_slash_scan_ids));

close all

for b = 1:length(brain_id_slash_scan_ids)

    brain_id_slash_scan_id = brain_id_slash_scan_ids{b};
    brain_id_underscore_scan_id = strrep(brain_id_slash_scan_id,'/','_');

   for r = 2 %2 %1:length(regions)

        region = regions{r};
        figc = figc + 3;

        filename = (['/data/pt_02101/results/mt_calibration_7t/',brain_id_underscore_scan_id,'_',variation,'_datapoints_',region,'.csv']);

        if exist(filename) == 2

            %%% get data
            dataimport = importdata(filename);
            datamat = dataimport.data;

            %%% define variables
           mt_data_orig = datamat(:,2);
           confound = datamat(:,1);
           mtsatloc = mt_data_orig;
           local_ft_brain = datamat(:,4); ft = local_ft_brain(1);
           alpha_sat = datamat(:,3);
           alpha_local = alpha_sat .* local_ft_brain;

            %%% threshold for modeling
            idx = find(mtsatloc > 0);

            alpha_sat_for_model = alpha_sat(idx);
            alpha_local_for_model = alpha_local(idx);
            mtsatloc_for_model = mtsatloc(idx);

            dv = mtsatloc; % .* ft.^2;
            dv_for_model = dv(idx);
            iv = (alpha_local_for_model - nominal_mt_pulse_to_calibrate_to_in_rad);

            %% fit model
            pars0 = [0, 0]; %%% starting point
            newfun = @(cpars,iv)(cpars(1) + iv * cpars(1) * cpars(2)); %%% cpars(1) = delta_corr, cpars(2) = A
            [BETA,R,J,COVB,MSE] = nlinfit(iv,dv_for_model,newfun,pars0); 
            prediction = newfun(BETA, iv); %

             %% plot 1, fit for all
            fh = figure(100);
            set(fh,'Position',[100 100 500 300]);
            hold on;
            set(0,'DefaultAxesTitleFontWeight','normal');

            aah = plot(alpha_local_for_model./nominal_mt_pulse_to_calibrate_to_in_rad,dv./BETA(1),'Color',colors_to_use(b,:),'LineWidth',1.2);
            aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
            plot(alpha_local_for_model./nominal_mt_pulse_to_calibrate_to_in_rad,prediction./BETA(1),'--','Color',colors_to_use(b,:),'LineWidth',1.2);


         end

   end

end
xlabel('\beta_{local} / \beta_{ref}')
ylabel('\delta_{MT}(\beta_{loc}) / \delta_{MT}(\beta_{ref})')
ylabel('MTsat(\beta_{loc}) / MTsat(\beta_{ref})')
legend({'brain 1','brain 2','brain 3','brain 4','brain 5'},'Location','NorthWestOutside')
hold on
aah = plot([1 1],[0 1],'k--'); aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
aah = plot([1 0],[1 1],'k--'); aah.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([0.2 1.2])
ylim([0 1.2])
set(gca, 'Color', 'white');
set(gca,'FontSize',12)
saveas(fh,[outpath,'/Figure_',variation,'_',region,'_example_fits_LJE.png']);



