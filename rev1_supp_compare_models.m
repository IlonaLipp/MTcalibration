% Supplementary figure showing the comparison of the new model intended for
% high MT pulse flip angles (Lipp, et al. 2022, 
% https://www.biorxiv.org/content/10.1101/2022.07.12.498197v1.abstract) 
% with the Helms model which is intended for lower MT pulse flip angles 
% (https://arxiv.org/abs/2104.14878)
function rev1_supp_compare_models

outpath = './supp';
[~,~] = mkdir(outpath);

brain_id_slash_scan_ids = {'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB'...
    '018_C_C_TOJO/mr/200922_Terra_7T_32Ch_WB'...
    '025_C_W_RAVEL_TAI_S/mr/200923_Terra_7T_32Ch_WB'...
    '026_C_W_OSCAR_TAI_S/mr/200921_Terra_7T_32Ch_WB'...
    '032_C_C_SAMSON_SL_TAC/mr/210206_Terra_7T_32Ch_WB'};

% careful, coords are from fsleyes, with 1 added because fsleyes starts
% to count at 0, so when checking in fsleyes, subtract a 1
caudate_roi_coords = {[24 28 21], [21 20 22], [20 24 21], [27 16 18], [24 32 20]}; %%% hoima tojo ravel oscar samson
splenium_roi_coords = {[21 26 22], [20 24 23], [17 27 20], [25 17 20], [22 27 25]}; %%% hoima tojo ravel oscar samson

nominal_mt_pulse_to_calibrate_to_in_rad = deg2rad(700);

variation = 'MTvariation';

% uncomment to see results in other brains
for scan_to_analyse = 2 %1:length(brain_id_slash_scan_ids)

    brain_id_slash_scan_id = brain_id_slash_scan_ids{scan_to_analyse};
    processed_data_folder = ['/data/pt_02101/preprocessed/',brain_id_slash_scan_id,'_',variation,'_V2'];

    %% BRAIN-SPECIFIC settings
    % set mt pulse strengths that were used
    if strcmp(brain_id_slash_scan_id,'008_C_W_HOIMA-1/mr/200624_Terra_7T_32Ch_WB')
        fa_variation_pulse_strengths_in_rad = deg2rad((200:100:700)'); %%% these are the MT amplitudes used
    elseif strcmp(brain_id_slash_scan_id,'007_C_C_NEGRA_ID/mr/200708_Terra_7T_32Ch_WB')...
            || strcmp(brain_id_slash_scan_id,'016_C_C_ROSIE/mr/200715_Terra_7T_32Ch_WB')
        fa_variation_pulse_strengths_in_rad = deg2rad((210:40:730)'); %%% these are the MT amplitudes used
    elseif strcmp(brain_id_slash_scan_id,'009_C_W_HOIMA-2/mr/201112_Terra_7T_32Ch_WB')
        fa_variation_pulse_strengths_in_rad = deg2rad((200:20:700)');
    else
        fa_variation_pulse_strengths_in_rad = deg2rad((200:20:760)'); %%% these are the MT amplitudes used
    end

    n_runs = length(fa_variation_pulse_strengths_in_rad);

    %% read data
    dataexist = true;
    for run = 1:n_runs
        filename = [processed_data_folder,'/MPM_calc/MTsat_real_TE0_WLS_processed_s_u_2p1_run',sprintf('%.02d',run),'_brain_masked.nii'];

        if exist(filename,"file")
            % read mtsat
            mtvol = spm_vol(filename);
            mtmap = spm_read_vols(mtvol);

            % concatenated MTsat maps at these pulse strengths
            mtmaps(:,:,:,run) = mtmap; %#ok<AGROW>
        else
            dataexist = false;
            break
        end

    end

    %% fit and plot data
    if dataexist

        % get b1+ map
        b1file = [processed_data_folder,'/B1mapping/B1map_AL_processed_2p1_mean.nii'];
        b1map = spm_read_vols(spm_vol(b1file));
        ft_brain = b1map / 100;

        test_models_LJE(mtmaps, ft_brain, fa_variation_pulse_strengths_in_rad, nominal_mt_pulse_to_calibrate_to_in_rad, caudate_roi_coords{scan_to_analyse}, splenium_roi_coords{scan_to_analyse});

    end
end
end

function test_models_LJE(mtmaps, fT, fa_variation_pulse_strengths_in_rad, nominal_mt_angle_in_rad, caudate_roi_coords, splenium_roi_coords)

roinames={'splenium','caudate'};
ROIs{1}=splenium_roi_coords;
ROIs{2}=caudate_roi_coords;

arange=[min(fT(fT>0)),max(fT(fT>0))]*rad2deg(nominal_mt_angle_in_rad);
ylimits=[0,15];

figure('Units','inches','Position',[1,5,16.5,4.6])
colors=get(gca,'colororder');

apredict=linspace(0,1e3); % deg

hold on
% show range of flip angles in sample
rectangle('Position',[arange(1),ylimits(1),diff(arange),diff(ylimits)],'FaceColor',[0.9*[1 1 1],0.5],'LineStyle','none')
plot(rad2deg(nominal_mt_angle_in_rad)*[1,1],ylimits,'k--','LineWidth',2)

for ROIidx=1:length(ROIs)
    ROI=ROIs{ROIidx};
    aobs=rad2deg(fa_variation_pulse_strengths_in_rad)*fT(ROI(1),ROI(2),ROI(3));
    MTobs=squeeze(mtmaps(ROI(1),ROI(2),ROI(3),:));
    plot(aobs,MTobs,'x','Color','k','MarkerSize',15,'LineWidth',2); %,'Color',colors(3+3*(ROIidx-1),:));
    text(500,MTobs(round(end/2))+0.5,roinames{ROIidx},"HorizontalAlignment","right","VerticalAlignment","bottom","FontSize",18)

    [~,helmspred]=helms(aobs,MTobs,apredict);
    plot(apredict,helmspred,"LineWidth",2,"Color",colors(1+2*(ROIidx-1),:));
    text(apredict(end)-10,helmspred(end),"Helms' model","HorizontalAlignment","right","VerticalAlignment","top","Color",colors(1+2*(ROIidx-1),:),"FontSize",18)

    [~,paperpred]=paper(aobs,MTobs,rad2deg(nominal_mt_angle_in_rad),apredict);
    plot(apredict,paperpred,"LineWidth",2,"Color",colors(2+2*(ROIidx-1),:));
    text(apredict(end)-10,paperpred(end),"Lipp's model","HorizontalAlignment","right","VerticalAlignment","bottom","Color",colors(2+2*(ROIidx-1),:),"FontSize",18)

end
xlabel('\beta_{loc} / degrees')
ylabel('MTsat / p.u.')
ylim(ylimits)
xlim([0,1000])
plot(xlim,[0,0],'k') % add black line for x axis as rectangles overlap it...

set(gca,"FontSize",18)

exportgraphics(gcf,"supp/fitcomparison.png",'Resolution',300)

end

function [fitted,prediction] = helms(beta,MTobs,beta0)

pars0 = [0, 0]; %%% starting point
fitfun = @(cpars,beta) beta.^2.*cpars(1).*(1 + beta * cpars(2)); % cpars(1) = A_0, cpars(2) = A_1
fitted = nlinfit(beta,MTobs,fitfun,pars0);
prediction = fitfun(fitted, beta0); %

end

function [fitted,prediction] = paper(beta,MTobs,betaref,beta0)

pars0 = [0, 0]; %%% starting point
fitfun = @(cpars,beta)(cpars(1) + beta * cpars(1) * cpars(2)); %%% cpars(1) = delta_corr, cpars(2) = A
fitted = nlinfit(beta-betaref,MTobs,fitfun,pars0);
prediction = fitfun(fitted, beta0-betaref);

end

